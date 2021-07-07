function reading_gtec_EKSO_mode(file_path,SUBJ)

%% include the filepath to the EKSO file you are provided and the suffix for saving the final results in SUBJ
close all;
rmpath(genpath('/path_to_EEGLab/eeglab14_1_1b'));
DATA=load(file_path);
data_accel=DATA.data(:,[9:11 15:17]);
for k=1:size(data_accel,2)
    data_val(k,:)=filtsig(data_accel(:,k),0,2e3,500,0,1.5);
    data_val(k,:)=data_val(k,:)./max(data_val(k,:));
end;
DATA_ref=mean(data_val(4:5,:),1)./max(mean(data_val(4:5,:),1));
t=linspace(0,length(DATA_ref)*(1/2000),length(DATA_ref));
plot(t,DATA_ref,'LineWidth',2);
hold on;
[point_a,time_a]=findpeaks(DATA_ref,t,'MinPeakDistance',1.5,'MinPeakHeight',0.8);
plot(time_a,point_a,'x','MarkerSize',10,'LineWidth',2);
addpath(genpath('/path_to_EEGLab/eeglab14_1_1b'));
rmpath(genpath('/path_to_EEGLab/eeglab14_1_1b/functions/octavefunc'));

%% segment the signal and process artifact rejection on each segmented trial
n=1;
for l=1:2:length(time_a)-2
    pos1=min(find(t>=time_a(l)));
    pos2=max(find(t<=time_a(l+2)));
    DATA_VAL{n}=DATA.data(pos1:pos2,45:108)';
    timet=linspace(0,size(DATA_VAL{n},2)*(1/2000),size(DATA_VAL{n},2));
    EEGdata{n}=define_structure_gtec(file_path,SUBJ,DATA_VAL{n},timet);
    EEGdata{n}=pop_resample(EEGdata{n},250);
    EEGdata{n} = pop_eegfiltnew(EEGdata{n},0.1,100,6600);
    [EEGdata{n}.icaweights,EEGdata{n}.icasphere]=runica(EEGdata{n}.data(:,:),'sphering','on','lrate',1e-5,'maxsteps',50); 
    %% ADJUST application
    EEGdata{n}.icawinv=inv(EEGdata{n}.icaweights*EEGdata{n}.icasphere);
    EEGdata{n}.icachansind=[1:1:64];
    rmpath(genpath('/path_to_EEGLab/eeglab14_1_1b/fieldtrip-20210507'));
    EEGdata{n}=cleanline('EEG',EEGdata{n},'Bandwidth',1,'ChanCompIndices',[1:EEGdata{n}.nbchan],'SignalType','Channels','ComputeSpectralPower','true','LineFrequencies',[50]);
    addpath(genpath('/home/jmm_vivobook_asus/DeepGaze_project/eeglab14_1_1b/fieldtrip-20210507'));
    EEG_no{n}=EEGdata{n};
    EEGdata{n}=pop_autobssemg(EEGdata{n},1,1,'bsscca',{'eigratio',1e6},'emg_psd',{'ratio',10,'fs',250,'range',[1 10]});
    EEGdata{n}.data= repmat(EEGdata{n}.data(:,:),1,1,60);
    EEGdata{n}.trials=60;
    [art_channels{n}]=ADJUST(EEGdata{n},'report.txt');
    EEGdata{n}=pop_subcomp(EEGdata{n},art_channels{n});
    EEGdata{n}.data=mean(EEGdata{n}.data(:,:,:),3);
    figure;
    EEG_no_spec(n,:,:)=pop_spectopo(EEG_no{n},1,[0 max(timet)*1000],'EEG','freq',[0.1 1 2],'freqrange',[0.1 100],'electrodes','off');
    figure;
    EEG_spec(n,:,:)=pop_spectopo(EEGdata{n},1,[0 max(timet)*1000],'EEG','freq',[0.1 1 2],'freqrange',[0.1 100],'electrodes','off');
    close all;
    n=n+1;
end;

rmpath(genpath('/path_to_EEGLab/eeglab14_1_1b'));
for i=1:size(EEG_spec,1)
    for l=1:size(EEG_spec,2)
        r_spec(i,l)=abs(thd(squeeze(EEG_spec(i,l,1:101))',linspace(0,100,101)',3,'psd'));
        r_no_spec(i,l)=abs(thd(squeeze(EEG_no_spec(i,l,1:101))',linspace(0,100,101)',3,'psd'));
    end;
end;

index=r_spec==Inf;
r_spec(index==1)=0;
index=r_spec==-Inf;
r_spec(index==1)=0;
index=r_no_spec==Inf;
r_no_spec(index==1)=0;
index=r_no_spec==-Inf;
r_no_spec(index==1)=0;
Val_no=nanmean(r_no_spec,2)';
Val_spec=nanmean(r_spec,2)';

%% plotting and statistical evaluation
notBoxPlot([Val_no Val_spec],[zeros([1,length(Val_no)]),ones([1,length(Val_spec)])],'style','sdline');
anova1([Val_no Val_spec],[zeros([1,length(Val_no)]),ones([1,length(Val_spec)])])
save(['res_gtec_clean_' SUBJ '.mat' ],'r_spec','r_no_spec','Val_no','Val_spec','EEG_spec','EEG_no_spec','EEGdata','EEG_no');
A=1;
