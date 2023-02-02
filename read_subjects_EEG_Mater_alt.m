function [accsal,EEGsal,ind,ind_bouts,EMG_sal]=read_subjects_EEG_Mater_alt(subj_ind,Run_ind,val_in,n_cycles,plot_sel,sel_eval_bouts,emg_sel,size_around,EMG)
close all
ind=0;
absolute_path_EEG=[pwd '/EEG'];
%absolute_path_EEG='/home/jmm_vivobook_asus/DeepGaze_project/Mater_EEG_EMG_dataset/EEG';
%% load the previous files containing the samples alignment
load('ACC_sync.mat');
load('EEG_gait_sync.mat');
load('syncBouts.mat');
if emg_sel==1
    EMG_data=EMG;
end;
%% read the corresponding EEG data 
EEGdata=load([absolute_path_EEG '/' subj_ind '/' Run_ind ' (with exo)/EEG_rec.mat']);
sel=str2num(subj_ind(end));
ind_bouts=0;
if sel==1
    DataFC=FC.Subj01;
    DataIC=IC.Subj01;
    acc=acc_sync.Subj01;
    Bouts=syncBouts.Subj01.EKSO.bouts;
    if emg_sel==1
       EMG_SO_Left=resample(EMG_data.Subj01.EKSO.EMG.SO.Left,250,1910);
       EMG_ST_Left=resample(EMG_data.Subj01.EKSO.EMG.ST.Left,250,1910);
       EMG_TA_Left=resample(EMG_data.Subj01.EKSO.EMG.TA.Left,250,1910);
       EMG_RF_Left=resample(EMG_data.Subj01.EKSO.EMG.RF.Left,250,1910);
       EMG_SO_Right=resample(EMG_data.Subj01.EKSO.EMG.SO.Right,250,1910);
       EMG_ST_Right=resample(EMG_data.Subj01.EKSO.EMG.ST.Right,250,1910);
       EMG_TA_Right=resample(EMG_data.Subj01.EKSO.EMG.TA.Right,250,1910);
       EMG_RF_Right=resample(EMG_data.Subj01.EKSO.EMG.RF.Right,250,1910);
    end;
elseif sel==2
    DataFC=FC.Subj02;
    DataIC=IC.Subj02;
    acc=acc_sync.Subj02;
    Bouts=syncBouts.Subj02.EKSO.bouts;
    if emg_sel==1
       EMG_SO_Left=resample(EMG_data.Subj02.EKSO.EMG.SO.Left,250,1910);
       EMG_ST_Left=resample(EMG_data.Subj02.EKSO.EMG.ST.Left,250,1910);
       EMG_TA_Left=resample(EMG_data.Subj02.EKSO.EMG.TA.Left,250,1910);
       EMG_RF_Left=resample(EMG_data.Subj02.EKSO.EMG.RF.Left,250,1910);
       EMG_SO_Right=resample(EMG_data.Subj02.EKSO.EMG.SO.Right,250,1910);
       EMG_ST_Right=resample(EMG_data.Subj02.EKSO.EMG.ST.Right,250,1910);
       EMG_TA_Right=resample(EMG_data.Subj02.EKSO.EMG.TA.Right,250,1910);
       EMG_RF_Right=resample(EMG_data.Subj02.EKSO.EMG.RF.Right,250,1910);
    end;
elseif sel==3
    DataFC=FC.Subj03;
    DataIC=IC.Subj03;
    acc=acc_sync.Subj03;
    Bouts=syncBouts.Subj03.EKSO.bouts;
    if emg_sel==1
       EMG_SO_Left=resample(EMG_data.Subj03.EKSO.EMG.SO.Left,250,1910);
       EMG_ST_Left=resample(EMG_data.Subj03.EKSO.EMG.ST.Left,250,1910);
       EMG_TA_Left=resample(EMG_data.Subj03.EKSO.EMG.TA.Left,250,1910);
       EMG_RF_Left=resample(EMG_data.Subj03.EKSO.EMG.RF.Left,250,1910);
       EMG_SO_Right=resample(EMG_data.Subj03.EKSO.EMG.SO.Right,250,1910);
       EMG_ST_Right=resample(EMG_data.Subj03.EKSO.EMG.ST.Right,250,1910);
       EMG_TA_Right=resample(EMG_data.Subj03.EKSO.EMG.TA.Right,250,1910);
       EMG_RF_Right=resample(EMG_data.Subj03.EKSO.EMG.RF.Right,250,1910);
    end;
elseif sel==4
    DataFC=FC.Subj04;
    DataIC=IC.Subj04;
    acc=acc_sync.Subj04;
    Bouts=syncBouts.Subj04.EKSO.bouts;
    if emg_sel==1
       EMG_SO_Left=resample(EMG_data.Subj04.EKSO.EMG.SO.Left,250,1910);
       EMG_ST_Left=resample(EMG_data.Subj04.EKSO.EMG.ST.Left,250,1910);
       EMG_TA_Left=resample(EMG_data.Subj04.EKSO.EMG.TA.Left,250,1910);
       EMG_RF_Left=resample(EMG_data.Subj04.EKSO.EMG.RF.Left,250,1910);
       EMG_SO_Right=resample(EMG_data.Subj04.EKSO.EMG.SO.Right,250,1910);
       EMG_ST_Right=resample(EMG_data.Subj04.EKSO.EMG.ST.Right,250,1910);
       EMG_TA_Right=resample(EMG_data.Subj04.EKSO.EMG.TA.Right,250,1910);
       EMG_RF_Right=resample(EMG_data.Subj04.EKSO.EMG.RF.Right,250,1910);
    end;
elseif sel==5
    DataFC=FC.Subj05;
    DataIC=IC.Subj05;
    acc=acc_sync.Subj05;
    Bouts=syncBouts.Subj05.EKSO.bouts;
    if emg_sel==1
       EMG_SO_Left=resample(EMG_data.Subj05.EKSO.EMG.SO.Left,250,1910);
       EMG_ST_Left=resample(EMG_data.Subj05.EKSO.EMG.ST.Left,250,1910);
       EMG_TA_Left=resample(EMG_data.Subj05.EKSO.EMG.TA.Left,250,1910);
       EMG_RF_Left=resample(EMG_data.Subj05.EKSO.EMG.RF.Left,250,1910);
       EMG_SO_Right=resample(EMG_data.Subj05.EKSO.EMG.SO.Right,250,1910);
       EMG_ST_Right=resample(EMG_data.Subj05.EKSO.EMG.ST.Right,250,1910);
       EMG_TA_Right=resample(EMG_data.Subj05.EKSO.EMG.TA.Right,250,1910);
       EMG_RF_Right=resample(EMG_data.Subj05.EKSO.EMG.RF.Right,250,1910);
    end;
else
    DataFC=FC.Subj06;
    DataIC=IC.Subj06;
    acc=acc_sync.Subj06;
    Bouts=syncBouts.Subj06.EKSO.bouts;
    if emg_sel==1
       EMG_SO_Left=resample(EMG_data.Subj06.EKSO.EMG.SO.Left,250,1910);
       EMG_ST_Left=resample(EMG_data.Subj06.EKSO.EMG.ST.Left,250,1910);
       EMG_TA_Left=resample(EMG_data.Subj06.EKSO.EMG.TA.Left,250,1910);
       EMG_RF_Left=resample(EMG_data.Subj06.EKSO.EMG.RF.Left,250,1910);
       EMG_SO_Right=resample(EMG_data.Subj06.EKSO.EMG.SO.Right,250,1910);
       EMG_ST_Right=resample(EMG_data.Subj06.EKSO.EMG.ST.Right,250,1910);
       EMG_TA_Right=resample(EMG_data.Subj06.EKSO.EMG.TA.Right,250,1910);
       EMG_RF_Right=resample(EMG_data.Subj06.EKSO.EMG.RF.Right,250,1910);
    end;
end;
sel_run=str2num(Run_ind(end));
%% check the attribute names length
names_val=fieldnames(DataFC);
names_val_other=fieldnames(DataIC);
if length(names_val{1})==5
    if sel_run==1
        DataFC=DataFC.Run01;
        if length(names_val{1})==length(names_val_other{1})
           DataIC=DataIC.Run01;
           acc=acc.Run01;
        else
           DataIC=DataIC.Run01_with_exo;
           acc=acc.Run01_with_exo;
        end;
    elseif sel_run==2
        DataFC=DataFC.Run02;
        DataIC=DataIC.Run02;
        acc=acc.Run02;
    elseif sel_run==3
        DataFC=DataFC.Run03;
        DataIC=DataIC.Run03;
        acc=acc.Run03;
    elseif sel_run==4
        DataFC=DataFC.Run04;
        DataIC=DataIC.Run04;
        acc=acc.Run04;
    end;
else
    if sel_run==1
        DataFC=DataFC.Run3;
        DataIC=DataIC.Run3;
        acc=acc.Run3;
    elseif sel_run==2
        DataFC=DataFC.Run4;
        DataIC=DataIC.Run4;
        acc=acc.Run4;
    elseif sel_run==3
        DataFC=DataFC.Run3;
        DataIC=DataIC.Run3;
        acc=acc.Run3;
    elseif sel_run==4
        DataFC=DataFC.Run4;
        DataIC=DataIC.Run4;
        acc=acc.Run4;
    end;
end;
if length(DataFC.Left)<=length(DataFC.Right)
   Datad=DataFC.Left;
else
   Datad=DataFC.Right;
end;
if val_in<length(Datad)-1
    t1l=DataFC.Left(val_in);
    t2l=DataFC.Left(val_in+n_cycles); %+n_cycles+1);
    t1r=DataFC.Right(val_in);
    t2r=DataFC.Right(val_in+n_cycles); %+n_cycles+1);
    t1hl=DataFC.Left(val_in);
    t2hl=DataIC.Left(val_in); %+n_cycles+1);
    t1hr=DataFC.Right(val_in);
    t2hr=DataIC.Right(val_in); %+n_cycles+1);
    if t1l>=t2l
        t2l=DataFC.Left(val_in+2);
    end;
    if t1r>=t2r
        t2r=DataFC.Right(val_in+2);
    end;
    if t1hl>=t2hl
        t2hl=DataIC.Left(val_in+1);
    end;
    if t1hr>=t2hr
        t2hr=DataIC.Right(val_in+1);
    end;
    %%% definition Bouts
    if sel_eval_bouts==1
        for l=1:size(Bouts,2)
            if (abs(t1l-Bouts(1,l))/250<=size_around || abs(t2l-Bouts(1,l))/250<=size_around || abs(t1r-Bouts(2,l))/250<=size_around || abs(t2r-Bouts(2,l))/250<=size_around)
                ind_bouts=1;
                break;
            end;
        end;
    end;
    accelL=acc.LShank(t1l:t2l,:);
    accelR=acc.RShank(t1r:t2r,:);
    accelLT=acc.LThigh(t1l:t2l,:);
    accelRT=acc.RThigh(t1r:t2r,:);
    EEGL=EEGdata.EEG_rec(2:33,t1l-5:t2l+5)';
    EEGR=EEGdata.EEG_rec(2:33,t1r-5:t2r+5)';
    [b,a]=butter(5,5/(250/2),'high');
    [d,c]=butter(5,120/(250/2),'low');
    if emg_sel==1
        EMG_SO_G_L=EMG_SO_Left(t1l-5:t2l+5);
        EMG_SO_G_R=EMG_SO_Right(t1r-5:t2r+5);
        EMG_ST_G_L=EMG_ST_Left(t1l-5:t2l+5);
        EMG_ST_G_R=EMG_ST_Right(t1r-5:t2r+5);
        EMG_TA_G_L=EMG_TA_Left(t1l-5:t2l+5);
        EMG_TA_G_R=EMG_TA_Right(t1r-5:t2r+5);
        EMG_RF_G_L=EMG_RF_Left(t1l-5:t2l+5);
        EMG_RF_G_R=EMG_RF_Right(t1r-5:t2r+5);
        EMG_SO_G_L(isnan(EMG_SO_G_L))=0;
        EMG_sal{1,1}=filter(b,a,EMG_SO_G_L.*1e6);
        EMG_sal{1,1}=filter(d,c,EMG_sal{1,1});
        EMG_sal{1,1}(isnan(EMG_sal{1,1}))=0;
        EMG_SO_G_R(isnan(EMG_SO_G_R))=0;
        EMG_sal{1,2}=filter(b,a,EMG_SO_G_R.*1e6);
        EMG_sal{1,2}=filter(d,c,EMG_sal{1,2});
        EMG_sal{1,2}(isnan(EMG_sal{1,2}))=0;
        EMG_ST_G_L(isnan(EMG_ST_G_L))=0;
        EMG_sal{2,1}=filter(b,a,EMG_ST_G_L.*1e6);
        EMG_sal{2,1}=filter(d,c,EMG_sal{2,1});
        EMG_ST_G_R(isnan(EMG_ST_G_R))=0;
        EMG_sal{2,2}=filter(b,a,EMG_ST_G_R.*1e6);
        EMG_sal{2,2}=filter(d,c,EMG_sal{2,2});
        EMG_TA_G_L(isnan(EMG_TA_G_L))=0;
        EMG_sal{3,1}=filter(b,a,EMG_TA_G_L.*1e6);
        EMG_sal{3,1}=filter(d,c,EMG_sal{3,1});
        EMG_TA_G_R(isnan(EMG_TA_G_R))=0;
        EMG_sal{3,2}=filter(b,a,EMG_TA_G_R.*1e6);
        EMG_sal{3,2}=filter(d,c,EMG_sal{3,2});
        EMG_RF_G_L(isnan(EMG_RF_G_L))=0;
        EMG_sal{4,1}=filter(b,a,EMG_RF_G_L.*1e6);
        EMG_sal{4,1}=filter(d,c,EMG_sal{4,1});
        EMG_RF_G_R(isnan(EMG_RF_G_R))=0;
        EMG_sal{4,2}=filter(b,a,EMG_RF_G_R.*1e6);
        EMG_sal{4,2}=filter(d,c,EMG_sal{4,2});
        EMG_sal{1,1}(isnan(EMG_sal{1,1}))=0;
        EMG_sal{1,2}(isnan(EMG_sal{1,2}))=0;
        EMG_sal{2,1}(isnan(EMG_sal{2,1}))=0;
        EMG_sal{2,2}(isnan(EMG_sal{2,2}))=0;
        EMG_sal{3,1}(isnan(EMG_sal{3,1}))=0;
        EMG_sal{3,2}(isnan(EMG_sal{3,2}))=0;
        EMG_sal{4,1}(isnan(EMG_sal{4,1}))=0;
        EMG_sal{4,2}(isnan(EMG_sal{4,2}))=0;
    else
        EMG_sal=[];
    end;
    EEGLh=EEGdata.EEG_rec(2:33,t1hl-5:t2hl+5)';
    EEGRh=EEGdata.EEG_rec(2:33,t1hr-5:t2hr+5)';
    if t2hl>=t1r
        t1r=DataFC.Right(val_in+1);
    end;
    if t2hr>=t1l
        t1l=DataFC.Left(val_in+1);
    end;
    EEGLp=EEGdata.EEG_rec(2:33,t2hl-5:t1r+5)';
    EEGRp=EEGdata.EEG_rec(2:33,t2hr-5:t1l+5)';
    t=linspace(0,size(EEGL,1)*(1/250),size(EEGL,1));
    tr=linspace(0,size(EEGR,1)*(1/250),size(EEGR,1));
    th=linspace(0,size(EEGLh,1)*(1/250),size(EEGLh,1));
    thr=linspace(0,size(EEGRh,1)*(1/250),size(EEGRh,1));
    tp=linspace(0,size(EEGLp,1)*(1/250),size(EEGLp,1));
    tpr=linspace(0,size(EEGRp,1)*(1/250),size(EEGRp,1));
    EEGsal{1,1}=t;
    EEGsal{1,2}=tr;
    EEGsal{4,1}=th;
    EEGsal{4,2}=thr;
    EEGsal{6,1}=tp;
    EEGsal{6,2}=tpr;
    EEGsal{2,1}=EEGL;
    EEGsal{2,2}=EEGR;
    EEGsal{3,1}=EEGLh;
    EEGsal{3,2}=EEGRh;
    EEGsal{5,1}=EEGLp;
    EEGsal{5,2}=EEGRp;
    accsal{1}=accelLT;
    accsal{2}=accelRT;
else
    ind=1;
    EEGsal=[];
    accsal=[];
    EMG_sal=[];
end;
if plot_sel==1
    %% Shank Signal
    figure;
    subplot(311),
    plot(t,accelL,'LineWidth',2.5),
    grid on;
    xlabel('time [s]');
    ylabel('[g or m/s^2]?')
    set(gca,'FontSize',12);
    legend('x','y','z');
    title('Shank L')
    subplot(312),
    plot(tr,accelR,'LineWidth',2.5);
    grid on;
    xlabel('time [s]');
    ylabel('[g or m/s^2]?')
    set(gca,'FontSize',12);
    legend('x','y','z');
    title('Shank R')
    subplot(313);
    plot(t,EEGL,'LineWidth',2.5);
    grid on;
    xlabel('time [s]');
    ylabel('[uv]')
    set(gca,'FontSize',12);
    legend('FP1','FP2','AF3','AF4','F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6', 'T7','C3','CZ','C4','T8','CP5','CP1','CP2','CP6','P7','P3','PZ','P4','P8','PO7','PO3','PO4','PO8','OZ');
    title('EEG');
    %% Thigh Signal
    figure;
    subplot(311),
    plot(t,accelLT,'LineWidth',2.5),
    grid on;
    xlabel('time [s]'),
    ylabel('[g or m/s^2]?')
    set(gca,'FontSize',12),
    legend('x','y','z'),
    title('Thigh L')
    subplot(312),
    plot(tr,accelRT,'LineWidth',2.5),
    grid on;
    xlabel('time [s]'),
    ylabel('[g or m/s^2]?')
    set(gca,'FontSize',12),
    legend('x','y','z'),
    title('Thigh R')
    subplot(313),
    plot(tr,EEGR,'LineWidth',2.5);
    grid on;
    ylabel('[uv]');
    xlabel('time [s]');
    set(gca,'FontSize',12);
    legend('FP1','FP2','AF3','AF4','F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6', 'T7','C3','CZ','C4','T8','CP5','CP1','CP2','CP6','P7','P3','PZ','P4','P8','PO7','PO3','PO4','PO8','OZ');
    title('EEG')
end;
