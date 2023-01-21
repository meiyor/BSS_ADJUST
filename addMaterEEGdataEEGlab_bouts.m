function addMaterEEGdataEEGlab_bouts(subj,runind,size_resamp,sel_emg)
%% add the EEGlab path here to employ EEG signal preprocessing if necessary %%
%addpath(genpath('/home/jmm_vivobook_asus/DeepGaze_project/eeglab14_1_1b'));%%%% add path eeglab
%
%rmpath('/home/jmm_vivobook_asus/DeepGaze_project/eeglab14_1_1b/field_trip'); %%%% remove field_trip from eeglab
ind=0;
ccount=1;
ckcount=1;
if sel_emg==1
    EMG=load('syncData.mat');
else
    EMG=[];
end;
while ind==0
     EEGL=eeg_emptyset();
     EEGR=eeg_emptyset();
     [accel{ccount},EEG{ccount},ind,indbout,EMG_sal{ccount}]=read_subjects_EEG_Mater_alt(subj,runind,ckcount,1,0,1,sel_emg,5,EMG.syncData); %% processing the window on 5 seconds around the Bout detection
     %ind
     if ind==1
         %save([subj '_' runind '_res_vals.mat'],'erspl','itcpl','freql','erspr','itcpr','freqr','tl','tr');
         break;
     end;
     EEGL_set{ccount}=EEGlabstructure_def(EEGL,EEG{ccount}{2,1}',EEG{ccount}{4,1},EEG{ccount}{6,1},EEG{ccount}{1,1},subj,runind);
     EEGR_set{ccount}=EEGlabstructure_def(EEGR,EEG{ccount}{2,2}',EEG{ccount}{4,2},EEG{ccount}{6,2},EEG{ccount}{1,2},subj,runind);
     
     if indbout==0
         EEGL_set{ccount} = pop_eegfiltnew(EEGL_set{ccount},0.1,100,8250);
         EEGR_set{ccount} = pop_eegfiltnew(EEGR_set{ccount},0.1,100,8250);
         for k=1:32
                    t_data=conv((1/2)*ones([1 2]),EEGL_set{ccount}.data(k,:));
                    t_data=detrend(t_data);
                    EEGL_set{ccount}.data(k,:)=t_data(2/2:EEGL_set{ccount}.pnts+(2/2-1));
                    t_data=conv((1/2)*ones([1 2]),EEGR_set{ccount}.data(k,:));
                    t_data=detrend(t_data);
                    EEGR_set{ccount}.data(k,:)=t_data(2/2:EEGR_set{ccount}.pnts+(2/2-1));
         end;
         
         baseline_channame_L={EEGL_set{ccount}.chanlocs(:).labels};
         baseline_chanloc_L=EEGL_set{ccount}.chanlocs;
         
         baseline_channame_R={EEGR_set{ccount}.chanlocs(:).labels};
         baseline_chanloc_R=EEGR_set{ccount}.chanlocs;
         
         EEGL_set{ccount}=pop_autobssemg(EEGL_set{ccount},1,1,'bsscca',{'eigratio',1e6},'emg_psd',{'ratio',10,'fs',250,'range',[1 10]});
         EEGR_set{ccount}=pop_autobssemg(EEGR_set{ccount},1,1,'bsscca',{'eigratio',1e6},'emg_psd',{'ratio',10,'fs',250,'range',[1 10]});
         EEGL_set{ccount}=clean_rawdata(EEGL_set{ccount},5,[0.25 0.75],0.85,-1,-1,-1);
         EEGR_set{ccount}=clean_rawdata(EEGR_set{ccount},5,[0.25 0.75],0.85,-1,-1,-1);
         
         %% organize data after ASR
         EEGchanL=32;
         EEGchanR=32;
         
         EEG_temp_L=EEGL_set{ccount};
         EEG_temp_R=EEGR_set{ccount};
         if EEGL_set{ccount}.nbchan~=32
               valdiff=setdiff(baseline_channame_L,{EEGL_set{ccount}.chanlocs(:).labels});
               inddiff=find(ismember(baseline_channame_L,valdiff)==1);
               EEG_temp_L.chanlocs=baseline_chanloc_L;
               c_Prev=1;
               c_Prev_orig=1;
               temp=zeros([EEGchanL,size(EEGL_set{ccount}.data,2),size(EEGL_set{ccount}.data,3)]);
               if length(inddiff)==1
                  temp(1:inddiff(1)-1,:,:)=EEGL_set{ccount}.data(1:inddiff(1)-1,:,:);
                  temp(inddiff(1)+1:end,:,:)=EEGL_set{ccount}.data(inddiff(1):end,:,:);
               else
                   for j=1:length(inddiff)
                       if j~=length(inddiff)
                           if j==1
                              temp(c_Prev:inddiff(j)-1,:,:)=EEGL_set{ccount}.data(c_Prev:inddiff(j)-1,:,:);
                           else
                               if c_Prev+(inddiff(j+1)-inddiff(j))-1>size(EEGL_set{ccount}.data,1)
                                    val_diff=((c_Prev+(inddiff(j+1)-inddiff(j))-1)-size(EEGL_set{ccount}.data,1));
                                    temp(inddiff(j)+1:inddiff(j+1)-1,:,:)=EEGL_set{ccount}.data((c_Prev+1:c_Prev+(inddiff(j+1)-inddiff(j))-1)-val_diff,:,:);
                               else
                                    temp(inddiff(j)+1:inddiff(j+1)-1,:,:)=EEGL_set{ccount}.data(c_Prev+1:c_Prev+(inddiff(j+1)-inddiff(j))-1,:,:);
                               end;
                           end;
                       else
                           if inddiff(j)~=EEGchanL
                             if size(temp(inddiff(j)+1:end,:,:),1)~=size(EEGL_set{ccount}.data(c_Prev:end,:,:),1)   
                                 temp(inddiff(j)+1:end,:,:)=EEGL_set{ccount}.data(c_Prev-1:end,:,:);
                             else
                                 temp(inddiff(j)+1:end,:,:)=EEGL_set{ccount}.data(c_Prev:end,:,:);
                             end
                           end;
                       end;
                       %temp(inddiff(j),:)=zeros([1,size(EEG{tm_sb}.data,2)]);
                       if j==1
                          if length(inddiff)>2
                             c_Prev=c_Prev+(inddiff(j+1)-inddiff(j));%inddiff(j);
                          else
                             c_Prev=inddiff(j+1)-1;
                          end;
                       elseif j~=length(inddiff)
                           c_Prev=c_Prev+(inddiff(j+1)-inddiff(j))-1;
                       end;
                   end;
          end;
          EEGL_set{ccount}.data=temp;
          EEGL_set{ccount}.nbchan=EEGchanL;
         end;

         if EEGR_set{ccount}.nbchan~=32
               valdiff=setdiff(baseline_channame_R,{EEGR_set{ccount}.chanlocs(:).labels});
               inddiff=find(ismember(baseline_channame_R,valdiff)==1);
               EEG_temp_R.chanlocs=baseline_chanloc_R;
               c_Prev=1;
               c_Prev_orig=1;
               temp=zeros([EEGchanR,size(EEGR_set{ccount}.data,2),size(EEGR_set{ccount}.data,3)]);
               if length(inddiff)==1
                  temp(1:inddiff(1)-1,:,:)=EEGR_set{ccount}.data(1:inddiff(1)-1,:,:);
                  temp(inddiff(1)+1:end,:,:)=EEGR_set{ccount}.data(inddiff(1):end,:,:);
               else
                   for j=1:length(inddiff)
                       if j~=length(inddiff)
                           if j==1
                              temp(c_Prev:inddiff(j)-1,:,:)=EEGR_set{ccount}.data(c_Prev:inddiff(j)-1,:,:);
                           else
                               if c_Prev+(inddiff(j+1)-inddiff(j))-1>size(EEGR_set{ccount}.data,1)
                                    val_diff=((c_Prev+(inddiff(j+1)-inddiff(j))-1)-size(EEGR_set{ccount}.data,1));
                                    temp(inddiff(j)+1:inddiff(j+1)-1,:,:)=EEGR_set{ccount}.data((c_Prev+1:c_Prev+(inddiff(j+1)-inddiff(j))-1)-val_diff,:,:);
                               else
                                    temp(inddiff(j)+1:inddiff(j+1)-1,:,:)=EEGR_set{ccount}.data(c_Prev+1:c_Prev+(inddiff(j+1)-inddiff(j))-1,:,:);
                               end;
                           end;
                       else
                           if inddiff(j)~=EEGchanR
                             if size(temp(inddiff(j)+1:end,:,:),1)~=size(EEGR_set{ccount}.data(c_Prev:end,:,:),1)   
                                 temp(inddiff(j)+1:end,:,:)=EEGR_set{ccount}.data(c_Prev-1:end,:,:);
                             else
                                 temp(inddiff(j)+1:end,:,:)=EEGR_set{ccount}.data(c_Prev:end,:,:);
                             end
                           end;
                       end;
                       %temp(inddiff(j),:)=zeros([1,size(EEG{tm_sb}.data,2)]);
                       if j==1
                          if length(inddiff)>2
                             c_Prev=c_Prev+(inddiff(j+1)-inddiff(j));%inddiff(j);
                          else
                             c_Prev=inddiff(j+1)-1;
                          end;
                       elseif j~=length(inddiff)
                           c_Prev=c_Prev+(inddiff(j+1)-inddiff(j))-1;
                       end;
                   end;
          end;
         EEGR_set{ccount}.data=temp;
         EEGR_set{ccount}.nbchan=EEGchanR;
         end;

         
         [EEGL_set{ccount}.icaweights,EEGL_set{ccount}.icasphere]=runica(EEGL_set{ccount}.data(:,:),'sphering','on','lrate',1e-5,'maxsteps',150,'reset_randomseed','off'); 
         [EEGR_set{ccount}.icaweights,EEGR_set{ccount}.icasphere]=runica(EEGR_set{ccount}.data(:,:),'sphering','on','lrate',1e-5,'maxsteps',150,'reset_randomseed','off'); 
         delete('report_left_n.txt');
         delete('report_right_n.txt');
         %% implement ADJUST here
         %% left
         EEGL_set{ccount}.icawinv=inv(EEGL_set{ccount}.icaweights*EEGL_set{ccount}.icasphere);
         EEGL_set{ccount}.data= repmat(EEGL_set{ccount}.data(:,:),1,1,25);
         EEGL_set{ccount}.trials=25;
         EEGL_set{ccount}.icachansind=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
         [art_channels_left]=ADJUST(EEGL_set{ccount},'report_left_n.txt');
         EEGL_set{ccount}=pop_subcomp(EEGL_set{ccount},art_channels_left);
         EEGL_set{ccount}.data=mean(EEGL_set{ccount}.data(:,:,:),3);
         %% right
         EEGR_set{ccount}.icawinv=inv(EEGR_set{ccount}.icaweights*EEGR_set{ccount}.icasphere);
         EEGR_set{ccount}.data= repmat(EEGR_set{ccount}.data(:,:),1,1,25);
         EEGR_set{ccount}.trials=25;
         EEGR_set{ccount}.icachansind=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
         delete('report_left_n.txt');
         delete('report_right_n.txt');
         [art_channels_right]=ADJUST(EEGR_set{ccount},'report_right_n.txt');
         EEGR_set{ccount}=pop_subcomp(EEGR_set{ccount},art_channels_right);
         EEGR_set{ccount}.data=mean(EEGR_set{ccount}.data(:,:,:),3);
         tl(ccount)=EEGL_set{ccount}.pnts*(1/250);
         tr(ccount)=EEGR_set{ccount}.pnts*(1/250);
         thl(ccount)=EEGL_set{ccount}.pntsh*(1/250);
         thr(ccount)=EEGR_set{ccount}.pntsh*(1/250);
         tpl(ccount)=EEGL_set{ccount}.pntsp*(1/250);
         tpr(ccount)=EEGR_set{ccount}.pntsp*(1/250);
         if tpl(ccount)>=10.5
            tpl(ccount)=0;
         end;
         if tpr(ccount)>=10.5
            tpr(ccount)=0;
         end;
         
    close all;
    if ~isreal(EEGL_set{ccount}.data(:,:))
      EEGL_set{ccount}.data(:,:)=real(EEGL_set{ccount}.data(:,:));
    end;
    if ~isreal(EEGR_set{ccount}.data(:,:))
      EEGR_set{ccount}.data(:,:)=real(EEGR_set{ccount}.data(:,:));
    end;
    
    if any(any(isinf(EEGL_set{ccount}.data(:,:))))
      EEGL_set{ccount}.data(:,:)=0;
    end;
    if any(any(isinf(EEGR_set{ccount}.data(:,:))))
      EEGR_set{ccount}.data(:,:)=0;
    end;     
    
    if any(any(isnan(EEGL_set{ccount}.data(:,:))))
      EEGL_set{ccount}.data(:,:)=0;
    end;
    if any(any(isnan(EEGR_set{ccount}.data(:,:))))
      EEGR_set{ccount}.data(:,:)=0;
    end;   
    
    if any(any((EEGL_set{ccount}.data(:,:)==0)))
         EEGL_set{ccount}.data(EEGL_set{ccount}.data(:,:)==0)=0.01;
    end;
    
    if any(any((EEGR_set{ccount}.data(:,:)==0)))
         EEGR_set{ccount}.data(EEGR_set{ccount}.data(:,:)==0)=0.01;
    end;
    
    for ch=1:32
              close all;
              [erspl{ccount,ch},itcpl{ccount,ch},powbasel{ccount,ch},timescl{ccount,ch},freql{ccount,ch}]=newtimef(EEGL_set{ccount}.data(ch,:),EEGL_set{ccount}.pnts,[EEGL_set{ccount}.xmin EEGL_set{ccount}.xmax]*1000,EEGL_set{ccount}.srate,0,'winsize',50,'nfreqs',1000,'freqs',[0 50],'padratio',32,'plotersp','off','plotitc','off');
              [erspr{ccount,ch},itcpr{ccount,ch},powbaser{ccount,ch},timescr{ccount,ch},freqr{ccount,ch}]=newtimef(EEGR_set{ccount}.data(ch,:),EEGR_set{ccount}.pnts,[EEGR_set{ccount}.xmin EEGR_set{ccount}.xmax]*1000,EEGR_set{ccount}.srate,0,'winsize',50,'nfreqs',1000,'freqs',[0 50],'padratio',32,'plotersp','off','plotitc','off');
              %[erspr{ccount,ch},itcpr{ccount,ch},powbaser{ccount,ch},timescr{ccount,ch},freqr{ccount,ch}]=newtimef(EEGR_set{ccount}.data(ch,:),EEGR_set{ccount}.pnts,[EEGR_set{ccount}.xmin EEGR_set{ccount}.xmax]*1000,EEGR_set{ccount}.srate,0,'itctype','phasecoher','winsize',40,'freqs',[0.1 20],'nfreqs',250,'padratio',16,'plotersp','off','plotitc','off');
              if sel_emg==1
                  [PL1{ccount,ch},fL1{ccount}]=pwelch(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),250/2,1000,250,'psd','onesided');
                  [PR1{ccount,ch},fR1{ccount}]=pwelch(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),250/2,1000,250,'psd','onesided');
                  [SL1{ccount,ch}]=xspectrogram(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),250/2,10,200,250,'psd');
                  [SR1{ccount,ch}]=xspectrogram(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),250/2,10,200,250,'psd');
                  SL1{ccount,ch}=imresize(SL1{ccount,ch},size_resamp);
                  SR1{ccount,ch}=imresize(SR1{ccount,ch},size_resamp);
                  [PL2_SO{ccount,ch},fL2_SO{ccount}]=pwelch(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{1,1}./max(EMG_sal{ccount}{1,1}),250/2,1000,250,'psd','onesided');
                  [PR2_SO{ccount,ch},fR2_SO{ccount}]=pwelch(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{1,2}./max(EMG_sal{ccount}{1,2}),250/2,1000,250,'psd','onesided');
                  [SL2_SO{ccount,ch}]=xspectrogram(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{1,1}./max(EMG_sal{ccount}{1,1}),250/2,10,200,250,'psd');
                  [SR2_SO{ccount,ch}]=xspectrogram(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{1,2}./max(EMG_sal{ccount}{1,2}),250/2,10,200,250,'psd');
                  SL2_SO{ccount,ch}=imresize(SL2_SO{ccount,ch},size_resamp);
                  SR2_SO{ccount,ch}=imresize(SR2_SO{ccount,ch},size_resamp);
                  [PL2_ST{ccount,ch},fL2_ST{ccount}]=pwelch(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{2,1}./max(EMG_sal{ccount}{2,1}),250/2,1000,250,'psd','onesided');
                  [PR2_ST{ccount,ch},fR2_ST{ccount}]=pwelch(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{2,2}./max(EMG_sal{ccount}{2,2}),250/2,1000,250,'psd','onesided');
                  [SL2_ST{ccount,ch}]=xspectrogram(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{2,1}./max(EMG_sal{ccount}{2,1}),250/2,10,200,250,'psd');
                  [SR2_ST{ccount,ch}]=xspectrogram(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{2,2}./max(EMG_sal{ccount}{2,2}),250/2,10,200,250,'psd');
                  SL2_ST{ccount,ch}=imresize(SL2_ST{ccount,ch},size_resamp);
                  SR2_ST{ccount,ch}=imresize(SR2_ST{ccount,ch},size_resamp);
                  [PL2_TA{ccount,ch},fL2_TA{ccount}]=pwelch(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{3,1}./max(EMG_sal{ccount}{3,1}),250/2,1000,250,'psd','onesided');
                  [PR2_TA{ccount,ch},fR2_TA{ccount}]=pwelch(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{3,2}./max(EMG_sal{ccount}{3,2}),250/2,1000,250,'psd','onesided');
                  [SL2_TA{ccount,ch}]=xspectrogram(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{3,1}./max(EMG_sal{ccount}{3,1}),250/2,10,200,250,'psd');
                  [SR2_TA{ccount,ch}]=xspectrogram(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{3,2}./max(EMG_sal{ccount}{3,2}),250/2,10,200,250,'psd');
                  SL2_TA{ccount,ch}=imresize(SL2_TA{ccount,ch},size_resamp);
                  SR2_TA{ccount,ch}=imresize(SR2_TA{ccount,ch},size_resamp);
                  [PL2_RF{ccount,ch},fL2_RF{ccount}]=pwelch(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{4,1}./max(EMG_sal{ccount}{4,1}),250/2,1000,250,'psd','onesided');
                  [PR2_RF{ccount,ch},fR2_RF{ccount}]=pwelch(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{4,2}./max(EMG_sal{ccount}{4,2}),250/2,1000,250,'psd','onesided');
                  [SL2_RF{ccount,ch}]=xspectrogram(squeeze(EEGL_set{ccount}.data(ch,:))./max(squeeze(EEGL_set{ccount}.data(ch,:))),EMG_sal{ccount}{4,1}./max(EMG_sal{ccount}{4,1}),250/2,10,200,250,'psd');
                  [SR2_RF{ccount,ch}]=xspectrogram(squeeze(EEGR_set{ccount}.data(ch,:))./max(squeeze(EEGR_set{ccount}.data(ch,:))),EMG_sal{ccount}{4,2}./max(EMG_sal{ccount}{4,2}),250/2,10,200,250,'psd');
                  SL2_RF{ccount,ch}=imresize(SL2_RF{ccount,ch},size_resamp);
                  SR2_RF{ccount,ch}=imresize(SR2_RF{ccount,ch},size_resamp);                  
              end;
              erspl{ccount,ch}=imresize(erspl{ccount,ch},size_resamp);
              erspr{ccount,ch}=imresize(erspr{ccount,ch},size_resamp);
              itcpl{ccount,ch}=imresize(itcpl{ccount,ch},size_resamp);
              itcpr{ccount,ch}=imresize(itcpr{ccount,ch},size_resamp);
         end;
         if sel_emg==1
              [PL3_SO{ccount},fL3_SO{ccount}]=pwelch(EMG_sal{ccount}{1,1}./max(EMG_sal{ccount}{1,1}),EMG_sal{ccount}{1,1}./max(EMG_sal{ccount}{1,1}),250/2,1000,250,'psd','onesided');
              [PR3_SO{ccount},fR3_SO{ccount}]=pwelch(EMG_sal{ccount}{1,2}./max(EMG_sal{ccount}{1,2}),EMG_sal{ccount}{1,2}./max(EMG_sal{ccount}{1,2}),250/2,1000,250,'psd','onesided');
              [SL3_SO{ccount}]=xspectrogram(EMG_sal{ccount}{1,1}./max(EMG_sal{ccount}{1,1}),EMG_sal{ccount}{1,1}./max(EMG_sal{ccount}{1,1}),250/2,10,200,250,'psd');
              [SR3_SO{ccount}]=xspectrogram(EMG_sal{ccount}{1,2}./max(EMG_sal{ccount}{1,2}),EMG_sal{ccount}{1,2}./max(EMG_sal{ccount}{1,2}),250/2,10,200,250,'psd');
              SL3_SO{ccount}=imresize(SL3_SO{ccount},size_resamp);
              SR3_SO{ccount}=imresize(SR3_SO{ccount},size_resamp);
              [PL3_ST{ccount},fL3_ST{ccount}]=pwelch(EMG_sal{ccount}{2,1}./max(EMG_sal{ccount}{2,1}),EMG_sal{ccount}{2,1}./max(EMG_sal{ccount}{2,1}),250/2,1000,250,'psd','onesided');
              [PR3_ST{ccount},fR3_ST{ccount}]=pwelch(EMG_sal{ccount}{2,2}./max(EMG_sal{ccount}{2,2}),EMG_sal{ccount}{2,2}./max(EMG_sal{ccount}{2,2}),250/2,1000,250,'psd','onesided');
              [SL3_ST{ccount}]=xspectrogram(EMG_sal{ccount}{2,1}./max(EMG_sal{ccount}{2,1}),EMG_sal{ccount}{2,1}./max(EMG_sal{ccount}{2,1}),250/2,10,200,250,'psd');
              [SR3_ST{ccount}]=xspectrogram(EMG_sal{ccount}{2,2}./max(EMG_sal{ccount}{2,2}),EMG_sal{ccount}{2,2}./max(EMG_sal{ccount}{2,2}),250/2,10,200,250,'psd');
              SL3_ST{ccount}=imresize(SL3_ST{ccount},size_resamp);
              SR3_ST{ccount}=imresize(SR3_ST{ccount},size_resamp);
              [PL3_TA{ccount},fL3_TA{ccount}]=pwelch(EMG_sal{ccount}{3,1}./max(EMG_sal{ccount}{3,1}),EMG_sal{ccount}{3,1}./max(EMG_sal{ccount}{3,1}),250/2,1000,250,'psd','onesided');
              [PR3_TA{ccount},fR3_TA{ccount}]=pwelch(EMG_sal{ccount}{3,2}./max(EMG_sal{ccount}{3,2}),EMG_sal{ccount}{3,2}./max(EMG_sal{ccount}{3,2}),250/2,1000,250,'psd','onesided');
              [SL3_TA{ccount}]=xspectrogram(EMG_sal{ccount}{3,1}./max(EMG_sal{ccount}{3,1}),EMG_sal{ccount}{3,1}./max(EMG_sal{ccount}{3,1}),250/2,10,200,250,'psd');
              [SR3_TA{ccount}]=xspectrogram(EMG_sal{ccount}{3,2}./max(EMG_sal{ccount}{3,2}),EMG_sal{ccount}{3,2}./max(EMG_sal{ccount}{3,2}),250/2,10,200,250,'psd');
              SL3_TA{ccount}=imresize(SL3_TA{ccount},size_resamp);
              SR3_TA{ccount}=imresize(SR3_TA{ccount},size_resamp);
              [PL3_RF{ccount},fL3_RF{ccount}]=pwelch(EMG_sal{ccount}{4,1}./max(EMG_sal{ccount}{4,1}),EMG_sal{ccount}{4,1}./max(EMG_sal{ccount}{4,1}),250/2,1000,250,'psd','onesided');
              [PR3_RF{ccount},fR3_RF{ccount}]=pwelch(EMG_sal{ccount}{4,2}./max(EMG_sal{ccount}{4,2}),EMG_sal{ccount}{4,2}./max(EMG_sal{ccount}{4,2}),250/2,1000,250,'psd','onesided');
              [SL3_RF{ccount}]=xspectrogram(EMG_sal{ccount}{4,1}./max(EMG_sal{ccount}{4,1}),EMG_sal{ccount}{4,1}./max(EMG_sal{ccount}{4,1}),250/2,10,200,250,'psd');
              [SR3_RF{ccount}]=xspectrogram(EMG_sal{ccount}{4,2}./max(EMG_sal{ccount}{4,2}),EMG_sal{ccount}{4,2}./max(EMG_sal{ccount}{4,2}),250/2,10,200,250,'psd');
              SL3_RF{ccount}=imresize(SL3_RF{ccount},size_resamp);
              SR3_RF{ccount}=imresize(SR3_RF{ccount},size_resamp);
              %% calculate coherence
              for ch=1:32
                  coh_SO_L{ccount,ch}=((PL2_SO{ccount,ch}))./sqrt(PL1{ccount,ch}.*PL3_SO{ccount});
                  coh_SO_R{ccount,ch}=((PR2_SO{ccount,ch}))./sqrt(PR1{ccount,ch}.*PR3_SO{ccount});
                  coh_ST_L{ccount,ch}=((PL2_ST{ccount,ch}))./sqrt(PL1{ccount,ch}.*PL3_ST{ccount});
                  coh_ST_R{ccount,ch}=((PR2_ST{ccount,ch}))./sqrt(PR1{ccount,ch}.*PR3_ST{ccount});
                  coh_TA_L{ccount,ch}=((PL2_TA{ccount,ch}))./sqrt(PL1{ccount,ch}.*PL3_TA{ccount});
                  coh_TA_R{ccount,ch}=((PR2_TA{ccount,ch}))./sqrt(PR1{ccount,ch}.*PR3_TA{ccount});
                  coh_RF_L{ccount,ch}=((PL2_RF{ccount,ch}))./sqrt(PL1{ccount,ch}.*PL3_RF{ccount});
                  coh_RF_R{ccount,ch}=((PR2_RF{ccount,ch}))./sqrt(PR1{ccount,ch}.*PR3_RF{ccount});
                  coh_SO_L{ccount,ch}=coh_SO_L{ccount,ch}./max(coh_SO_L{ccount,ch});
                  coh_SO_R{ccount,ch}=coh_SO_R{ccount,ch}./max(coh_SO_R{ccount,ch});
                  coh_ST_L{ccount,ch}=coh_ST_L{ccount,ch}./max(coh_ST_L{ccount,ch});
                  coh_ST_R{ccount,ch}=coh_ST_R{ccount,ch}./max(coh_ST_R{ccount,ch});
                  coh_TA_L{ccount,ch}=coh_TA_L{ccount,ch}./max(coh_TA_L{ccount,ch});
                  coh_TA_R{ccount,ch}=coh_TA_R{ccount,ch}./max(coh_TA_R{ccount,ch});
                  coh_RF_L{ccount,ch}=coh_RF_L{ccount,ch}./max(coh_RF_L{ccount,ch});
                  coh_RF_R{ccount,ch}=coh_RF_R{ccount,ch}./max(coh_RF_R{ccount,ch});
                  coh_SO_Ls{ccount,ch}=(abs(SL2_SO{ccount,ch}).^2)./((SL1{ccount,ch})*(SL3_SO{ccount}));
                  coh_SO_Rs{ccount,ch}=(abs(SR2_SO{ccount,ch}).^2)./((SR1{ccount,ch})*(SR3_SO{ccount}));
                  coh_ST_Ls{ccount,ch}=(abs(SL2_ST{ccount,ch}).^2)./((SL1{ccount,ch})*(SL3_ST{ccount}));
                  coh_ST_Rs{ccount,ch}=(abs(SR2_ST{ccount,ch}).^2)./((SR1{ccount,ch})*(SR3_ST{ccount}));
                  coh_TA_Ls{ccount,ch}=(abs(SL2_TA{ccount,ch}).^2)./((SL1{ccount,ch})*(SL3_TA{ccount}));
                  coh_TA_Rs{ccount,ch}=(abs(SR2_TA{ccount,ch}).^2)./((SR1{ccount,ch})*(SR3_TA{ccount}));
                  coh_RF_Ls{ccount,ch}=(abs(SL2_RF{ccount,ch}).^2)./((SL1{ccount,ch})*(SL3_RF{ccount}));
                  coh_RF_Rs{ccount,ch}=(abs(SR2_RF{ccount,ch}).^2)./((SR1{ccount,ch})*(SR3_RF{ccount}));
              end;
         end;
         ccount=ccount+1;
     else
         accel{ccount}=[];
         EEG{ccount}=[];
         EEGL_set{ccount}=[];
         EEGR_set{ccount}=[];
     end;
     ckcount=ckcount+1;
     ckcount
end;
save([subj '_' runind '_res_vals.mat'],'erspl','itcpl','freql','erspr','itcpr','freqr','tl','tr','accel','thl','thr','tpl','tpr');
save([subj '_' runind '_res_coherence.mat'],'coh_SO_L','coh_SO_R','coh_ST_L','coh_ST_R','coh_TA_L','coh_TA_R','coh_RF_L','coh_RF_R','coh_SO_Ls','coh_SO_Rs','coh_ST_Ls','coh_ST_Rs','coh_TA_Ls','coh_TA_Rs','coh_RF_Ls','coh_RF_Rs','PL1','PR1','PL2_SO','PR2_SO','PL2_ST','PR2_ST','PL2_RF','PR2_RF','PL3_SO','PR3_SO','PL3_ST','PR3_ST','PL3_TA','PR3_TA','PL3_RF','PR3_RF');
