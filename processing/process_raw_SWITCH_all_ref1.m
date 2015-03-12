

 files=dir('CR_*.vhdr');
 files={files.name};
% 
% files_s=dir('*_Prt2_*.vhdr');
% files_s={files_s.name};
% 
% files = [files_p files_s];
% files{13}=[];
% 
 for ffil = 1:length(files)
   
    path_n='I:\GSK\ActCh\P10_sc_disect';
  % EIT_fname = 'CR_EIT009_T_EIT_100uA_6kHz_20s_200ms_Prt1_ref1_NMBA.vhdr';
   EIT_fname = files{ffil};
   %if ffil<=length(files_p)
        Prt_name='Prt_1.txt';
 %  else
  %      Prt_name='Prt_2.txt';
  % end

  map_p=[[ 3 4 7 8 9 12 13 14 15 16 19 21 22 23 26 27 28]-1]; % Montage
   map_=[ 3 4 7 8 9 12 13 14 15 16 19 21 22 23 26 27 28]; 
   
   bad_ch=[1,5,10,16];
   
   ring = [7 28	12	23	21	16	26	9	13	4	3	15	8	27	14	22];
     
   EP_cutoff = 2000;  % cutoff frequency for EPs (low-pass freq.)
   dZ_BW = 2000;      % Bandwidth for demodulation
   N_butter_EP = 3; % Butterworth order for filtering EPs
   N_butter_dZ = 3;  % Butterworth order for filtering dZs
   N_butter_notch = 3; % Butterworth order for notch filtering    
   Noise_thres = 1000;  % Noise threshold (in uV)
   Filter_50Hz = false; % 50 Hz notch filter on EPs
   
   T_window = 0.1; % time to take around trigger
 
%% find all prt switches
   chdir(path_n);    
   HDR= sopen([EIT_fname(1:end-5) '.eeg']);
   Fs=HDR.SampleRate;
   
   
   SWITCH=HDR.EVENT.POS;

    for i=1:length(SWITCH)
       if ~strcmp([HDR.EVENT.Desc{i}], ['S 13']);
          SWITCH(i)=0;
       end
    end

    SWITCH(SWITCH==0)=[];
    SWITCH=SWITCH(3:end-4);
% Load Prt and check that switches are ok
 disp(['Loading protocol file ',Prt_name, ' ...']);
 Prt = dlmread(Prt_name,'\t');
 Prt_size = size(Prt,1);
 
 if (Prt_size+1==size(SWITCH,1))
     disp('All data seemed to be there');
 elseif (Prt_size+1>size(SWITCH,1))
     disp('Protocol size is bigger than the number of switches in the file');
     Prt_size=size(SWITCH,1)-1;
 else
     disp('Protocol size is smaller than the number of switches in the file!'); %%, assuming repeated protocol');
     %%write stuff to handle this
 end

 EEG = pop_loadbv('',EIT_fname,[SWITCH(1)+Fs SWITCH(1)+10*Fs]);
 EEG.data=EEG.data(map_p,:);
 V_inj = detrend(double(EEG.data(map_==Prt(1,2),:)'),'constant'); 
 
 NFFT = 2^nextpow2(length(V_inj)); % Next power of 2 from length of y
Y = fft(V_inj,NFFT)/length(V_inj);
f = Fs/2*linspace(0,1,NFFT/2+1);
w_ind=2*abs(Y(1:NFFT/2+1));
 %plot(f,w_ind) 
 
 [~,maxw] = max((w_ind));
 Fc = f(maxw);
  disp(sprintf('****** Detected carrier frequency: Fc = %i Hz ******',Fc));
 
 %%
EIT = [];  
 
strtt = EIT_fname(1:end-4);
strtt(strtt=='_')= ' ';
        
 hhplots = figure('Position',[10,50,1900,950],'PaperPositionMode','auto','name',strtt);
                  
           
                  
       
            
  
for iPair = 1:Prt_size
      
      disp(sprintf('Processing protocol line %02i / %02i',iPair,Prt_size));
      
      EEG = pop_loadbv( '',EIT_fname,[SWITCH(iPair) SWITCH(iPair+1)]);
      
      EEG.data=EEG.data(map_p,:);
      
      T_trig=cell2mat({EEG.event.latency})';
 
      for i=1:length(T_trig)
         if strcmp([EEG.event(1,i).type], ['S 12']);
           T_trig(i)=0;
         end
      end
      
      T_trig(T_trig==0)=[];
      T_trig=T_trig(4:end-3);
      
      T_stim = mean(T_trig(2:end)-T_trig(1:end-1))/Fs;
      disp(['Siimulation every ' num2str(round(T_stim*1000)) ' ms']);
      
 
  
 %%
 
    Data = double(EEG.data');
    
    if Fc<=300
        dZ_BW = 50;
        EP_cutoff = 100;
    end
    
 
    % Low-pass filter to retrieve EP's
    [b,a] = butter(N_butter_EP,EP_cutoff/(Fs/2),'low');
    X_ep = filtfilt(b,a,Data);
    if Filter_50Hz
        [b,a] = iirnotch(50/(Fs/2),(50/(Fs/2))/35);
        X_ep = filtfilt(b,a,X_ep);
    end
    
    % Band-pass filter and demodulate
    [b,a] = butter(N_butter_dZ,(Fc+dZ_BW*[-1,1])/(Fs/2));
    X_dz = filtfilt(b,a,Data);
    X_dz = hilbert(X_dz);
    A_dz   = abs(X_dz);              % amplitude
%     if Filter_50Hz
          [b,a] = iirnotch(50/(Fs/2),(50/(Fs/2))/25);
          A_dz = filtfilt(b,a,A_dz);
%      end
    
    % Number of channels, triggers, bins
    N_chan = size(Data,2);
    N_trig = length(T_trig);
    N_bin = round(T_window*Fs);
    w = (1:N_bin)-round(N_bin*0.5);    % window around trigger
    T = 1e3*w/Fs;
    t0 = find(T>-10 & T<-5);
    
    % Get results of each channel in matrices
    EP   = cell(1,N_chan);
    dZ   = cell(1,N_chan);
    BV = cell(1,N_chan);
    for iChan = 1:N_chan
        EP{iChan} = zeros(N_bin,N_trig);
        dZ{iChan} = zeros(N_bin,N_trig);
        BV{iChan} = zeros(1,N_trig);
        
        % iterate over triggers
        for jTrig = 1:N_trig
            ival = T_trig(jTrig)+w;
            EP{iChan}(:,jTrig) = detrend(X_ep(ival,iChan),'constant');
            dZ{iChan}(:,jTrig) = detrend(A_dz(ival,iChan),'constant');
            BV{iChan}(jTrig) = mean(A_dz(ival,iChan));
        end
    end

    % Find trials with low noise
    Trial_sel = [];
    for iChan = setdiff(1:N_chan,[Prt(iPair,:)])
        Trial_sel = [Trial_sel; BV{iChan}>1e2 ...
                    & max(abs(dZ{iChan}))<Noise_thres];
    end
    good_chan = sum(Trial_sel,2)>size(Trial_sel,2)/2;
% $$$     disp(sprintf('%i good channels',nnz(good_chan)))
    Trial_sel = sum(Trial_sel(good_chan,:))>nnz(good_chan)/2;
    N_trial = nnz(Trial_sel);
    Trial_sel = Trial_sel;
    disp(sprintf('Removed %i trials (out of %i) for excessive noise',N_trig-N_trial,N_trig))
    
    % Compute average / std. error for EP, dZ, dPhi 
    avg_EP = zeros(N_bin,N_chan);
    avg_dZ_rel = zeros(N_bin,N_chan);
    avg_dZ_abs = zeros(N_bin,N_chan);
    BV0=zeros(length(ring),1);
    
    [~,mm]=max(Data(:,1));
    invert = Data(mm,:)<0;
    
    for iChan = 1:N_chan
        ep = EP{iChan}(:,Trial_sel);
        dz = dZ{iChan}(:,Trial_sel);
        bv = repmat(BV{iChan}(Trial_sel),N_bin,1);        
        
        avg_EP(:,iChan) = mean(ep,2)-mean(reshape(ep(t0,:),[],1));
        avg_dZ_abs(:,iChan) = mean(dz,2)-mean(reshape(dz(t0,:),[],1));
        avg_dZ_rel(:,iChan) = 100*mean(dz./bv,2);            
        avg_dZ_rel(:,iChan) = avg_dZ_rel(:,iChan)-mean(avg_dZ_rel(t0,iChan));     
        if (~invert(iChan))
            BV0(ring==map_(iChan))=mean(mean(bv,2));
        else
            BV0(ring==map_(iChan))=-mean(mean(bv,2));
        end
        
    end

    
    
    EIT{iPair}.Pair = Prt(iPair,:);
    EIT{iPair}.EP_avg=avg_EP;
    EIT{iPair}.dZ_avg=avg_dZ_abs;
    EIT{iPair}.dZ=dZ;
    EIT{iPair}.BV0=BV0;
    
% Plot
    
%     if ~isempty(Plot_chan)
%         % create plot
%         figure('Position',[680,50,560,950],'PaperPositionMode','auto')
%         subplot(3,1,1)
%         plot(T,detrend(avg_EP(:,Plot_chan)))
%        % ylim([-500,500])
%         %xlim([-5,50])
%         grid on
%         ylabel('EP (uV)')
%         subplot(3,1,2)
%         plot(T,detrend(avg_dZ_abs(:,Plot_chan)))
%         %ylim([-30,30])
%         %xlim([-5,50])
%         grid on
%         ylabel('dZ (uV)')
%         
%         strtt = EIT_fname(1:end-4);
%         strtt(strtt=='_')= ' ';
%         
%         title(sprintf('EIT %s iPair %02i (%i, %i) [Fc: %i]',...
%                       strtt,iPair,Prt(iPair,1),...
%                       Prt(iPair,2),Fc))
%         subplot(3,1,3)
%         plot(T,detrend(avg_dZ_rel(:,Plot_chan)))
%         ylim([-0.2,0.2])
%         %xlim([-5,50])
%         grid on    
%         ylabel('dZ (%)')
%         xlabel('Time (ms)')
% %        print(gcf,sprintf('png2/EIT_%s_Ipair_%02i.png',EIT_serial,iPair),...
% %               '-dpng','-r300')
%        % close(gcf)
%     end

        [~,k]=ismember(EIT{iPair}.Pair,ring);
        Plot_chan = setdiff([1:length(map_)],[bad_ch k]);
        
        subplot(round(sqrt(Prt_size))+1,round(Prt_size/round(sqrt(Prt_size))),iPair);
      
        
        plot(T,detrend(avg_dZ_abs(:,Plot_chan)))
        title(['Pair=' num2str(iPair)]);
        %ylim([-30,30])
        %xlim([-5,50])
        grid on
        ylabel('dZ (uV)')
        
        EIT{iPair}.Pair_0 = k;
        
        strtt = EIT_fname(1:end-4);
        strtt(strtt=='_')= ' ';
        
 
    drawnow;
     
end
      
%% plot boundary voltages for checking the electrode

figure('name',sprintf('EIT %s iPair %02i (%i, %i) [Fc: %i]',...
                      strtt,iPair,Prt(iPair,1),...
                      Prt(iPair,2),Fc),'Position',[10,50,1900,950],'PaperPositionMode','auto'); 
for iPair=1:Prt_size  
        
         subplot(round(sqrt(Prt_size))+1,round(sqrt(Prt_size)),iPair);
       
        bar(EIT{iPair}.BV0);
        [~,k]=ismember(EIT{iPair}.Pair,ring);
        title(['Pair  ' num2str(k)]);
end

%end
save([EIT_fname(1:end-5) '_bw4.mat'],'EIT','-v7.3');
      
 end
 
%% 

% 
% dZ=[];
% Prt_0=[];
% n=1;
% for iPair=1:Prt_size
%     
%     for i = [2:11,13:size(EIT{iPair}.dZ,2)] 
%         if max(abs(EIT{iPair}.dZ_avg(:,i)))<300 && std(EIT{iPair}.dZ_avg(1:3000,i)) <10
%             dZ(:,n)=EIT{iPair}.dZ_avg(:,i);
%            
%             g=1:16;
%             Prt_0(n,:) = [g(ring==EIT{iPair}.Pair(1)), g(ring==EIT{iPair}.Pair(2)), ...
%                 g(ring==map_(i)) ];
%             n=n+1;
%         end
%     end
% end
%     
%      plot(dZ);
%     
%     save(['for_image_bw4_no_big_' EIT_fname(1:end-5) '.mat'],'dZ','Prt_0','-v7.3');
%     
%    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      