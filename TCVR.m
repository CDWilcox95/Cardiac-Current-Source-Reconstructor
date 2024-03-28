close all; clear all;


% Set Program Guide
new_ECG_filter=false;
new_cond_recon=false;
new_BSM_interp=false;
new_H_find=false;
new_fwdmap_gen=false;
new_G_map=false;

num_recon_src=1;   % Number of dipole sources to reconstruct


%% Load in Data
load("C:\Users\cdwil\Desktop\EIT_EKG_Projects\Cardiac_Imaging\Data\Subject63_22_06_29\sbj63_152kHz_prone_0mA_22_06_29_16_27_43.mat");   % Only ECG Data (No EIT Data -- 0 Amp Currents applied)
load("EKg_noEIT.mat", 'Vekg');  Vekg_0Amp=Vekg;
load("L2ekg_noEIT.mat");    L2ekg_0Amp=ecg_wave;

% load("C:\Users\cdwil\Desktop\EIT_EKG_Projects\Cardiac_Imaging\Data\Subject63_22_06_29\sbj63_93kHz_prone_22_06_29_16_14_55.mat");    % ECG Data with EIT Data taken at 93 kHz

if new_ECG_filter
    FilterEKGData(frame_ECG);
end
EKG_save_file="AdaptiveFilter_ECG_subj63.mat";  load(EKG_save_file);    [~, num_samples]=size(Vekg);
ecgwave_save_file="AdaptiveFilter_ecgwave_sbj63.mat";   load(ecgwave_save_file);
Vekg=(10^-3).*Vekg;         %Convert from mV



I=(1/current_amp).*cur_pattern(:,1:31); [L,numCurr]=size(I);
[Vreal, Vimag]=CorrectVoltageData(frame_voltage);   Veit=Vreal+1i.*Vimag;


%% Set Subject Paramaters

% Subject 63 Paramaters
Circum=circumference;
if ~exist("vert_gap")
    vert_gap=7;
end

% Data Collection Frequency
global Fs;
Fs = 864.0553;                                      %Sampling Frequency (1/s)

%Conversion of Liters to mL
LtomL=10^-3;

%Set Data Structure "model_info" to be used in functions & Set Forward
%Model of Subject Domain

model_info=GetSubjectParamaters(Circum, vert_gap);model_info.R=10^-2 *model_info.R;   model_info.H=10^-2 *model_info.H; 
fmdl=GetReconMesh(model_info.R,model_info.H, new_fwdmap_gen, 'c'); model_info.FEM_Mesh=fmdl;   elec_pts=fmdl.nodes(fmdl.electrode_idx,:);   





%% Reconstruct Conductivity Using Todler 
if new_cond_recon
    [gamma_real_vec, gamma_imag_vec, sigma_b]=fTodler_act5(frame_voltage, cur_pattern, false);
    save("subject63_sigma0.mat", 'sigma_b');   save('sbj63_sigma_r.mat', 'gamma_real_vec');
    save('sbj63_sigma_i.mat', 'gamma_imag_vec');
end

load('sbj63_sigma_r.mat'); load('subject63_sigma0.mat');    


%% Set Analysis Frames ( [s0, sN]= EIT Frames |   [s0ekg, sNekg]= ECG Frames )
s0=1725;    sN=2000;
s0ekg=s0*32;    sNekg=sN*32;




%% Determine Frames Corresponding to num_cycles Complete Cardiac Cycle(s)
num_cycles=5;
Pt=model_info.FEM_Mesh.nodes;  
[mappedWaves, first_good_slide, approx_cycle_length]=MapEKG2Color(ecg_wave', model_info.num_elec, 1, length(ecg_wave));
[first_ecg_idx, end_ecg_idx, int_cycles]=isolateCardiacCycle(mappedWaves,num_cycles);  num_ecg_slides=end_ecg_idx-first_ecg_idx; cardiac_length=num_ecg_slides/Fs;  
if num_cycles>1
    [~, R_slide_idx]=findpeaks(ecg_wave(first_ecg_idx:end_ecg_idx), 'MinPeakDistance',550);
end
sigma0=sigma_b;
frames_range=[first_ecg_idx, end_ecg_idx, int_cycles'];
cnt=1;
for t=frames_range(1):frames_range(2)
    [wave_color, ekg_wave_color]=EKG2Color(mappedWaves(t));
    plot((cnt:cnt+1)./Fs, ecg_wave(t:t+1), ekg_wave_color);hold on;
    cnt=cnt+1;
end

title("ECG Voltage Difference of Electrodes 16 and 25");
xlabel("Time (s)"); ylabel("Voltage (mV)");
grid on;

% save("ECGframes_wEIT.mat", 'frames_range'); save("ECGmapframes_wEIT.mat", 'mappedWaves');
save("ECGframes_0Amp.mat", 'frames_range'); save("ECGmapframes_0Amp.mat", 'mappedWaves');  

%% Interpolate Data for ECG Body Surface Map 
if new_BSM_interp
    Vq(:,:, first_ecg_idx:end_ecg_idx)=InterpolateVekgBSM(model_info, Vekg(:, first_ecg_idx:end_ecg_idx), true);  save("interpV.mat", 'Vq', '-v7.3');
end
load("interpV.mat");


%% Generate Heart Region from blood fraction information
if new_H_find
    [Bh, H_CH, LL_CH, RL_CH, wiggers_data]=VentillationPerfusionIdx(model_info.R, gammaR, cardiac_length, ecg_wave, first_ecg_idx, end_ecg_idx, int_cycles, end_qrs, sigma_b, true);  segment_info.H=Bh;  segment_info.LLch=LL_CH;    segment_info.RLch=RL_CH;
    save("VQBall.mat", 'Bh');  save("VQHch.mat", 'H_CH'); save("VQLLch.mat", 'LL_CH');    save("VQRLch.mat", 'RL_CH');    save("segmentation_info.mat", 'segment_info');
    save("WiggersData.mat", 'wiggers_data');
    Bh(1)=Bh(1)+0.02; Bh(2)=Bh(2)+0.02;
    [Qc, multi_Qc]=FindBestQ0([], Bh, elec_pts, model_info, Vekg(:, event_frames), "FitFwdSims","FEM", true, num_recon_src, "sbj_data_");      save("BestFitQc_multi"+num2str(num_recon_src)+".mat", 'multi_Qc'); save("BestFitQ0c"+num2str(num_recon_src)+".mat", 'Qc');
end

load("VQBall.mat"); load("VQHch.mat");  load("VQLLch.mat"); load("VQRLch.mat"); load("segmentation_info.mat");  load("WiggersData.mat");


load("BestFitQc_multi"+num2str(num_recon_src)+".mat");
load("BestFitQ0c"+num2str(num_recon_src)+".mat");



event_frames=first_ecg_idx:end_ecg_idx;
num_recon_src=1;


%% Generate Forward Map G
if new_G_map
    G=EKGFwdMap(model_info, Qc, 1, "FEM");        save("FEMMap_simGc.mat", 'G');
end
load("FEMMap_simGc.mat");


%% Compute Reconstructions of TCV



cnt=1;

for s=frames_range(1):frames_range(2)
    Vs=Vekg(:,s);
    Vs_0Amp=Vekg_0Amp(:,s);

    sigma0=real(sigma_b(1));    sigma_g=0.2;
    m_0Amp=sigma_g*((G'*G)\G')*Vs_0Amp;
    m=sigma0*((G'*G)\G')*Vs;

    Vcsim=(1/(sigma0)).*G*m; SimulatedDataC(:,s)=Vcsim;

    Mc_0Amp(s,:)=m;
    Mc(s,:)=m;
   sim2data_corr_const(s)=GetCorrelationValue(Vcsim, Vs);

    cnt=cnt+1;
end
save("TCV_0Amp.mat", 'Mc_0Amp');
% save("TCV_wEIT.mat", 'Mc');

ec_bfcc=CC_RelativeError(SimulatedDataC(:,event_frames), Vekg(:,event_frames))




%% Plot Total Cardiac Vector

TCVicons=["+", "*", "o", "hexagram", "diamond"];

figure;
tiledlayout(3,2);
ax1=nexttile(1, [2,1]);
title("Path of Total Cardiac Vector During Cardiac Cycle Using EIT Data");
xlabel("X (Amp meters)");    ylabel("Y (Amp meters)");    zlabel("Z (Amp meters)");    view([45, 45]);
ax1.FontSize=16;
% legend('','Location', 'eastoutside')

grid on;hold on;

ax2=nexttile(2, [2,1]);
title("Path of Total Cardiac Vector During Cardiac Cycle W/Out EIT Data")
xlabel("X (Amp meters)");    ylabel("Y (Amp meters)");    zlabel("Z (Amp meters)");    view([45, 45]);
ax2.FontSize=16;
% legend('','Location', 'eastoutside')

grid on;hold on;


ax3=nexttile(5);
title("Lead 2 ECG W/ EIT Data");
xlabel("Time (s)"); ylabel("Voltage (V)");
ax3.FontSize=12;
grid on; hold on;

ax4=nexttile(6);
title("Lead 2 ECG W/Out EIT Data");
xlabel("Time (s)"); ylabel("Voltage (V)");
ax4.FontSize=12;
grid on; hold on;

load("TCV_0Amp.mat");
load("ECGframes_0Amp.mat", 'frames_range'); int_cycles=frames_range(3:end);
load("ECGmapframes_0Amp.mat", 'mappedWaves');
cnt=1;  end_cycle=int_cycles(1);    
for t=frames_range(1):frames_range(2)
    if t<=end_cycle
        wave_marker=TCVicons(cnt);
    else
        if cnt==length(int_cycles)
            end_cycle=end_ecg_idx;
        else
            cnt=cnt+1;
            end_cycle=int_cycles(cnt);
        end
    end
    [wave_color, ekg_wave_color, color]=EKG2Color(mappedWaves(t)); 
    TCV_marker=color+wave_marker;

    Mc_0Amp_plt=OrderSrcVectors(num_recon_src, Mc_0Amp(t,:));
    for k=1:num_recon_src
        plot3(ax2, Mc_0Amp_plt(k,1), Mc_0Amp_plt(k,2), Mc_0Amp_plt(k,3), TCV_marker, 'MarkerSize', 7); hold on;
    end

    plot(ax4, (t:t+1)./Fs, L2ekg_0Amp(t:t+1), ekg_wave_color);hold on;
end

load("TCV_wEIT.mat");
load("L2ekg_wEIT.mat");
load("ECGframes_wEIT.mat", 'frames_range'); int_cycles=frames_range(3:end);
load("ECGmapframes_wEIT.mat", 'mappedWaves');
cnt=1;  end_cycle=int_cycles(1);    
for t=frames_range(1):frames_range(2)
    if t<=end_cycle
        wave_marker=TCVicons(cnt);
    else
        if cnt==length(int_cycles)
            end_cycle=end_ecg_idx;
        else
            cnt=cnt+1;
            end_cycle=int_cycles(cnt);
        end
    end

    [wave_color, ekg_wave_color, color]=EKG2Color(mappedWaves(t));
    TCV_marker=color+wave_marker;

    Mc_plt=OrderSrcVectors(num_recon_src, Mc(t,:));
    for k=1:num_recon_src
        plot3(ax1, Mc_plt(k,1), Mc_plt(k,2), Mc_plt(k,3), TCV_marker, 'MarkerSize', 7); hold on;
    end

    plot(ax3, (t:t+1)./Fs, ecg_wave(t:t+1), ekg_wave_color);hold on;
end

grid on; hold on;