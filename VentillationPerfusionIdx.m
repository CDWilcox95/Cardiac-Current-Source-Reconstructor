function [Bh, Heart_CH, Llung_CH, Rlung_CH, wiggers_data]=VentillationPerfusionIdx(R, sigma, cardiac_length, ecg_wave, start_slide, end_slide, int_cycles, diff_idx, sigma_b, plot_region)
    %% Input:  
    %   sigma=Reconstructed Conductivity
    %   T=Mesh Reconstruction Performed On
    % Assumptions:
    %   Data Collected at 27*32 sampling frequency so each slide represents ~1ms

    %%Collect Maximum and Minimum Conductivities over all voxels at each
    %%slide

    sigma_R=real(sigma);
    load joshmap3D_circle;

    cd EKGEIT_V2;
    Fs=864.0553;

    [M,numSlides]=size(sigma_R);
    f_a=zeros(M, numSlides); f_b=f_a;
    
    [Vol,Pt]=GetVoxelVolumeJTMesh3D(R); cd ../;
    num_elems=length(Vol); mesh_idx=1:num_elems;    top_layer=1:num_elems/2;    bot_layer=num_elems/2+1:num_elems;


    ElemNums=zeros(length(Vol),1);
        
    
    cd ../;
    [s0, sN]=isolateRespiratoryCycle(sigma_b);    
    
    
    
    
    dt_b=cardiac_length;
    dt_a=4.4;

    
    [S0, SN]=GetStartEndFrames(start_slide, int_cycles, end_slide);
    if ~isempty(int_cycles)
        num_cycles=length(int_cycles)+1;
    else
        num_cycles=1;
    end
    

    %% Calculate Fraction of Air (f_a)  and Fraction of Blood (f_b)

    [sigma_cc_max, sigma_cc_min]=ComputeExtCycleCond(sigma_R, S0, SN);
   
    cnt=1;  i0=S0(1);   iN=SN(1);   f_b=zeros(M, end_slide);   sigma_ref=zeros(size(sigma_R));
    for s=start_slide:end_slide
        if s>= i0 && s<=iN
            ref_idx=diff_idx(1);
        else
            cnt=cnt+1;
            i0=S0(cnt); iN=SN(cnt);
            ref_idx=diff_idx(1);
        end
        sigma_ref(:,s)=sigma_R(:,ref_idx);

    end

    f_b(:, 1:end_slide)=ComputeVolumeFractions(sigma_R, min(sigma_cc_min), max(sigma_cc_max), S0, SN, 'b');



    S0rr=[1569, 4537,  8130];   SNrr=[4536, 8129, 12331];
    [sigma_rr_max, sigma_rr_min]=ComputeExtCycleCond(sigma_R, S0rr, SNrr);

    f_a=ComputeVolumeFractions(sigma_R, sigma_rr_min, sigma_rr_max, S0rr, SNrr, 'a');

    %% Calculate Reference Frames
    % Max Blood Volume in Heart should be at R Slide--Max Ventricular
    % Depolarization.  
    % fb_ref=zeros(num_elems,numSlides);  sigma_ref=zeros(num_elems, numSlides);
    cnt=1;  i0=S0(1);   iN=SN(1); 
    for s=start_slide:end_slide
        if s>= i0 && s<=iN
            ref_idx=diff_idx(1);
        else
            cnt=cnt+1;
            i0=S0(cnt); iN=SN(cnt);
            ref_idx=diff_idx(1);
        end
        fb_ref(:, s)=f_b(:, ref_idx);
    end
    

    %% Isolate Voxel Signals
    elem0=62;       % Heuristic guess for heart element
%     elem0=317;
    elem0=[267,892];       % Heuristic guess for heart element
    elem0=267;       % Heuristic guess for heart element

    search_num=60;
    adj_elems_bot=AdjacencyMatrix(Jash, elem0, search_num);
    adj_elems_top=AdjacencyMatrix(Jash, 892, search_num);
    adj_elems=union(adj_elems_top, adj_elems_bot);

    LLelem0=[177, 720];
    % RLelem0=195;
    RLelem0=[218,773];


    heart_elemscorr=GetHeartRegionByCorrelationVals(sigma(:, start_slide:end_slide)-real(sigma_b(start_slide:end_slide)), mesh_idx, elem0, 0.98, true);
    heart_elemscorr=mesh_idx(heart_elemscorr);

    plt_heart_elems=zeros(M,1);
    plt_heart_elems(heart_elemscorr)=1;  plt_heart_elems(elem0)=5;
    % cd EKGEIT_V2; plt_onto_JoshMesh3D(plt_heart_elems, plt_heart_elems, [],[], "Binary Plot Lung Elements", "Binary Plot Heart Elements", false); cd ../;

    %% Lung Segmentation 
    lung_adj_elems_LTOP=AdjacencyMatrix(Jash, LLelem0(1), 80);
    lung_adj_elems_LBOT=AdjacencyMatrix(Jash, LLelem0(2), 80);

    lung_adj_elems_RTOP=AdjacencyMatrix(Jash, RLelem0(1), 80);
    lung_adj_elems_RBOT=AdjacencyMatrix(Jash, RLelem0(2), 80);


    lung_adj_elems_left=union(lung_adj_elems_LTOP, lung_adj_elems_LBOT);
    lung_adj_elems_right=union(lung_adj_elems_RTOP, lung_adj_elems_RBOT);



    Llung_elems1=GetHeartRegionByCorrelationVals(sigma(:, start_slide:end_slide), lung_adj_elems_left, LLelem0(1), 0.7, false);
    Llung_elems2=GetHeartRegionByCorrelationVals(sigma(:, start_slide:end_slide), lung_adj_elems_left, LLelem0(2), 0.7, false);
    Llung_elems=union(mesh_idx(Llung_elems1), mesh_idx(Llung_elems2));


    Rlung_elems1=GetHeartRegionByCorrelationVals(sigma(:, start_slide:end_slide), lung_adj_elems_right, RLelem0(1), 0.7, false);
    Rlung_elems2=GetHeartRegionByCorrelationVals(sigma(:, start_slide:end_slide), lung_adj_elems_right, RLelem0(2), 0.7, false);
    Rlung_elems=union(mesh_idx(Rlung_elems1), mesh_idx(Rlung_elems2));


    plt_adj_elems=zeros(M,1);    plt_lung_elems=zeros(M,1);
    plt_adj_elems(adj_elems)=1;  
    plt_lung_elems(Llung_elems)=2;  plt_lung_elems(LLelem0)=5;  
    plt_lung_elems(Rlung_elems)=3;  plt_lung_elems(RLelem0)=5;

    corr_lung_elems=union(Rlung_elems, Llung_elems);


    lung_elems=union(mesh_idx(Llung_elems), mesh_idx(Rlung_elems));


    %% Segment Heart Elements

    heart_elems=GetLungElems(sigma, f_b, 0.65, S0, SN, mesh_idx);
    
    %Remove any overlapping elements or elements too close to boundary
    if ~isempty(intersect(lung_elems, heart_elems))
        lung_elems(lung_elems==intersect(lung_elems, heart_elems))=[];
    end
    lung_elems=RemoveBoundaryElems(lung_elems, 0.8);
    [Llung_elems, Rlung_elems]=SeperateLungRegions(lung_elems);

    %% Plot Segmented Regions onto Joshua Tree Mesh
    plt_heart_elems=zeros(M,1); plt_lung_elems=zeros(M,1);
    plt_lung_elems(Llung_elems)=2;  plt_lung_elems(LLelem0)=5;  plt_lung_elems(Rlung_elems)=3;  plt_lung_elems(RLelem0)=5;
    plt_heart_elems(heart_elemscorr)=1;  plt_heart_elems(elem0)=5;
    % cd EKGEIT_V2; plt_onto_JoshMesh3D(plt_lung_elems, plt_heart_elems, [],[], "Binary Plot Lung Elements", "Binary Plot Heart Elements", false); cd ../;



    % save("Heart_Elems.mat", 'heart_elems');   save("LeftLung_Elems.mat", 'Llung_elems');    save("RightLung_elems.mat", 'Rlung_elems');





    %% Plot mean Conductivity and Volume Fractions
    figure;tiledlayout(3,2);
    
    nexttile(1);
    xaxis=start_slide:end_slide;
    plot(xaxis, mean(f_b(heart_elems,xaxis),1),  'r-');
    title("Volume Fraction of Blood of Heart Element");
    xlabel("Frame Index"); ylabel("Volume Fraction");


    nexttile(2);
    plot(xaxis, ecg_wave(xaxis), 'k-');
    title("Lead 2 ECG");
    xlabel("Frame Index"); ylabel("Voltage (mV)");

    nexttile(3);
    Heart_cond=mean(sigma(heart_elems,xaxis),1);
    plot(xaxis, Heart_cond, 'c-');
    title("Conductivity of Heart Region");
    xlabel("Frame Index"); ylabel("Conductivity (mS/M)");


    nexttile(4);
    Llung_cond=mean(sigma(Llung_elems,xaxis),1);
    plot(xaxis, Llung_cond./norm(Llung_cond,2), 'b-', xaxis, Heart_cond./norm(Heart_cond), 'r-');
    title("Conductivity of Left Lung vs. Heart Conductivity");
    xlabel("Frame Index"); ylabel("Conductivity (mS/M)");
    legend("Avg Left Lung Conductivity", "Avg Heart Conductivity");


    nexttile(5);
    Rlung_cond=mean(sigma(Rlung_elems, xaxis),1);
    plot(xaxis, Rlung_cond./norm(Rlung_cond,2), 'b-', xaxis, sigma(elem0, xaxis)./norm(sigma(elem0, xaxis)), 'r-');
    title("Conductivity of Right Lung vs. Heart Conductivity");
    xlabel("Frame Index"); ylabel("Conductivity (mS/M)");
    legend("Avg Right Lung Conductivity", "Avg Heart Conductivity");

    nexttile(6);    plt_scale=max(abs(ecg_wave(xaxis)));
    plt_ecg=ecg_wave(xaxis);    plt_ecg=ecg_wave(xaxis)./plt_scale; 
    plt_heart_cond=Heart_cond;   plt_heart_cond=plt_heart_cond; plt_heart_cond=(1/max(abs(plt_heart_cond))).*plt_heart_cond;
    plt_lung_cond=Llung_cond;   plt_lung_cond=plt_lung_cond+(plt_ecg(1)-plt_lung_cond(1));      
    plot(xaxis, plt_ecg, 'k-', xaxis, plt_heart_cond, 'r-', xaxis, plt_lung_cond, 'b-');
    title("ECG - Avg Heart Conductivity - Avg Left Lung Conductivity");
    xlabel("Frame Index");  
    legend("ECG", "Avg Heart Conductivity", "Avg Left Lung Conductivity");



    figure;
    tiledlayout(3,1);
    xaxis=start_slide:end_slide;    vox_vol=Vol(1);    m3tomL=10^6;

    nexttile(1);
    plot(xaxis./Fs, sum(f_b(heart_elems,xaxis).*vox_vol*m3tomL,1),  'r-');
    title("Volume of Blood During Cardiac Cycles");
    xlabel("Time (s)"); ylabel("Volume (mL)");
    grid on;

    nexttile(2);
    plot(xaxis./Fs, mean(sigma(heart_elems,xaxis),1),  'b-');
    title("Average Conductivity During Cardiac Cycles");
    xlabel("Time (s) "); ylabel("Conductivity (S/M)");
    grid on;


    nexttile(3);
    plot(xaxis./Fs, ecg_wave(xaxis), 'k-');
    title("Lead 2 ECG");
    xlabel("Time (s)"); ylabel("Voltage (mV)");
    grid on;


    %% Compute and Plot Wiggers Diagram for Volume Flow 

    numHelems=length(heart_elems);   tidalH_fvol=f_b(heart_elems,xaxis);   
    tidalLL_vol=sum(f_b(Llung_elems, xaxis),1)*vox_vol*m3tomL;  tidalRL_vol=sum(f_b(Rlung_elems, xaxis),1)*vox_vol*m3tomL;

    global_vol=sum(f_b(:,xaxis),1)*vox_vol*m3tomL;
    wiggers_data=zeros(1, end_slide);

    wiggers_data(xaxis)=sum(tidalH_fvol, 1)*vox_vol*m3tomL;    
    % wiggers_data(xaxis)=wiggers_data(xaxis)-min(wiggers_data(xaxis));

    figure;
    tiledlayout(2,2);
    nexttile(1);
    title("Volume of Blood in  H");  xlabel("Time (s)"); ylabel("Volume (mL)");   grid on;  hold on;
    plot(xaxis./Fs, wiggers_data(xaxis), 'r-');

    nexttile(2);
    title("Volume of Blood in Left Lung");  xlabel("Time (s)"); ylabel("Volume (mL)");  grid on; hold on;
    plot(xaxis./Fs, tidalLL_vol-min(tidalLL_vol), 'b-', xaxis./Fs, wiggers_data(xaxis), 'r-');

    nexttile(3);
    title("Lead 2 ECG");    xlabel("Time (s)"); ylabel("Voltage (V)");  grid on; hold on;
    plot(xaxis./Fs, ecg_wave(xaxis), 'k-');


    nexttile(4);
    title("Volume of Blood in Right Lung");  xlabel("Time (s)"); ylabel("Volume (mL)");  grid on; hold on;
    plot(xaxis./Fs, tidalRL_vol-min(tidalRL_vol), 'b-',xaxis./Fs, wiggers_data(xaxis), 'r-');






    
    LLelem0=177;

    




    lung_cond=mean(sigma(Llung_elems,start_slide:end_slide), 1);    heart_cond=mean(sigma(heart_elems,start_slide:end_slide),1);

    plot(start_slide:end_slide, lung_cond./norm(lung_cond,2), 'g-', start_slide:end_slide, heart_cond./norm(heart_cond,2), 'r-');

    plt_adj_elems=zeros(M,1);               plt_heart_elems=zeros(M,1); plt_lung_elems=zeros(M,1);
    plt_adj_elems(adj_elems)=1;    plt_heart_elems(heart_elems)=1;  plt_adj_elems(elem0)=5; plt_heart_elems(elem0)=5;
    plt_lung_elems(Llung_elems)=2;  plt_lung_elems(LLelem0)=5;  
    plt_lung_elems(Rlung_elems)=3;  plt_lung_elems(RLelem0)=5;
    plt_lung_elems(heart_elems)=10; plt_lung_elems(elem0)=5;
    % cd EKGEIT_V2; plt_onto_JoshMesh3D(plt_lung_elems, plt_heart_elems, [],[], "Binary Plot Adjacent Elements", "Binary Plot Heart Elements", false); cd ../;
    % cd EKGEIT_V2; plt_onto_JoshMesh3D(plt_lung_elems, plt_heart_elems, [],[], "Binary Plot Adjacent Elements", "Binary Plot Heart Elements", false); cd ../;


    organ_elems=union(union(Llung_elems, Rlung_elems), heart_elems);


    %% Plot onto segmented heart
    fb_ref=zeros(num_elems,numSlides);
    cnt=1;  i0=S0(1);   iN=SN(1); 
    for s=start_slide:end_slide
        if s>= i0 && s<=iN
            ref_idx=diff_idx(cnt);
        else
            cnt=cnt+1;
            i0=S0(cnt); iN=SN(cnt);
            ref_idx=diff_idx(cnt);
        end
        fb_ref(:, s)=f_b(:, ref_idx);
    end
    close all;

        peak_ecg_frames=[6421, 6488,6520, 6552,6733]; % Sbj 63
    Pwave=peak_ecg_frames(1);
    Qwave=peak_ecg_frames(2);
    Rwave=peak_ecg_frames(3);
    Swave=peak_ecg_frames(4);
    Twave=peak_ecg_frames(5);
    
    plt_fb_seg=zeros(M,numSlides);
    plt_fb_seg(:,start_slide:end_slide)=f_b(:,start_slide:end_slide);
    % plt_fb_seg(Llung_elems, start_slide:32:end_slide)=f_b(Llung_elems,start_slide:32:end_slide)-fb_ref(Llung_elems, start_slide:32:end_slide);
    % plt_fb_seg(Rlung_elems, start_slide:32:end_slide)=f_b(Rlung_elems,start_slide:32:end_slide)-fb_ref(Rlung_elems, start_slide:32:end_slide);


    % cd EKGEIT_V2;    plt_onto_JoshMesh3D( f_b(:, start_slide:32:end_slide),[],  ecg_wave(start_slide:32:end_slide),wiggers_data(start_slide:32:end_slide),  "fb over Omega", "fb over H",false);   cd ../;   Llung_CH=GetLungRegion(Llung_elems, Pt);    Rlung_CH=GetLungRegion(Rlung_elems,Pt);     Heart_CH=GetLungRegion(heart_elems,Pt);
    cd EKGEIT_V2;    plt_onto_JoshMesh3D( plt_fb_seg(:, Pwave),[Pwave-start_slide, max(max(f_b(:,start_slide:end_slide))), min(min(f_b(:,start_slide:end_slide)))],  ecg_wave(start_slide:end_slide),wiggers_data(start_slide:end_slide),  "fb over Omega", "fb over H",false); 
    plt_onto_JoshMesh3D( plt_fb_seg(:, Qwave),[Qwave-start_slide, max(max(f_b(:,start_slide:end_slide))), min(min(f_b(:,start_slide:end_slide)))],  ecg_wave(start_slide:end_slide),wiggers_data(start_slide:end_slide),  "fb over Omega", "fb over H",false); 
    plt_onto_JoshMesh3D( plt_fb_seg(:, Rwave),[Rwave-start_slide,  max(max(f_b(:,start_slide:end_slide))), min(min(f_b(:,start_slide:end_slide)))],  ecg_wave(start_slide:end_slide),wiggers_data(start_slide:end_slide),  "fb over Omega", "fb over H",false); 
    plt_onto_JoshMesh3D( plt_fb_seg(:, Swave),[Swave-start_slide,  max(max(f_b(:,start_slide:end_slide))), min(min(f_b(:,start_slide:end_slide)))],  ecg_wave(start_slide:end_slide),wiggers_data(start_slide:end_slide),  "fb over Omega", "fb over H",false); 
    plt_onto_JoshMesh3D( plt_fb_seg(:, Twave),[Twave-start_slide,  max(max(f_b(:,start_slide:end_slide))), min(min(f_b(:,start_slide:end_slide)))],  ecg_wave(start_slide:end_slide),wiggers_data(start_slide:end_slide),  "fb over Omega", "fb over H",false); 

    cd ../;
    [c0, r0]=GetHeartPtsFromTodler(mesh_idx(heart_elems), Pt);
    Bh=[c0, r0];

    L=16;
    theta_l=(2*pi/L).*[1:L];
    GenerateSourceOnSphere(r0, c0, 1000); hold on;
    plot3(c0(1), c0(2), c0(3), 'r+','MarkerSize', 10); hold on;
    t=linspace(0,2*pi, 1000); rplt=0.159;   H=0.15;
    plot3(rplt.*cos(t), rplt.*sin(t), zeros(size(t)), 'k-');   hold on;
    plot3(rplt.*cos(t), rplt.*sin(t), H.*ones(size(t)), 'k-'); hold on;
    plot3(rplt.*cos(theta_l), rplt.*sin(theta_l), zeros(L,1), 'co', 'MarkerSize', 10); hold on;
    plot3(rplt.*cos(theta_l), rplt.*sin(theta_l), H.*ones(L,1), 'co', 'MarkerSize', 10); hold on;
    for i=1:L
       text(rplt.*cos(theta_l(i)), rplt.*sin(theta_l(i)), 0, num2str(i)); hold on;
       text(rplt.*cos(theta_l(i)), rplt.*sin(theta_l(i)), H, num2str(i+L)); hold on;
    end
    xlabel("X");    ylabel("Y");     zlabel("Z");
    
    
    
    %% Compute Ventillation and Perfusion Index
    if num_cycles>1

        num_CCcycles=num_cycles;   num_Rcycles=1;
        V=zeros(M, num_Rcycles);
        Q=zeros(M, num_CCcycles);



        i0=start_slide; iN=int_cycles(1);   Delta_t=(iN-i0)/Fs;
        for i=1:num_CCcycles
            for p=1:M
                [max_fb_val, max_fb_idx]=max(f_b(p, i0:iN));
                [min_fb_val, min_fb_idx]=min(f_b(p, i0:iN));

                Q(p,i)=(max_fb_val-min_fb_val)/Delta_t;

            end

            i0=S0(i);    iN=SN(i); Delta_t=(iN-i0)/Fs;
            if i==1
                plt_Q(:,i)=Q(:,i);
            else
                plt_Q(:,i)=Q(:,i)-Q(:,1);
            end
        end

        cd EKGEIT_V2; plt_onto_JoshMesh3D(Q, [], ecg_wave(start_slide:end_slide), "Perfusion", "Ventillation", false); cd ../;
    end
end