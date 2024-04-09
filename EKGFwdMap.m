function G=EKGFwdMap(model_info, Q0, sigma_dist, DOMAIN_CONST)
%% Input:
% R=Cylinder Radius
% H=Cylinder Height
% Q0=Current Source Point
% fmdl=Mesh Information
% sigma0=Reconstructed Conductivity
% L=Number of Electrodes

%% Construct Forward Map

% load("sigmaOmega_simdata.mat"); 

R=model_info.R; H=model_info.H; 

plt_orientation=false;

e1=[1,0,0]'; e2=[0,1,0]'; e3=[0,0,1]';  E=[e1, e2, e3];         %Elementary Basis
[K,~]=size(Q0);
el_area=0.0254^2;

switch DOMAIN_CONST
    case "Cylinder"     
        fmdl=model_info.FEM_Mesh;   L=model_info.num_elec;
        Pt=fmdl.nodes;  elec_pts=Pt(fmdl.electrode_idx,:);
        fprintf("Constructing Cylinder Forward Map \n");
        line_length=fprintf("Number of Columns Computed:  %d/%d", 0, K*3);

        electrode_nodes=fmdl.electrode_idx;

        col_cnt=1;
        for k=1:K
            Qk=Q0(k,:);
            for i=1:3
                 [~, ~, ~, Ve]=CalculateEKG_SolnCyl(R, fmdl, E(:,i)', 1, Qk,H, 15, true);
                Ve=Ve-sum(Ve)/L;                    %Normalize Voltages
                G(:,col_cnt)=Ve;

                fprintf(repmat('\b',1,line_length));
                line_length=fprintf("Number of Elements Integrated Over:  %d/%d",col_cnt, K*3);
                col_cnt=col_cnt+1;
            end

        end
        
        fprintf("\n Cylinder Forward Map Computed! \n");
    case "FEM"

        fmdl=GetReconMesh(R,H, false, model_info.shape);


        Pt=fmdl.nodes; 
        electrode_nodes=fmdl.electrode_idx;
        T=fmdl.elems;                           % Mesh Elements       (T)
        N=length(Pt); 

        C=fmdl.coef;    vol_tetr=fmdl.volumes;


        
        if plt_orientation
            theta_lH=(2*pi/16)*13; xH=(model_info.R/2)*cos(theta_lH);   yH=(model_info.R/2)*sin(theta_lH); zH=model_info.H/2;
            x0=[xH, yH, zH];

            figure;
            t=linspace(0,2*pi, 1000);
            x=R*cos(t);    y=R.*sin(t);
            plot3(x, y, 0.04.*ones(size(x)), 'k-'); hold on;
            plot3(x, y, 0.11.*ones(size(x)), 'k-'); hold on;

            plot3(Q0(1), Q0(2), Q0(3), 'r+', 'MarkerSize', 10); hold on;

            plot3(x0(1), x0(2), x0(3), 'g+', 'MarkerSize', 10); hold on;
            for l=1:length(electrode_nodes)
                el=Pt(electrode_nodes(l),:);
                if l<=length(electrode_nodes)/2
                    zl=0.04;
                else
                    zl=0.11;
                end
                text(el(1), el(2), zl, num2str(l)); hold on;
            end
            xlabel("X (m)"); ylabel("Y (m)"); zlabel("Z (m)"); 

        end

        %Forward Map Setting
        EVAL_NODES=electrode_nodes;
        G=zeros(length(EVAL_NODES),3);

        if length(sigma_dist)==1

            rhs_setting="const";
            sigma_Omega=ones(length(T),1);
            sigma0=1;
        elseif length(sigma_dist)==4
            rhs_setting="var";
            sigma_Omega=SetCondBodySims(R,H, fmdl, sigma_dist(1), sigma_dist(2), sigma_dist(3), sigma_dist(4), false);
            
            % true_sigma=sigmaOmega_t;
            elem_ind=1:length(T);
            fmdl.heart_elems=elem_ind(sigma_Omega==sigma_dist(2));
            sigma0=sigma_dist(2);
        else
            rhs_setting="var";
            sigma_Omega=sigma_dist;
            Bh=model_info.seg_info.H;
            Helems=GetHeartElemsFromSpecPt(Bh(1:3), Bh(4), fmdl);

            fmdl.heart_elems=Helems;
            sigma0=mean(sigma_Omega(Helems));
            sigma_Omega(Helems)=sigma0;
        end         


        mesh_data.elem_labels=fmdl.mat_idx;
        mesh_data.Pts=Pt;   mesh_data.radius=R; mesh_data.height=R; mesh_data.total_elems=T;
        switch model_info.shape
            case 'c'
                S=Create3DStiffMatrix(C, vol_tetr, sigma_Omega, T, Pt, "FWD_MAP");  save("cyl_FWD_MAPFEMStiffMatrix.mat", 'S');
                load("cyl_FWD_MAPFEMStiffMatrix.mat");

            case 'e'
                S=Create3DStiffMatrix(C, vol_tetr, sigma_Omega, T, Pt, "FWD_MAP");  save("ell_FWD_MAPFEMStiffMatrix.mat", 'S');
                load("ell_FWD_MAPFEMStiffMatrix.mat");

        end
        % load("FWD_MAPFEMStiffMatrix.mat");


        g=CreateUniqueVector("pt_elec", electrode_nodes, fmdl, Pt, T);


        col_cnt=1;
        for k=1:K
            Qk=Q0(k,:);

            for i=1:3

                b=CalculateRHSvec(rhs_setting, T, Pt, Qk, sigma_dist, C, vol_tetr, fmdl, E(:,i)');

                %% Compute Potentials Over Entire Mesh
                S(1:N, N+1)=g;    S(N+1,1:N)=g';    b(N+1)=0;
                U_ki=S\b;  V_ki=U_ki(EVAL_NODES);
                
                u0_mesh=Computeu0OnMesh(fmdl, sigma0, Qk, E(:,i)');

                U0=ComputeSolnVoltages(el_area, fmdl, U_ki, u0_mesh);

                U0inf=ComputeSolnVoltages(el_area, fmdl, [], u0_mesh);

                G(:,col_cnt)=U0-sum(U0)/model_info.num_elec;

                % G(:,col_cnt)=G0(:,col_cnt)+V_ki;   G(:,col_cnt)=G(:,col_cnt)-sum(G(:,col_cnt))/length(EVAL_NODES);
                col_cnt=col_cnt+1;
            end
        end


    case "Sbj"
        fmdl=GetReconMesh(R,H, false, model_info.shape);


        Pt=fmdl.nodes;
        electrode_nodes=fmdl.electrode_idx;
        T=fmdl.elems;                           % Mesh Elements       (T)
        N=length(Pt);

        C=fmdl.coef;    vol_tetr=fmdl.volumes;


        %Forward Map Setting
        EVAL_NODES=electrode_nodes;
        G=zeros(length(EVAL_NODES),3);
        idx=model_info.frame;

        if length(sigma_dist)==1

            rhs_setting="const";
            sigma_Omega=ones(length(T),1);
            sigma0=1;
        elseif length(sigma_dist)==4
            rhs_setting="var";
            sigma_Omega=SetCondBodySims(R,H, fmdl, sigma_dist(1), sigma_dist(2), sigma_dist(3), sigma_dist(4));

            % true_sigma=sigmaOmega_t;
            elem_ind=1:length(T);
            fmdl.heart_elems=elem_ind(sigma_Omega==sigma_dist(2));
            sigma0=sigma_dist(2);
        else
            rhs_setting="var";
            sigma_Omega=sigma_dist;
            Bh=model_info.seg_info.H;
            Helems=GetHeartElemsFromSpecPt(Bh(1:3), Bh(4), fmdl);

            fmdl.heart_elems=Helems;
            sigma0=mean(sigma_Omega(Helems));
            sigma_Omega(Helems)=sigma0;
        end


        mesh_data.elem_labels=fmdl.mat_idx;
        mesh_data.Pts=Pt;   mesh_data.radius=R; mesh_data.height=R; mesh_data.total_elems=T;
        switch model_info.shape
            case 'c'
                
                % S=Create3DStiffMatrix(C, vol_tetr, sigma_Omega, T, Pt, "FWD_MAP");  save("cyl_FWD_MAPFEMStiffMatrix.mat", 'S');
                % load("cyl_FWD_MAPFEMStiffMatrix.mat");
                load("C:\Users\cdwil\Desktop\EIT_EKG_Projects\Cardiac_Imaging\varCondStiffMat\St_"+num2str(idx)+".mat");
            case 'e'
                S=Create3DStiffMatrix(C, vol_tetr, sigma_Omega, T, Pt, "FWD_MAP");  save("ell_FWD_MAPFEMStiffMatrix.mat", 'S');
                load("ell_FWD_MAPFEMStiffMatrix.mat");

        end
        % load("FWD_MAPFEMStiffMatrix.mat");


        g=CreateUniqueVector("pt_elec", electrode_nodes, fmdl, Pt, T);


        col_cnt=1;
        for k=1:K
            Qk=Q0(k,:);

            for i=1:3

                b=CalculateRHSvec(rhs_setting, T, Pt, Qk, sigma_dist, C, vol_tetr, fmdl, E(:,i)');

                %% Compute Potentials Over Entire Mesh
                S(1:N, N+1)=g;    S(N+1,1:N)=g';    b(N+1)=0;
                U_ki=S\b;  V_ki=U_ki(EVAL_NODES);

                u0_mesh=Computeu0OnMesh(fmdl, sigma0, Qk, E(:,i)');

                U0=ComputeSolnVoltages(el_area, fmdl, U_ki, u0_mesh);
                G(:,col_cnt)=U0-sum(U0)/model_info.num_elec;

                % G(:,col_cnt)=G0(:,col_cnt)+V_ki;   G(:,col_cnt)=G(:,col_cnt)-sum(G(:,col_cnt))/length(EVAL_NODES);
                col_cnt=col_cnt+1;
            end
        end

    case "Infinite"
        
        fmdl=model_info.FEM_Mesh;   L=model_info.num_elec;
        fmdl=GetReconMesh(R,H, false, 'c');



        Pt=fmdl.nodes;  elec_pts=Pt(fmdl.electrode_idx,:);

        cnt=1;
        for k=1:K
            Qk=Q0(k,:);
            for i=1:3
                u0_mesh=Computeu0OnMesh(fmdl, 1, Qk, E(:,i)');
                V(:,cnt)=ComputeSolnVoltages(el_area, fmdl, [], u0_mesh);
                cnt=cnt+1;
            end
        end
        G=V;
        save("InfMap_G.mat", 'G');

end

load("u0_vals.mat");    load("fem_solns.mat");
%Orthogonalize Columns

end