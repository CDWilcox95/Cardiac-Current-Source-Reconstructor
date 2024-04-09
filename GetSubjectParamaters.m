function model_info=GetSubjectParamaters(Circum, vert_gap)
    unit=1;

    R=Circum*2.54/(2*pi);                                 %cm
    H=vert_gap+8;            H=R;                           %Add 4cm above and below electrode planes
    Z1=4;   Z2=Z1+vert_gap;    Z=[Z1, Z2];              %Z-Planes Where Electrodes Are Placed
    L=32;                                               %Number of Total Electrodes
    Fs = 864.0553;                                      %Sampling Frequency (1/s)
    
    % Load Thorax Mesh
    load("Thorax_FEMmodel.mat");            % Thorax Data
    load("FEM_ThoraxMesh_elec_nodes.mat");  % Electrode Nodes     (el)
    

    Pt=fmdl.nodes;  Pt(:,1:2)=R.*Pt(:,1:2); Pt(:,3)=H.*Pt(:,3); fmdl.nodes=Pt;  
    [C, vol_tetr]=CalculateMeshCoef(fmdl.elems, Pt);             

    fmdl.mesh_coef=C;   fmdl.volumes=vol_tetr;
    fmdl.electrode_nodes=electrode_nodes;
    model_info.R=R; model_info.H=H; model_info.elec_planes=Z;   model_info.num_elec=L;  
    model_info.sample_rate=Fs;  

    model_info.elec_width=2.54; model_info.elec_height=3.4; model_info.FEM_Mesh=fmdl;


    shift=2*pi/(L/2);
    theta_l=(2*pi/(L/2)).*[1:L/2]'-shift;
    elec_pts=[R.*cos(theta_l), R.*sin(theta_l), Z1.*ones(L/2,1); R.*cos(theta_l), R.*sin(theta_l), Z2.*ones(L/2,1)];

    model_info.elec_pts=elec_pts;

    model_info.elec_planes=Z;
    model_info.shape='c';
end