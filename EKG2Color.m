function [wave_color, ekg_color, color]=EKG2Color(EKG_Wave)

    
    if EKG_Wave=="p"
        color="r";
        wave_color="ro";
        ekg_color="r-";
    elseif EKG_Wave=="qrs"
        color="b";
        wave_color="bo";
        ekg_color="b-";
    elseif EKG_Wave=="t"
        color="k";
        wave_color="ko";
        ekg_color="k-";
    else
        fprintf("Region is not Labeled.  Plotting in default (blue)\n");
        wave_color="yo";
        ekg_color="y-";
    end




end