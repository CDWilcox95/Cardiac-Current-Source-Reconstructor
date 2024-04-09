function [Vr, Vi]=CorrectVoltageData(frame_voltage)
    [~, L, num_frames]=size(frame_voltage);
    num_curs=L-1;
    V=zeros(L,num_curs, num_frames);

    for s=1:num_frames
        volt_frame=frame_voltage(:,:,s);
        volt_frame=volt_frame(1:num_curs, :)';
        volt_frame=volt_frame-repmat(sum(volt_frame)/L, L, 1);
        V(:,:, s)=volt_frame;

    end

    Vr=real(V); Vi=imag(V);
end