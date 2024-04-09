function e=CC_RelativeError(U, V)
[L, num_frames]=size(U);

e=0;
num=0;  den=0;
for s=1:num_frames
    for l=1:L
        num=num+abs(U(l,s)-V(l,s))^2;
        den=den+abs(V(l,s))^2;
    end
end

e=sqrt(num/den);

end