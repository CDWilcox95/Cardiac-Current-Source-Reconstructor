function M=OrderSrcVectors(num_src, M)

M0=zeros(num_src,3);

for k=1:num_src

    M0(k,:)=M(3*k-2:3*k);
end

M=M0;

end