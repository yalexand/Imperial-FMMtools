function y = dilate_labels( x )

y = zeros(size(x));

L = length(x);

z = watershed(bwdist(0~=x));

    for k=1:L         
        LX=x(k);
        if 0~=LX
            LZ=z(k);
            y(z==LZ)=LX;
        end
    end    
            
end

