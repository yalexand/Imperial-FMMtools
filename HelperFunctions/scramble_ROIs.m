function scrambled = scramble_ROIs( sgm )

    Lb1 = bwlabel(sgm);
    Lb2 = bwlabel(~sgm);

    N1 = max(Lb1);
    N2 = max(Lb2);

    i1 = randperm(N1);
    i2 = randperm(N2);

    z = []; % destination    
    
    for k=1:N1
        z = [z ones(1,numel(find(Lb1==i1(k))))]; % drawback - always start from ROI
        if k<=N2
            z = [z zeros(1,numel(find(Lb2==i2(k))))];
        end    
    end

    if length(z)<length(sgm)
        z = [z zeros(1,numel(find(Lb2==i2(N2))))];
    end

    scrambled = z;
   
end

