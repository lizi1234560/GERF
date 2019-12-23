function [Acc, JQ_Acc, trainsam, z] = RF_CART(FC, Chrom, N)

    EMCI = -ones(1, 37);
    NC = ones(1, 36);
    labels = [EMCI, NC];
    n = 0.6 * 73;
    x = randperm(73);
    y = x(1:n);          
    z = x(n+1:73);      
    
    test_num = length(z);
    sam_num = 40;               
    fc_num = 57;

    output = zeros(N, test_num);
    trainsam = nan(N, sam_num);   

    for i = 1:N
        sample_num = randsample(length(y), sam_num, 'true'); 
        trainsam(i, :) = y(sample_num);                      
        train_lab = labels(trainsam(i, :));                  

        temp = find(Chrom(i,:) == 1);      
        fc(i, :) = randsample(length(temp), fc_num);
        fc_temp = fc(i, :);                   

        train_set = FC(trainsam(i, :), fc_temp);   
        tree = fitctree(train_set, train_lab);   

        test_set = FC(z, fc_temp);               
        test_lab = labels(z);                    
        output(i, :) = predict(tree, test_set);  
        Acc(i) = sum(test_lab == output(i, :))/length(test_lab);

    end
    
    Res_sum = sum(output);
    Result = sign(Res_sum);
    R = sum(Result == test_lab);
    ['×¼È·ÂÊ£º' num2str(R/length(z))]
    JQ_Acc = R/length(z);
    
end