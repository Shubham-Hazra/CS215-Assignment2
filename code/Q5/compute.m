function coordinates = compute(digits_train,labels_train)
    count = 0;
    digits= zeros(28,28,1); %Matrix to store the data for each digit
    for i=1:60000
        if labels_train(i) == 0  %For the digit 0
            count = count+1;
            digits(:,:,count) = digits_train(:,:,i);
        end
    end
    digits= cast(digits,'double'); %Casting to make the integers into floating point 
    data_matrix = reshape(digits,28*28,count); %Reshaping to get a 784x1 coulumn vector for each image
    mean_vector = sum(data_matrix,2)/count; %Calculating the mean vector
    mean_matrix = ones(28*28,count); %Mean matrix i.e. each column is the mean vector
        for i = 1:count
            mean_matrix(:,i) = mean_vector;
        end
    cov_matrix = (data_matrix-mean_matrix)*((data_matrix-mean_matrix)')/(count-1); %Calculating the covariance matrix
    H = (data_matrix-mean_matrix)/sqrt(count-1); %Defining H matrix so as to calculate eigen vectors using H*H' which gives only real values
    [V,D] = eig(H*H'); %Calculating this way because other ways might give complex eigen vectors due to numerical errors 
    [~,ind] = sort(abs(diag(D)));
    coordinates = sort(ind(1:84,1)); %The coordinates of the top 84 eigenvalues 
end