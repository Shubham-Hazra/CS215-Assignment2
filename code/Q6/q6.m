data = ones(19200,1);
folder = 'C:\Users\shubh\OneDrive\Desktop\CS215 Assignment2\data\data_fruit' ;
files = [folder filesep '\*.png'] ;
files = dir(files);
for k=1:16
    fullFileName = fullfile(folder, files(k).name);
    idx=k;
    im = files(idx).name;
    img = imread(fullfile(folder,im));
    data(:,k) = reshape(img,19200,1);
end
mean_vector = sum(data,2)/16;
mean_matrix = ones(19200,16);
    for i = 1:16
        mean_matrix(:,i) = mean_vector;
    end
cov_matrix = (data - mean_matrix)*((data-mean_matrix)')/15;
H = (data-mean_matrix)/sqrt(15);
[V,~] = eigs(H*H',4);
subplot(3,4, [2 3])
image(rescale(reshape(mean_vector,80,80,3),0,1))
title 'Mean Image'
for i=1:4
subplot(3,4,[3+2*i 4+2*i])
image(rescale(reshape(V(:,i),80,80,3),0,1))
end
subplot(1,1,1)
[V,D] = eigs(H*H',10);
[d,ind] = sort(abs(diag(D)));
plot(d,'-')
xlabel 'number of eigenvalues'
ylabel 'Eigenvalues'
title 'Top 10 Eigenvalues'
