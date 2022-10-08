load('mnist.mat');
axis equal
count = 0;
digits= zeros(28,28,1);
for i=1:60000
    if labels_train(i) == 0
        count = count+1;
        digits(:,:,count) = digits_train(:,:,i);
    end
end
digits= cast(digits,'double');
data_matrix = reshape(digits,28*28,count);
mean_vector = sum(data_matrix,2)/count;
mean_matrix = ones(28*28,count);
    for i = 1:count
        mean_matrix(:,i) = mean_vector;
    end
cov_matrix = (data_matrix-mean_matrix)*((data_matrix-mean_matrix).')/(count-1);
[V,D] = eig(cov_matrix);
[d,ind] = sort(abs(diag(D)),'descend');
principle_eigenval= d(1,1);
principle_eigenvector = V(:,ind(1,1));
D_sorted = D(ind,ind);
V_sorted = V(:,ind);
V_compressed = V(:,1:84);
Projection = ((data_matrix-mean_matrix).')*V_compressed;
Reprojection = (Projection * V_compressed.').';
imagesc(reshape(digits(:,:,1),28,28));
imagesc(reshape(Reprojection(:,1),28,28));

