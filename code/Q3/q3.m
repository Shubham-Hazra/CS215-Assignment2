load('points2D_set1.mat');
data_matrix(1,:) = x';
data_matrix(2,:) = y';
[count, ~] = size(x);
mean_vector = sum(data_matrix,2)/count;
mean_matrix = ones(2,count);
    for i = 1:count
        mean_matrix(:,i) = mean_vector;
    end
cov_matrix = (data_matrix-mean_matrix)*((data_matrix-mean_matrix)')/(count-1);
H = (data_matrix-mean_matrix)/sqrt(count-1);
[V,D] = eig(H*H');
[d,ind] = sort(abs(diag(D)),'descend');
principle_eigenval= d(1,1);
principle_eigenvector = V(:,ind(1,1));
scatter(x,y);
x = linspace(0,1,1000);
slope = principle_eigenvector(2,1)/principle_eigenvector(1,1);
c = mean_vector(2,1) - slope*mean_vector(1,1);
y = slope.*x + c;
hold on
plot(x,y,'r-',LineWidth=3)
xlabel 'X values'
ylabel 'Y values'
title 'Points2D set1'
legend('','Calculated Linear relation b/w x and y')
hold off
clear
load('points2D_set2.mat');
data_matrix(1,:) = x';
data_matrix(2,:) = y';
[count, ~] = size(x);
mean_vector = sum(data_matrix,2)/count;
mean_matrix = ones(2,count);
    for i = 1:count
        mean_matrix(:,i) = mean_vector;
    end
cov_matrix = (data_matrix-mean_matrix)*((data_matrix-mean_matrix)')/(count-1);
H = (data_matrix-mean_matrix)/sqrt(count-1);
[V,D] = eig(H*H');
[d,ind] = sort(abs(diag(D)),'descend');
principle_eigenval= d(1,1);
principle_eigenvector = V(:,ind(1,1));
scatter(x,y);
x = linspace(-2,2,1000);
slope = principle_eigenvector(2,1)/principle_eigenvector(1,1);
c = mean_vector(2,1) - slope*mean_vector(1,1);
y = slope.*x + c;
hold on
plot(x,y,'r-',LineWidth=3)
xlabel 'X values'
ylabel 'Y values'
title 'Points2D set2'
legend('','Calculated Linear relation b/w x and y')
hold off
clear


