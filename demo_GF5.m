clc
clear

%% load data
load('data/GF5_Shanghai.mat')
Y = im2double(data); % Y - noisy data
[M,N,B] = size(Y); 

%% noisy image
mkdir('results\noisy')
for i=1:B
    imwrite(Y(:,:,i), ['results\noisy\',num2str(i),'.jpg'])
end

%% BALMF
Rank = 3;
start =  tic;
[DY,U,V,model1] = BALMF(Y, Rank);
exe_time = toc(start);

mkdir('results\BALMF')
for i=1:B
    imwrite(DY(:,:,i), ['results\BALMF\',num2str(i),'.jpg'])
end
Y_BALMF = uint8(255*DY);
save('results\BALMF.mat','Y_BALMF')

%% Pseudo-color images 
subplot(1,2,1), rsshow(Y_BALMF(:,:,[152,96,43])), title('BALMF')
subplot(1,2,2), rsshow(Y(:,:,[152,96,43])), title('Noisy Image')
