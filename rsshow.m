function A = rsshow(I, scale)
% Remote sensing image enhancement for visualization.
%
% Usage:
% display an image: rsshow(I, 0.05)
% write an image: A = rsshow(I, 0.05); imwrite(A, 'output.jpg')
%
% If the code is used in your scientific research, please cite the paper.
% [1] Shuang Xu, Xiangyong Cao, Jiangjun Peng, Qiao Ke, Cong Ma and Deyu
% Meng. Hyperspectral Image Denoising by Asymmetric Noise Modeling. IEEE
% TGRS, 2023.

if nargin==1
    scale=0.005;
end
I = double(I);

if ismatrix(I)

    q = quantile(I(:),[scale, 1-scale]);
    [low, high] = deal(q(1),q(2));
    I(I>high) = high;
    I(I<low) = low;
    I = (I-low)/(high-low);
else
  
    for i=1:size(I,3)
        temp = I(:,:,i);
        q = quantile(temp(:),[scale, 1-scale]);
        [low, high] = deal(q(1),q(2));
        temp(temp>high) = high;
        temp(temp<low) = low;
        temp = (temp-low)/(high-low);
        I(:,:,i) = temp;
    end
end

if nargout==1
    A=I;
else
    imshow(I)
end
end