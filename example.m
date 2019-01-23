% This is an usage example of the Interactive Image Segmentation using 
% Label Propagation through Complex Networks method.
%
% BREVE, Fabricio Aparecido. Interactive Image Segmentation using Label
% Propagation through Complex Networks. Expert System With Applications, 
% v. 123, p.18 – 33, 2019.
%
% by Fabricio Breve - 23/01/2019
%
% Loading example image
img = imread('ralph.jpg');
% Loading scribbles image
imgslab = imread('ralph-scribble.png');
% Loading ground-truth image
gt = imread('ralph-gt.png');
% Calling method to segment the image with default parameters
disp('Running image segmentation...');
tStart = tic;
[owner, pot] = lpcn(img, imgslab);
tElapsed = toc(tStart);
% Converting output to image
imgres = own2img(owner,img,0);
% Calculating error rate
y = imgeval(imgres, gt, imgslab);
% Displaying numerical results
fprintf('Error Rate: %0.4f - Execution Time: %0.4fs\n',y,tElapsed);
% Applying mask to the original image
imgres3 = repmat(imgres,1,1,3)./255;
imgseg = img .* imgres3;
imgseg(imgres3==0)=255;
% Showing segmentation results
figure('units','normalized','outerposition',[0 0 1 1])
imgslabmix = imread('ralph-scribblemix.jpg');
subplot(2,3,1), imshow(img), title('Original image'),
subplot(2,3,2), imshow(imgslab), title('Scribbles (user-input)'),
subplot(2,3,3), imshow(imgslabmix), title('Scribbles overlaying original image'),
subplot(2,3,4), imshow(gt), title('Ground-truth'),
subplot(2,3,5), imshow(imgres), title('Segmentation results (mask)'),
subplot(2,3,6), imshow(imgseg); title('Segmentation results (extracted foreground)'),