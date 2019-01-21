% This is an usage example of the Interactive Image Segmentation using 
% Label Propagation through Complex Networks method.
%
% BREVE, Fabricio Aparecido. Interactive Image Segmentation using Label
% Propagation through Complex Networks. Expert System With Applications, 
% v. 123, p.18 � 33, 2019.
%
% by Fabricio Breve - 21/01/2019
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
subplot(2,2,1), imshow(img),
subplot(2,2,2), imshow(imgslab),
subplot(2,2,3), imshow(imgres),
subplot(2,2,4), imshow(imgseg);