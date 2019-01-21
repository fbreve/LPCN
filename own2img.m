%Uso: own2img(owner,img,1) - exibe imagem
%     imgres = own2img(owner,img,0); - não exibe imagem
function imgres = own2img(owner,img,exhibitimg,imgslab)
if (nargin < 3) || isempty(exhibitimg)
    exhibitimg = 1;
end
dim = size(img);
if (nargin <4) || isempty(imgslab)
    owner = owner-min(owner);
    owner = owner./max(owner)*255;    
else
    colors = unique(imgslab);
    colors = [0; colors(colors~=64 & colors~=128)];
    owner = colors(owner);
end
imgres = uint8(reshape(owner,dim(1),dim(2)));
if exhibitimg==1
    imshow(imgres,gray(256))
end
end