% Avaliação de imagens segmentadas do Microsoft GrabCut dataset
% Uso: error = imgeval(imgres, gt, imgslab)
% imgres - imagem resultante da segmentação
% gt - ground truth
% imgslab - imagem 
function error = imgeval(imgres, gt, imgslab)
    totunlpix = sum(sum(imgslab==128 & gt~=128));
    toterrpix = sum(sum(abs(double(imgres)-double(gt))>1 & imgslab==128 & gt~=128));
    error = toterrpix/totunlpix;
end