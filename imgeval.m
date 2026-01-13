% Avaliação de imagens segmentadas do Microsoft GrabCut dataset
% Uso: error = imgeval(imgres, gt, imgslab)
%   imgres  - imagem resultante da segmentação (máscara 0/255 ou similar)
%   gt      - ground truth
%   imgslab - imagem de scribbles (128 = região avaliada)
function error = imgeval(imgres, gt, imgslab)

    % total de pixels a serem avaliados:
    % gt ~= 128  -> ignora pixels marcados como indecisos no GT
    % imgslab==128 -> só onde o usuário realmente quis segmentar
    totunlpix = sum(sum(imgslab == 128 & gt ~= 128));

    % pixels em que a diferença entre resultado e GT é > 1 (tolerância mínima),
    % restrito à mesma máscara de avaliação
    toterrpix = sum(sum( ...
        abs(double(imgres) - double(gt)) > 1 & ...
        imgslab == 128 & gt ~= 128));

    error = toterrpix / totunlpix;
end
