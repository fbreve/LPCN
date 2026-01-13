% Uso: own2img(owner,img,1)          - exibe imagem
%      imgres = own2img(owner,img,0) - não exibe imagem
function imgres = own2img(owner, img, exhibitimg, imgslab)

if (nargin < 3) || isempty(exhibitimg)
    exhibitimg = 1;
end

dim = size(img);

if (nargin < 4) || isempty(imgslab)
    % Normaliza rótulos de owner para [0, 255]
    owner = owner - min(owner);
    maxv  = max(owner);
    if maxv > 0
        owner = owner ./ maxv * 255;
    else
        owner = zeros(size(owner));
    end
else
    % Usa as mesmas cores da imagem de scribbles
    colors = unique(imgslab);
    % remove 64 (fundo) e 128 (não rotulado)
    colors = [0; colors(colors ~= 64 & colors ~= 128)];
    % owner é índice em 'colors'
    owner = colors(owner);
end

imgres = uint8(reshape(owner, dim(1), dim(2)));

if exhibitimg == 1
    imshow(imgres, gray(256));
end

end
