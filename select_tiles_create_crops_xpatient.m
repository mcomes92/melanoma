clear all; close all; clc
%% load tiles images, annotations saved by QuPath with info about each tile of WSI (e.g., number of cell deection)
nome_dir_tiles = 'D:\digital pathology\melanoma project\Ti2les\';
nome_dir_annotations = 'D:\digital pathology\melanoma project\annotation results\';
names_pazienti = dir(nome_dir_tiles);
names_annotations = dir(nome_dir_annotations);
addpath(genpath([nome_dir_tiles]))
addpath(genpath([nome_dir_annotations]))
%% load code for image normalization with Macenko's method (the code is available at the following link: https://github.com/mitkovetta/staining-normalization)
addpath(genpath('D:\digital pathology\melanoma project\staining-normalization-master\'))
%% set thresholds
th = 0.9; % thresh for cell detection
normal = true;
th_crop = 0.25;% thresh for crop backgroud
%% create folders to save selected files and crops
nome_dir_tiles_new = ['D:\digital pathology\melanoma project\Ti2les_selected_',num2str(th),'_normal\'];
nome_dir_crops_new = ['D:\digital pathology\melanoma project\CropTi2les2_selected_',num2str(th),'_normal\'];

mkdir(nome_dir_crops_new)
mkdir(nome_dir_tiles_new)
numtotTiles = [];
pztot = [];
labtot = [];
%% define selected tiles and create crops for each patient
for np = 1:length(names_pazienti)-2
    disp(names_pazienti(2+np).name)
    nome_dir = [nome_dir_tiles,names_pazienti(2+np).name,'\'];
    addpath(genpath(nome_dir))
    images = dir(fullfile(nome_dir,'*tif'));

    %% find annotation related to the WSI of the np-th patient
    A = importdata([nome_dir_annotations,names_pazienti(2+np).name,'.svs Annotations.txt']);
    % A.data contains the following columns:
    % Centroid X µm	Centroid Y µm	Num Detections	Area µm^2	Perimeter µm
    %% select squared tiles of the same dimensions
    id_im = [];
    for im = 1:length(images)
        if size(imread(images(im).name),1)==size(imread(images(im).name),2) %prendo solo tiles quadrate
            id_im = [id_im im];
        end
    end
    if isempty(id_im)
        id_im = 1:length(images);
    end
    num_tiles1 = [];
    for jj = 1:length(id_im)
        if not(isempty(str2num(images(id_im(jj)).name(end-6)))) && not(isempty(str2num(images(id_im(jj)).name(end-5))))
            num_tile =str2num(images(id_im(jj)).name(end-6:end-4));
        elseif not(isempty(str2num(images(id_im(jj)).name(end-5)))) && isempty(str2num(images(id_im(jj)).name(end-6)))
            num_tile =str2num(images(id_im(jj)).name(end-5:end-4));
        elseif  not(isempty(str2num(images(id_im(jj)).name(end-4)))) && isempty(str2num(images(id_im(jj)).name(end-5)))
            num_tile =str2num(images(id_im(jj)).name(end-4));
        end
        num_tiles1 = [num_tiles1 num_tile];
    end
    Adata2 = A.data(num_tiles1,:);
    % create associations between tiles and annotations
    for k = 1:size(A.data,1)
        if not(isempty(str2num(A.textdata{k+1,2}(end-3:end))))
            tiles(k) =str2num(A.textdata{k+1,2}(end-3:end));
        elseif not(isempty(str2num(A.textdata{k+1,2}(end-2:end))))
            tiles(k) =str2num(A.textdata{k+1,2}(end-2:end));
        elseif  not(isempty(str2num(A.textdata{k+1,2}(end-1:end))))
            tiles(k) =str2num(A.textdata{k+1,2}(end-1:end));
        end
    end
    tiles2 = tiles(num_tiles1);
    %% retain only tiles whose cell density exceed the 90th percentile of the distribution of the number of cell detenctions
    detect_quant = quantile(Adata2(:,3),th);
    id_im2 = find(Adata2(:,3)>=detect_quant);
    Adata3 = Adata2(id_im2,:);
    tiles3 = tiles2(id_im2);
    numtotTiles = [numtotTiles length(tiles3)];
    %% save the retained tiles
    nome_dir_new = [nome_dir_tiles_new,names_pazienti(2+np).name,'\'];
    mkdir(nome_dir_new)
    nome_dir_new2 = [nome_dir_crops_new,names_pazienti(2+np).name,'\'];
    mkdir(nome_dir_new2)
    for jj = 1:length(images)
        if not(isempty(str2num(images(jj).name(end-6)))) && not(isempty(str2num(images(jj).name(end-5))))
            num_tile =str2num(images(jj).name(end-6:end-4));
        elseif not(isempty(str2num(images(jj).name(end-5)))) && isempty(str2num(images(jj).name(end-6)))
            num_tile =str2num(images(jj).name(end-5:end-4));
        elseif  not(isempty(str2num(images(jj).name(end-4)))) && isempty(str2num(images(jj).name(end-5)))
            num_tile =str2num(images(jj).name(end-4));
        end
        %% create random crops per tiles (of dimensions equal to 1/4 of the tile dimensions)
        if not(isempty(find(tiles3==num_tile))) % le immagini sono già rgb

            Ncrop = 50;
            inputSize = size(imread(images(jj).name));
            targetSize = [round(size(imread(images(jj).name),1)/4) round(size(imread(images(jj).name),2)/4)];
            for cc = 1:Ncrop
                rect = randomWindow2d(inputSize,targetSize);
                rectXYWH = [rect.XLimits(1) rect.YLimits(1) ...
                    diff(rect.XLimits)+1 diff(rect.YLimits)+1];
                J = imcrop(imread(images(jj).name),rectXYWH);
                % select tiles with less than 25% background pixels
                % (Luma >170)
                if (length(find(rgb2gray(J)>170))/(targetSize(1)*targetSize(2)))<th_crop
                    [images(jj).name cc]
                    imwrite(J,[nome_dir_new2,'tile',num2str(num_tile),'_crop',num2str(cc),'.png'])
                end
                clear rect J
            end

            imwrite(normalizeStaining(imread(images(jj).name)),[nome_dir_new,'tile',num2str(num_tile),'.png'])
        end
    end

    clearvars -except normal labtot numtotTiles pztot np th  th_crop size_pred nome_dir_tiles nome_dir_annotations nome_dir_tiles_new nome_dir_crops_new names_pazienti names_annotations
end





