%% feature extraction by transfer learning
clear all; close all; clc;
net_name = 'inceptionv3'; % choose net
net= inceptionv3;
layer='avg_pool'; % choose layer
th = 0.9;
mkdir(['D:\digital pathology\melanoma project\extracted_data\',net_name,'\feat_crops\'])
%% load the folder containing crops per patient
name_dir_ROI_files = ['D:\digital pathology\melanoma project\CropTi2les2_selected_',num2str(th),'_normal\'];
ROI_files = dir(name_dir_ROI_files );
addpath(name_dir_ROI_files )
all_ROIs = ROI_files([ROI_files(:).isdir]);
num_ROIs = numel(all_ROIs(3:end));
N = num_ROIs; % number of patients

features = [];
id = [];
id_count = [];
id_unique = [];
id_count_unique = [];
for k = 1 : N
    disp(k)
    path_N=[name_dir_ROI_files,ROI_files(2+k).name,'\'];
    addpath(genpath(path_N))
    img_dir = dir([path_N,'*png']); % load all the crops per patient
    %% resize images as required to be input of the network (function my_read_inceptionv3.m)
    imds = imageDatastore(path_N,'ReadFcn',@(x) my_read_inceptionv3(x));
    I=readall(imds);
    figure(1); imshow(I);
    N_roi=length(I);
    %% extract features for each crop and concatenate features for each patient
    features = [features;activations(net,imds,layer,'OutputAs','rows')];
    %% create arrays with the name of crops and the index of patient of belonging
    id = [id; repmat(num2str(ROI_files(2+k).name(5:9)),N_roi,1) ];
    id_count = [id_count; repmat(k,N_roi,1)];
    id_unique = [id_unique; num2str(ROI_files(2+k).name(5:9))];
    id_count_unique = [id_count_unique; k];
    clear imds I
end

save(['D:\digital pathology\melanoma project\extracted_data\',net_name,'\feat_crops\feat_crops',num2str(th),'_',net_name,'_avgpool_norm.mat'],'features','id','id_count','id_unique','id_count_unique')
writematrix([1:size(features,2);features],['D:\digital pathology\melanoma project\extracted_data\',net_name,'\feat_crops\feat_crops',num2str(th),'_',net_name,'_avgpool_norm.csv'])
