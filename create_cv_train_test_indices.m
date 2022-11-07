%% create indices of test sets in 5-fold cross-validation for 5 rounds
clear all; close all; clc
numRounds = 5;
numFolds = 5;
net_name = 'inceptionv3';
th = 0.9;
load(['D:\digital pathology\melanoma project\extracted_data\',net_name,'\feat_crops\feat_crops',num2str(th),'_',net_name,'_avgpool_norm.mat'],'features','id','id_count','id_unique','id_count_unique')
load(['D:\digital pathology\melanoma project\patient_label\labels_patients.mat'],'gt_unique','gt')
id_count_unique1 =id_count_unique;
gt_unique1 = gt_unique; % patients' labels
id_unique1 = str2num(convertCharsToStrings(id_unique));
features1 = features;

id_tot = [id_unique1];
gt_tot = [gt_unique1];
id_count_tot = [id_count_unique1];
features_tot = [features1];
matToSave = zeros(length(id_tot),numRounds);
for round=1:numRounds
cv = cvpartition(gt_tot,'KFold',numFolds,'Stratify',true);
for i = 1:numFolds
    matToSave(cv.test(i),round)= i;
end
  clear cv
end

%% create training/test splitting 
name_dir = ['D:\digital pathology\melanoma project\extracted_data\',net_name,'\feat_crops\training_test_split\'];
mkdir(name_dir)

id_tot = [id];
gt_tot = [gt];
id_count_tot = [id_count];
features_tot = [features];

for round=1:numRounds
    for i = 1:numFolds
        [i round]
        id0 = find(matToSave(:,round)==i);
        idTest = [];
        for l = 1:length(id0)
            idTest = [idTest;find(id_count_tot==id0(l))];
            
        end
            gtTest = gt_tot(idTest);
            featuresTest = features_tot(idTest,:);
            matrixTest = [idTest,gtTest,featuresTest];
            writematrix([[1:size(matrixTest,2);matrixTest]],[name_dir,'test_fold',num2str(i),'_round',num2str(round),'.csv'])
            idTrain = setdiff([1:length(id_count_tot)],idTest)';
            gtTrain = gt_tot(idTrain);
            featuresTrain = features_tot(idTrain,:);
            matrixTrain = [idTrain,gtTrain,featuresTrain];
            writematrix([[1:size(matrixTrain,2);matrixTrain]],[name_dir,'training_fold',num2str(i),'_round',num2str(round),'.csv'])
           
    end
end

%% create unique indices for test (it will be used for the vote score thresholding procedure)
for round=1:numRounds
    for i = 1:numFolds
        [i round]
        id0 = find(matToSave(:,round)==i);
        idTest = [];
        idTestunique = [];
        for l = 1:length(id0)
            idTest = [idTest;find(id_count_tot==id0(l))];
            idTestunique = [idTestunique;repmat(id0(l),length(find(find(id_count_tot==id0(l)))),1)];
            
        end
         save([name_dir,'test_fold',num2str(i),'_round',num2str(round),'.mat'],'idTestunique','idTest','-append')
    end
end