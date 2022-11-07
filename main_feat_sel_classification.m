clear all; close all; clc
net_name = 'inceptionv3';
dir_db = ['D:\digital pathology\melanoma project\extracted_data\',net_name,'\feat_crops\training_test_split\'];

NumRounds = 5; % number of rounds
NumFolds = 5; % number of folds
th_q = 0.75; % threshold for vote score thresholding


for nr = 1: NumRounds
    for nf = 1: Numfolds
        %% load train set
        A = load([dir_db,'training_fold',num2str(nf),'_round',num2str(nr),'.csv']);
        A(1,:)=[];
        feat_train0 = normalize(A(:,3:end));
        GT_train = A(:,2);
        id_train = A(:,1);
        %% load test set
        B = load([dir_db,'test_fold',num2str(nf),'_round',num2str(nr),'.csv']);
        B(1,:) = [];
        load([dir_db,'test_fold',num2str(nf),'_round',num2str(nr),'.mat'],'idTestunique')
        idpz = idTestunique;
        idpzunique = unique(idTestunique);
        feat_test0 = normalize(B(:,3:end));
        GT_test = B(:,2);
        id_test = B(:,1);
        %% feature selection: features whose AUC value was greater than 0.6 were retained
        id_auc = [];
        for k = 1:size(feat_train0,2)
            [~,~,~,auc] = perfcurve(label_tot,feat_train0(:,k),1);  
            if auc>=0.6
                id_auc = [id_auc k];   
            end
        end
        feat_train = feat_train(:,id_auc);
        feat_test = feat_test(:,id_auc);
        clear id_auc
        %% principal component analysis
        [coeff,scoreTrain,latent,tsquared,explained,mu] = pca(feat_train);
        sum_explained = 0;
        idpca = 0;
        while sum_explained < 80
            idpca = idpca + 1;
            sum_explained = sum_explained + explained(idpca);
        end
     
        feat_train2 = scoreTrain(:,1:idpca);
        feat_test2 = (feat_test-mu)*coeff(:,1:idpca);
       
        %% svm classifier
       
        classifier = fitcsvm(feat_train2,GT_train,'Standardize', true, 'KernelFunction','RBF','Kernelscale','auto'); %
        SVMModel =  fitPosterior(classifier,feat_train2,GT_train);
        [yhat0,scores] = predict(SVMModel, feat_test2);
        
        %% vote score thresholding
         for j = 1:length(idpzunique)
            dbnew(j).idnew = find(idTestunique==idpzunique(j));
            dbnew(j).labelnew = unique(GT_test(dbnew(j).idnew));
            dbnew(j).scoresnew = scores(dbnew(j).idnew,2);
            scores_new(j) = quantile(dbnew(j).scoresnew,th_q); 
         end
        %% scores on tiles
         scores_totxslide(idpzunique,nr) = scores_new;
        %% scores on crops
         scores_totxcrop(id_test,nr) = scores(:,2);
        clearvars -except scores_totxslide scores_totxcrop nr nf NumRounds NumFolds th_q
    end
  
end


