% Implemented by Shobeir Fakhraei for experimental evaluation in the following paper:
%
%     @article{fakharei2014network,
%     title={Network-Based Drug-Target Interaction Prediction with Probabilistic Soft Logic},
%     author={Fakhraei, Shobeir and Huang, Bert and Raschid, Louiqa and Getoor, Lise},
%     journal={IEEE/ACM Transactions on Computational Biology and Bioinformatics},
%     year={2014},
%     }
%
% The implementation is based on the following paper:
%
%     @article{perlman2011combining,
%     title={Combining drug and gene similarity measures for drug-target elucidation},
%     author={Perlman, Liat and Gottlieb, Assaf and Atias, Nir and Ruppin, Eytan and Sharan, Roded},
%     journal={Journal of computational biology},
%     year={2011},
%     }

% This file contains the cross validation experiment

clear all
addpath(genpath('./'));
addpath(genpath('../'));

% Making similarities global
global DrugSim_ATCHierDrugsCommonSimilarityMat DrugSim_chemicalDrugsCommonSimilarityMat DrugSim_ligandJaccardDrugsCommonSimilarityMat DrugSim_newCMapJaccardDrugsCommonSimilarityMat DrugSim_pSideEffectDrugsCommonSimilarityMat;
global TargetSim_distTargetsCommonSimilarityMat TargetSim_GOTargetsCommonSimilarityMat TargetSim_seqTargetsCommonSimilarityMat;

fprintf('\nLoading the data...\n');
load Perlman_Data.mat;

% Setting the parameters
nDrugs = 315; % No of drugs
nTargets = 250; % No of targets
num_folds = 10; % No of folds for cross validation
indices =  PSLFolds(:,3); % Setting the folds to the one used with PSL

fprintf('Starting Cross-validation!\n');
for currentFold = 1:num_folds,
    
    % removing the held out interactions
    Interactions_Matrix_Fold =  Interactions_Matrix;
    
    for i=1:nDrugs,
        for j=1:nTargets,
            if indices((i-1)*nTargets+j)==currentFold
                Interactions_Matrix_Fold(i,j)=0;
            end
        end
    end
    
    % Creating the dataset
    fprintf('\nStarting to creat the dataset for fold %i:\n',currentFold);
    Perlman_Dataset_Fold = generate_dataset( nDrugs, nTargets, Interactions_Matrix_Fold, Interactions_Matrix);
    
    % Splitting the data
    test = (indices == currentFold);
    train = ~test;
    Labels = Perlman_Dataset_Fold(:,18);
    FoldData = Perlman_Dataset_Fold(:,3:17);
    Fold_Labels(currentFold,1:length(Labels(test))) = Labels(test);
    
    % Training
    [Model,~,~] = train_model(FoldData, train, Labels);
    
    % Testing
    Fold_Predictions(currentFold,1:length(Labels(test))) = test_model(Model, FoldData, test);
    
end

% Report the results

% AUROC
for i = 1:num_folds
    [xROC,yROC,tROC,AUROC(i)]=perfcurve(Fold_Labels(i,:),Fold_Predictions(i,:),1);
end

% AUPR
for i = 1:num_folds
    [xPR,yPR,tPR,AUPR(i)] = perfcurve(Fold_Labels(i,:),Fold_Predictions(i,:),1, 'xCrit', 'reca', 'yCrit', 'prec');
end

% Displaying the results
fprintf('\nFinal Results:\n');
fprintf('AUROC: %i +/- %i \n',mean(AUROC),std(AUROC));
fprintf('AUPR: %i +/- %i \n',mean(AUPR),std(AUPR));



