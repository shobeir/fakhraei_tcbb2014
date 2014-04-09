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

% This file contains the new interactions experiment
% We extracted the new interactions list from Drugbank on March 2014.

clear all
addpath(genpath('./'));
addpath(genpath('../'));

% Making similarities global
global DrugSim_ATCHierDrugsCommonSimilarityMat DrugSim_chemicalDrugsCommonSimilarityMat DrugSim_ligandJaccardDrugsCommonSimilarityMat DrugSim_newCMapJaccardDrugsCommonSimilarityMat DrugSim_pSideEffectDrugsCommonSimilarityMat;
global TargetSim_distTargetsCommonSimilarityMat TargetSim_GOTargetsCommonSimilarityMat TargetSim_seqTargetsCommonSimilarityMat;

fprintf('\nLoading the data...\n');
load Perlman_Data.mat;
load New_Interactions.mat;

% Setting the parameters
nDrugs = 315; % No of drugs
nTargets = 250; % No of targets
num_folds = 10; % No of folds for cross validation
indices =  NewInteractionFolds(:,3); % Setting the folds to the same one used with PSL

% Creating the dataset
fprintf('\nStarting to creat the dataset... \n');
Perlman_Dataset = generate_dataset( nDrugs, nTargets, Interactions_Matrix, Interactions_Matrix);

% Adding the new interactions lables
for i=1:length(new_Interactions),
    Perlman_Dataset(find((Perlman_Dataset(:,1)==new_Interactions(i,1))&(Perlman_Dataset(:,2)==new_Interactions(i,2))),18)=1;
end

% Splitting feature values and labels
Labels = Perlman_Dataset(:,18);
Data = Perlman_Dataset(:,3:17);

fprintf('\nTesting the new interactions with 10 differnt negative samples... \n');
% Testing the new interactions with ten different negative samples
for currentFold = 1:num_folds,
    
    % Splitting the data
    test = ((indices == currentFold) | (indices == -1));
    train = ~test;
    Fold_Labels(currentFold,1:length(Labels(test))) = Labels(test);
    
    % Training
    [Model,~,~] = train_model(Data, train, Labels);
    
    % Testing
    Fold_Predictions(currentFold,1:length(Labels(test))) = test_model(Model, Data, test);
    
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
