function [ Perlman_Dataset ] = generate_dataset( nDrugs, nTargets, Interactions_Matrix, Interactions_Matrix_Full)

global DrugSim_ATCHierDrugsCommonSimilarityMat DrugSim_chemicalDrugsCommonSimilarityMat DrugSim_ligandJaccardDrugsCommonSimilarityMat DrugSim_newCMapJaccardDrugsCommonSimilarityMat DrugSim_pSideEffectDrugsCommonSimilarityMat;
global TargetSim_distTargetsCommonSimilarityMat TargetSim_GOTargetsCommonSimilarityMat TargetSim_seqTargetsCommonSimilarityMat;

Perlman_Dataset = zeros(nDrugs*nTargets,18);

fprintf('Will spin for %i iterations:\n',nDrugs*nTargets);

for i=1:nDrugs,
    for j=1:nTargets,
        
        Perlman_Dataset((i-1)*nTargets+j,1)=i;
        Perlman_Dataset((i-1)*nTargets+j,2)=j;
        
        Perlman_Dataset((i-1)*nTargets+j,3)=calcualte_feature(i,j,DrugSim_ATCHierDrugsCommonSimilarityMat,TargetSim_GOTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,4)=calcualte_feature(i,j,DrugSim_ATCHierDrugsCommonSimilarityMat,TargetSim_distTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,5)=calcualte_feature(i,j,DrugSim_ATCHierDrugsCommonSimilarityMat,TargetSim_seqTargetsCommonSimilarityMat,Interactions_Matrix);
        
        Perlman_Dataset((i-1)*nTargets+j,6)=calcualte_feature(i,j,DrugSim_chemicalDrugsCommonSimilarityMat,TargetSim_GOTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,7)=calcualte_feature(i,j,DrugSim_chemicalDrugsCommonSimilarityMat,TargetSim_distTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,8)=calcualte_feature(i,j,DrugSim_chemicalDrugsCommonSimilarityMat,TargetSim_seqTargetsCommonSimilarityMat,Interactions_Matrix);
        
        Perlman_Dataset((i-1)*nTargets+j,9)=calcualte_feature(i,j,DrugSim_ligandJaccardDrugsCommonSimilarityMat,TargetSim_GOTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,10)=calcualte_feature(i,j,DrugSim_ligandJaccardDrugsCommonSimilarityMat,TargetSim_distTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,11)=calcualte_feature(i,j,DrugSim_ligandJaccardDrugsCommonSimilarityMat,TargetSim_seqTargetsCommonSimilarityMat,Interactions_Matrix);
        
        Perlman_Dataset((i-1)*nTargets+j,12)=calcualte_feature(i,j,DrugSim_newCMapJaccardDrugsCommonSimilarityMat,TargetSim_GOTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,13)=calcualte_feature(i,j,DrugSim_newCMapJaccardDrugsCommonSimilarityMat,TargetSim_distTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,14)=calcualte_feature(i,j,DrugSim_newCMapJaccardDrugsCommonSimilarityMat,TargetSim_seqTargetsCommonSimilarityMat,Interactions_Matrix);
        
        Perlman_Dataset((i-1)*nTargets+j,15)=calcualte_feature(i,j,DrugSim_pSideEffectDrugsCommonSimilarityMat,TargetSim_GOTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,16)=calcualte_feature(i,j,DrugSim_pSideEffectDrugsCommonSimilarityMat,TargetSim_distTargetsCommonSimilarityMat,Interactions_Matrix);
        Perlman_Dataset((i-1)*nTargets+j,17)=calcualte_feature(i,j,DrugSim_pSideEffectDrugsCommonSimilarityMat,TargetSim_seqTargetsCommonSimilarityMat,Interactions_Matrix);
        
        Perlman_Dataset((i-1)*nTargets+j,18)=Interactions_Matrix_Full(i,j);
        
        % Printing the status
        if rem(((i-1)*nTargets)+j,1000)==0
            fprintf('.');
        end
    end
end

end

