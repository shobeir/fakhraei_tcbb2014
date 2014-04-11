function y = calculate_feature(Drug,Target,DrugSimMat,TargetSimMat,InteractionsMat)

    r=0.4; % Set based on the value stated in the paper
    InteractionsMat(Drug,Target)=0;
    y=max(max(((DrugSimMat(Drug,:).^(r))'*(TargetSimMat(Target,:).^(1-r))).*InteractionsMat));

end
