function [b,dev,stats] = train_model( X, train, Labels )

    % Fitting the Logisitc Regression with Cross Validation
    [b,dev,stats] = glmfit(X(train,:),Labels(train),'binomial','logit'); % Logistic regression

end

