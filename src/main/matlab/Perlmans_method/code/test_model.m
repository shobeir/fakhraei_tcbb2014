function [ Predictions ] = test_model( b, X, test )

    Predictions = glmval(b,X(test,:),'logit')';

end

