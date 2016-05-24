mgForward <-
function(genD, vectorsMEM, perm=100, alpha=0.05) {
    
    X <- as.matrix(vectorsMEM)
    genD <- as.matrix(genD)
    Number_Predictors <- ncol(X)
    n <- nrow(X)
    X <- apply(X, 2, scale)
    result <- mgRDA(genD, X, full=FALSE)
    F_Observed <- result$F
    
    # global test with all predictors
    Prob_Global <- 1/perm
    for (i in 1:(perm-1)) {
        X_permuted <- X[sample(n,replace=FALSE), ]
        result <- mgRDA(genD, X_permuted, full=FALSE)
        F_Random <- result$F
        if (F_Random >= F_Observed) {
            Prob_Global <- Prob_Global +1/perm
        }
    }
    
    Variables_In_Model <- NA
    Rsq_Final_Model <- NA
    
    if (Prob_Global < alpha) {    
        # calculate contribution of each predictor separately as they're orthgonal
        F_Ind_X <- mat.or.vec(Number_Predictors,1)
        for (i in 1:Number_Predictors) {
            F_Ind_X[i] <- mgRDA(genD, X[, i], full=FALSE)$F
        }
        # start selection
        Variables_In_Model <- as.matrix(which.max(F_Ind_X))
        F_Ind_X[Variables_In_Model[1]] <- NA
        found <- FALSE
        while (found == FALSE) { 
            # contrast current model with the largest variable contribution not entered in the model
            candidate_model <- cbind(X[, Variables_In_Model],X[, which.max(F_Ind_X)])
            result_candidate <- mgRDA(genD, candidate_model, full=FALSE)
            F_candidate_Obs <- result_candidate[2]
            # test whether the entered variable improves fit
            Prob_F <- 1/perm
            for (i in 1:(perm-1)) {
                candidate_predictor <- which.max(F_Ind_X)
                candidate_model <- cbind(X[,Variables_In_Model],X[sample(n,replace=FALSE),candidate_predictor])
                result_candidate <- mgRDA(genD, candidate_model, full=FALSE);
                F_candidate_Rnd <- result_candidate$F
                if (F_candidate_Rnd >= F_candidate_Obs) {
                    Prob_F <- Prob_F + 1/perm}
            }
            if (Prob_F > alpha) {
                found <- TRUE
            } else {
                F_Ind_X[candidate_predictor] <- NA
                Variables_In_Model <- append(Variables_In_Model,candidate_predictor)
            }
          
        }
        result <- mgRDA(genD, X[,Variables_In_Model], full=FALSE)
        Rsq_Final_Model<-result$RsqAdj
    }
    return(list(GlobalP=Prob_Global,
                selectedRsqAdj=Rsq_Final_Model,
                selectedMEM=Variables_In_Model))
}
