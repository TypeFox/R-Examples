opt.splitval <-
    function(optFUN="opt1D",testset="equal",scaling=TRUE,...){
        extra.vars <- list(...)
        if(testset[1]=="equal"){
            testset <- rep(TRUE,nrow(extra.vars$penalized))
            testset[sample(1:nrow(extra.vars$penalized),round(nrow(extra.vars$penalized)/2))] <- FALSE
            testset <- which(testset)
        }
        trainingset <- (1:nrow(extra.vars$penalized))[-testset]
        ##Define training and test data
        penalized.training <- extra.vars$penalized[trainingset,]
        if("unpenalized" %in% names(extra.vars)){
            unpenalized.training <- extra.vars$unpenalized[trainingset,]
        }
        pheno.training <- extra.vars$response[trainingset]
        penalized.test <- extra.vars$penalized[testset,]
        if("unpenalized" %in% names(extra.vars)){
            unpenalized.test <- extra.vars$unpenalized[testset,]
        }
        pheno.test <- extra.vars$response[testset]
        if(scaling){
            penalized.training <- scale(penalized.training)
            ##use scaling factors determined in training set to scale the test set:
            penalized.test <- sweep(penalized.test,2,attr(penalized.training,"scaled:center"))
            penalized.test <- sweep(penalized.test,2,attr(penalized.training,"scaled:scale"),FUN="/")
            penalized.training <- data.frame(penalized.training)
        }
        ##Assign training samples to extra.vars for regression:
        extra.vars$response <- pheno.training
        extra.vars$penalized <- penalized.training
        if("unpenalized" %in% names(extra.vars)){
            extra.vars$unpenalized <- unpenalized.training
        }
        ##regression call:
        output <- do.call(get(optFUN), args = extra.vars)
        ##coefficients:
        cc <- output[which.max(output[, "cvl"]),]
        cc <- cc[-na.omit(match(c("L1","L2","cvl","convergence","fncalls"),names(cc)))]
        ##does the model have an intercept?
        if(grepl("(Intercept)",names(cc)[1],fixed=TRUE)){  
            intercept <- cc[1]
            cc <- cc[-1]
        }
        ##predictions in test set:
        if("unpenalized" %in% names(extra.vars)){
            dat.test <- cbind(unpenalized.test,penalized.test)
            preds.test <- as.matrix(dat.test) %*% cc
        }else{
            preds.test <- as.matrix(penalized.test) %*% cc
        }
        output <- preds.test[,1]
        if(exists("intercept")) output <- output + intercept
        names(output) <- rownames(preds.test)
        return(output)
    }


