train.doctor <- function(F_,L_,training.samples, validating.samples, considered.features, maximum.features.num, 
                         balance=TRUE,return_linear.models=TRUE, report.fitting.failure=FALSE){
    ##****	liner modeling: MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    ## In this function, a class is analized by cross-validating a portion of samples (gamma) which are passed as X_.
    ## In consequence, num.of.features.score, error.log, mislabeling.record and feature.scores are computed.
    ## Given the feature values, linear and logistic coefficients
    ## returns the probablity of being in the target class and other results.
    ##____________________________________________________________________________________________________________________________________
    ## INPUT: 	F_ is the full features matrix.
    ##			L_ should be the vector of labels ordered according to rows of F. The lables are 0 or +1 for negative and positive instances.
    ##			training.samples are the names of the rows which corespond training samples.  
    ## 			validating.samples are the names of the rows which corespond validating samples.  
    ##____________________________________________________________________________________________________________________________________
    X_ <- F_[training.samples,considered.features]			
    ## X_ is the feature matrix for training samples. Some features might have been omited.
    contribution.to.feature.scores <- c() ## Training this doctor should later on, add to total feature scores according to this vector.
    contribution.to.feature.scores[considered.features] <- 0
    
    contribution.to.features.frequency <- c()	## Training this doctor should later on, add to total frequency of features 
    ## (number of apprearance) according to this vector.
    contribution.to.features.frequency[considered.features] <- 0
    contribution.to.mislabeling.record <- c()	## Training this doctor should later on, add to mislabeling.record.
    contribution.to.mislabeling.record[union(training.samples, validating.samples)] <- 0
    ##
    f.measure.log <-c()	## A vector keeping track of f.meatures for all tried number of featurs.
    best.f.measure <- 0					
    best.feature.num <- 0
    best.logistic.scores <-c()
    best.logistic.scores[validating.samples] <- 0
    best.logistic.cof <- 0

    Y_ <- L_[rownames(X_)] ## We wish we could fit a w such that: X_W=Y_
    ## Input check (partial):
    if(length(unique(Y_))==1)
        stop("Input labels are all identical, i.e. Y=const")

    ## balancing:
    if (balance){	
        XY <- compute.balanced(F_=X_,L_=Y_)
        X_ <- XY$F_
        Y_ <- XY$L_
    }##End if (balance).
    ##message("compute.balanced()"); a() ## CHECK POINT for testing
    
    ## LARC
    lars.result <- try(lars(x=X_, y=Y_, trace=FALSE,use.Gram=FALSE),silent=!report.fitting.failure)
    if(inherits(lars.result, "try-error"))
        stop("Lars failed.")
    ## "When the number m of variables is very large, i.e. larger than number of samples, 
    ## then you may not want LARS to precompute the Gram matrix. Default is use.Gram=TRUE",
    ## from  lars manual.
    ##for(i in 1:30){p<- predict(object=lars.result,F_[ ,considered.features]); plot(Re(p$fit[,i]),ylim=c(-2,2),main=i)}
    ##message("lars()"); a()	## CHECK POINT for testing
    global.errors <- (1/(2*dim(X_)[1])) *lars.result$RSS[-1] + (lars.result$lambda * rowSums(abs(coef(lars.result)))[-1])	
    ## This vector is needed for Bolasso.
    ## The global error including the classification error and the regularization penalty.
    ## ||WX-Y||_2 + lambda * ||W||_1
    ## Bach's methodology will need its minimum  over the best value for the regularization parameter.
    ## The following might also be needed somewhere in this function.
    ##best.global.error.index <- which.min(global.error[allowed.num.of.features])[1]		
    ##features.with.best.global.error <- selected.features.sequence[[best.global.error.index]]

    selected.features.sequence <- list()	## The ith entry will be the names are features when i features are allowed. 		
    all.linear.scores <- predict(object=lars.result,F_[ ,considered.features])$fit
    all.cofs.lars <- coef(lars.result)
    allowed.num.of.features <- 1:(min(maximum.features.num, dim(all.cofs.lars)[2]))	## Lasso for these number of features will be studied.
    for (num.of.features in allowed.num.of.features ){		##num.of.featuresnum.of.featuresnum.of.featuresnum.of.featuresnum.of.feat
        
        ##message(num.of.features); a()
        ##ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
        ## linear coefficients:
        linear.cofs <- c()
        linear.cofs[considered.features] <- 0
        cof.lars <- try(all.cofs.lars[num.of.features+1, ],silent=!report.fitting.failure)
        if(inherits(cof.lars, "try-error")){
            if (report.fitting.failure)
                warning("Linear model could not be fitted for ",num.of.features," features.")
            next
        }##End if.
        linear.cofs[names(cof.lars)] <- cof.lars	
        ## length of cof may not be appropriate!! 
        ## This is because some features might be ignored in the 
        ## training sample selection stage. 			
        ## selected features:
        selected.feature.names <- names(which(linear.cofs!=0))
        selected.features.sequence[[num.of.features]] <- selected.feature.names
        
        ## Keeping track of features,
        contribution.to.feature.scores[selected.feature.names] <- contribution.to.feature.scores[selected.feature.names] + 1/length(selected.feature.names)
        ## score_{i+1} = score_i + 1/(number of features)
        contribution.to.features.frequency[selected.feature.names] <- contribution.to.features.frequency[selected.feature.names] + 1
        ##In above line, we add to the score of currently live features.
        ##ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

        ## Computing score from linear model 
        linear.scores <- all.linear.scores[ ,num.of.features+1]

        ## fitting sigmoid (logistic regression):	
        compute.logistic.score.result <- try(compute.logistic.score(F_=F_,L_=L_,
                                                                    considered.features=considered.features,
                                                                    training.samples=training.samples,												
                                                                    validating.samples=validating.samples,
                                                                    linear.scores=linear.scores,
                                                                    report.fitting.failure=report.fitting.failure)
                                             ,silent=!report.fitting.failure)
        if(inherits(compute.logistic.score.result, "try-error"))
            stop("compute.logistic.score() failed.")	
        ## fitting a Logistic Regression Model might not be possible => error.
        ## It may not be able to fit the model if the sores are too far from the labels.
        ## This can heppen because the scores are computed for the validating patients and maybe far from reality. 	
        logistic.cofs <- compute.logistic.score.result$logistic.cofs			
        logistic.scores <- compute.logistic.score.result$logistic.scores
        ##message("c.l.s()"); plot(linear.scores,ylim=c(-3,3),col="green",main=num.of.features); points(L_); points(logistic.scores,col="red"); a(); 	
        ## CHECK POINT for testing						


        ## validating: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        ## to validate the model and comparing the errors.
        ## message("v()-start"); plot(L_[validating.samples],ylim=c(-2,2),main=num.of.features); points(logistic.scores[validating.samples],col='red');a()
        validation <- doctor.validate(true.labels=L_[validating.samples],predictions=logistic.scores[validating.samples])
        f.measure <- validation$f.measure
        mislabeled <- validation$mislabeled

        ## Updating best scores:
        if ( f.measure > best.f.measure ){	## This f.measure is better than others,
            best.logistic.scores <- logistic.scores[validating.samples]
            best.logistic.cof <- logistic.cofs
            best.feature.num <- num.of.features
        }##End if.
        
        ## Tracking mislabeled:
        if(length(mislabeled) > 0 )
            contribution.to.mislabeling.record[mislabeled] <- contribution.to.mislabeling.record[mislabeled]+1

        f.measure.log[num.of.features] <- f.measure
        ##message("validate()"); message(mislabeled);a()
        ## CHECK POINT for testing						
        ##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv		
    }##End for (num.of.features.num.of.featuresnum.of.featuresnum.of.featuresnum.of.featuresnum.of.featuresnum.of.features.

    ## accuracy:	
    accuracy <- f.measure.log[best.feature.num]		## the best f.measure for this doctor.

    ## selecting features based on global error
    best.global.error.index <- which.min(global.errors[allowed.num.of.features])
    features.with.best.global.error <- selected.features.sequence[[best.global.error.index]] 

    ## The doctor is trained.
    if(return_linear.models){
        lars.result_for_doctor=lars.result
    } else {
        lars.result_for_doctor= "return_linear.models is set to FALSE"
    }##End elseif.
    doctor <- list(linear.models=lars.result_for_doctor, 
                   best.number.of.features=best.feature.num,
                   probabilities=best.logistic.scores,
                   accuracy = accuracy,
                   best.logistic.cof=best.logistic.cof,
                   contribution.to.feature.scores= contribution.to.feature.scores,
                   contribution.to.features.frequency=contribution.to.features.frequency,
                   training.samples=training.samples,
                   validating.samples=validating.samples,
                   precision=validation$precision,
                   recall =validation$recall,
                   selected.features.sequence=selected.features.sequence,
                   global.errors=global.errors,
                   features.with.best.global.error=features.with.best.global.error)

    return(doctor)
    
    ##MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
}##End train.doctor <- function.	

