FeaLect <- function(F, L, maximum.features.num=dim(F)[2], total.num.of.models, gamma =3/4,
                    persistence = 1000,talk=FALSE, minimum.class.size=2,report.fitting.failure=FALSE,
                    return_linear.models=TRUE,balance=TRUE, replace=TRUE,plot.scores=TRUE){
    ##FeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLect
    ## In the following loop,  a class is analized by cross-validating a portion of samples (gamma).
    ## In consequence, num.of.features.score, error.log, mislabeling.record and feature.scores are computed.
    ##____________________________________________________________________________________________________________________________________
    ## INPUT: 	F is the features matrix.
    ##			L should be the vector of lebels ordered according to rows of F. The lables are 0 or +1 for negative and positive instances.
    ##			persistence:	maximum.num.of.tries.for_randomly.choosing.samples (=1000 as default),
    ##							If we try this many times and the optained lables are all the same,
    ##							we give up (maybe the WHOLE lables are the same) with an error message.
    ##			gamma:	=3/4 by default. gamma fraction of samples are considered as training ones.
    ##			return_linear.models:	They are memory intensive, so for total.num.of.models >>1000 we may decide to ignore them.
    ##			minimum.class.size:		The size of both positive and negative classes after sampling should be greater than this threshould.
    ##			replace:				If true, the random subset will be sampled by replacment.
    ##____________________________________________________________________________________________________________________________________

    ## Input check:
    ##message("input.check.FeaLect()"); a()	## CHECK POINT for testing
    FL <- input.check.FeaLect(F_=F,L_=L,maximum.features.num=maximum.features.num,gamma=gamma)
    F_ <- FL$F_		## More readable code by renaming F to F_.
    L_ <- FL$L_
    updated.maximum.features.num <- FL$maximum.features.num
    all.feature.names <- colnames(F_)
    ##total.num.of.features <- length(all.feature.names)
    total.input.samples <- rownames(F_)
    
    ## Timing:
    time.start <-Sys.time()
    if(talk) message("***********************************************")
    if(talk) message("Scoring ",dim(F_)[2]," features using ",dim(F_)[1], " samples.")
    if(talk) message(" - started at: ", time.start)


    ## Initializing output obojects:			
    mislabeling.record <- c(); mislabeling.record[total.input.samples] <- 0
    feature.scores <-c() 				## It keeps track of the frequency of features.
    feature.scores[all.feature.names] <- 0
    features.frequency <- c()
    features.frequency[all.feature.names] <- 0;
    num.of.features.score <- c()
    num.of.features.score[all.feature.names] <- 0
    ## The score will be added if for a number of features the error is minimized.
    best.feature.num <- c()				## best.feature.num[sampling.index] will be the number with best f.measure.
    best.features.intersection <-c()	
    feature.scores.history <- data.frame()	
    ## feature.scores.history[5,] will be feature.scores at the time of fifth samling. 
    probability.history <- data.frame()	## the row probability.history["PAT8934",] will contain the probablity of beeing in class_F.
    
    doctors <- list()		## All information that is produced by fitting a linear model is stored here.
    ## In this code, the output of check() will resemble the knowledge of a doctor.

    features.with.best.global.error <- list()
    ## features.with.best.global.error[[4]] contains the set features which provide the 
    ## best global error over all possible values of bound.

    global.error.history <- data.frame()
    
    

    ## This is the MAIN LOOP.	LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    for (sampling.index in 1:total.num.of.models){	## Every time we choose some samples for training, we give an index to them.  
        ## counter and flags
        number.of.global.tries.for_this.samling.index <- 0
        enough.variation <- TRUE

        ## Loop
        repeat{ ##until you get no error
            number.of.global.tries.for_this.samling.index <- number.of.global.tries.for_this.samling.index +1
            if(number.of.global.tries.for_this.samling.index >= persistence ){
                stop(paste("FeaLect() was not successful at sampling.index= ", sampling.index,
                           "\n Meybe because the positive and negative classes in the sampled data are too unbalanced."))
            } ##End if.		
            
            try_sampling.index <- try({ 
                if(talk) message(" - sampling.index: ",sampling.index)						
                
##### selecting a random subset of samples:----
                ##message("random.subset()-starts"); a()	## CHECK POINT for testing
                XY <- try(random.subset(F_=F_, L_=L_,gamma=gamma, persistence=persistence, minimum.class.size=minimum.class.size,replace=replace))
                if(inherits(XY, "try-error")){
                    enough.variation=FALSE	
                    stop("random.subset() was not successful.") 
                }##End if.
                ## random subset is selected.
                X_ <- XY$X_		## a random submatrix of the feature matrix that contains training samples.
                ## feature names:
                features.names.for_this.sample <- colnames(X_)
                ## Training and validating samples:
                training.samples <- unique(rownames(X_)) 
                validating.samples.indices <- which(is.na(match(total.input.samples,training.samples))) ##the rest
                validating.samples <- total.input.samples[validating.samples.indices]
                ##message("random.subset()-ends"); a()	## CHECK POINT for testing

###	Training th edoctor: MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
                ##message("train.doctor()-start"); a()	## CHECK POINT for testing				
                doctor <- try(train.doctor(F_=F_,L_=L_, 
                                           considered.features=features.names.for_this.sample,
                                           balance=balance,return_linear.models=return_linear.models,
                                           training.samples=training.samples,validating.samples=validating.samples,
                                           maximum.features.num=updated.maximum.features.num,
                                           report.fitting.failure=report.fitting.failure),
                              silent=!report.fitting.failure)
                if(inherits(doctor,"try-error")){
                    if(talk) message("----------------again")
                    next
                }##End if.
                doctors[[sampling.index]] <- doctor							
                ##message("train.doctor()-end"); a()	## CHECK POINT for testing				


                ## learning from the doctor's experience
                ctfs <- doctor$contribution.to.feature.scores
                feature.scores[names(ctfs)] <- feature.scores[names(ctfs)] + ctfs
                feature.scores.history <- rbind(feature.scores.history, feature.scores/sampling.index)				
                ctff <- doctor$contribution.to.feature.frequency
                features.frequency[names(ctff)] <- features.frequency[names(ctff)] + ctff
                num.of.features.score[doctor$best.number.of.features]	<- num.of.features.score[doctor$best.number.of.features] + 1					

                ## global error needed by Bolasso
                features.with.best.global.error[[sampling.index]] <- doctor$features.with.best.global.error
                ## Storing the results of this run of the loop.
                ## At right, we want the indices of all minimum points of error.log[sampling.index,]
                ## global.error.history <- rbind(global.error.history, t(data.frame(global.error)) ) 
                ##The above does not work because of variable length of global.error
                log.scores <- sort(log(feature.scores),decreasing=FALSE)	
                color.temp <- match(names(log.scores),names(feature.scores))
                if(plot.scores)				
                    plot(log.scores,main=sampling.index, col=color.temp,pch=19)
                ##a()	## CHECK POINT for testing		
                
                
                
                ##MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
                ## Avareging error:
                ##most.frequent.num.of.features <- which.max(num.of.features.score)
                ##message(message("Most frequent number of features:", which.max(num.of.features.score)))
                ##average.error <- mean(t(error.log[ ,most.frequent.num.of.features]))
                ##message(message("Average errror is :", round(average.error,2)))
                ## best features according to Bach:
                if(sampling.index==1){
                    best.features.intersection <- features.with.best.global.error[[sampling.index]]
                }else{	
                    best.features.intersection <- intersect(features.with.best.global.error[[sampling.index]], best.features.intersection)
                }##End elseif(sampling.index==1).		
            })##End try_sampling.index <- try.	
            if(!enough.variation)
                stop("Variation in labels is not enough, maybe all lables are the same!")	
            ## It does not mean to continue sampling any more.
            if(!inherits(try_sampling.index,"try-error"))
                break	## That's fine for this sampling index, we do not need to repeat.
        }##End repeat
    }##End for (sampling.index.
    ##LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL	

    ## Timing:
    time.taken <- Sys.time()-time.start
    if(talk) message("****************************************************")
    if(talk) message("validation ended at: ", Sys.time(), "   taking:   ", time.taken)
    if(talk) message("****************************************************")
    
    ## Outputting:
    return(list(
        log.scores=log.scores, 
        feature.matrix=F_,
        labels=L_,
        total.num.of.models=total.num.of.models,
        maximum.features.num=maximum.features.num,
        feature.scores.history=feature.scores.history,
        num.of.features.score=num.of.features.score,
        best.feature.num=best.feature.num,
        mislabeling.record=mislabeling.record, 
        doctors=doctors, 
        best.features.intersection=best.features.intersection, 
        features.with.best.global.error=features.with.best.global.error,
        time.taken=time.taken))
    ##FeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLectFeaLec	
}##End FeaLect <- function().	

