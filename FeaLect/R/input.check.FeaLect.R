input.check.FeaLect <- function(F_, L_, maximum.features.num,gamma){
    ##__________________________________________________________________________________________
    ## Checks if the input are desirable for FeLect() function.
    ## INPUT: 	F is the features matrix.
    ##		Y should be the vector of lebels ordered according to rows of F.
    ##                 The lables are -1 or +1 for negative and positive instances.
    ##		maximum.num.of.tries.for_randomly.choosing.samples (=1000 as default),
    ##		  If we try this many times and the optained lables are all the same,
    ##		  we give up (maybe the WHOLE lables are the same) with an error message.
    ##_________________________________________________________________________________________________________________________________________
    
    ## Input check:
    if(class(F_)!="matrix" | class(L_)!="numeric")
        stop(paste("Input error! F should be a matrix, L should be a numeric vector."))

    if(dim(F_)[1]!=length(L_))
        stop(paste("Number of rows of F (input feature matrix) and length of L (vector of labels) should be the same."))		
    
    ## feature names:
    if (is.null(colnames(F_))){
        colnames(F_) <- 1:dim(F_)[2]
    }##End if.
    inds <- which(is.na(colnames(F_)))
    colnames(F_)[inds] <- paste("NA",inds,sep="")

    ## sample names:
    if (is.null(rownames(F_)) | is.null(names(L_))){
        rownames(F_) <- 1:dim(F_)[1]
        names(L_) <- 1:length(L_)
    }##End if.	
    
    ## upper bound on the number of features
    if (maximum.features.num > gamma*length(unique(rownames(F_))) ){
        maximum.features.num <- round(gamma * length(unique(rownames(F_))),0)
        warning("\nmaximum.features.num was more than gamma times the number of instances (number of distinct rows of F_). \n",
                "It was automatically reduced to: ",maximum.features.num, 
                " (gamma times the number of instances,) \n", 
                "because Lasso would not select more than this number of features by any bound value for lambda. \n",
                "The technical reason is that the linear equation will have exact solution with error zero!\n")
    }##End if.

    ## upper bound on the number of features
    if (prod(rownames(F_) == names(L_))==0 ){
        warning(" rownames(F_) != rownames(L_) !!! ",
                "\n The names of row of the labels (L_) should be exatcly the same as row names of features matrix (F_). \n")
    }##End if.
    
    return(list(F_=F_,L_=L_,maximum.features.num=maximum.features.num))
}##End input.check.FeaLect <- function.

