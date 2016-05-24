RPParallel <-
function(XTrain #n by p trining data matrix
       , YTrain #n vector of the classes of the trining samples
       , XVal #n.val by p validation data matrix
       , YVal #n.val vector of the classes of the validation samples
       , XTest  #n.test by p test data matrix
       , d      #dimension to project into
       , B1 = 100 #number of blocks
       , B2 = 100 #block size
       , base = "LDA" # base classifier, eg "knn","LDA","QDA", or "Other"
       , projmethod = "Haar"
       , estmethod = "resub"
       , k = c(3,5) # possible k
       , cores = 2 # number of computer cores available
       , splitsample = FALSE #split sample Yes/No
       , ... )
{
    if (splitsample == FALSE){
        RP.out <- simplify2array(mclapply(rep(1,B1), function(x){return(RPChoose(XTrain, YTrain, XTest, d, B2, base, k, projmethod, estmethod))}, mc.cores = cores))
    }
    if (splitsample == TRUE){
        n.val <- length(YVal)
        RP.out <- simplify2array(mclapply(rep(1,B1), function(x){return(RPChooseSS(XTrain, YTrain, XVal, YVal, XTest, d, B2, base, k, projmethod))}, mc.cores = cores))
    }
    return (RP.out)
}
