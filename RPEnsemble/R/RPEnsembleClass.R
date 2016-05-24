RPEnsembleClass <-
function(RP.out # the result of a call to RPEnsembleOut
        ,n  # training sample size
        ,n.val #validation set size if samplesplit = TRUE
        ,n.test #test sample size
        ,p1 #(estimate of) prior probability
        ,splitsample = FALSE #split sample Yes/No
        ,alpha  #voting cutoff alpha
        ,... )
    {
    if (splitsample == FALSE){
        Test.Class <- RP.out[n + 1:n.test, ]
    }
    if (splitsample == TRUE){
        Test.Class <- RP.out[n.val + 1:n.test, ]
    }
    vote <- rowMeans(Test.Class, na.rm = TRUE)
    Class  <- 1 + as.numeric(vote > alpha)
    return(Class)
}
