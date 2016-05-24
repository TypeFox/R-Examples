`scoresFromFit` <- function(icFIT, scores, dqfunc=NULL){
    A<-icFIT$A
    k<-dim(A)[[2]]
    n<-dim(A)[[1]]
    if (is.null(dqfunc) & scores=="general") stop("must provide dqfunc for general scores")
    phat<-c(0,icFIT$pf)
    Shat<- 1-cumsum(phat)

    ## ckstar gives scores associated with each element of the pf vector
    if (scores=="logrank1"){
        Stilde<- exp( - c(0,cumsum(phat[2:(k+1)]/Shat[1:k])) )
        ckstar<- (1/phat[2:(k+1)])*( Shat[1:k]* log(Stilde[1:k]) - Shat[2:(k+1)]* log(Stilde[2:(k+1)]) )
    } else if (scores=="logrank2"){
        ckstar<- (1/phat[2:(k+1)])*c( 
            Shat[1:(k-1)]* log( Shat[1:(k-1)]) - Shat[2:(k)]* log( Shat[2:(k)]),
            Shat[k]* log( Shat[k]))
    } else if (scores=="wmw"){
        ckstar<-Shat[1:k] + Shat[2:(k+1)] - 1
    } else if (scores=="normal"){
        ckstar<- (-1/phat[2:(k+1)])*( 
            dnorm( qnorm( Shat[1:k])) - dnorm( qnorm( Shat[2:(k+1)] )))
    } else if (scores=="general"){
        if (is.null(dqfunc)) stop("must provide dqfunc for general scores")
        ckstar<- (-1/phat[2:(k+1)])*c( 
            dqfunc( Shat[1:k]) - dqfunc(Shat[2:(k+1)]))
    }

    p<-phat[2:(k+1)]
    tempfunc<-function(Arow){
        sum(Arow*p*ckstar)/sum(Arow*p)
    }     
    ## cc is the vector of ci values for each individual
    ## it is just the weighted average of the ckstar values for the possible intervals 
    ## (for the ith individual it is intervals with A[i,j]=1) weighted by the pf values    
    cc<- apply(A,1,tempfunc)
    cc
}
