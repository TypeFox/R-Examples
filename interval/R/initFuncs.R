initcomputeMLE<-function(L,R,Lin=NULL,Rin=NULL,A=NULL,max.inner=10,max.outer=1000,tol=1e-09){
    n<-length(L)
    ## computeMLE was developed for bivariate interval censored data
    ## but we can trick it into doing univariate by making the second
    ## variable all the same, that is why columns 3 and 4 of R are all 0 and 1
    R<-matrix(c(L,R,rep(0,n),rep(1,n)),n,4)
    ## put Lin and Rin into matrix form for computeMLE, 1=TRUE, 0=FALSE
    ## again the 3 and 4th column represent the second variable that are all 0 and 1
    if (is.null(Lin) & is.null(Rin)){
        B<-c(0,1)
    } else if (length(Lin)==1 & length(Rin)==1){
        if (Lin){ bL<-1
        } else bL<-0
        if (Rin){ bR<-1
        } else bR<-0
        B<-c(bL,bR)
    } else if (length(Lin)==n & length(Rin)==n){
        B<-matrix(c(Lin,Rin,rep(0,n),rep(1,n)),n,4)
    } 
    mle<-computeMLE(R,B,max.inner=max.inner,max.outer=max.outer,tol=tol)
    ## computeMLE outputs a rects matrix with coordinates of bivariate 
    ## rectangles where the masses of the MLE reside
    ## ignore 3rd and 4th columns
    intmap<-t(mle$rects[,1:2])
    ## take boundary and convert it to Lin and Rin attributes
    k<- length(mle$p)
    if (is.vector(mle$bounds)){
        intmapLin<-rep(mle$bounds[1]==1,k)
        intmapRin<-rep(mle$bounds[2]==1,k)
    } else {
        intmapLin<- mle$bounds[,1]==1
        intmapRin<- mle$bounds[,2]==1
    }
    attr(intmap,"LRin")<-matrix(c(intmapLin,intmapRin),byrow=TRUE,nrow=2)

    out<-list(bounds=mle$bounds,conv=mle$conv,llh=mle$llh,intmap=intmap,pf=mle$p)
    class(out)<-"icfit"
    out
}

initEMICM<-function(L=NULL,R=NULL,Lin=NULL,Rin=NULL,A=NULL,maxiter=1000,tol=1e-7){
    ## since EMICM allows inputing the A matrix, we use that form
    ## except we need to transpose our input A matrix from n x m
    ## to m x n, where n is the number of observations
    EMICM(t(A),maxiter=maxiter,tol=tol)
}
