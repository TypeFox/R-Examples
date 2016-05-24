#####################################################
# Author: Chong Wu
# Email: wuxx0845@umn.edu
# The time consuming part is written in C
# Intermidate function, we mask them for the users
######################################################
MiSPUR <- function(Y, X, X2, cov = NULL, model=c("gaussian","binomial"), pow=c(2:8, Inf), n.perm=1000){
    
    model = match.arg(model)
    #pow=c(2:8, Inf)
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    #### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        r<-Y-mean(Y)
        U<-as.vector(t(X) %*% r)
    } else {
        tdat1 <-data.frame(trait=Y, cov)
        fit1 <-glm(trait~.,family=model,data=tdat1)
        pis <-fitted.values(fit1)
        r<-Y - pis
        
        U<-t(X) %*% r
    }
    
    ##observed statistics
    Ts=rep(NA,length(pow))
    for (j in 1:length(pow)){
        if (pow[j]<Inf) Ts[j] = sum(U^pow[j]) else Ts[j] = max(abs(U))
    }
    
    ################################################
    # For the X2, the second variable           ####
    ################################################
    #### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        r2<-Y-mean(Y)
        U2<-as.vector(t(X2) %*% r2)
    } else {
        # the null model is the same, so we just use pis.
        r2<-Y - pis
        U2<-t(X2) %*% r2
    }
    
    ##observed statistics
    Ts2=rep(NA,length(pow))
    for (j in 1:length(pow)){
        if (pow[j]<Inf) Ts2[j] = sum(U2^pow[j]) else Ts2[j] = max(abs(U2))
    }
    
    ########################################
    pow[pow==Inf] = 0 # pass 0 as infitiy
    TsC = MiSPUC(t(X),as.matrix(r), t(X2),as.matrix(r2),as.matrix(pow),n.perm)
    
    ###########################################################
    ## residual permutation for the first variable X        ###
    ###########################################################
    pPerm0 = rep(NA,length(pow))
    
    T0s = TsC$X1
    for (j in 1:length(pow)){

        pPerm0[j] =  sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm
        P0s = ( (n.perm-rank(abs(T0s[,j]))) + 1 )/ n.perm
        if (j==1) minp0 = P0s else minp0[which(minp0>P0s)] = P0s[which(minp0>P0s)]
    }
    
    
    #    cat("P0s caculated","\n")
    Paspu<-(sum(minp0<=min(pPerm0))+1)/(n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    
    pow[pow==0] = Inf
    
    names(Ts) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs) = names(Ts)
    
    Unweighted = list(Ts = Ts, pvs = pvs)

    ############################################################
    ## residual permutation for the second variable X        ###
    ############################################################
    pPerm02 = rep(NA,length(pow))

    T02s = TsC$X2
    for (j in 1:length(pow)){
        
        pPerm02[j] =  sum(abs(Ts2[j])<=abs(T02s[,j])) / n.perm
        P0s = ( (n.perm-rank(abs(T02s[,j]))) + 1 )/(n.perm)
        if (j==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
    }
    
    
    #    cat("P0s caculated","\n")
    Paspu2<-(sum(minp0<=min(pPerm02))+1)/(n.perm+1)
    pvs2 <- c(pPerm02, Paspu2)
    
    Ts2 <- c(Ts2, min(pPerm02))
    
    pow[pow==0] = Inf
    
    names(Ts2) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs2) = names(Ts2)
    Weighted = list(Ts = Ts2, pvs = pvs2)
    
    #############################################################
    ## residual permutation to summarize two variables        ###
    #############################################################
    pPerm03 = rep(NA,2*length(pow))
    
    for (j in 1:(2*length(pow))){
        
        if(j <= length(pow)) {
            pPerm03[j] = sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm
            P0s = ( (n.perm-rank(abs(T0s[,j]))) + 1 )/(n.perm)
            if (j==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
        } else {
            j1 = j - length(pow)
            pPerm03[j] = sum(abs(Ts2[j1])<=abs(T02s[,j1])) / n.perm
            P0s = ( (n.perm-rank(abs(T02s[,j1]))) + 1 ) / n.perm
            minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
            
        }
    }
    Paspu3<-(sum(minp0<=min(pPerm03))+1)/(n.perm+1)

    aMiSPU = list(Ts =min(pPerm03) ,pvalue = Paspu3)
    
    res= list(Unweighted = Unweighted,Weighted = Weighted,aMiSPU= aMiSPU)
    res

}














