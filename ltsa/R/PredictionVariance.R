`PredictionVariance` <-
function(r, maxLead=1, DLQ=TRUE){
    if (abs(r[1]) <= 0) 
        stop(paste("error - lag zero autocovariance not positive!  r[1] =", r[1], ", must be <1"))
    if (DLQ) {
        R<-r[-1]/r[1]
        out<-DLAcfToAR(R, useC=FALSE)
        sigsq <- out[nrow(out),3] #innovation variance
        phi<-out[,1]
        Q<-max(length(phi), maxLead-1, length(r))
        if (Q > 0 && maxLead > 1)
                psi<-c(1,ARMAtoMA(ar=phi, lag.max=Q))[1:maxLead]
        else
                psi<-1
        VL<-r[1]*sigsq*cumsum(psi^2)
    }
    else {
        n<-length(r)-maxLead
        if (n < 1) stop ("error: length(r) needs to be > maxLead!")
        gk<-t(matrix(r[n + 1 + outer(1:maxLead,1:n,"-")],ncol=maxLead,byrow=TRUE))
        GI<-TrenchInverse(toeplitz(r[1:n]))
        gkGI<-crossprod(t(gk),GI)
        VL<-numeric(maxLead)
        for (j in 1:maxLead)
            VL[j]<-r[1] - sum(gkGI[j,]*gk[j,])
        }
        if (any(VL<=0)) stop("error: ACVF has not pd sequence!")
        VL
    }

