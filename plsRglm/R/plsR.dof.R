plsR.dof <- function(modplsR,naive=FALSE){
temp.object <- vector("list",0)
dof.object <- vector("list",0)
temp.object$RSS <- as.vector(modplsR$RSS)
temp.object$Yhat <- matrix(mean(modplsR$dataY),nrow=modplsR$nr,ncol=1)
for(i in 1:modplsR$computed_nt){
temp.object$Yhat <- cbind(temp.object$Yhat,attr(modplsR$RepY, "scaled:center") + attr(modplsR$RepY, "scaled:scale") * modplsR$tt[,1:i,drop=FALSE] %*% modplsR$CoeffC[1:i])
}
colnames(temp.object$Yhat) <- paste("Nt_",0:modplsR$computed_nt,sep="")
if(!naive){
temp.object$TT <- sweep(modplsR$tt,MARGIN=2,FUN="/",sqrt(diag(crossprod(modplsR$tt))))

require(plsdof)
pls.doftemp <- function (pls.object, n, y, K, m, DoF.max) 
{
    TT <- pls.object$TT
    Yhat <- pls.object$Yhat[, 2:(m + 1), drop=FALSE]
    TK = matrix(, m, m)
    KY <- plsdof::krylov(K, K %*% y, m)
    lambda <- eigen(K)$values
    tr.K <- vector(length = m)
    for (i in 1:m) {
        tr.K[i] <- sum(lambda^i)
    }
    BB = t(TT) %*% KY
    BB[row(BB) > col(BB)] = 0
    b <- t(TT) %*% y
    DoF = vector(length = m)
    Binv <- backsolve(BB, diag(m))
    tkt <- rep(0, m)
    ykv <- rep(0, m)
    KjT <- array(dim = c(m, n, m))
    dummy <- TT
    for (i in 1:m) {
        dummy <- K %*% dummy
        KjT[i, , ] <- dummy
    }
    trace.term = rep(0, m)
    for (i in 1:m) {
        Binvi <- Binv[1:i, 1:i, drop = FALSE]
        ci <- Binvi %*% b[1:i]
        Vi <- TT[, 1:i, drop = FALSE] %*% t(Binvi)
        trace.term[i] <- sum(ci * tr.K[1:i])
        ri <- y - Yhat[, i]
        for (j in 1:i) {
            KjTj = KjT[j, , ]
            if(is.null(dim(KjTj))){KjTj <- matrix(KjTj,ncol=1)}
            tkt[i] <- tkt[i] + ci[j] * plsdof::tr(t(TT[, 1:i, drop = FALSE]) %*% 
                KjTj[, 1:i, drop = FALSE])
            ri <- K %*% ri
            ykv[i] <- ykv[i] + sum(ri * Vi[, j])
        }
    }
    DoF <- trace.term + 1:m - tkt + ykv
    DoF[DoF > DoF.max] = DoF.max
    sigmahat = sqrt(pls.object$RSS[-1]/(n - DoF))
    return(list(DoF = DoF, sigmahat = sigmahat))
}

dof.object <- pls.doftemp(temp.object,n=modplsR$nr,y=modplsR$dataY,K=modplsR$ExpliX%*%t(modplsR$ExpliX),m=modplsR$computed_nt,DoF.max=min(modplsR$nr-1,modplsR$nc+1)-1)
temp.object$sigmahat <- c(sqrt(temp.object$RSS[1]/(modplsR$nr-1)), sqrt(temp.object$RSS[-1]/(modplsR$nr-dof.object$DoF)))
temp.object$DoF <- c(0, dof.object$DoF) + 1
} else
{
temp.object$DoF <- 1:(modplsR$computed_nt+1)
temp.object$sigmahat <- sqrt(temp.object$RSS/(modplsR$nr - temp.object$DoF))
}
temp.object$yhat <- vector("numeric",length=modplsR$computed_nt+1)
for(i in 1:(modplsR$computed_nt+1)){
temp.object$yhat[i] <- sum((temp.object$Yhat[,i])^2)
}
return(list(DoF=temp.object$DoF, sigmahat=temp.object$sigmahat, Yhat=temp.object$Yhat, yhat=temp.object$yhat, RSS=temp.object$RSS))
}
