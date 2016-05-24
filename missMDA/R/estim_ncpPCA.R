estim_ncpPCA <- function(X,ncp.min=0,ncp.max=5,method=c("Regularized","EM"),scale=TRUE,method.cv=c("gcv","loo","Kfold"),nbsim=100,pNA=0.05,threshold=1e-4){

## method = "em" or "Regularized"
## method.cv = "loo" (for leave-one-out) or "Kfold" (a percentage of pNA missing values is added and nbsim are done)
method <- match.arg(method,c("Regularized","EM","em","regularized"),several.ok=T)[1]
method.cv <- match.arg(method.cv,c("gcv","loo","Kfold","GCV","kfold","LOO"),several.ok=T)[1]

method <- tolower(method)
method.cv <- tolower(method.cv)
auxi = NULL
for (j in 1:ncol(X)) if (!is.numeric(X[,j])) auxi = c(auxi,colnames(X)[j])
if (!is.null(auxi)) stop(paste("\nThe following variables are not quantitative: ", auxi))
ncp.max <- min(ncp.max,ncol(X)-1,nrow(X)-2)
res <- NULL

if (method.cv=="gcv") {
p=ncol(X)
n=nrow(X)
if (is.null(ncp.max)) ncp.max <- ncol(X)-1
ncp.max <- min(nrow(X)-2,ncol(X)-1,ncp.max)
crit <- NULL
    if (ncp.min == 0) crit = mean((X - rep(colMeans(X, na.rm = TRUE), each = nrow(X)))^2, na.rm = TRUE)
    for (q in max(ncp.min, 1):ncp.max) {
#        res.pca = PCA(imputePCA(X,scale=scale,ncp=q,method=method,maxiter=1000)$completeObs, scale = scale, graph = FALSE, ncp = max(q, 2))
#        rec = reconst(res.pca, ncp = q)
        rec = imputePCA(X,scale=scale,ncp=q,method=method,maxiter=1000)$fittedX
        crit = c(crit, mean(((n * p - sum(is.na(X))) * (X - rec)/((n-1) * p - sum(is.na(X)) - q * (n + p - q-1)))^2, na.rm = T))
	}
  if (any(diff(crit)>0)) { ncp = which(diff(crit)>0)[1]
  } else ncp <- which.min(crit)
 names(crit) <- c(ncp.min:ncp.max)
 return(list(ncp = as.integer(ncp+ncp.min-1),criterion=crit))
}

if (method.cv=="loo"){
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
 for (nbaxes in ncp.min:ncp.max){
   Xhat <- X
   for (i in 1:nrow(X)){
    for (j in 1:ncol(X)){
     if (!is.na(X[i,j])){
      XNA <- as.matrix(X)
      XNA[i,j] <- NA
      if (nbaxes==0) Xhat[i,j] <- mean(XNA[,j],na.rm=TRUE)
      else Xhat[i,j] <- imputePCA(XNA,ncp=nbaxes,threshold=threshold,method=method,scale=scale)$completeObs[i,j]
    }
   }
    setTxtProgressBar(pb, round((((1:length(ncp.min:ncp.max))[which(nbaxes==(ncp.min:ncp.max))]-1)*nrow(X)+i)/(length(ncp.min:ncp.max)*nrow(X))*100))  
  }
  res <- c(res,mean((Xhat-X)^2,na.rm=TRUE))
 }
  close(pb)
 names(res) <- c(ncp.min:ncp.max)
 result = list(ncp = as.integer(which.min(res)+ncp.min-1),criterion=res)
}

if (method.cv=="kfold"){
  res <- matrix(NA,ncp.max-ncp.min+1,nbsim)
  pb <- txtProgressBar(min=1/nbsim*100, max=100,style=3)
  for (sim in 1:nbsim){
   XNA <- as.matrix(X)
   XNA[sample(1:(nrow(XNA)*ncol(XNA)),round(pNA*nrow(XNA)*ncol(XNA),0))] <- NA
   for (nbaxes in ncp.min:ncp.max){
    if (nbaxes==0) {
       Xhat <- XNA
       for (j in 1:ncol(X)) Xhat[,j] <- replace(XNA[,j],is.na(XNA[, j]),mean(XNA[,j],na.rm=TRUE))
    } else Xhat <- imputePCA(XNA,ncp=nbaxes,threshold=threshold,method=method,scale=scale)$completeObs  
   res[nbaxes-ncp.min+1,sim] <- sum((Xhat-X)^2,na.rm=TRUE)
  }
   setTxtProgressBar(pb, sim/nbsim*100)
 }
  close(pb)
 resu <- apply(res,1,mean)
 names(resu) <- c(ncp.min:ncp.max)
result <- list(ncp = as.integer(which.min(resu)+ncp.min-1),criterion=resu)
}
return(result)
}
