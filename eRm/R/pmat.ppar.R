`pmat.ppar` <-
function(object)
# computes a list of expected probabilities for objects of class "ppar" for each NA-subgroup
# without category!
{

X <- object$X
mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
mt_ind <- rep(1:length(mt_vek),mt_vek)

rp <- rowSums(X,na.rm=TRUE)
maxrp <- sum(mt_vek)
TFrow <- ((rp==maxrp) | (rp==0))

pmat.l <- lapply(object$thetapar, function(theta1) {                       #runs over missing structures
             theta <- theta1
             p.list <- tapply(object$betapar,mt_ind,function(beta.i) {     #matrices of expected prob as list (over items)
                     beta.i <- c(0,beta.i)
                     ind.h <- 0:(length(beta.i)-1)
                     theta.h <- ind.h %*% t(theta)
                     tb <- exp(theta.h+beta.i)
                     denom <- colSums(tb)
                     pi.mat <- apply(tb,1,function(y) {y/denom})
                     return(pi.mat)
                   })
    
    p.list0 <- lapply(p.list,function(pl) {rbind(pl)[,-1]})               #delete 0th category
    pmat <- matrix(unlist(p.list0),nrow=length(theta1))      #save as matrix
    return(pmat)
  }) 

#----------item-category labels----------
cnames <- substr(names(object$betapar),6,40)
for (i in 1:length(pmat.l)) dimnames(pmat.l[[i]]) <- list(names(object$thetapar[[i]]),cnames)
#-----------end labels-------   

if (length(object$pers.ex) > 0) {    
      X <- object$X[-object$pers.ex,]                                        #list with raw scores
      X01 <- object$X01[-object$pers.ex,]     
  } else {
      X <- object$X
      X01 <- object$X01
  }

NApos <- tapply(1:length(object$gmemb),object$gmemb,function(ind) {   #positions for NA replacement
                       xvec <- X01[ind[1],]
                       which(is.na(xvec))
                       })

pmat <- NULL
for (i in 1:length(pmat.l)) {
       pmat.l[[i]][,NApos[[i]]] <- NA            #insert NA's
       pmat <- rbind(pmat,pmat.l[[i]])
       }

#-------------- reorder the p-matrix ---------------      
ind.orig <- as.vector(unlist(tapply(1:length(object$gmemb), object$gmemb, function(ind) {ind})))
pmat.orig.list <- by(pmat, ind.orig, function(ii) return(ii))
pmat.orig <- as.matrix(do.call(rbind, pmat.orig.list))      #final P-matrix (corresponding to X)
rownames(pmat.orig) <- rownames(X)

return(pmat.orig)
}

