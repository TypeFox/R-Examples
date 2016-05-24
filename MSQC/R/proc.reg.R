proc.reg <-
function(x,alpha=0.0027,...){
p <- ncol(x) # number of quality characteristics
m <- nrow(x) # sample number
Xmv <- colMeans(x) # mean vector
S <- cov(x) # covariance matrix
Sinv <- solve(S)
s3<-matrix(0,p-1,p-1)
 LPL <- UPL <- matrix(0,nrow = p,ncol = 1)
         
for(i in 1:p){
             ifelse(nrow(s3)>1,s3 <- (Sinv[-i,-i]),s3 <- matrix(Sinv[-i,-i]))
             LPL[i,1] <- Xmv[i] - sqrt((det(s3) * (qchisq((1 - alpha),df = p))) / det(Sinv))
             UPL[i,1] <- Xmv[i] + sqrt((det(s3) * (qchisq((1 - alpha),df = p))) / det(Sinv))
             }
 return(list ("LPL" = LPL,"UPL" = UPL))
          }
