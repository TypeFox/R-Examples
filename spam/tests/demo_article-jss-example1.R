# This is file ../spam/tests/demo_article-jss-example1.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     




# JSS article:
#     "spam: A Sparse Matrix R Package with Emphasis on
#            MCMC Methods for Gaussian Markov Random Fields"
#
# Compared to the R code given in the article, here we give:
# - improved formatting
# - more comments
# - the R code to construct the figures



# SETUP:
library("spam")
spam.options(structurebased=TRUE)
data("UKDriverDeaths")

y <- sqrt(c(UKDriverDeaths))       # square root counts

n <- length(y)                     # n=192
m <- 12                            # We want to predict for one season.
nm <- n+m                          # Total length of s and t


priorshape <-  c(4, 1, 1)          # alpha's, as in Rue & Held (2005)
priorinvscale <- c(4, 0.1, 0.0005) # beta's 

# Construct the individual block precisions
# (based on unit precision parameters kappa, denoted with k):

# Qsy, Qty are trivial:
Qsy <- diag.spam(n)
pad(Qsy) <- c(n+m, n)  # previously:  dim(Qsy) <- c(n+m, n)

Qty <- Qsy

Qst <- spam(0, nm, nm)
Qst[cbind(1:n, 1:n)] <- rep(1, n)


# The form of Qss is given by (Rue and Held equation 3.59).
# Qss can be constructed with a loop:
Qss <- spam(0, nm, nm)
for (i in 0:(nm-m)) {
    Qss[i+1:m,i+1:m] <- Qss[i+1:m, i+1:m] + matrix(1,m,m)
#    Qss[i+1:m,i+1:m] <- Qss[i+1:m, i+1:m]+1  # previously...
}

# Note that for the final version we need:
# Qss <- k_s * Qss + k_y * diag.spam(nm)  




# The form of Qtt is given by (Rue and Held equation 3.40).
# Similar approaches to construct Qtt:

Qtt <- spam(0,nm,nm)
Qtt[cbind(1:(nm-1),2:nm)] <- -c(2,rep(4,nm-3),2)
Qtt[cbind(1:(nm-2),3:nm)] <- rep(1,nm-2)
Qtt <- Qtt + t( Qtt)
diag(Qtt) <- c(1,5,rep(6,nm-4),5,1)



# Create temporary kappa and precision matrix to illustrate
# adjacency matrix and ordering.
k <- c(1,1,1)
Qst_yk <- rbind(cbind(k[2]*Qss + k[1]*diag.spam(nm), k[1]*Qst),
       	        cbind(k[1]*Qst, k[3]*Qtt + k[1]*diag.spam(nm)))
                
struct <- chol(Qst_yk)

        

# Note that we do not provide the exactly the same ordering 
# algorithms. Hence, the following is sightly different than
# Figure RH4.2.
cholQst_yk <- chol(Qst_yk,pivot="RCM")
P <- ordering(cholQst_yk)
display(Qst_yk)
display(Qst_yk[P,P])



# Recall:
# k=( kappa_y, kappa_s, kappa_t)'

# Gibbs sampler
ngibbs <- 100   # In the original version is 500!
burnin <- 10    # > 0
totalg <- ngibbs+burnin
set.seed(14)

# Initialize parameters:
spost <- tpost <- array(0, c(totalg, nm))
kpost <- array(0, c(totalg, 3)) 

# Starting values:
kpost[1,] <- c(.5,28,500)
tpost[1,] <- 40

# calculation of a few variables:
postshape <- priorshape + c(	n/2, (n+1)/2, (n+m-2)/2) 


for (ig in 2:totalg) {
    
  Q <- rbind(cbind(kpost[ig-1,2]*Qss + kpost[ig-1,1]*Qst, 
                   kpost[ig-1,1]*Qst),
             cbind(kpost[ig-1,1]*Qst,  
                   kpost[ig-1,3]*Qtt + kpost[ig-1,1]*Qst))
  
  b <- c(kpost[ig-1,1]*Qsy %*% y, kpost[ig-1,1]*Qsy %*% y)
  
  tmp <- rmvnorm.canonical(1, b, Q, Lstruct=struct) 
  
      
  spost[ig,] <- tmp[1:nm]		 

  tpost[ig,] <- tmp[1:nm+nm]


  tmp <- y-spost[ig,1:n]-tpost[ig,1:n]
  
  postinvscale <- priorinvscale + # prior contribution
    c( sum( tmp^2)/2,     # Qyy_st is the identity
      t(spost[ig,]) %*% (Qss %*% spost[ig,])/2,
      t(tpost[ig,]) %*% (Qtt %*% tpost[ig,])/2)


  kpost[ig,] <- rgamma(3, postshape, postinvscale)	

  if( (ig%%10)==0) cat('.')

}



# Eliminate burn-in:
kpost <- kpost[-c(1:burnin),]
spost <- spost[-c(1:burnin),]
tpost <- tpost[-c(1:burnin),]

postquant <- apply(spost+tpost, 2, quantile,c(.025,.975))
postmean  <- apply(spost+tpost, 2, mean)
postmedi  <- apply(spost+tpost, 2, median)

if (F){

par(mfcol=c(1,1),mai=c(.6,.8,.01,.01))

plot( y^2, ylim=c(800,2900),xlim=c(0,nm),ylab="Counts")
#lines( postmean^2, col=2)
lines( postmedi^2, col=2)
matlines( t(postquant)^2, col=4,lty=1)

legend("topright",legend=c("Posterior median", "Quantiles of posterior sample",
                    "Quantiles of predictive distribution"),
       bty="n",col=c(2,4,3),lty=1)




# Constructing a predictive distribution:
ypred <- rnorm( ngibbs*nm, c(spost+tpost),sd=rep( 1/sqrt(kpost[,1]), nm)) 
dim(ypred) <- c(ngibbs,nm)
postpredquant <- apply(ypred, 2, quantile,c(.025,.975))
matlines( t(postpredquant)^2, col=3,lty=1)
points(y^2)
dev.off() 

kpostmedian <- apply(kpost,2,median)

par(mfcol=c(1,3),mai=c(.65,.65,.01,.01),cex=.85,mgp=c(2.6,1,0))

matplot( log( kpost), lty=1, type="l",xlab="Index")
abline(h=log(kpostmedian),col=3)
acf( kpost[,3],ylab=expression(kappa[t]))
plot(kpost[,2:3],ylab=expression(kappa[t]),xlab=expression(kappa[s]),cex=.8)
abline(h=kpostmedian[3],v=kpostmedian[2],col=3)
dev.off()


allkappas <- rbind(apply(kpost,2,mean),
                   apply(kpost,2,median),
                   apply(1/kpost,2,mean),
                   apply(1/kpost,2,median))
colnames(allkappas) <- c("kappa_y", "kappa_s", "kappa_t")
rownames(allkappas) <- c("Prec (mean)", "Prec (median)",
                         "Var (mean)", "Var (median) ")
print(allkappas,4)

png("example1_m1.png",width=300,height=300)
par(mai=c(.5,.5,.05,.05))
display(Qst_yk)
dev.off()

png("example1_m2.png",width=300,height=300)
par(mai=c(.5,.5,.05,.05))
display(struct)
dev.off()


summary(kpost)


}
