
## crimeClust_bayes
##==============================================================================
#' Bayesian model-based partially-supervised clustering for crime series identification
#'
#'  @param crimeID n-vector of criminal IDs for the n crimes in the dataset.
#'    For unsolved crimes, the value should be \code{NA}.
#'  @param spatial (n x 2) matrix of spatial locations, represent missing locations
#'    with \code{NA}
#'  @param t1 earliest possible time for crime
#'  @param t2 latest possible time for crime. Crime occurred between \code{t1}
#'    and \code{t2}.
#'  @param Xcat (n x q) matrix of categorical crime features.  Each column is a
#'    variable, such as mode of entry.  The different factors (window, door, etc)
#'    should be coded as integers 1,2,\dots,m.
#'  @param Xnorm (n x p) matrix of continuous crime features.  
#'  @param maxcriminals maximum number of clusters in the model.
#'  @param iters Number of MCMC samples to generate.
#'  @param burn	Number of MCMC samples to discard as burn-in.
#'  @param plot (logical) Should plots be produced during run.
#'  @param update	Number of MCMC iterations between graphical displays.
#'  @param seed seed for random number generation 
#'  @param use_space  (logical) should the spatial locations be used in clustering?
#'  @param use_time  (logical) should the event times be used in clustering?
#'  @param use_cats (logical) should the categorical crime features be used in
#'    clustering?  
#'  @return (list) p.equal is the (n x n) matrix of probabilities that each pair of
#'    crimes are committed by the same criminal.
#'    
#'    if \code{plot=TRUE}, then progress plots are produced.
#'  @seealso \code{\link{bayesPairs}}
#'  @references
#'  Reich, B. J. and Porter, M. D. (2015), Partially supervised spatiotemporal
#'    clustering for burglary crime series identification.
#'    \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}.
#'    178:2, 465--480.
#'  \url{http://www4.stat.ncsu.edu/~reich/papers/CrimeClust.pdf}
#'  @author Brian J. Reich
##   (\emph{brian_reich@@ncsu.edu})
#'  @export
#'  @examples
#'  # Toy dataset with 12 crimes and three criminals. 
#'  
#'  # Make IDs: Criminal 1 committed crimes 1-4, etc.
#'  id <- c(1,1,1,1,
#'          2,2,2,2,
#'                  3,3,3,3)
#'                  
#'  # spatial locations of the crimes:
#'  s <- c(0.8,0.9,1.1,1.2,  
#'         1.8,1.9,2.1,2.2,
#'         2.8,2.9,3.1,3.2)
#'  s <- cbind(0,s)
#'  
#'  # Categorical crime features, say mode of entry (1=door, 2=other) and 
#'  # type of residence (1=apartment, 2=other)
#'  Mode <- c(1,1,1,1,  #Different distribution by criminal
#'            1,2,1,2,
#'            2,2,2,2)
#'  Type <- c(1,2,1,2,  #Same distribution for all criminals
#'            1,2,1,2,
#'            1,2,1,2)
#'  Xcat <- cbind(Mode,Type)
#'  
#'  # Times of the crimes
#'  t <- c(1,2,3,4,
#'         2,3,4,5,
#'         3,4,5,6)
#'  
#'  # Now let's pretend we don't know the criminal for crimes 1, 4, 6, 8, and 12.
#'  id <- c(NA,1,1,NA,2,NA,2,NA,3,3,3,NA)
#'  
#'  # Fit the model (nb: use much larger iters and burn on real problem)
#'  fit <- crimeClust_bayes(crimeID=id, spatial=s, t1=t,t2=t, Xcat=Xcat, 
#'                    maxcriminals=12,iters=500,burn=100,update=100)
#'                    
#'  # Plot the posterior probability matrix that each pair of crimes was 
#'  # committed by the same criminal:
#'  if(require(fields,quietly=TRUE)){
#'  fields::image.plot(1:12,1:12,fit$p.equal,
#'             xlab="Crime",ylab="Crime",
#'             main="Probability crimes are from the same criminal")
#'  }
#'  
#'  # Extract the crimes with the largest posterior probability
#'  bayesPairs(fit$p.equal)
#'  bayesProb(fit$p.equal[1,])
#'  
##  Check Xnorm, this was originally set up really to hold time. But needs extending 
##   to general. use_time should also include use_norm, etc. 
##==============================================================================
crimeClust_bayes<-function(crimeID,spatial,t1,t2,Xcat,Xnorm,
                           maxcriminals=1000,iters=10000,burn=5000,
                           plot=TRUE,update=100,seed=NULL,
                           use_space=TRUE,use_time=TRUE,use_cats=TRUE){
  
#----------- Helper Functions ------------------
  ## taken from MCMCpack package
  rdirichlet <- function (n, alpha) {
      l <- length(alpha)
      x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
      sm <- x %*% rep(1, l)
      return(x/as.vector(sm))
  }
  
  ddir<-function(y,df,p){
     alpha<-df*p
     lll<-lgamma(sum(alpha))+sum(alpha*log(y)) - sum(lgamma(alpha))
  lll}
  
  count<-function(j,y){sum(y==j)}
  
  get.counts<-function(y,M){
         sapply(1:M,count,y=y)
  }
  sum.by.g<-function(y,g,M){
     unig<-sort(unique(g))
     sss<-rep(0,M)
     sss[unig]<-tapply(y,g,sum)
  sss}
  
  count.g<-function(g,M){
     unig<-sort(unique(g))
     sss<-rep(0,M)
     sss[unig]<-table(g)
  sss}
  
  #Draw samples from a truncated normal:
  rtruncnorm<-function(n,mu,sigma,lower,upper){
     lp<-pnorm(lower,mu,sigma)
     up<-pnorm(upper,mu,sigma)
     y<-qnorm(runif(n,lp,up),mu,sigma)
  y}
  
  ddir2<-function(probs,D){
    k<-length(probs)
    lgamma(k*D)-k*lgamma(D)+(D-1)*sum(log(probs))
  }

#-------------------------------------------------  
  
  
  
  if(!is.null(seed)) set.seed(seed)
  if(missing(Xnorm)) Xnorm = NULL
  if(!missing(t1)){  
    if(missing(t2)) t2 = t1     # if time is not interval censored
    #-- If time is Date-time (Class POSIXct), then convert to fractional days since first crime
    origen = min(t1,t2) 
    if("POSIXct" %in% class(t1)) t1 = as.numeric(difftime(t1,origen,units='days'))
    if("POSIXct" %in% class(t2)) t2 = as.numeric(difftime(t2,origen,units='days'))  
    Xnorm = cbind(t1,Xnorm)      # first column of Xnorm is time
  } else use_time = FALSE
  if(missing(spatial)) use_space = FALSE
  if(missing(Xcat))    use_cats = FALSE
  else Xcat = apply(Xcat,2,function(x) as.integer(factor(x)))   # convert factors to integers
  
  M <- min(maxcriminals,length(crimeID))
  g<-crimeID
  maxg<-max(g,na.rm=T)
  miss<-is.na(g)
  g[miss]<-sample(1:M,sum(miss),replace=T)

  n<-length(g)
  p<-ncol(Xnorm)
  q<-ncol(Xcat)

  miss.norm<-is.na(Xnorm)
  miss.cat<-is.na(Xcat)
  Xcat[miss.cat]<-1
  if(p>1){for(j in 2:p){
    Xnorm[miss.norm[,j],j]<-mean(miss.norm[,j],na.rm=T)
  }}

  ncats<-apply(Xcat,2,max)
  maxcat<-max(ncats)
  pcat<-matrix(0,maxcat,q)
  for(j in 1:q){
    for(l in 1:maxcat){
     pcat[l,j]<-1+sum(Xcat[,j]==l)
    }
    pcat[,j]<-pcat[,j]/sum(pcat[,j])
  }

  mus<-cbind(runif(M,min(spatial[,1],na.rm=T),max(spatial[,1],na.rm=T)),
             runif(M,min(spatial[,2],na.rm=T),max(spatial[,2],na.rm=T)))
  mncat<-array(0,c(maxcat,q,M))
  for(l in 1:q){mncat[1:ncats[l],l,]<-1/ncats[l]}
  mu<-matrix(0,M,p)
  tau1<-tau2<-rep(0,p)
  for(l in 1:p){
    mu[,p]<-mean(Xnorm[,l],na.rm=T)
    tau1[l]<-tau2[l]<-1/var(Xnorm[,l],na.rm=T)
  }
  taus1<-taus2<-1/mean(var(spatial,na.rm=T))
  theta<-colMeans(mu)
  thetas<-colMeans(mus)
  df<-rep(1,q)
  probs<-rep(1,M)
  D<-100

  keep.df<-matrix(0,iters,q)
  keep.D<-rep(0,iters)
  keep.sd1<-keep.sd2<-keep.theta<-matrix(0,iters,p)
  keep.sds<-matrix(0,iters,2)
  missing<-(1:n)[miss]
  missing_s<-(1:n)[is.na(spatial[,1])]
  missing_t<-(1:n)[t1!=t2]
  n.missing_s<-length(missing_s)
  n.missing_t<-length(missing_t)
  keep.s<-keep.t<-NULL


  p.equal<-0

  if(n.missing_s>0){spatial[missing_s,]<-0}
  for(i in 1:iters){

    ############   Missing Spatial locations   ##########
    sss<-1/sqrt(taus1)
    if(use_space & n.missing_s>0){
       j<-missing_s
       spatial[j,1]<-rnorm(n.missing_s,mus[g[j],1],sss)
       spatial[j,2]<-rnorm(n.missing_s,mus[g[j],2],sss)
    }

    ############   Censored times   ##########
    if(use_time & n.missing_t>0){
       j<-missing_t
       new<-rtruncnorm(n.missing_t,mu[g[j],1],1/sqrt(tau1[1]),t1[j],t2[j])
       new[abs(new)==Inf]<-NA
       Xnorm[j,1]<-ifelse(is.na(new),Xnorm[j,1],new)
    }


    ############   Missing cat vars   ##########
    if(use_cats){for(l in 1:q){
      t<-1:ncats[l]
      for(j in 1:n){if(miss.cat[j,l]){
        Xcat[j,l]<-sample(t,1,prob=mncat[t,l,g[j]])
      }}
    }}

    ############   Missing cont vars   ##########
    if(p>1){for(l in 2:p){
      j<-which(miss.norm[,l])
      if(length(j)>0){
        Xnorm[j,l]<-rnorm(1,mu[g[j],l],1/sqrt(tau1[l]))
      }
    }}


    ############   MISSING GROUP LABELS   ##########
    oldg<-g
    for(j in missing){
       R<-log(probs)
       if(use_space){
         R<-R-0.5*taus1*((spatial[j,1]-mus[,1])^2+(spatial[j,2]-mus[,2])^2)
       }
       if(use_time){
         R<-R-0.5*tau1[1]*(Xnorm[j,1]-mu[,1])^2
       }
       ppp<-exp(R-max(R,na.rm=T))
       ppp[is.na(ppp)]<-0
       cang<-g[j]
       if(sum(ppp)>0){
         cang<-sample(1:M,1,prob=ppp)
       }
       R<-0
       if(use_cats){for(ccc in 1:q){
        R<-R+log(mncat[Xcat[j,q],1,cang])-
             log(mncat[Xcat[j,q],1,g[j]])
       }}
       if(runif(1)<exp(R)){g[j]<-cang}
    }

    ccc<-get.counts(g,M)
    probs<-rdirichlet(1,D+ccc)

    canD<-exp(rnorm(1,log(D),.05))
    R<-dnorm(log(canD),0,10,log=T)-
       dnorm(log(D),0,10,log=T)+
       ddir2(probs,canD)-
       ddir2(probs,D)
    if(!is.na(exp(R))){if(runif(1)<exp(R)){
        D<-canD
    }}


    ############   Update spatial model   ##########
    if(use_space){
     VVV<-taus1*count.g(g,M)+taus2
     MMM<-taus1*sum.by.g(spatial[,1],g,M)+taus2*thetas[1]
     mus[,1]<-rnorm(M,MMM/VVV,1/sqrt(VVV))
     MMM<-taus1*sum.by.g(spatial[,2],g,M)+taus2*thetas[2]
     mus[,2]<-rnorm(M,MMM/VVV,1/sqrt(VVV))

     VVV<-taus2*M+.0001
     MMM<-taus2*colSums(mus)+c(-76.6,39.3)*.0001
     thetas<-rnorm(2,MMM/VVV,1/sqrt(VVV))

     taus1<-rgamma(1,2*n/2+.1,sum((spatial-mus[g,])^2)/2+.1)
     SS<-sum((mus[,1]-thetas[1])^2+(mus[,2]-thetas[2])^2)
     taus2<-rgamma(1,2*M/2+.1,SS/2+.1)
    }

    ############   Continuous predictors   ##########
    if(use_time){for(l in 1:p){
      VVV<-tau1[l]*count.g(g,M)+tau2[l]
      MMM<-tau1[l]*sum.by.g(Xnorm[,l],g,M)+tau2[l]*theta[l]
      mu[,l]<-rnorm(M,MMM/VVV,1/sqrt(VVV))

      VVV<-tau2[l]*M+.0001
      MMM<-tau2[l]*sum(mu[,l])+0*.0001
      theta[l]<-rnorm(1,MMM/VVV,1/sqrt(VVV))

      tau1[l]<-rgamma(1,n/2+.1,sum((Xnorm[,l]-mu[g,l])^2)/2+.1)
      tau2[l]<-rgamma(1,M/2+.1,sum((mu[,l]-theta[l])^2)/2+.1)
   }}


   ############   Categorical predictors   ##########
   if(use_cats){
    eps<-0.0001
    for(l in 1:q){
      t<-1:ncats[l]
      for(j in 1:M){
        DDD<-df[l]*pcat[t,l]+
             get.counts(Xcat[g==j,l],ncats[l])
        duh<-rdirichlet(1,DDD[t])
        duh<-ifelse(duh>1-eps,1-eps,duh)
        duh<-ifelse(duh<eps,eps,duh)
        mncat[t,l,j]<-duh/sum(duh)
      }
      candf<-exp(rnorm(1,log(df[l]),0.1))
      R<-dnorm(log(candf),0,10,log=T)-
         dnorm(log(df[l]),0,10,log=T)
      for(j in 1:M){
        R<-R+ddir(mncat[t,l,j],candf,pcat[t,l])-
             ddir(mncat[t,l,j],df[l],pcat[t,l])
      }
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        df[l]<-candf}}
    }
   }


   if(i%%update==0){
     plot(g,main=paste("Iteration",i,"of",iters),
           xlab="Crime number",ylab="Criminal ID",
           pch = ifelse(miss,1,19),
           ylim=c(1,M),las=1)
   }

   keep.df[i,]<-df
   keep.sd1[i,]<-1/sqrt(tau1)
   keep.sd2[i,]<-1/sqrt(tau2)
   keep.sds[i,]<-1/sqrt(c(taus1,taus2))
   keep.D[i]<-D
   if(i>burn){
     #ddd<-round(rdist(g))
     ddd <- outer(g,g,'-') # remove dependence on fields
     p.equal<-p.equal+(ddd==0)/(iters-burn)
   }
  }

list(p.equal=p.equal,
   D=keep.D,
   df=keep.df,
   sd1=keep.sd1,
   sd2=keep.sd2,
   sds=keep.sds,
   theta=keep.theta,
   s.miss=keep.s,
   t.censored=keep.t,
   missing_s=missing_s,missing_t=missing_t,crimeID=crimeID)
}



  
## bayesPairs, bayesProb
##==============================================================================
#' Extracts the crimes with the largest probability of being linked. 
#'
#' Extracts the crimes (from \code{\link{crimeClust_bayes}}) with the largest 
#'  probability of being linked. 
##  Inputs:
#'  @param p.equal the posterior probability matrix produced by 
#'    \code{\link{crimeClust_bayes}}
#'  @param prob a column (or row) of the posterior probability matrix produced by 
#'    \code{\link{crimeClust_bayes}}
#'  @param drop only return crimes with a posterior linkage probability that 
#'   exceeds drop. Set to NA to return all results. 
#'  @details This is a helper function to easily extract the crimes with a high 
#'   probability of being linked from the output of \code{\link{crimeClust_bayes}}.
#'   \code{bayesPairs} searches the full posterior probability matrix and
#'   \code{bayesProb} only searches a particular column (or row).
##  Outputs:
#'  @return  data.frame of the indices of crimes with estimated posterior 
#'   probabilities, ordered from largest to smallest
#'  @seealso \code{\link{crimeClust_bayes}}
#'  @export
#' @name bayesPairs
NULL

#' @rdname bayesPairs
#' @export
bayesPairs <- function(p.equal,drop=0){
  pp = p.equal
  diag(pp) = NA
  low.tri = lower.tri(pp)
  prob = pp[low.tri]
  row = row(pp)[low.tri]
  col = col(pp)[low.tri]
  flip = (row>col)
  a = data.frame(i1=ifelse(flip,col,row),i2=ifelse(flip,row,col),prob)
  a = a[order(-a$prob),]
  if(is.numeric(drop)) a = subset(a,prob>drop)  
  return(a)
}

#' @rdname bayesPairs
#' @export
bayesProb <- function(prob,drop=0){
  ord = order(prob,decreasing=TRUE,na.last=NA)
  a = data.frame(index=ord,prob=prob[ord])
  if(is.numeric(drop)) a = subset(a,prob>drop)  
  return(a)
}







