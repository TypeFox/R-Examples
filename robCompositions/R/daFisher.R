#' Discriminant analysis by Fisher Rule.
#' 
#' Discriminant analysis by Fishers rule using CoDa methods.
#' 
#' The Fisher rule leads only to linear boundaries. However, this method allows
#' for dimension reduction and thus for a better visualization of the
#' separation boundaries. For the Fisher discriminant rule (Fisher, 1938; Rao,
#' 1948) the assumption of normal distribution of the groups is not explicitly
#' required, although the method looses its optimality in case of deviations
#' from normality.
#' 
#' The classical Fisher discriminant rule is invariant to ilr and clr
#' transformations. The robust rule is invariant to ilr transformations if
#' affine equivariant robust estimators of location and covariance are taken.
#' 
#' Robustification is done (method \dQuote{robust}) by estimating the
#' columnwise means and the covariance by the Minimum Covariance Estimator.
#' 
#' @aliases daFisher print.daFisher
#' @param x a matrix or data frame containing the explanatory variables
#' (training set)
#' @param grp grouping variable: a factor specifying the class for each
#' observation.
#' @param coda TRUE, when the underlying data are compositions.
#' @param method \dQuote{classical} or \dQuote{robust}
#' @param plotScore TRUE, if the scores should be plotted automatically.
#' @param ... additional arguments for the print method passed through
#' @importFrom e1071 matchClasses
#' @return an object of class \dQuote{daFisher} including the following
#' elements \item{B }{Between variance of the groups} \item{W }{Within variance
#' of the groups} \item{loadings}{loadings} \item{scores}{fisher scores} \item{mc}{table indicating misclassifications} \item{mcrate}{misclassification rate}  \item{coda}{coda}
#' @author Peter Filzmoser, Matthias Templ.
#' @seealso \code{\link[rrcov]{Linda}}
#' @references Filzmoser, P. and Hron, K. and Templ, M. (2012) 
#' Discriminant analysis for compositional data and robust parameter estimation. 
#' \emph{Computational Statistics}, Vol. 27(4), pp. 585-604, 2012.
#' 
#' Fisher, R. A. (1938) The statistical utiliziation of multiple measurements.
#' \emph{Annals of Eugenics}, 8:376-386.
#' 
#' Rao, C.R. (1948) The utilization of multiple measurements in problems of
#' biological classification. \emph{Journal of the Royal Statistical Society},
#' Series B, 10:159-203.
#' @keywords multivariate
#' @export
#' @import rrcov
#' @examples
#' ## toy data (non-compositional)
#' require(MASS)
#' x1 <- mvrnorm(20,c(0,0,0),diag(3))
#' x2 <- mvrnorm(30,c(3,0,0),diag(3))
#' x3 <- mvrnorm(40,c(0,3,0),diag(3))
#' X <- rbind(x1,x2,x3)
#' grp=c(rep(1,20),rep(2,30),rep(3,40))
#' 
#' #par(mfrow=c(1,2))
#' d1 <- daFisher(X,grp=grp,method="classical",coda=FALSE)
#' d2 <- daFisher(X,grp=grp,method="robust",coda=FALSE)
#' d2
#' summary(d2)
#' predict(d2)
#' 
#' ## example with olive data:
#'\dontrun{
#' data(olives, package = "classifly")
#' # exclude zeros (alternatively impute them if 
#' # the detection limit is known using impRZilr())
#' ind <- which(olives==0, arr.ind = TRUE)[,1]
#' olives <- olives[-ind, ]
#' x <- olives[, 4:10]
#' grp <- olives$Region # 3 groups
#' res <- daFisher(x,grp)
#' res
#' summary(res)
#' predict(res)
#' res <- daFisher(x, grp, plotScore = TRUE)
#' res <- daFisher(x, grp, method = "robust")
#' res
#' summary(res)
#' predict(res)
#' res <- daFisher(x,grp, plotScore = TRUE, method = "robust")
#' 
#' # 9 regions
#' grp <- olives$Area
#' res <- daFisher(x, grp, plotScore = TRUE)
#' res
#' summary(res)
#' predict(res)
#' }
daFisher <- function(x, grp, coda=TRUE,
                   method = "classical", 
                   plotScore = FALSE){
  ## some checks
  if(class(x) == "data.frame") x <- as.matrix(x)
  ## Fisher LDA:
  if(length(grp) != dim(x)[1]){
    stop(paste("grp must be of length", dim(x)[1]))
  }
  if(dim(x)[2] < 1){
    stop("matrix or data.frame expected.")
  }
  if(coda){
    x <- isomLR(x)
  }
  n <- nrow(x)
  p <- ncol(x)
  glev <- unique(grp)
  g <- length(glev)
  pj <- rep(NA,g)
  meanj <- matrix(NA,nrow=p,ncol=g)
  cv <- list()
  for (j in 1:g){
    pj[j] <- sum(grp==glev[j])/n
    if(method == "classical"){
      meanj[,j] <- apply(x[grp==glev[j],],2,mean)
      cv[[j]] <- cov(x[grp==glev[j],])
    } else {
      robcov <- covMcd(x[grp==glev[j],])
      meanj[,j] <- robcov$center
      cv[[j]] <- robcov$cov
      #   else {
      #   #  require(rrcov)
      #     res <- by(x,factor(grp),CovMcd)
      #     muil <- lapply(res,getCenter)
      #     sigil <- lapply(res,getCov)
      #   }
    }
  }
  meanov <- t(t(meanj)*pj)
  B <- matrix(0,p,p)
  W <- matrix(0,p,p)
  for (j in 1:g){
    B <- B+pj[j]*((meanj-meanov)%*%t(meanj-meanov))
#    W <- W+pj[j]*cov(x[grp==glev[j],])
    W <- W+pj[j]*cv[[j]]
  }
  l <- min(g-1,p) # use this number of components
  #V=matrix(Re(eigen(solve(W)%*%B)$vec)[,1:l],ncol=l)
  #V=t(t(V)/(sqrt(diag(t(V)%*%W%*%V))))
  
  # besser:
  B.svd <- svd(B)
  B12 <- B.svd$u[,1:l]%*%diag(sqrt(B.svd$d[1:l]))%*%t(B.svd$u[,1:l])
  Bm12 <- B.svd$u[,1:l]%*%diag(1/sqrt(B.svd$d[1:l]))%*%t(B.svd$u[,1:l])
  K <- eigen(B12%*%solve(W)%*%B12)
  Vs <- Bm12%*%K$vec[,1:l]
  V <- t(t(Vs)/(sqrt(diag(t(Vs)%*%W%*%Vs))))
  
  
  # Fisher scores
  fs=matrix(NA,nrow=n,ncol=g)
  for (j in 1:g){
    xc <- scale(x,meanj[,j],scale=FALSE)
    xproj <- xc%*%V
    fs[,j] <- sqrt(apply(xproj^2,1,sum)-2*log(pj[j]))
  }
  
  ## predition:
  grppred <- apply(fs, 1, which.min)
  
  ## misclassification rate:
  mc <- table(grp, grppred)
  mc <- mc[, matchClasses(mc, method = "exact")]
  rate <- 1 - sum(diag(mc)) / sum(mc)
  
  ## plot scores (if TRUE)
  if(plotScore){
    proj <- xc %*%V [,1:2]
    proj <- data.frame(proj)
    proj$grp <- as.factor(grp)
    proj$grppred <- as.factor(grppred)
    firstscores <- NULL
    secondscores <- NULL
    colnames(proj) <- c("firstscores", "secondscores","grp", "grppred")
    gg <- ggplot(proj, aes(firstscores, secondscores, colour = grp, shape = grppred)) 
    gg <- gg + geom_point()
    gg <- gg + xlab("first fisher scores") + ylab("second fisher scores")
    print(gg)
#    plot(, col=grp, pch=grppred, 
#         xlab="first fisher scores", ylab="second fisher scores")
  }
  
  res <- list(B = B, 
              W = W,
              loadings = V,
              scores = fs,#classification=postgroup, 
            #  mu=muil, 
            #  sigma=sigil,
              mc = mc,
              mcrate =  rate,
              coda=coda)
  class(res) <- "daFisher"
  res
}

# daFisher <- function(x,grp,coda=TRUE,method="classical",plotScore=FALSE)
# {
#   if(class(x)=="data.frame") x <- as.matrix(x)
#   # Fisher LDA:
#   if(length(grp) != dim(x)[1]){
#   	stop(paste("grp must be of length",dim(x)[1]))
#   }
#   if(dim(x)[2] < 1){
#   	stop("matrix or data.frame expected.")
#   }
#   if(coda){
#   	x <- isomLR(x)
#   }
#   
#   	p <- ncol(x)
#   	ni <- table(grp)
#   	ng <- length(ni)
#   	n <- sum(ni)
#   	pi <- ni/n
#   if (method=="classical"){
#     muil <- by(x,factor(grp),colMeans)
#     sigil <- by(x,factor(grp),cov)
#   }
#   else {
#   #  require(rrcov)
#     res <- by(x,factor(grp),CovMcd)
#     muil <- lapply(res,getCenter)
#     sigil <- lapply(res,getCov)
#   }
#   
#   	mui <- matrix(unlist(muil),ng,p,byrow=TRUE)
#   	mu <- pi%*%mui
#   	hlp <- diag(sqrt(pi))%*%(mui-rep(1,ng)%*%mu)
#   	B <- t(hlp)%*%hlp
#   	sigi <- array(unlist(sigil),dim=c(p,p,ng))
#   	W <- apply(sigi*array(sort(rep(pi,p*p)),dim=c(p,p,ng)),c(1,2),sum)
#   	adir <- matrix(as.numeric(eigen(solve(W)%*%B)$vec),ncol=p)
#   	adirs <- t(t(adir)/(sqrt(diag(t(adir)%*%W%*%adir))))
#   	scores=x%*%adirs
#   if(plotScore){
#     pl <- as.numeric(factor(grp))
#     plot(scores[,1:2],col=pl, pch=pl, cex=1.5, xlab="Scores 1", ylab="Scores 2", cex.lab=1.2)
#     legend("topright", legend=levels(factor(grp)), pch=unique(pl), col=unique(pl), cex=1.3)
#   }
#   #    postgroup <- apply(scores, 1, which.min)
#   #	print(postgroup)
#   	res <- list(B=B,W=W,loadings=adir,scores=scores,#classification=postgroup, 
#   			mu=muil, sigma=sigil,
#   			coda=coda)
#   	class(res) <- "daFisher"
#   	
#   #	## fill in for class lda
#   #	g <- as.factor(grp)
#   #	lev <- lev1 <- levels(g)
#   #	counts <- as.vector(table(g))
#   #	prior <- counts/n
#   #	prior <- prior[counts > 0]
#   #
#   #if(method == "moment") fac <- 1/(n-ng) else fac <- 1/n
#   #X <- sqrt(fac) * (x - group.means[g,  ]) %*% scaling
#   #X.s <- svd(X, nu = 0)
#   #X <-  sqrt(nu/(nu-2)*(1 + p/nu)/n * w) * (x - group.means[g,  ]) %*% scaling
#   #X.s <- svd(X, nu = 0)
#   #	cl <- match.call()
#   #	cl[[1L]] <- as.name("daFisher")
#   #	
#   #	res <- structure(list(prior = prior, counts = counts, means = mui,
#   #					scaling = hlp, lev = lev, svd = hlp,
#   #					N = n, call = cl, B=B, W=W, loadings=adir, coda=coda),
#   #			class = "lda")
#   #
#   #	z1=z[grp=="arabica",]
#   #	z2=z[grp=="blended",]
#   #	n1=nrow(z1)
#   #	n2=nrow(z2)
#   #	n=n1+n2
#   #	p1=n1/n
#   #	p2=n2/n
#   #	m1=apply(z1,2,mean)
#   #	m2=apply(z2,2,mean)
#   #	S1=cov(z1)
#   #	S2=cov(z2)
#   #	Sp=((n1-1)/(n1-1+n2-1))*S1+((n2-1)/(n1-1+n2-1))*S2
#   #	Sp1=solve(Sp)
#   #	yLDA=as.numeric(t(m1-m2)%*%Sp1%*%t(z)-as.numeric(1/2*t(m1-m2)%*%Sp1%*%(m1+m2)))-log(p2/p1)	
#   #	plot(z, pch=21, bg=ifelse(grp=="arabica","red","blue"))#bg=ifelse(yLDA<0,"red","blue"))
#   #	y1=seq(from=min(z[,1])-1.5,to=max(z[,1])+1.9,by=0.05)
#   #	y2=seq(from=min(z[,2]),to=max(z[,2])+0.2,by=0.05)
#   #	y1a=rep(y1,length(y2))
#   #	y2a=sort(rep(y2,length(y1)))
#   #	ya=cbind(y1a,y2a)
#   #	yaLDA=as.numeric(t(m1-m2)%*%Sp1%*%t(ya)-
#   #					as.numeric(1/2*t(m1-m2)%*%Sp1%*%(m1+m2)))-log(p2/p1)
#   #	
#   #	boundLDA=abs(yaLDA)<0.05
#   #	lines(lowess(y1a[boundLDA],y2a[boundLDA]),col=gray(0.6),lwd=1.5,lty=1)
#   
#   invisible(res)
# }


#' @rdname daFisher
#' @method print daFisher
#' @export
print.daFisher <- function(x,...){
  cat("--------------------------------------")
  cat("\nResults from Fishers discriminant analysis, coda ==", x$coda)
  cat("\n- Variance between the classes: \n")
  print(x$B)
  cat("\n- Variance within the classes: \n")
  print(x$W)
  cat("\n- Loadings matrix: \n")
  print(x$load)
  cat("--------------------------------------\n")
}

#' @rdname daFisher
#' @method predict daFisher
#' @param object object of class \dQuote{daFisher}
#' @export
predict.daFisher <- function(object, ...){
  grppred <- apply(object$scores, 1, which.min)
  return(grppred)
}

#' @rdname daFisher
#' @method summary daFisher
#' @export
summary.daFisher <- function(object, ...){
  cat("--------------------------------------")
  cat("\nMisclassification rate from Fishers discriminant analysis, coda ==", object$coda)
  cat("\n")
  print(object$mcrate)
  cat("\n--------------------------------------")
  cat("\nMisclassifications from Fishers discriminant analysis, coda ==", object$coda)
  cat("\n")
  print(object$mc)
  cat("\n--------------------------------------\n")
}

