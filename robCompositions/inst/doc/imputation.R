## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----load package, echo=FALSE, results='hide'----------------------------
library(robCompositions)
constSum <- function(x, const=1){
	x / rowSums(x) * const
}

## ----echo=FALSE, results='hide'------------------------------------------
##require(compositions)
genData <- function(n=1000, out=0.05, 
            Sigma=5*c(1,1)%*%t(c(1,1))+0.05*c(1,-1)%*%t(c(1,-1))){
    ## Gruppe ohne Ausreisser:
    z <- mvrnorm(n, mu=c(0,2), Sigma=Sigma)
    N <- dim(z)[1]
    n1 <- N - floor(n*out)
    n2 <- N - floor(2*n*out)
    if(out > 0){
      z[(n1+1):N, ] <- mvrnorm(floor(n*out), mu=c(0,6), Sigma=Sigma) ## erste Ausreissergruppe (Euclidean+Aitchison)
    }
    z <- isomLRinv(z) #ilr.inv(z)
    sum=runif(n1,0,1)  #rnorm(n1,10,1)
    z[1:n1, ] <- z[1:n1,] * sum
    if(out > 0){ 
      sum1=runif(floor(2*n*out),13,17) #rnorm(n2-n1,15,1)
      z[(n2+1):N, ] <- z[(n2+1):N, ] * sum1
      z[(n1+1):n2, ] <- z[(n1+1):n2, ] * 10
    }
    ## generate missings
    zmiss <- z
    s <- c(0.2, 0.1, 0.05, 0.05)
    for(i in 1:ncol(z)){
      zmiss[sample(1:n2, floor(s[i]*n2)), i] <- NA #1:index
    }
    list(zmiss=zmiss, z2=z, good=n2)
}

## ----echo=FALSE, results='hide'------------------------------------------
genData <- function(n=1000, out=0.05, 
            Sigma=1*c(1,1)%*%t(c(1,1))+0.05*c(1,-1)%*%t(c(1,-1))){
    ## Gruppe ohne Ausreisser:
    z <- mvrnorm(n, mu=c(0,0), Sigma=Sigma)
    N <- dim(z)[1]
    n1 <- N - floor(n*out)
    n2 <- N - floor(2*n*out)
    if(out > 0){
      z[(n1+1):N, ] <- mvrnorm(floor(n*out), mu=c(0,6), Sigma=Sigma) ## erste Ausreissergruppe (Euclidean+Aitchison)
    }
    z <- isomLRinv(z) #ilr.inv(z)
    sum=runif(n1,0,1)  #rnorm(n1,10,1)
    z[1:n1, ] <- z[1:n1,] * sum
    if(out > 0){ 
      sum1=runif(floor(2*n*out),13,17) #rnorm(n2-n1,15,1)
      z[(n2+1):N, ] <- z[(n2+1):N, ] * sum1
      z[(n1+1):n2, ] <- z[(n1+1):n2, ] * 10
    }
    ## generate missings
    zmiss <- z
    s <- c(0.2, 0.1, 0.05, 0.05)
    for(i in 1:ncol(z)){
      zmiss[sample(1:n2, floor(s[i]*n2)), i] <- NA #1:index
    }
    list(zmiss=zmiss, z2=z, good=n2)
}

## ----seed, echo=FALSE----------------------------------------------------
set.seed(123)
library(MASS)

## ----new data, echo=FALSE------------------------------------------------
x <- genData(100)

## ----plot.acomp, echo=FALSE----------------------------------------------
plot.acomp <- 
function (x, ..., labels = colnames(X), cn = colnames(X), aspanel = FALSE,
    id = FALSE, idlabs = NULL, idcol = 2, center = FALSE, scale = FALSE,
    pca = FALSE, col.pca = par("col"), margin = "acomp", add = FALSE,
    triangle = !add, col = par("col"),
    cexT=1.5, 
    adj=-1    
    )
{
    col <- unclass(col)
    X <- oneOrDataset(x)
    oX <- X
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
    if (NCOL(X) > 3) {
        if (margin == "rcomp")
            infkt <- function(x, y, ...) {
                plot.acomp(rcompmargin(X, d = c(gsi.mapfrom01(x),
                  gsi.mapfrom01(y)), pos = 1)[, c(3, 2, 1)],
                  ..., aspanel = TRUE, center = center, scale = scale,
                  col = col)
            }
        else if (margin == "acomp") {
            infkt <- function(x, y, ...) {
                plot.acomp(acompmargin(X, d = c(gsi.mapfrom01(x),
                  gsi.mapfrom01(y)), pos = 1)[, c(3, 2, 1)],
                  ..., aspanel = TRUE, center = center, scale = scale,
                  col = col)
            }
        }
        else {
            if (!is.numeric(margin))
                margin <- match(margin, colnames(X))
            fest <- X[, margin, drop = FALSE]
            X <- X[, -margin]
            infkt <- function(x, y, ...) {
                plot.acomp(acomp(cbind(X[, c(gsi.mapfrom01(y),
                  gsi.mapfrom01(x))], fest)), ..., aspanel = TRUE,
                  center = center, scale = scale, col = col)
            }
        }
        nn <- NCOL(X)
        if (add)
            gsi.add2pairs(sapply(1:NCOL(X), gsi.mapin01), infkt,
                ...)
        else gsi.pairs(sapply(1:NCOL(X), gsi.mapin01), labels = labels,
            panel = infkt, ...)
    }
    else {
        if (is.null(cn)) {
            cn <- c(expression(x[1]), expression(x[2]), expression(x[3]))
        }
        if (aspanel) {
            usr <- par("usr")
            on.exit(par(usr))
            par(usr = c(0, 1, 0, 1), pty = "s")
            lines(x = c(0, c60, 1, 0), y = c(0, s60, 0, 0))
            text(0, 0.2, cn[1], pos = 4, offset = 0.01, xpd = TRUE, cex=2)
            text(1, 0.2, cn[2], pos = 2, offset = 0.01, xpd = TRUE, cex=2)
            text(0.5, s60, cn[3], pos = 3, offset = 0.01, xpd = TRUE, cex=2)
        }
        else {
            if (!add) {
                usr <- par("pty")
                on.exit(par(usr))
                par(pty = "s")
                plot(x = c(0, c60, 1, 0), y = c(0, s60, 0, 0),
                  xlim = c(0, 1), ylim = c(0, 1), type = "n",
                  xlab = "", ylab = "", axes = FALSE)
                gsi.plots[[dev.cur()]] <<- NULL
            }
            if (triangle) {
                segments(x0 = c(0, 1, c60), y0 = c(0, 0, s60),
                  x1 = c(1, c60, 0), y1 = c(0, s60, 0))
                mtext(cn[1], side = 1, adj = 0, padj=adj, line = 1.5, cex=cexT)
                mtext(cn[2], side = 1, adj = 1, padj=adj , line = 1.5, cex=cexT)
                text(0.5, s60 * 1.03, cn[3], pos = 3, offset = 0.01,
                  xpd = TRUE, cex=cexT)
            }
        }
        X <- acomp(X, c(1, 2, 3))
        Y <- scale.acomp(X, center = center, scale = scale)
        gsi.setCoorInfo(mean = if (center)
            -mean(acomp(X))
        else acomp(c(1, 1, 1)), scale = if (scale)
            1/msd(X)
        else 1)
        x <- Y[, 2] + Y[, 3] * c60
        y <- Y[, 3] * s60
        points(x, y, ..., col = col)
    }
    return(invisible(NULL))
}

## ----knn-----------------------------------------------------------------
library(robCompositions) 
packageDescription("robCompositions")$Version
xImp <- impKNNa(x$zmiss, k=6)

## ----class---------------------------------------------------------------
class(xImp)

## ----printSummary--------------------------------------------------------
methods(class = "imp")
xImp

## ----imp-----------------------------------------------------------------
xImp1 <- impCoda(x$zmiss, method='lm')
xImp2 <- impCoda(x$zmiss, method='ltsReg')

## ----da, echo=FALSE------------------------------------------------------
cda <- function(xOrig, xImp, w){
  da <- function(x,y){
	d <- 0
	p <- length(x)
	for(i in 1:(p-1)){
		for(j in (i+1):p){
			d <- d + (log(x[i]/x[j]) - log(y[i]/y[j]))^2
		}
	}
	d=d/p
	sqrt(d)
  }
  das <- 0
  for(i in 1:nrow(xOrig)){
	das <- das + da(x=xOrig[i,], y=xImp[i,])
  }
  das/w
}

## ----variations, echo=FALSE----------------------------------------------
v1 <- variation(constSum(x$z2[1:95,]), robust=FALSE)
v2 <- variation(constSum(xImp1$xImp[1:95,]), robust=FALSE)
v22 <- variation(constSum(xImp2$xImp[1:95,]), robust=FALSE)
variations1 <- sum(abs(v1[upper.tri(v1, diag=FALSE)] - v2[upper.tri(v2, diag=FALSE)]), na.rm=TRUE)/length(c(upper.tri(v2, diag=FALSE))) 
variations2 <- sum(abs(v1[upper.tri(v1, diag=FALSE)] - v22[upper.tri(v22, diag=FALSE)]), na.rm=TRUE)/length(c(upper.tri(v22, diag=FALSE))) 

## ----erg-variations, echo=FALSE------------------------------------------
paste("RDA: iterative lm approach:", round(cda(x$z2, xImp1$xImp, xImp1$w),3))
paste("RDA: iterative ltsReg approach:", round(cda(x$z2, xImp2$xImp, xImp2$w), 3))
paste("DV: iterative lm approach:", round(variations1, 3))
paste("DV: iterative ltsReg approach:", round(variations2, 3))

## ----bootstrap-old, echo=FALSE-------------------------------------------
bootimp <- function(x, R = 100, method = "lm") {
     d <- dim(x)[2]
     n <- dim(x)[1]
     thetaM <- matrix(NA, ncol = 2, nrow = d)
     xs <- theta_hat1 <- matrix(0, nrow = n, ncol = d)
     med <- matrix(NA, ncol = d, nrow = R)
     for (i in 1:R) {
                ind <- sample(1:n, replace = TRUE)
                forbid <- c(91:100)
                induse <- ind[apply(outer(ind,forbid,"!=")==FALSE,1,sum)==0]
                s1 <- x[ind, ]
             simp <- impCoda(s1, method=method)$xImp
             med[i, ] <- apply(simp[induse,], 2, geometricmean)
     }
     thetaHat <- apply(med, 2, mean)
     for (i in 1:d) {
         thetaM[i, 1] <- quantile(med[, i], 0.025)
         thetaM[i, 2] <- quantile(med[, i], 0.975)
     }
     res <- list(geometricmean = thetaHat, ci = thetaM)
     res
}

## ----bootstrap-new, echo=FALSE-------------------------------------------
bootimp <- function(x, R = 100, method = "lm") {
     d <- dim(x)[2]
     n <- dim(x)[1]
     thetaM <- matrix(NA, ncol = 2, nrow = d)
     xs <- theta_hat1 <- matrix(0, nrow = n, ncol = d)
     med <- matrix(NA, ncol = d, nrow = R)
     for (i in 1:R) {
                ind <- sample(1:n, replace = TRUE)
                forbid <- c(91:100)
                induse <- ind[apply(outer(ind,forbid,"!=")==FALSE,1,sum)==0]
                s1 <- x[ind, ]
             simp <- impCoda(s1, method=method)$xImp
             med[i, ] <- apply(simp[induse,], 2, geometricmean)
     }
     invisible(as.data.frame(constSum(med)))
	 ##thetaHatMed <- apply(med, 2, mean)
     ##for (i in 1:d) {
     ##    thetaM[i, 1] <- quantile(med[, i], 0.025)
     ##    thetaM[i, 2] <- quantile(med[, i], 0.975)
     ##}
     ##res <- list(geometricmean = thetaHat, ci = thetaM)
     ##res
}
bootimpEM <- function(x, R = 100) {
	d <- dim(x)[2]
	n <- dim(x)[1]
	thetaM <- matrix(NA, ncol = 2, nrow = d)
	xs <- theta_hat1 <- matrix(0, nrow = n, ncol = d)
	med <- matrix(NA, ncol = d, nrow = R)
	for (i in 1:R) {
                ind <- sample(1:n, replace = TRUE)
                forbid <- c(91:100)
                induse <- ind[apply(outer(ind,forbid,"!=")==FALSE,1,sum)==0]
                s1 <- x[ind, ]
	## und das hier produziert Unsinn:
		s <- prelim.norm(s1) 
		thetahat <- em.norm(s, showits=FALSE)   
		simp <- imp.norm(s, thetahat, s)[[1]] 
		print(simp)
		med[i, ] <- apply(simp[induse,], 2, geometricmean)
	}
	invisible(as.data.frame(constSum(med)))
}
geometricmean <- function (x) {
	if (any(na.omit(x == 0)))
		0
	else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
}
bootimpGM <- function(x, R = 100) {
	d <- dim(x)[2]
	n <- dim(x)[1]
	thetaM <- matrix(NA, ncol = 2, nrow = d)
	xs <- theta_hat1 <- matrix(0, nrow = n, ncol = d)
	med <- matrix(NA, ncol = d, nrow = R)
	for (i in 1:R) {
                ind <- sample(1:n, replace = TRUE)
                forbid <- c(91:100)
                induse <- ind[apply(outer(ind,forbid,"!=")==FALSE,1,sum)==0]
                s1 <- x[ind, ]
		w <- is.na(s1)
		gm <- apply(s1, 2, function(x) {
					geometricmean(x[complete.cases(x)])
				})
		xmean <- x
		for(j in 1:ncol(x)){
			xmean[w[,j], j] <- gm[j]
		}
		med[i, ] <- apply(xmean[induse,], 2, geometricmean)
	}
	invisible(as.data.frame(constSum(med)))
}
bootimpM <- function(x, R = 100) {
        d <- dim(x)[2]
        n <- dim(x)[1]
        thetaM <- matrix(NA, ncol = 2, nrow = d)
        xs <- theta_hat1 <- matrix(0, nrow = n, ncol = d)
        med <- matrix(NA, ncol = d, nrow = R)
        for (i in 1:R) {
		ind <- sample(1:n, replace = TRUE)	
		forbid <- c(91:100)
		induse <- ind[apply(outer(ind,forbid,"!=")==FALSE,1,sum)==0]
                s1 <- x[ind, ]
                simp <- impute(s1, what="mean")
                med[i, ] <- apply(simp[induse,], 2, geometricmean)
        }
        invisible(as.data.frame(constSum(med)))
}


## ----bootstat-erg--------------------------------------------------------
R <- 5
bootimp(x$z2, R=R)

