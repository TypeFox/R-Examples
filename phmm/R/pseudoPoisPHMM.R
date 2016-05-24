#' Pseudo poisson data for fitting PHMM via GLMM
#' 
#' Function for generating a pseudo Poisson data set which can be used to fit a
#' PHMM using GLMM software. This follows the mixed-model extension Whitehead
#' (1980), who described how to fit Cox (fixed effects) models with GLM
#' software.
#' 
#' 
#' @param x an object of class \code{phmm}.
#' @return A \code{data.frame} with columns:
#' @returnItem time the event time;
#' @returnItem N the number at risk at time \code{time};
#' @returnItem m the number at risk (in the same cluster with same covariates)
#' at time \code{time};
#' @returnItem cluster the integer cluster indicator;
#' @returnItem N the number at risk at time \code{time};
#' @returnItem fixedeffectscovariates denoted \code{z1}, \code{z2}, etc.;
#' @returnItem randomeffectscovariates denoted \code{w1}, \code{w2}, etc.;
#' @returnItem linear.predictors the linear predictors from the \code{phmm} fit
#' (excluding the cumulative hazard estimates.
#' @seealso \code{\link{phmm}}, \code{\link{traceHat}}
#' @references Whitehead, J. (1980). Fitting Cox's Regression Model to Survival
#' Data using GLIM. Journal of the Royal Statistical Society. Series C, Applied
#' statistics, 29(3). 268-.
#' @keywords survival
#' @examples
#' \dontrun{
#' n <- 50      # total sample size
#' nclust <- 5  # number of clusters
#' clusters <- rep(1:nclust,each=n/nclust)
#' beta0 <- c(1,2)
#' set.seed(13)
#' #generate phmm data set
#' Z <- cbind(Z1=sample(0:1,n,replace=TRUE),
#'            Z2=sample(0:1,n,replace=TRUE),
#'            Z3=sample(0:1,n,replace=TRUE))
#' b <- cbind(rep(rnorm(nclust),each=n/nclust),rep(rnorm(nclust),each=n/nclust))
#' Wb <- matrix(0,n,2)
#' for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
#' Wb <- apply(Wb,1,sum)
#' T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
#' C <- runif(n,0,1)
#' time <- ifelse(T<C,T,C)
#' event <- ifelse(T<=C,1,0)
#' mean(event)
#' phmmd <- data.frame(Z)
#' phmmd$cluster <- clusters
#' phmmd$time <- time
#' phmmd$event <- event
#' 
#' fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (-1 + Z1 + Z2 | cluster), 
#'    phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
#'    NINIT = 10, MAXSTEP = 100, CONVERG=90)
#' 
#' # Same data can be fit with lmer,
#' # though the correlation structures are different.
#' poisphmmd <- pseudoPoisPHMM(fit.phmm)
#' 
#' library(lme4)
#' fit.lmer <- lmer(m~-1+as.factor(time)+z1+z2+
#'   (-1+w1+w2|cluster)+offset(log(N)), 
#'   as.data.frame(as(poisphmmd, "matrix")), family=poisson)
#' 
#' fixef(fit.lmer)[c("z1","z2")]
#' fit.phmm$coef
#' 
#' VarCorr(fit.lmer)$cluster
#' fit.phmm$Sigma
#' 
#' logLik(fit.lmer)
#' fit.phmm$loglik
#' 
#' traceHat(fit.phmm)
#' }
pseudoPoisPHMM <- function (x) UseMethod("pseudoPoisPHMM")
pseudoPoisPHMM.phmm <- function(x){	
	dd <- cBind(x$cluster, x$Z, x$W)
	group <- apply(dd,1,paste,collapse="XX")
	groups <- unique(group)
	dd <- cBind(x$Y[, 1], x$Y[, 2], x$linear.predictors, dd)
	colnames(dd) <- c("time", "delta", "linear.predictors",
		"cluster",
		paste("z", 1:x$nfixed, sep=''),
		paste("w", 1:x$nrandom, sep=''))
	ddext <- lapply(sort(unique(dd[dd[,"delta"]==1,"time"])), function(t){
		tdd <- lapply(1:length(groups), function(i){
			unlist(c(t,sum(dd[group==groups[i],"time"]>=t),
			sum(dd[group==groups[i]&dd[,"delta"]==1,"time"]==t),
			dd[group==groups[i]&!duplicated(group==groups[i]),
				c("cluster",
				paste("z", 1:x$nfixed, sep=''),
				paste("w", 1:x$nrandom, sep=''),
				"linear.predictors")]))
			})
		})
	cnames <- c('time','N','m',"cluster",
				paste("z", 1:x$nfixed, sep=''),
				paste("w", 1:x$nrandom, sep=''),
				"linear.predictors")
	ddext <- Matrix(unlist(ddext), byrow=TRUE, ncol=length(cnames), sparse=TRUE)
	colnames(ddext) <- cnames
	times <- sort(unique(dd[dd[,"delta"]==1,"time"]))
	timematrix <- Matrix(0,nrow(ddext),ncol=length(times))
	colnames(timematrix) <- paste("t",1:length(times),sep='')
	for(i in 1:length(times) ){		
		timematrix[ddext[,"time"]==times[i],paste("t",i,sep='')] <- 1
	}
	ddext <- cBind(ddext, timematrix)
	ddext <- ddext[ddext[,'N']!=0,]
	ddext <- ddext[order(ddext[,'cluster']),]
	return(ddext)
}

pseudoPoisPHMM.coxph <- function(x){	
	group <- apply(x$x,1,paste,collapse="XX")
	groups <- unique(group)
	dd <- cBind(x$y[, 1], x$y[, 2], x$x%*%x$coef, x$x)
	colnames(dd) <- c("time", "delta", "linear.predictors",
		paste("z", 1:ncol(x$x), sep=''))
	ddext <- lapply(sort(unique(dd[dd[,"delta"]==1,"time"])), function(t){
		tdd <- lapply(1:length(groups), function(i){
			unlist(c(t,sum(dd[group==groups[i],"time"]>=t),
			sum(dd[group==groups[i]&dd[,"delta"]==1,"time"]==t),
			dd[group==groups[i]&!duplicated(group==groups[i]),
				c(paste("z", 1:ncol(x$x), sep=''),
				  "linear.predictors")]))
			})
		})
	cnames <- c('time','N','m',
				paste("z", 1:ncol(x$x), sep=''),
				"linear.predictors")
	ddext <- Matrix(unlist(ddext), byrow=TRUE, ncol=length(cnames), sparse=TRUE)
	colnames(ddext) <- cnames
	times <- sort(unique(dd[dd[,"delta"]==1,"time"]))
	timematrix <- Matrix(0,nrow(ddext),ncol=length(times))
	colnames(timematrix) <- paste("t",1:length(times),sep='')
	for(i in 1:length(times) ){		
		timematrix[ddext[,"time"]==times[i],paste("t",i,sep='')] <- 1
	}
	ddext <- cBind(ddext, timematrix)
	ddext <- ddext[ddext[,'N']!=0,]
	# bh <- basehaz(x, centered = FALSE)
	# lambda <- bh$hazard - c(0,bh$hazard[1:(length(bh$hazard)-1)])
	# alpha <- log(ddext[,paste("t",1:length(times),sep='')]%*%lambda)
	# ddext <- cBind(ddext, alpha)
	return(ddext)
}