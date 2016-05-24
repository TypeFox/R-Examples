#' MHalgoGen is a function to use mcmc with Metropolis-Hastings algorithm
#' @title Monte-Carlo Markov-chain with Metropolis-Hastings algorithm
#' @author Marc Girondot
#' @return A mcmcComposite object with all characteristics of the model and mcmc run
#' @param n.iter Number of iterations for each chain
#' @param parameters A data.frame with priors; see description and examples
#' @param likelihood The function that returns -ln likelihood using data and parameters
#' @param n.chains Number of chains
#' @param n.adapt Number of iteration to stabilize likelihood
#' @param thin Interval for thinning likelihoods
#' @param trace Or FALSE or period to show progress
#' @param intermediate Or NULL of period to save intermediate result
#' @param filename Name of file in which intermediate results are saved
#' @param previous The content of the file in which intermediate results are saved
#' @param ... Parameters to be transmitted to likelihood function
#' @description The parameters must be stored in a data.frame with named rows for each parameter with the following columns:\cr
#' \itemize{
#'   \item Density. The density function name, example \code{dnorm}, \code{dlnorm}, \code{dunif}
#'   \item Prior1. The first parameter to send to the \code{Density} function
#'   \item Prior2. The second parameter to send to the \code{Density} function
#'   \item SDProp. The standard error from new proposition value of this parameter
#'   \item Min. The minimum value for this parameter
#'   \item Max. The maximum value for this parameter
#'   \item Init. The initial value for this parameter
#' }
#' This script has been deeply modified from a MCMC script provided by Olivier Martin (INRA, Paris-Grignon).
#' @family mcmcComposite functions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' require(coda)
#' x <- rnorm(30, 10, 2)
#' dnormx <- function(x, par) return(-sum(dnorm(x, mean=par['mean'], sd=par['sd'], log=TRUE)))
#' parameters_mcmc <- data.frame(Density=c('dnorm', 'dlnorm'), 
#' Prior1=c(10, 0.5), Prior2=c(2, 0.5), SDProp=c(0.35, 0.2), 
#' Min=c(-3, 0), Max=c(100, 10), Init=c(10, 2), stringsAsFactors = FALSE, 
#' row.names=c('mean', 'sd'))
#' mcmc_run <- MHalgoGen(n.iter=10000, parameters=parameters_mcmc, x=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' plot(mcmc_run, xlim=c(0, 20))
#' plot(mcmc_run, xlim=c(0, 10), parameters="sd")
#' mcmcforcoda <- as.mcmc(mcmc_run)
#' #' heidel.diag(mcmcforcoda)
#' raftery.diag(mcmcforcoda)
#' autocorr.diag(mcmcforcoda)
#' acf(mcmcforcoda[[1]][,"mean"], lag.max=20, bty="n", las=1)
#' acf(mcmcforcoda[[1]][,"sd"], lag.max=20, bty="n", las=1)
#' batchSE(mcmcforcoda, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmcforcoda)$statistics[,"Time-series SE"]
#' summary(mcmc_run)
#' as.parameters(mcmc_run)
#' lastp <- as.parameters(mcmc_run, index="last")
#' parameters_mcmc[,"Init"] <- lastp
#' # The n.adapt set to 1 is used to not record the first set of parameters
#' # then it is not duplicated (as it is also the last one for 
#' # the object mcmc_run)
#' mcmc_run2 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, x=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=1, thin=1, trace=1)
#' mcmc_run3 <- merge(mcmc_run, mcmc_run2)
#' ####### no adaptation, n.adapt must be 0
#' parameters_mcmc[,"Init"] <- c(mean(x), sd(x))
#' mcmc_run3 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, x=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=0, thin=1, trace=1)
#' }
#' @export

# Algo Metropolis-Hastings
# ------------------------

MHalgoGen<-function(likelihood=stop("A likelihood function must be supplied"), 
                    parameters=stop("Priors  must be supplied"), ..., 
                    n.iter=10000, n.chains = 1, n.adapt = 100, thin=30, trace=FALSE, 
                    intermediate=NULL, filename="intermediate.Rdata",
                    previous=NULL)

{
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("coda package is necessary for this function")
  }
  
# likelihood=NULL; parameters=NULL; n.iter=10000; n.chains = 1; n.adapt = 100; thin=30; trace=FALSE; intermediate=NULL; filename="intermediate.Rdata"; previous=NULL
# datax <- list()
  
  if (is.null(previous)) {
    nbvar <- dim(parameters)[1]
    deb_kk <- 1
    deb_i <- 2
    res <- as.list(NULL)
    deb_varp2 <- matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
    deb_varp<-matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
    colnames(deb_varp)<-c(rownames(parameters), "Ln L")
    colnames(deb_varp2)<-c(rownames(parameters), "Ln L")
    t <- as.character(trace)
    cpt_trace <- 0
    res<-as.list(NULL)
    resL<-as.list(NULL)
    datax <- list(...)
  } else {
    n.iter <- previous$n.iter
    parameters <- previous$parameters
    # Initialisation
    nbvar<-dim(parameters)[1]
    datax <- previous$datax
    likelihood <- previous$likelihood
    n.chains <- previous$n.chains
    n.adapt <- previous$n.adapt
    thin <- previous$thin
    trace <- previous$trace
    intermediate <- previous$intermediate
    filename <- previous$filename
    deb_kk <- previous$chain
    deb_i <- previous$iter
    res <- previous$res
    deb_varp2 <- previous$varp2
    deb_varp <- previous$varp
    cpt <- previous$cpt
    cpt_trace <- previous$cpt_trace
    resL <- previous$resL
    sdg <- previous$sdg
    MaxL <- previous$MaxL
    t <- as.character(previous$trace)
  }
  
  pt <- NULL
  if (t=="TRUE") {pt <- 1;tf <- TRUE}
  if (t=="FALSE") {pt <- 0;tf <- FALSE}
  if (is.null(pt)) {
    tf <- TRUE
    pt <- floor((n.adapt+n.iter)/trace)
  }
  

for (kk in deb_kk:n.chains) {



varp <- deb_varp
varp2 <- deb_varp2

if (is.null(previous)) {

deb_varp2 <- matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
deb_varp <- matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
colnames(deb_varp)<-c(rownames(parameters), "Ln L")
colnames(deb_varp2)<-c(rownames(parameters), "Ln L")


varp[1, 1:nbvar] <- as.numeric(parameters[1:nbvar, 'Init'])

varp[1, "Ln L"]<- -do.call(likelihood, list(data=datax, x=varp[1, 1:nbvar]))
cpt<-1
varp2[cpt, 1:nbvar] <- varp[1, 1:nbvar]
varp2[cpt, "Ln L"] <- varp[1, "Ln L"]
cpt <- 2


if (trace) {
	cat(paste("Chain ", kk, ": [", 1, "] ",as.numeric(varp[1, nbvar+1]), "\n", sep=""))
} else {
	cat(paste("Chain ", kk, "\n", sep=""))
}

MaxL<-varp[1,]

sdg=NULL


# 18/1/2013
# for(i in 1:nbvar) sdg<-c(sdg, as.numeric(parameters[i, 'SDProp']))
sdg<-c(sdg, as.numeric(parameters[1:nbvar, 'SDProp']))

# previous <- NULL
}


if (is.data.frame(parameters[,2:3])) {
  Prior <- as.matrix(parameters[,2:3])
  Limites <- as.matrix(parameters[,5:6])
} else {
  Prior<-matrix(as.numeric(parameters[,2:3]), ncol=2)
  Limites<-matrix(as.numeric(parameters[,5:6]), ncol=2)
}

dfun<-parameters[,"Density"]


# Itérations
for (i in deb_i:(n.adapt+n.iter)) {
  
  # est-ce que je sauve où j'en suis
  if (!is.null(intermediate))
    if (i %% intermediate==0) {
      itr <- list(chain=kk, iter=i, varp2=varp2, res=res, 
                  n.iter=n.iter,
                  parameters=parameters,
                  datax=datax,
                  likelihood=likelihood,
                  n.chains=n.chains,
                  n.adapt=n.adapt,
                  thin=thin,
                  trace=trace,
                  intermediate=intermediate,
                  filename=filename, 
                  varp=varp, 
                  cpt=cpt,
                  cpt_trace=cpt_trace, 
                  resL=resL, 
                  sdg=sdg, MaxL=MaxL)
      save(itr, file=filename)
    }
  
  
  ###### Pas clair si previous
    deb_i <- 2
	newvarp<-varp[i-1, 1:nbvar]
	for (j in 1:nbvar) {	
		# Nouveau paramètre
		propvarp<-newvarp
		propvarp[j]<-propvarp[j]+rnorm(1, mean=0, sd=sdg[j])

		if (propvarp[j]<=Limites[j,2] && propvarp[j]>=Limites[j,1]) 
			{
			logratio <- get(dfun[j])(propvarp[j],Prior[j,1],Prior[j,2],log=T)+
					-do.call(likelihood, list(data=datax, x=propvarp))-
		        	(get(dfun[j])(newvarp[j],Prior[j,1],Prior[j,2],log=T)+varp[i-1, "Ln L"])
			alpha<-min(c(1,exp(logratio)))
			# 15/2/2015 Pour éviter des erreurs
			if (!is.finite(alpha)) alpha <- -1
			if (runif(1, min=0, max=1)<=alpha) {newvarp<-propvarp} 
			}
	}		
	varp[i, 1:nbvar]<-newvarp
	varp[i, "Ln L"] <- -do.call(likelihood, list(data=datax, x=newvarp))

	if (MaxL["Ln L"]<varp[i, "Ln L"]) {MaxL<-varp[i,]}
  

# 6/10/2012: Je stocke tout	
#	if ((i>n.adapt) && ((i%%thin)==0)) {
		varp2[cpt, 1:nbvar]<-newvarp
		varp2[cpt, "Ln L"]<-varp[i, "Ln L"]
		cpt<-cpt+1
#	}
	if (tf) {
    cpt_trace <- cpt_trace+1
    if (cpt_trace>=pt) {
	    cat(paste("Chain ", kk, ": [", i, "] ", as.numeric(varp[i, "Ln L"]), "\n", sep=""))
      cpt_trace <- 0
    }
	}
		
		
}

lp <- getFromNamespace("as.mcmc", ns="coda")(varp2[(n.adapt+1):(cpt-1), 1:nbvar])
lp <- getFromNamespace("mcmc", ns="coda")(data=lp, start=n.adapt+1, end=n.iter+n.adapt, thin=thin)

res<-c(res, list(lp))
resL <- c(resL, list(varp2[(n.adapt+1):(cpt-1), "Ln L"]))

}

names(res) <- 1:n.chains

res <- getFromNamespace("as.mcmc.list", ns="coda")(res)

cat("Best likelihood for: \n")
for (j in 1:nbvar) {cat(paste(names(MaxL[j]), "=", MaxL[j], "\n"))}



out <- (list(resultMCMC=res, resultLnL=resL, parametersMCMC=list(parameters=parameters, n.iter=n.iter, n.chains=n.chains, n.adapt=n.adapt, thin=thin)))
class(out) <- "mcmcComposite"
return(out)

}
