#' Produce quantile summaries of model posterior samples
#'
#' Inputs MCMC samples for model parameters and constructs c(2.5\%,50\%,97.5\%) quantile summaries.
#'
#' @param model.output  A list vector of objects returned by MCMC sampling functions.  e.g. mmCplusDpPost for \code{option = "mmcar"}.
#' @param dosemat An \code{P x (T+1)} \code{matrix} object that maps \code{subjects} to treatment dosages.  The first column should be an
#'			intercept column (filled with 1's).
#' @param Nfixed Number of total fixed effects, both time-based and nuisance.
#' @param Nrandom Number of total random effects, both time-based and nuisance, all grouped by subject.
#' @param Nsubject Number of unique subjects (on which repeated measures are observed).
#' @param typet A numeric vector of length equal to the number of treatments that contains the base distribution for each treatment.
#'		\code{1 = "car"}, \code{2 = "mvn"}, \code{3 = "ind"}.
#' @return A list object containing quantile summaries for all sampled model parameters.
#'     \item{deviance.summary}{vector of length 3 summarizing quantiles for model deviance.}
#'     \item{beta.summary}{\code{Nfixed x 3} quantile summaries of model fixed effects.}
#'     \item{alpha.summary}{quantile summary of model global intercept parameter.}
#'     \item{theta.summary}{list object of length \code{Nrandom}, each cell containing a \code{n x 3} matrix of by-subject random effect parameter quantile summaries.}
#'     \item{lambda.summary}{\code{Nrandom^2 x 3} quantile summaries of by-polynomial order precision parameters used in base distributions.}
#'     \item{lambda.mean}{\code{Nrandom x Nrandom} posterior means of by-polynomial order precision parameters used in base distributions.}
#'     \item{alphacar.summary}{\code{numcar x 3} quantile summaries of proper CAR strength of correlation parameters for CAR base distribution on subject-dose random effects, where numcar <= nty treatments.}
#'     \item{taucar.summary}{\code{numcar x 3} quantile summaries of proper CAR precision parameters for CAR base distribution on subject-dose random effects, where numcar <= nty treatments.}
#'     \item{dosetrt.summary}{list object of length \code{nty}, each cell containing a \code{Nsubject*(Nrandom*numt[m]) x 3} matrix of quantile summaries for subject-dose random effects.}
#'     \item{dosetrt.mean}{list object of length \code{nty}, each cell containing a \code{Nsubject x (Nrandom*numt[m])} matrix of posterior mean values for subject-dose random effects.}
#'     \item{pind.summary}{list object of length \code{numind}, each cell contains a \code{numt[m] x 3} matrix of quantile summaries for precision values under IND base distribution, where numind <= nty treatments.}
#'     \item{pmvn.summary}{list object of length \code{nummvn}, each cell containing quantile summaries for precision parameters used for MVN base distribution on subject-dose random effects, where nummvn <= nty treatments.}
#'     \item{pmvn.summary}{list object of length \code{nummvn}, each cell containing posterior mean for \code{numt[m] x numt[m]} matrix of precision parameters used for MVN base distribution on subject-dose random effects, where nummvn <= nty treatments.}
#'     \item{doseint.summary}{list object of length \code{Nrandom}, each cell containing a \code{Nsubject x 3} matrix of quantile summaries for the intercept parameter in the subject-dose random effects.}
#'     \item{taue.summary}{quantile summary for model error precision parameter.}
#'     \item{M.summary}{quantile summary for number of DP posterior clusters formed.}
#'     \item{Dbar}{Model fit statistics.}
#'     \item{pD}{Model fit statistics.}
#'     \item{pV}{Model fit statistics.}
#'     \item{DIC}{Model fit statistics.}
#'     \item{lpml}{Model fit statistics.}
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases ddp_quantiles 
#' @export
ddp_quantiles = function(model.output, dosemat, Nfixed, Nrandom, Nsubject, typet) 
{
  ##
  ## Input checks
  ##

  stopifnot(ncol(model.output$Beta) == Nfixed)
  stopifnot(ncol(model.output$Theta) == (Nrandom*Nsubject)) ## same dimensions and loading structure as B under MM models
  stopifnot(length(model.output$DoseEffects) == (length(typet)) + 1)

  ##
  ## compose model fit statistics
  ##

  pV		= var(model.output$Deviance)/2
  DIC 		= model.output$devres[1] 
  Dbar 		= model.output$devres[2] 
  Dhat 		= model.output$devres[3] 
  pD 		= model.output$devres[4]
  lpml		= model.output$lpml

  ##
  ## compose summaries for DDP effects
  ##

  ## Extract dimensions
  nty		<- length(model.output$DoseEffects) - 1 ## "-1" for the q intercept anova effects
  ## extract size for each treatment type
  poscar	<- which(typet == 1); numcar = length(poscar); 
  posmvn	<- which(typet == 2); nummvn = length(posmvn); countmvn = 1;
  posind	<- which(typet == 3); numind = length(posind); countind = 1;

  ## DoseEffects summary and mean matrix, the latter for plotting.
  dosetrt.summary	<- vector("list", nty)
  dosetrt.mean		<- vector("list", nty)
  numt			<- vector("numeric", nty)

  ## Precision parameters for base distribution treatments.
  if( nummvn > 0 )
  {
	pmvn.summary  	<- vector("list", nummvn) 
	pmvn.mean	<- vector("list", nummvn)
  }
  if( numind > 0 ) tauind.summary <- vector("list", numind)
  if( numcar > 0 ) 
  {
	Alphacar 		<- model.output$CAR_Q[[1]] ## nkeep x numcar matrix
	Taucar			<- model.output$CAR_Q[[2]] ## nkeep x numcar matrix
	mcmc.data		= cbind(Alphacar, Taucar)
	mcmc.mean		= colMeans(mcmc.data)
  	mcmc.low		= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
 	mcmc.high		= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
	mcmc.summary		= t(rbind(mcmc.low,mcmc.mean,mcmc.high)) ## transpose to 3 columns for c("mcmc.low","mcmc.mean","mcmc.high") and rows ordered as Deviance, Beta, .. 
  	alphacar.summary	= matrix(mcmc.summary[1:numcar,], numcar, 3) ## numcar x 3
	taucar.summary		= matrix(mcmc.summary[(numcar+1):(2*numcar),], numcar, 3) ## numcar x 3
  }
	
  start.labs 		= 2 
  for(m in 1:nty)
  {
	 ## DoseEffects summary and mean matrix, the latter for plotting.
	 numt[m]			= ncol(model.output$DoseEffects[[m]]) / (Nsubject*Nrandom)
	 mcmc.data			= model.output$DoseEffects[[m]]
	 mcmc.mean			= colMeans(mcmc.data)
  	 mcmc.low			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
 	 mcmc.high			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
         ## for dosetrt.summary[[m]], subject is slowest moving index, followed by polynomial order, followed by numt[m]
	 ## <--> Nsubject*(Nrandom*numt[m]) x 1, index speed order: subject, order, dose
  	 dosetrt.summary[[m]]		= t(rbind(mcmc.low,mcmc.mean,mcmc.high)) ## transpose to 3 columns for c("mcmc.low","mcmc.mean","mcmc.high")  
 	 dosetrt.mean[[m]]		= matrix(dosetrt.summary[[m]][,2], (Nrandom*numt[m]), Nsubject, byrow = FALSE) ## column 2 captures the means
	 dosetrt.mean[[m]]		= t(dosetrt.mean[[m]]) ## Nsubject x (Nrandom*numt[m]) - will take row i to construct delta_i for plotting heat map
	 if( !is.null(colnames(dosemat)) ) ## assign colnames to dosetrt.mean[[m]] to anticipate plotting
	 {
		end.labs			= start.labs + numt[m] - 1
		doselabs			= colnames(dosemat)[start.labs:end.labs]
		start.labs			= end.labs + 1
	 	colnames(dosetrt.mean[[m]])	= paste(rep(1:Nrandom, each = numt[m]), rep(doselabs, Nrandom), sep=".")
	 }else{ ## no colnames for dosemat
		colnames(dosetrt.mean[[m]])	= paste(rep(1:Nrandom, each = numt[m]), rep(1:numt[m], Nrandom), sep=".") ## 1.1, 1.2, ...., 1.numt[m]
	 }

   	## Precision parameters for base distribution treatments.
	if(typet[m] == 2) ## mvn
	{
		mcmc.data			= model.output$Pmvn[[countmvn]]
		mcmc.mean			= colMeans(mcmc.data)
  	 	mcmc.low			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
 	 	mcmc.high			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
		pmvn.summary[[countmvn]]	= t(rbind(mcmc.low,mcmc.mean,mcmc.high)) 
		pmvn.mean[[countmvn]]		= matrix(pmvn.summary[[countmvn]][,2], numt[m], numt[m], byrow = FALSE) ## loading by column reverses vectorization in DDP.cpp
		countmvn			= countmvn + 1
	}else{ 	
		if( typet[m] == 3 ) ## typet[m] == 3, ind
		{
			mcmc.data			= model.output$Tauind[[countind]]
			mcmc.mean			= colMeans(mcmc.data)
  	 		mcmc.low			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
 	 		mcmc.high			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
			tauind.summary[[countind]]	= t(rbind(mcmc.low,mcmc.mean,mcmc.high)) ## numt[m] x 3
			countind			= countind + 1
		}
	}
  }

  ##
  ## compose summaries for remaining model parameters
  ##

  mcmc.data		= cbind(model.output$Deviance,model.output$Beta,model.output$Alpha,model.output$Theta,model.output$Lambda,
						model.output$Taue,model.output$M,model.output$DoseEffects[[nty+1]]) #nrow(mcmc.data) = n.keep
  mcmc.mean		= colMeans(mcmc.data)
  mcmc.low		= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
  mcmc.high		= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
  rm(mcmc.data) #free up memory
  mcmc.summary		= t(rbind(mcmc.low,mcmc.mean,mcmc.high)) ## transpose to 3 columns for c("mcmc.low","mcmc.mean","mcmc.high") and rows ordered as Deviance, Beta, .. 
  deviance.summary	= mcmc.summary[1,]; 
  beta.summary		= mcmc.summary[2:(Nfixed+1),]
  alpha.summary		= mcmc.summary[(Nfixed+2),]
  theta.summary	  	= vector(mode = "list", length = Nrandom)
  for(i in 1:Nrandom)
  {
	theta.summary[[i]] = mcmc.summary[(Nfixed+3+(i-1)*Nsubject):(Nfixed+2+i*Nsubject),]
  }
  Nrandom.2			= Nrandom^2
  lambda.summary		= mcmc.summary[(Nfixed+3+Nrandom*Nsubject):(Nfixed+2+Nrandom*Nsubject+Nrandom.2),]
  taue.summary			= mcmc.summary[(Nfixed+3+Nrandom*Nsubject+Nrandom.2),]
  M.summary			= mcmc.summary[(Nfixed+4+Nrandom*Nsubject+Nrandom.2),]
  doseint.summary	  	= vector(mode = "list", length = Nrandom)
  for(i in 1:Nrandom)
  {
	doseint.summary[[i]] = mcmc.summary[(Nfixed+5+Nrandom*Nsubject+Nrandom.2+(i-1)*Nsubject):(Nfixed+4+Nrandom*Nsubject+Nrandom.2+i*Nsubject),]
  }
  ## create lambda.mean (Nrandom x Nrandom)
  lambda.mean			= matrix(lambda.summary[,2], Nrandom, Nrandom, byrow = FALSE) ## reverse by column stacking into a vector		

  ##
  ## return results
  ##

  res		<- list(deviance.summary = deviance.summary, beta.summary = beta.summary, alpha.summary = alpha.summary, theta.summary = theta.summary, 
					taue.summary = taue.summary, M.summary = M.summary, lambda.summary = lambda.summary, lambda.mean = lambda.mean,
					dosetrt.summary = dosetrt.summary, dosetrt.mean = dosetrt.mean, 
					doseint.summary = doseint.summary, pV = pV, DIC = DIC, Dbar = Dbar, Dhat = Dhat, pD = pD, lpml = lpml)
  if( numcar > 0 )
  {
	res$alphacar.summary		<- alphacar.summary
	res$taucar.summary		<- taucar.summary
  }
  if( numind > 0 )
  {
	res$tauind.summary		<- tauind.summary
  }
  if( nummvn > 0 )
  {
	res$pmvn.summary		<- pmvn.summary
	res$pmvn.mean			<- pmvn.mean
  }

  return(res)

} ## end function ddp_quantiles.R