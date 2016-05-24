#' Produce quantile summaries of model posterior samples
#'
#' Inputs MCMC samples for model parameters and constructs c(2.5\%,50\%,97.5\%) quantile summaries.
#'
#' @param model.output  An output object of class within \code{c("dpgrow", "dpgrowmm")}
#' @param Nfixed Number of total fixed effects, both time-based and nuisance.
#' @param Nrandom Number of total random effects, both time-based and nuisance, all grouped by subject.
#' @param Nsubject Number of unique subjects (on which repeated measures are observed).
#' @param Nsubj.aff Number of subjects, \code{P.aff}, receiving multiple membership effects
#' @param Nmv Number of multivariate MM effects.  Defaults to 1 for univariate MM if left blank.
#' @param Nsession Number of multiple membership effects for each entry in "Nmv".  May be left blank for univariate MM.
#' @return A list object containing quantile summaries for all sampled model parameters.
#'     \item{deviance.summary}{vector of length 3 summarizing quantiles for model deviance.}
#'     \item{beta.summary}{\code{Nfixed x 3} quantile summaries of model fixed effects.}
#'     \item{alpha.summary}{quantile summary of model global intercept parameter.}
#'     \item{bmat.summary}{list object of length \code{Nrandom}, each cell containing a \code{Nsubject x 3} matrix of by-subject parameter quantile summaries.}
#'     \item{tauu.summary}{\code{Nmv x 3} quantile summary for prior precision parameters employed for multiple membership random effects. An \code{nty x 3} matrix in the case of \code{nty} multiple membership effect terms.}
#'     \item{taue.summary}{quantile summary for model error precision parameter.}
#'     \item{taub}{\code{Nrandom x 3} quantile summaries for subject effect precision parameters.}
#'     \item{u.summary}{\code{S*Nmv x 3} quantile summaries for multiple membership random effect parameters. A list of such matrices in the case of \code{nty} multiple membership effect terms.}
#'     \item{mm.summary}{\code{Nsubj.aff x 3} quantile summaries derived from multiplying the affected subject weight matrix by the multiple membership random effects.}
#'     \item{M.summary}{quantile summary for number of DP posterior clusters formed.}
#'     \item{Dbar}{Model fit statistics.}
#'     \item{pD}{Model fit statistics.}
#'     \item{pV}{Model fit statistics.}
#'     \item{DIC}{Model fit statistics.}
#'     \item{lpml}{Model fit statistics.}
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases summary_quantiles quantiles
#' @export
summary_quantiles = function(model.output, Nfixed, Nrandom, Nsubject, Nsubj.aff = NULL, Nmv = 1, Nsession = NULL) ## Nsubj.aff for capturing MM mapping of U only under single MM term.
{
  stopifnot(ncol(model.output$Beta) == Nfixed)
  stopifnot(ncol(model.output$B) == (Nrandom*Nsubject))

  if( (Nmv > 1) & is.null(Nsession) ) stop("\nMust input 'Nsession' in the case of multivariate MM.\n")

  if( !is.null(model.output$U) ) ## DP + MM
  {
	## if( ncol(model.output$Tauu) > 1 )
	if( is.list(model.output$U) ) ## multiple MM terms
	{
		nty	= ncol(model.output$Tauu)
		stopifnot(length(model.output$U) == nty)
		numt	= vector("numeric",nty)
		for( i in 1:nty )
		{
			numt[i] = ncol(model.output$U[[i]])
		}
		U 		= do.call("cbind", model.output$U) ## creates a global matrix cbinding each U[[i]]
		stopifnot(ncol(U) == sum(numt))
		MM 		= NULL ## will allow this to pass through cbind
	}else{ ## single term MM - either univariate or multivariate
		U		= model.output$U
	}
	mcmc.data		= cbind(model.output$Deviance,model.output$Beta,model.output$Alpha,model.output$B,U,model.output$M,model.output$Tauu,
					model.output$Taue,model.output$Taub,model.output$MM) #nrow(mcmc.data) = n.keep
  }else{
	if( !is.null(model.output$M) ) ## LGM or DP
	{
		mcmc.data		= cbind(model.output$Deviance,model.output$Beta,model.output$Alpha,model.output$B,
						model.output$Taue,model.output$Taub,model.output$M) #nrow(mcmc.data) = n.keep
	}else{ ## LGM
		mcmc.data		= cbind(model.output$Deviance,model.output$Beta,model.output$Alpha,model.output$B,
						model.output$Taue,model.output$Taub) #nrow(mcmc.data) = n.keep
	}
  }

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
  ## create summary matrix of sampled parater quantiles
  ## 
  mcmc.mean		= colMeans(mcmc.data)
  mcmc.low		= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
  mcmc.high		= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
  rm(mcmc.data) #free up memory
  mcmc.summary		= t(rbind(mcmc.low,mcmc.mean,mcmc.high)) ## transpose to 3 columns for c("mcmc.low","mcmc.mean","mcmc.high") and rows ordered as Deviance, Beta, .. 
  deviance.summary	= mcmc.summary[1,]; 
  beta.summary		= mcmc.summary[2:(Nfixed+1),]
  alpha.summary		= mcmc.summary[(Nfixed+2),]
  bmat.summary	  	= vector(mode = "list", length = Nrandom)
  for(i in 1:Nrandom)
  {
	bmat.summary[[i]] = mcmc.summary[(Nfixed+3+(i-1)*Nsubject):(Nfixed+2+i*Nsubject),]
  }
  stopifnot(nrow(bmat.summary[[1]]) == Nsubject)

  ## fork for model-unique parameters
  if( !is.null(model.output$U) ) ## DP + MM
  { 
	## if( ncol(model.output$Tauu) > 1 )
	if( is.list(model.output$U) )
	{
		u.summary 	= vector(mode = "list", length = nty)
		n.cumsessi	= 0
		for(i in 1:nty)
		{
			u.summary[[i]] 	= mcmc.summary[(Nfixed+3+Nrandom*Nsubject+n.cumsessi):(Nfixed+2+Nrandom*Nsubject+n.cumsessi+numt[i]),]
			n.cumsessi 	= n.cumsessi + numt[i]
		}
		stopifnot(nrow(u.summary[[nty]]) == numt[nty])
		Nsession 		= sum(numt)
		stopifnot(n.cumsessi == Nsession)
		M.summary		= mcmc.summary[(Nfixed+3+Nrandom*Nsubject+Nsession),]
		tauu.summary		= mcmc.summary[(Nfixed+4+Nrandom*Nsubject+Nsession):(Nfixed+3+Nrandom*Nsubject+Nsession+nty),]
		stopifnot(nrow(tauu.summary) == nty)
		taue.summary		= mcmc.summary[(Nfixed+4+Nrandom*Nsubject+Nsession+nty),]
		taub.summary		= mcmc.summary[(Nfixed+5+Nrandom*Nsubject+Nsession+nty):(Nfixed+4+Nrandom*Nsubject+Nsession+nty+Nrandom),]
		mm.summary		= NULL ## will show up in summary.results under "mm.summary", but with value = NULL
		stopifnot(length(u.summary) == nty)  #dimension check
		## return summary objects for dpgrowmult
		return(list(deviance.summary = deviance.summary, beta.summary = beta.summary, alpha.summary = alpha.summary, bmat.summary = bmat.summary, 
					M.summary = M.summary, u.summary = u.summary, tauu.summary = tauu.summary,  
					taue.summary = taue.summary, taub.summary = taub.summary,
					mm.summary = mm.summary, pV = pV, DIC = DIC, Dbar = Dbar, Dhat = Dhat, pD = pD, lpml = lpml))
	}else{ ## only one MM term
		u.summary			= mcmc.summary[(Nfixed+3+Nrandom*Nsubject):(Nfixed+2+Nrandom*Nsubject+Nmv*Nsession),]
		rownames(u.summary)		= paste(rep(paste("eff",1:Nsession,sep="_"),times=Nmv),rep(paste("order",1:Nmv,sep="_"),each=Nsession),sep=".")
  		M.summary			= mcmc.summary[(Nfixed+3+Nrandom*Nsubject+Nmv*Nsession),]
  		tauu.summary			= mcmc.summary[(Nfixed+4+Nrandom*Nsubject+Nmv*Nsession):(Nfixed+3+Nrandom*Nsubject+Nmv*Nsession+Nmv),] ## vector of length Nmv
  		taue.summary			= mcmc.summary[(Nfixed+4+Nrandom*Nsubject+Nmv*Nsession+Nmv),]
  		taub.summary			= mcmc.summary[(Nfixed+5+Nrandom*Nsubject+Nmv*Nsession+Nmv):(Nfixed+4+Nrandom*Nsubject+Nmv*Nsession+Nmv+Nrandom),]
		mm.summary			= mcmc.summary[(Nfixed+5+Nrandom*Nsubject+Nmv*Nsession+Nmv+Nrandom):(Nfixed+4+Nrandom*Nsubject+Nmv*Nsession+Nmv+Nrandom+Nmv*Nsubj.aff),]
		rownames(mm.summary)		= paste(rep(paste("subj",1:Nsubj.aff,sep="_"),times=Nmv),rep(paste("order",1:Nmv,sep="_"),each=Nsubj.aff),sep=".")
		stopifnot(nrow(u.summary) == Nmv*Nsession)  #dimension check
		if( Nmv > 1 ) ## multivariate MM
		{
			## rhotauu
			mcmc.data			<- model.output$Rhotauu
			mcmc.mean			= colMeans(mcmc.data)
  			mcmc.low			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
  			mcmc.high			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
  			rm(mcmc.data) #free up memory
 		 	rhotauu.summary			= t(rbind(mcmc.low,mcmc.mean,mcmc.high))
			colnames(tauu.summary)		= paste("order",1:Nmv,sep="_")
		}else{ ## univariate MM
			rhotauu.summary			= NULL	
		} ## end conditional statement on whether dpgrowmm object is univariate or multivariate
		## return summary objects for dpgrowmm
		return(list(deviance.summary = deviance.summary, beta.summary = beta.summary, alpha.summary = alpha.summary, bmat.summary = bmat.summary, 
				M.summary = M.summary, u.summary = u.summary, tauu.summary = tauu.summary,  rhotauu.summary = rhotauu.summary,
				taue.summary = taue.summary, taub.summary = taub.summary,
				mm.summary = mm.summary, pV = pV, DIC = DIC, Dbar = Dbar, Dhat = Dhat, pD = pD, lpml = lpml))
	}
	
  }else{ ## DP or LGM
  	taue.summary		= mcmc.summary[(Nfixed+3+Nrandom*Nsubject),]
  	taub.summary		= mcmc.summary[(Nfixed+4+Nrandom*Nsubject):(Nfixed+3+Nrandom*Nsubject+Nrandom),]
  	if( !is.null(model.output$M) ) ## DP
  	{
    		M.summary		= mcmc.summary[(Nfixed+4+Nrandom*Nsubject+Nrandom),]
		return(list(deviance.summary = deviance.summary, beta.summary = beta.summary, alpha.summary = alpha.summary, bmat.summary = bmat.summary, 
					M.summary = M.summary, taue.summary = taue.summary, taub.summary = taub.summary,
					pV = pV, DIC = DIC, Dbar = Dbar, Dhat = Dhat, pD = pD, 
					dev8 = model.output$devresp, logcpo = model.output$logcpo, lpml = lpml))
  	}else{ ## LGM
		return(list(deviance.summary = deviance.summary, beta.summary = beta.summary, alpha.summary = alpha.summary, bmat.summary = bmat.summary, 
					taue.summary = taue.summary, taub.summary = taub.summary,
					pV = pV, DIC = DIC, DIC3 = model.output$devres3[1], Dbar = Dbar, Dbar3 = model.output$devres3[2], Dhat = Dhat, 
					pD = pD, logcpo = model.output$logcpo, lpml = lpml))
	}
  } ## end fork for returning model-unique parameters

} ## end function summary_quantiles.R

  

