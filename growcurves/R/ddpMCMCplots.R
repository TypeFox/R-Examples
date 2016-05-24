#' generate plots of posterior samples under \code{ddpgrow} model
#'
#' Constructs plots for subject effects, theta, with credible intervals as well as selected trace plots.
#' Returns a list of objects of class \code{ggplot}.
#' 
#' @param subjecti.u A vector of length \code{P}, number of unique subjects, containing unique set of user input values for \code{subject}.
#' @param labt An vector object (of the same length as \code{typet}) providing user names for each treatment.  The names are
#'		used in plot objects.
#' @param typet A numeric vector of length equal to the number of treatments that specifies prior option for each treatment.  Options must be one of: 1 = "car",
#' 		2 = "mvn", 3 = "ind".
#' @param numt A numeric vector with same length as \code{typet} where each entry counts the number of doses for that treatment.
#' @param theta.summary A list object of \code{q} elements, each containing an \code{P x 3} matrix of c(2.5\%,50\%,97.5\%) quantile summaries
#'		for each subject of the applicable subject random effect parameter.  \code{P} = number of subjects, \code{q} = number of random effect parameters, per subject.
#' @param lambda.mean A \code{q x q} matrix of mean values of the polyomial order covariance matrix, \code{Lambda}, returned from \code{ddpgrow}.
#' @param pmvn.mean A list object of length equal to the number of treatments with \code{"mvn" \%in\% typetreat}.  Each list object contains an \code{numt[m] x numt[m]}
#'			matrix of mean elements of the "mvn" treatments covariance matrices, \code{Pmvn}, returned from \code{ddpgrow}.
#' @param taucar.summary A \code{numcar x 3} numeric matrix of c(2.5\%,50\%,97.5\%) quantile summaries for the scale parameter for \code{numcar} treatments where \code{"car" \%in\% typetreat}.
#' @param alphacar.summary A \code{numcar x 3} numeric matrix of c(2.5\%,50\%,97.5\%) quantile summaries for the strength parameter for \code{numcar} treatments where \code{"car" \%in\% typetreat}.
#' @param Taucar An \code{iter.keep x numcar} matrix of posterior samples capturing the CAR precision parameter for each treatment where \code{"car" \%in\% typetreat}.
#' @param Alphacar An \code{iter.keep x numcar} matrix of posterior samples capturing the CAR strength parameter for each treatment where \code{"car" \%in\% typetreat}.
#' @param tauind.summary A list object of length equal to \code{numcar}, the number of treatment with \code{"in" \%in\% typetreat}.
#'				Each list element contains an \code{numt[m] x 3} matrix of c(2.5\%,50\%,97.5\%) quantile summaries for the dosage scale parameters associated
#'				to that treatment \code{m}.
#' @param Tauind A list object of length \code{numind}, the number treatments where \code{"car" \%in\% typetreat} with each element
#'		holding an \code{nkeep x numt[m]} matrix of sampled by-dose precision parameters for treatment \code{m}.
#' @param M The \code{iter.keep x 1} matrix of posterior samples for the parameter capturing the number of clusters formed under the DP prior on the client effects.
#' @param Taue \code{iter.keep x 1} matrix of posterior samples capturing the precision parameter for the model error term.
#' @param Deviance \code{iter.keep x 1} matrix of posterior samples for the model deviance.
#' @return A list of plot objects of class \code{ggplot2} including:
#'     	\item{p.theta}{stacked plots of b0,...,b(q-1) - vertical lines for 
#'     	each client span 2.5\% - 97.5\% values with mean noted.} 
#'     	\item{p.M}{MCMC trace plot of M, number of clusters.} 
#'     	\item{p.acar}{MCMC trace plots of alphacar, the CAR strength parameter. 
#'     	Plot is faceted for more than one set of treatments under \code{"car"}.} 
#'     	\item{p.tcar}{MCMC trace plots of alphacar, the CAR precision parameter. 
#'     	Plot is faceted for more than one set of treatments under \code{"car"}.} 
#'     	\item{p.taue}{MCMC trace plots of tau.e.} 
#'     	\item{p.dev}{MCMC trace plots of deviance.} 
#'	    \item{p.lam}{Heatmap (tiled) plot of posterior mean covariance, 
#'	    Lambda, of polynomial orders of random effects, Delta_i.}
#'	    \item{p.mvn}{Heatmap plots of posterior mean covariance, 
#'	    Pmvn, of treatment dosages under \code{"mvn"} base distribution. 
#'	    Plots are faceted by treatment.}
#'	    \item{p.aband}{95\% credible bands for alphacar CAR strength parameters 
#'	    associated to treatments under \code{"car"} base distribution.}
#'	    \item{p.tband}{95\% credible bands for taucar CAR precision parameters 
#'	    associated to treatments under \code{"car"} base distribution.}
#'	    \item{p.iband}{95\% credible bands for tauind dose precision parameters 
#'	    for \code{"ind"} base distribution, faceted by treatment (if more than one).}
#' @seealso \code{\link{ddpgrow}}, \code{\link{dpgrowmm}}, \code{\link{dpgrow}}, \code{\link{dpgrowmult}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{ddpgrow}}
#' @aliases ddpMCMCplots ddpplots
#' @export
ddpMCMCplots = function(subjecti.u, labt = NULL, typet, numt, theta.summary, lambda.mean, pmvn.mean = NULL, taucar.summary = NULL,
			alphacar.summary = NULL, Taucar = NULL, Alphacar = NULL,  tauind.summary = NULL, Tauind = NULL, M, Taue, Deviance)
{
   ##
   ## some dimensions
   ##	

   iter.keep	= length(Taue)
   Nrandom	= length(theta.summary)  ## theta.summary is a list object of length Nrandom

   ##
   ## produce plots
   ##

   ## theta_1, .., theta_nr - by subject - with credible ranges
   ## data.frame
   dat.long		  	= vector(mode = "list", length = Nrandom)
   for(i in 1:Nrandom)
   {
	## note: using unique values of subject labels of length P, not case formatted of length N.
	dat			= as.data.frame(cbind(subjecti.u,theta.summary[[i]])) ## using unique values of subject.input for plotting
	names(dat)		= c("subject","low","mean","high")
	dat.long[[i]]		= melt(dat,id="subject",measure = c("low","mean","high"))
	records			= nrow(dat.long[[i]])
	labs			= rep(paste("theta",i-1,sep=""),records)
	dat.long[[i]]		= as.data.frame(cbind(dat.long[[i]],labs))
	names(dat.long[[i]])[4] = "type"
   } ## end loop i to create long data.frame with label identifier for each bi
      	datT			= as.data.frame(do.call("rbind",dat.long))
	datT$type		= factor(datT$type)
	rm(dat); rm(dat.long)

   ## plot
   p.theta 	= ggplot(data=datT,aes(x=subject, y=value , group=factor(subject)))
   l.1 		= geom_line(colour = "steelblue4")
   ## l.2	= stat_summary(fun.y = mean, geom = "point", size = 3, shape = 16, colour = "pink")
   l.2		= geom_smooth(aes(group=1),method = "loess", span = .2, size = 1.2, colour = "pink")
   f 		= facet_wrap(~type, scales="free_y",ncol=1)
   axis		= labs(x = "Subject", y = "Effect Size")
   options 	= labs(title="Subject Effects, Theta")
   p.theta 	= p.theta + l.1 + l.2 + f + axis + options 

   ## tau.e - MCMC trace
   dattaue		= as.data.frame(cbind(1:iter.keep,Taue))
   names(dattaue) 	= c("iteration","value")
   p.taue		= ggplot(data=dattaue,aes(x=iteration,y=value))
   l			= geom_line()
   axis			= labs(x = "Iterations", y = "Sampled Value")
   options		= labs(title="MCMC Trace plot for tau.e")
   p.taue		= p.taue + l + axis + options

   ## deviance - MCMC trace
   datdev		= as.data.frame(cbind(1:iter.keep,Deviance))
   names(datdev) 	= c("iteration","value")
   p.dev		= ggplot(data=datdev,aes(x=iteration,y=value))
   l			= geom_line()
   axis			= labs(x = "Iterations", y = "Sampled Value")
   options		= labs(title="MCMC Trace plot for Deviance")
   p.dev		= p.dev + l + axis + options

   ## M - MCMC trace
   datM		= as.data.frame(cbind(1:iter.keep,M))
   names(datM) 	= c("iteration","value")
   p.M		= ggplot(data=datM,aes(x=iteration,y=value))
   l		= geom_line()
   axis		= labs(x = "Iterations", y = "Sampled Value")
   options	= labs(title="MCMC Trace plot for M")
   p.M		= p.M + l + axis + options

   ## Heatmap of covariance matrix, Lambda
   covmat	= lambda.mean
   dat		= melt(covmat)
   names(dat)	= c("dose","column","effects")
   p.lam	= ggplot(data=dat,aes(x=dose, y = column, fill = effects))
   l		= geom_tile()
   yaxis 	= ylab("Dose")
   xaxis 	= xlab("Dose")
   options	= labs(title="Heatmap for Base Dist. Columns Covariance Matrix, Lambda.")
   p.lam	= p.lam + l + xaxis + yaxis + scale_fill_gradient(low = "white", high = "steelblue") + options


  ## Heatmap of covariance matrix, Pmvn
  if( any(typet == 2) ) ## if(!is.null(pmvn.mean))
  {
	## dimensions for number of mvn effect terms and associated names
  	nt		= length(typet[typet == 2]) ## = length(pmvn.mean)
	locmvn		= which(typet == 2)
	dat		= vector("list",nt)
	if( !is.null(labt) )
	{
		namesmvn	= labt[locmvn]
	}else{
		namesmvn	= paste("trt",locmvn,sep="_")
	}
	## construct plot data.frame objecct
 	for(m in 1:nt)
  	{
		numt.m		= numt[locmvn[m]]
		covmat.m	= pmvn.mean[[m]]
		dat[[m]] 	= melt(covmat.m)
		dat[[m]]$trt	= namesmvn[[m]]
  	}
  	dat		= do.call("rbind", dat)
  	names(dat)	= c("dose","column","effects","trt")	

  	## output plots
  	## tile geom 
  	p.mvn		= ggplot(data=dat,aes(x=dose, y = column, fill = effects))
  	l		= geom_tile()
  	yaxis 		= ylab("Dose")
  	xaxis 		= xlab("Dose")
  	f		= facet_wrap(~trt, scales="free_y", ncol=min(nt,5) )
	options		= labs(title="Heatmap for Covariance Matrices under MVN Prior.")
  	p.mvn		= p.mvn + l + xaxis + yaxis + f + scale_fill_gradient(low = "white", high = "steelblue") + options
  }else{
	p.mvn <- NULL
  }

  if( any(typet == 1) ) ## !is.null(alphacar.summary); nrow(Alphacar) > 0 
   {
	##
	## Trace plots (faceted by treatment) of Alphacar and Taucar
	##

	## extract treatment names for CAR treatments
	loccar		= which( typet == 1 )
	if( !is.null(labt) )
	{
		namescar	= labt[loccar]
	}else{ ## is.null(labt)
		namescar	= paste("trt",loccar,sep="_")
	}

	## Trace plot Alphacar
	numcar		= ncol(Alphacar) ## nrow(alphacar.summary)
   	## Alphacar - MCMC trace
   	datA		= as.data.frame(cbind(1:iter.keep,Alphacar))
	names(datA)	= c("iteration",namescar)
	datA		= melt(datA,id="iteration")
	names(datA) 	= c("iteration","label","value")
   	p.acar		= ggplot(data=datA,aes(x=iteration,y=value))
   	l		= geom_line()
	f 		= facet_wrap(~label, scales="free_y",ncol=1)
   	axis		= labs(x = "Iterations", y = "Sampled Value")
   	options		= labs(title="MCMC Trace plot for CAR strength parameter, Alphacar")
   	p.acar		= p.acar + l + f + axis + options

	## Taucar - MCMC trace
   	datt		= as.data.frame(cbind(1:iter.keep,Taucar))
   	names(datt)	= c("iteration",namescar)
	datt		= melt(datt,id="iteration")
	names(datt) 	= c("iteration","label","value")
   	p.tcar		= ggplot(data=datt,aes(x=iteration,y=value))
   	l		= geom_line()
	f 		= facet_wrap(~label, scales="free_y",ncol=1)
   	axis		= labs(x = "Iterations", y = "Sampled Value")
   	options		= labs(title="MCMC Trace plot for CAR precision parameter, Taucar")
   	p.tcar		= p.tcar + l + f + axis + options

	##
	## Credible Bands (by treatment) of Alphacar and Taucar
	##

	## alphacar
	alphacar.summary	= matrix(alphacar.summary, numcar, 3)
	dat			= data.frame(namescar,alphacar.summary) ## using unique values of subject.input for plotting
	names(dat)		= c("treatment","low","mean","high")
	dat.long		= melt(dat,id="treatment",measure = c("low","mean","high"))

	## plot
   	p.aband		= ggplot(data=dat.long,aes(x=treatment, y=value , group=factor(treatment)))
   	l.1 			= geom_line(colour = "steelblue4")
   	l.2			= geom_smooth(aes(group=1),method = "loess", span = .8, size = 1.2, colour = "pink")
   	axis			= labs(x = "Treatment", y = "95% Credible Interval")
   	options 		= labs(title="CAR Strength Parameters, alpha")
   	p.aband 		= p.aband + l.1 + l.2 + axis + options 

	## taucar
	taucar.summary		= matrix(taucar.summary, numcar, 3)
	dat			= data.frame(namescar,taucar.summary) ## using unique values of subject.input for plotting
	names(dat)		= c("treatment","low","mean","high")
	dat.long		= melt(dat,id="treatment",measure = c("low","mean","high"))

	## plot
   	p.tband		= ggplot(data=dat.long,aes(x=treatment, y=value , group=factor(treatment)))
   	l.1 		= geom_line(colour = "steelblue4")
   	l.2		= geom_smooth(aes(group=1),method = "loess", span = .8, size = 1.2, colour = "pink")
   	axis		= labs(x = "Treatment", y = "95% Credible Interval")
   	options 	= labs(title="CAR Precision Parameters, tau")
   	p.tband 	= p.tband + l.1 + l.2 + axis + options 

   }else{ 
		p.acar <- p.tcar <- p.aband <- p.tband <- NULL
   }

   if( !is.null(tauind.summary) ) ## plot trace and confidence intervals for precision parameter under independent prior
   {
	##
	## Credible Bands (by treatment) of Tauind
	##
	
	## dimensions for number of "ind" effect terms and associated names
  	nt		= length(typet[typet == 3]) ## = length(tauind.summary), which is a list object since each dose.m gets it's own precision parameter
	locind		= which(typet == 3)
	dat		= vector("list",nt)
	if( !is.null(labt) )
	{
		namesind	= labt[locind]
	}else{
		namesind	= paste("trt",locind,sep="_")
	}

	## construct plot data.frame objecct
 	for(m in 1:nt)
  	{
		numt.m		= numt[locind[m]]
		dat[[m]] 	= data.frame(1:numt.m,tauind.summary[[m]],namesind[[m]])
  	}
  	dat			= do.call("rbind", dat)
  	names(dat)		= c("dose","low","mean","high","trt")
	dat.long		= melt(dat,id=c("trt","dose"), measure = c("low","mean","high"))

	## plot
   	p.iband		= ggplot(data=dat.long,aes(x=dose, y=value , group=factor(dose)))
   	l.1 		= geom_line(colour = "steelblue4")
   	l.2		= geom_smooth(aes(group=1),method = "loess", span = .8, size = 1.2, colour = "pink")
	f 		= facet_wrap(~trt, scales="free_y",ncol=1)
   	axis		= labs(x = "Dose Precision Parameter", y = "95% Credible Interval")
   	options 	= labs(title="Precision Parameter under Independent Gaussian base, tauind")
   	p.iband 	= p.iband + l.1 + l.2 + f + axis + options 		

   }else{
   	p.iband <- NULL
   }

   ##
   ## return results
   ##

   res		<- list(p.theta = p.theta, p.taue = p.taue, p.dev = p.dev, p.M = p.M, p.lam = p.lam, p.acar = p.acar, p.tcar = p.tcar, p.aband = p.aband, p.tband = p.tband,
			p.iband = p.iband, p.mvn = p.mvn)
   res		<- res[!sapply(res, is.null)]
   return(res)

   value <- iteration <- subject <- dose <- column <- trt <- treatment <- effects <- NULL; 
   rm(value); rm(iteration);

} ## end function ddpMCMCplots