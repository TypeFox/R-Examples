#' Plot comparison of Mean Effects for Any Two Treatments 
#'
#' Produces a set of box plots of the posterior distributions of the fixed effects mean difference between any 2 chosen treatments
#' across a chosen subset models (priors) at each of a selection of time points.  Each box plot spans the 95% credible interval.
#'
#' @param run.objects A list object where each entry is either a \code{dpgrow}, \code{dgrowmm} or \code{dpgrowmult} object
#'	returned from a model run.  The list length is equal to the number of models on which a treatment comparison is desired.
#' @param run.models A character vector supplying the names of models to compare of the same length as the number of elements in \code{run.objects}.
#' @param trt.labs A vector of exactly 2-elements that contains the labels (used in modeling) for the two treatments desired to compare.
#' @param time.points A numeric vector of data time points on which to compare the fixed effects treatment means.
#' @param y.label An optional character scalar for label for the y axis.
#' @param time.labels An optional character vector of labels for \code{time.points} to use in plotting.
#' @param n.thin Gap between successive MCMC sampling iterations on fixed effects parameters to use for generating distribution for each treatment-time.  
#'		Defaults to \code{n.thin = 10}.  Higher values produces faster plot generation.
#' @return A list object containing quantile summaries for all sampled model parameters.
#'     \item{dat.trt}{A \code{data.frame} object used to generate the fixed effects mean comparison for two chosen treatments.
#'			Fields are titled, \code{c("Mu_diff","models","time")}.}
#'     \item{p.trt}{A \code{ggplot2} object of box plots split by time point.}
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}, \code{\link{dpgrowmult}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases trtplot
#' @export

trtplot <- function(run.objects, run.models, trt.labs, time.points, y.label = NULL, time.labels = NULL, n.thin = 10)
{
  ## check inputs
  if( !is.list(run.objects) ) stop("\n'run.objects' must be a list object.\n")
  if(length(run.objects) != length(run.models)) stop("\nNumber of objects must be equal to number of models.\n")
  if( length(trt.labs) > 2 ) 
  {
	trt.labs = trt.labs[1:2]
	warning("\nOnly the first two entries in 'trt.labs' will be used for plotting\n")
  }
  if(length(time.points) > 4 ) warning("\nPlot quality will decrease for more than 4 time points.\n")
  if( !is.null(time.labels) ) stopifnot(length(time.points) == length(time.labels))

  ### save fixed effects posterior samples.
  models	= run.models
  num.mods	= length(models)
  num.time	= length(time.points)

  ## Extract design matrix - Assuming same data used to generate all models
  map.trt 	= summary(run.objects[[1]])$summary.results$map.trt ## data.frame
  map.trt	= unique(map.trt)
  if( length(setdiff(map.trt$label.input,trt.labs)) > 0 ) stop("\nTreatment labels must be subset of those used for running the models.\n")
  cntl.trt	= subset(map.trt, (label.input %in% trt.labs) & (label.new == 0) )
  if( nrow(cntl.trt) > 0  )
  {
	trt.cntl = cntl.trt$label.input
	pos.cntl = which(trt.labs == trt.cntl)
	pos.trt  = which(trt.labs != trt.cntl)
  }else{
	pos.cntl = NULL
  }

  ## set-up run 
  Mu.dat	= vector(mode = "list",num.mods*num.time)

  for(m in 1:num.mods) ## over models
  {
	X.template.m		= as.data.frame(summary(run.objects[[m]])$summary.results$X)
	Z			= as.data.frame(summary(run.objects[[m]])$summary.results$Z)
	n.random		= ncol(Z)
	n.fix_degree		= summary(run.objects[[m]])$summary.results$n.fix_degree
	Beta			<- samples(run.objects[[m]])$Beta
	Beta			= Beta[,-(1:n.fix_degree)]
	Beta			= as.matrix(Beta)
	M			= nrow(Beta) ## number of MCMC iterations
	keep			= floor(M/n.thin)
	samps			= as.integer(seq(1,M, length.out = keep)) ## set of thinned MCMC iterations to sample
	for( j in 1:num.time ) ## over time points
	{
			X									= X.template.m
			X									= subset(X,X[,1] == time.points[j])
			X									= X[,-(1:n.fix_degree)]
			## for on whether one of the treatments is the control
			if( is.null(pos.cntl) )
  			{
				## note: fixing the time point and trt, covariate values are the same across subjects.
				## trt1
				selpos.trt1					= grep(trt.labs[1],names(X))
				trt1.name					= as.name(names(X)[selpos.trt1[1]]) ## isolate the treatment covariate in {0,1}
				xtrt.1						= X[,c(selpos.trt1)]
				xtrt.1						= subset(xtrt.1,(eval(trt1.name) == 1))[1,] ## all subjects have same covariate valuse
				## trt2
  				selpos.trt2					= grep(trt.labs[2],names(X))
				trt2.name					= as.name(names(X)[selpos.trt2[1]])
				xtrt.2						= X[,c(selpos.trt2)]
				xtrt.2						= subset(xtrt.2,(eval(trt2.name) == 1))[1,]
				## compute Mu.f (fixed effects treatment values); time covariates don't vary across treatments
				Mutrt.1						= Beta[samps,c(selpos.trt1)] %*% t(as.matrix(xtrt.1))
				Mutrt.2						= Beta[samps,c(selpos.trt2)] %*% t(as.matrix(xtrt.2))
				Mu.diff						= Mutrt.2 - Mutrt.1
			}else{ ## one of the treatments is a base treatment
				## pos.trt
				selpos.trt					= grep(trt.labs[pos.trt],names(X))
				trt.name					= as.name(names(X)[selpos.trt[1]]) ## isolate the treatment covariate in {0,1}
				xtrt						= X[,c(selpos.trt)]
				xtrt						= subset(xtrt,(eval(trt.name) == 1))[1,]
				## compute Mu.f (fixed effects treatment values). Trt effect for control is exactly 0.
				Mu.diff						= Beta[samps,c(selpos.trt)] %*% t(as.matrix(xtrt))
			} ## end fork on whether one of the comparison treatments is the control
			## create data.frame for each model
			Mudwlab							= data.frame(Mu.diff,models[m],time.points[j])
			names(Mudwlab)						= c("Mu_diff","models","time")
			Mu.dat[[(m-1)*num.time + j]]				= Mudwlab	
 	} ## end loop j over measurement wave times
  } ## end loop m over models  


  ## Mu.dat
  Mu.dat			= do.call("rbind",Mu.dat)
  if( !(is.null(time.labels)) )
  {
  	Mu.dat$time		= factor(Mu.dat$time,labels = time.labels)
  }else{
	Mu.dat$time		= factor(Mu.dat$time)
  }

  ## plot Y-axis label
  if( !is.null(y.label) )
  {
	y.lab	<- y.label
  }else{
  	if( is.null(pos.cntl) )
  	{
		y.lab	<- paste("TRT_",trt.labs[2]," - TRT_",trt.labs[1]," Fixed Mean Difference",sep="")
  	}else{
		y.lab 	<- paste("TRT_",trt.labs[pos.trt]," -  TRT_",trt.labs[pos.cntl],"  Fixed Mean Difference",sep="")
  	}
  }
  ## distribution boxplots
  f <- function(x) {
    r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  p	<- ggplot(data=Mu.dat,aes(x = models, y = Mu_diff, fill = models))
  l	<- stat_summary(fun.data = f, geom="boxplot")
  f	<- facet_wrap(~time,scales = "fixed")
  axis	<- labs(x = "Model", y = y.lab)
  if ( length(unique(Mu.dat$models)) >= 3 )
  {
  	p	<- p + l + f + axis + scale_fill_brewer(palette="OrRd")  ## color scale supports bw printing
  }else{ ## need >= 3 categories to use OrRd palette
	p	<- p + l + f + axis
  }
  dev.new()
  print(p)

  return(invisible(list(p.trt = p, dat.trt = Mu.dat)))

  Mu_diff <- models <- time <- label.input <- label.new <- NULL; rm(Mu_diff); rm(models); rm(time); rm(label.input); rm(label.new);

} ## end function trtplot