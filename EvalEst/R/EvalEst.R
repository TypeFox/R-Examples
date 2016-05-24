############################################################################
# Comments here need to be cleaned up to reflect re-org
# Functions in the next group are mainly for evaluating estimation techniques.

# The first group are for generating simulations (ie- generate multiple
#   stochastic simulations of a model using simulate.)
#   These are good for looking at the stochastic properties of a model,
#   but mostly these are intended
#   only for verification purposes since other functions also generate 
#   simulations and it is usually more efficient to regenerate by setting
#   the RNG and seed than it is to save data.
#   The main function in this group is MonteCarloSimulations().
#   The main object class in this group is "simulation"

# The second group are for analysing the convergence of estimators. This
#  is extended to functions of estimates (such as model roots). These
#  functions evaluate a single given estimation method with multiple
#  simulated data sets from a given "true" model.
#  The main function in this group is EstEval().
#  The main object classes in this group are
#    "EstEval"
#     c("rootsEstEval","EstEval")
#     c("TSmodelEstEval","EstEval")
#     c("TSestModelEstEval","EstEval")
#     c("coefEstEval","EstEval")           and
#     c("rootsEstEval","EstEval")


# The third group applies multiple estimation techniques to a given data set.
#  This is primarily a utility for other functions 
#  The main (only) function in this group is estimateModels().
#  It returns an object of class c("estimatedModels")


# The fourth group looks at the forecasts and the covariance of forecasts 
#   for multiple horizons. 
#  The simplest case horizonForecasts() which calculates the forecast for 
#    different horizons at all periods in the sample. This is primarily
#    a utility for calculating the forecast error.
#  The next case is is estimatorsHorizonForecastsWRTdata() which 
#   is an extention of horizonForecasts(). It takes specified data 
#   and estimation techniques and calculates forecasts from the estimated 
#   models.
#  The generic function forecastCov() which considers mulitple  
#   models and calculates the cov of forecasts relative to a given data set.
#   It takes a list of models (+ trend, zero)  and calculates the cov 
#   of predictions. It uses forecastCovSingleModel()
#  The next case, forecastCovEstimatorsWRTdata() uses a list of estimation
#   methods to estimate a list of models and calculate the cov of predictions 
#   relative to one given data set.

#   The next case forecastCovWRTtrue() takes a list of models (+ trend,
#   zero)  and calculates the cov of forecasts relative to data sets 
#   simulated with a true.model.
#   The next case, forecastCovEstimatorsWRTtrue simulates data and 
#   uses a list of estimation methods to estimate a list of models, then
#   calculates the cov of predictions relative to other simulated data set.
#  The main object classes in this group are
#     c("estimatorsHorizonForecastsWRTdata") # ? "horizonForecasts")
#     "horizonForecasts"
#     c("multiModelHorizonForecasts","horizonForecasts")
#     "forecastCov"
#     c("forecastCovWRTdata", "forecastCov")
#     c("forecastCovWRTtrue", "forecastCov")
#     c("forecastCovEstimatorsWRTdata",  "forecastCov")
#     c("forecastCovEstimatorsWRTtrue",  "forecastCov")


# The fifth group are some experimental estimation techniques.

############################################################################
#
#       methods for MonteCarloSimulations  <<<<<<<<<<
#
############################################################################

generateSSmodel <- function(m,n,p, stable=FALSE)
 {#randomly generate an innov state space model. Discard models with largest root
  # greater than 1 (if stable=F) or equal to or greater than 1 if stable=T.
  repeat
    {FF <- matrix(runif(n^2, min=-1,max=1),n,n)
     if (m!=0) G <- matrix(runif(n*m, min=-1,max=1),n,m)
     else G <- NULL
     H <- matrix(runif(n*p, min=-1,max=1),p,n)
     K <- matrix(runif(n*p, min=-1,max=1),n,p)
     model <- SS(F.=FF, G=G,H=H,K=K)
     if (stable) {if (max(Mod(roots(model))) <  1.0) break()}
     else        {if (max(Mod(roots(model))) <= 1.0) break()}
    }
  model
 }


MonteCarloSimulations <- function(model, simulation.args=NULL, 
           replications=100, rng=NULL, quiet=FALSE, ...) UseMethod("MonteCarloSimulations")

MonteCarloSimulations.default <- function (model, simulation.args=NULL, 
 		replications=100, rng=NULL, quiet=FALSE, ...){
        #  (... further arguments, currently disregarded)
 	if (is.null(rng)) rng <- setRNG()
	else {
		old.rng <- setRNG(rng)
		on.exit(setRNG(old.rng))
	}
	arglist <- append(list(model), simulation.args)
	r <- do.call("simulate", arglist)
	if (! is.matrix(r)) stop("simulate(model) must return a matrix.")
        result <- array(NA, c(dim(r), replications))
        tfr <- tframe(r)
        result[, , 1] <- r
        if (1 < replications) 
            for (i in 2:replications)
	      result[, , i] <- do.call("simulate", arglist)
	# default does not work on array seriesNames(result) <- seriesNames(r)
        if (length(seriesNames(r)) != dim(result)[2])
         stop("length of names (",length(seriesNames(r)),
	      ") does not match number of series(",dim(result)[2],").")
        attr(result,"seriesNames") <- seriesNames(r)
	result <- tframed(result, tfr)
	invisible(classed(list(simulations = result, model = model, 
        rng = rng,  simulation.args = simulation.args, 
        description = "data generated by MonteCarloSimulations.default"), 
        c("MonteCarloSimulations")))
}


MonteCarloSimulations.TSestModel <- function(model, simulation.args=NULL, 
           replications=100, rng=NULL, quiet=FALSE, ...)
  {if (is.null(simulation.args$sd) & is.null(simulation.args$Cov)) 
     simulation.args$Cov <- model$estimates$cov
   if (is.null(simulation.args$input)) simulation.args$input <- inputData(model)
   MonteCarloSimulations(TSmodel(model), simulation.args=simulation.args, 
           replications=replications, rng=rng, quiet=quiet, ...)
  }

MonteCarloSimulations.EstEval <- function(model, simulation.args=NULL,
            replications=100, rng=getRNG(model),  quiet=FALSE, ...)
     MonteCarloSimulations(TSmodel(model), simulation.args=simulation.args,
         replications=replications, rng=rng,  quiet=quiet, ...)
# this looks like a candidate for NextMethod

MonteCarloSimulations.MonteCarloSimulations <- function(model, simulation.args=NULL,
            replications=100, rng=getRNG(model),  quiet=FALSE, ...)
     MonteCarloSimulations(TSmodel(model), simulation.args=simulation.args,
         replications=replications, rng=rng,  quiet=quiet, ...)
# this looks like a candidate for NextMethod


MonteCarloSimulations.TSmodel <- function(model, simulation.args=NULL,
          replications=100, rng=NULL, quiet=FALSE, ...)
#	  Spawn=if (exists(".SPAWN")) .SPAWN else FALSE, ...)
{ 
 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
 arglist <- append(list(model), simulation.args)
# if (Spawn)
#  {if (!quiet)cat("Spawning processes to calculate ", replications, " replications.\n")
#   assign("sim.forloop.n", replications, where = 1)
#   assign("sim.forloop.result", list(NULL), where = 1)
# #  assign("sim.forloop.model", model, where = 1)
#   assign("sim.forloop.arglist", arglist, where = 1)
#   on.exit(remove(c("sim.forloop.i", "sim.forloop.n", "sim.forloop.result",
#       "sim.forloop.arglist"),where = 1))
#   For(sim.forloop.i = 1:sim.forloop.n, sim.forloop.result[[sim.forloop.i]] <- 
#       do.call("simulate",  sim.forloop.arglist),
#       first=options(warn=-1), sync = TRUE)
#   result <- array(NA, c(dim(sim.forloop.result[[1]]$output),replications))
#   tfr <- tframe(sim.forloop.result[[1]]$output)
#   for (i in 1:replications) result[,,i] <- sim.forloop.result[[i]]$output
#  }
# else {
   #r <- simulate(model, list(...))$output
   r <- do.call("simulate", arglist)$output
   result <- array(NA, c(dim(r),replications))
   tfr <- tframe(r)
   result[,,1] <- r
   if (1 < replications)
     for (i in 2:replications) 
        result[,,i] <- outputData(do.call("simulate", arglist))
#  }
# default does not work on array seriesNames(result) <- seriesNamesOutput(model)
if (length(seriesNamesOutput(model)) != dim(result)[2])
 stop("length of names (",length(seriesNamesOutput(model)),
      ") does not match number of series(",dim(result)[2],").")
attr(result,"seriesNames") <- seriesNamesOutput(model)

result <- tframed(result, tfr)  # my more general multidimensional ts
invisible(classed( # constructor MonteCarloSimulations
         list(simulations=result, model=model, rng=rng, simulation.args=simulation.args,
              description = "data generated by MonteCarloSimulations.TSmodel"),
   c("MonteCarloSimulations") ))
}

is.MonteCarloSimulations <- function(obj) 
   {inherits(obj,"MonteCarloSimulations")}

print.MonteCarloSimulations <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
 cat("Simulation with RNG ", x$rng, " from model:\n")
 print(x$model)
 invisible(x)
}

nseriesOutput.MonteCarloSimulations <- function(x)
   {dim(x$simulations)[2]}

nseriesInput.MonteCarloSimulations <- function(x)
  {nseriesInput(x$simulation.args$data)}

tframe.MonteCarloSimulations <- function(x) tframe(x$simulations)

Tobs.MonteCarloSimulations <- function(x) Tobs(tframe(x))

seriesNamesOutput.MonteCarloSimulations <- function(x)
   {dimnames(x$simulations)[[2]]}

seriesNamesInput.MonteCarloSimulations <- function(x)
  {seriesNamesInput(x$simulation.args$data)}

testEqual.MonteCarloSimulations <- function(obj1, obj2, fuzz=1e-16)
 {if (length(obj1$result) != length(obj2$result)) r <- FALSE
  else  r <- all(fuzz > abs(obj1$simulations - obj2$simulations))
  r
 }


summary.MonteCarloSimulations <- function(object,
        series=NULL, periods=1:3, ...)
 {#  (... further arguments, currently disregarded)
  stats <- NULL
  if (!is.null(series))
    {if (dim(object$simulations)[3] <20) 
        warning("SD calculation is not very good with so few simulations.")
     names <- seriesNamesOutput(object)
     if(is.null(names)) names <- seriesNamesOutput(object$model)
     if (!is.numeric(series)) series <-match(series, names)
     names <- names[series]
     mn<-apply(object$simulations[periods,series,,drop=FALSE],c(1,2),mean)
     sd<-apply(object$simulations[periods,series,,drop=FALSE],c(1,2),var)^0.5
     stats <- rbind(mn,sd) 
     dimnames(stats)<- list(c(paste("mean period", periods), 
                          paste("S.D. period",periods)), names)
    }
  classed(list(  # constructor summary.MonteCarloSimulations
     description=object$description,
     sampleT=dim(object$simulations)[1],
     p=      dim(object$simulations)[2], 
     simulations=dim(object$simulations)[3], 
     summary.stats=stats,
     rng=getRNG(object)), 
  "summary.MonteCarloSimulations")
}


print.summary.MonteCarloSimulations <- function(x, digits=options()$digits, ...)
 {#  (... further arguments, currently disregarded)
  cat("Object class MonteCarloSimulations\n")
  cat(x$description, "\n")
  cat("Tobs=",x$sampleT, "variables=", x$p,"simulations=",x$simulations,"\n")
  cat("rng= ", x$rng, "\n")
  if (!is.null(x$summary.stats))   print(x$summary.stats, digits=digits)
  invisible(x)
 }


tfplot.MonteCarloSimulations <- function(x, 
    tf=tframe(x$simulations), start=tfstart(tf), end=tfend(tf),
    series=seq((dim(x$simulations)[2])), 
    select.simulations=seq(dim(x$simulations)[3]),
    graphs.per.page=5, mar=par()$mar, ...)
  {#  (... further arguments, currently disregarded)
   names <- seriesNames(x$simulations)
   tf.p <- tframe(x$simulations) # actual may be different (but not in default)
   Ngraph <- min(length(series), graphs.per.page)
   old.par <-par(mfcol = c(Ngraph, 1), mar= mar, no.readonly=TRUE) #c(5.1,6.1,4.1,2.1))
   on.exit(par(old.par))
   #zz<- matrix(NA, dim(sim)[1], length(x$simulations))
   if (!is.numeric(series)) series <- match(series, names)
   for(i in series) 
        {zz <- (x$simulations)[,i,select.simulations]
         tframe(zz) <- tf.p
	 seriesNames(zz) <- NULL #otherwise name length does not match in tfwindow.default(zz)
         tfOnePlot(zz,start=start,end=end, ylab=names[i]) #tsplot
         if(i == series[1])  title(main = "Monte Carlo Simulations")
        }
   invisible()
}



distribution <- function(obj, ...)UseMethod("distribution")

#distribution.TSdata <- function(obj, bandwidth=0.2, 
#        select.inputs = seq(length= nseriesInput(obj)),
#        select.outputs= seq(length=nseriesOutput(obj)), ...)
distribution.TSdata <- function(obj, ..., bandwidth=0.2, 
        select.inputs = seq(length= nseriesInput(obj)),
        select.outputs= seq(length=nseriesOutput(obj)))
  {#  (... further objects, currently disregarded)
   if (0 !=  nseriesInput(obj))
      distribution( inputData(obj), bandwidth=bandwidth, series=select.inputs)
   if (0 != nseriesOutput(obj))
      distribution(outputData(obj), bandwidth=bandwidth, series=select.outputs)
   invisible(obj)
  }

# this should be something like distribution.tframed if distribution becomes
#  generic in EstEval.  ALSO SEE distribution.factorsEstEval RE ...
#distribution.default <- function(obj, bandwidth=0.2, series=NULL, ...)
distribution.default <- function(obj, ..., bandwidth=0.2, series=NULL)
  {#  (... further objects, currently disregarded)
   # obj should be a ts matrix (perhaps this should be a tf method).
   # If series is NULL then all series are ploted.
   # note that this graphic can be fairly misleading:
   #    distribution(runif(1000))  should be uniform
   names <- seriesNames(obj)
   if (!is.matrix(obj) ) obj <- matrix(obj, length(obj), 1)
   if(!is.null(series)) 
     {obj   <-   obj[,series, drop=FALSE]
      names <- names[series]
     }
   par(mfcol=c(ncol(obj),1))
   for ( i in 1:ncol(obj))
      {if      (exists("density")) rd <- density(obj[,i], bw= bandwidth)
       else if (exists("ksmooth") & !is.R()) rd <- ksmooth(obj[,i], bandwidth=bandwidth) 
       else     stop("Neither ksmooth nor density are available.")
       plot(rd, type="l", ylab="density", ylim=c(0, max(rd$y)), xlab=names[i] )
      }
   invisible()
  }


distribution.MonteCarloSimulations <- function(obj,
     series=seq(dim(obj$simulations)[2]),
     x.sections=TRUE, periods=1:3, graphs.per.page=5, ...)
  {#  (... further arguments, currently disregarded)
  if (dim(obj$simulations)[3] <20) 
     warning("This is not very good with so few simulations.")
   names <- seriesNamesOutput(obj)
   if(is.null(names)) names <- seriesNamesOutput(obj$model)
   if (!is.numeric(series)) series <- match(series, names)
   names <- names[series]
   Ngraph <- min(length(series), graphs.per.page)
   if (x.sections)
       {data <- obj$simulations[periods, series,,drop=FALSE]
        old.par <-par(mfrow =c(Ngraph, length(periods)),
	              mar=c(5.1,6.1,4.1,2.1), no.readonly=TRUE)
       }
   else 
       {old.par <-par(mfrow =c(Ngraph, 1), mar=c(5.1,6.1,4.1,2.1), no.readonly=TRUE)
        mn <- apply(obj$simulations[, series,,drop=FALSE], c(1,2), mean)
        sd <- apply(obj$simulations[, series,,drop=FALSE], c(1,2), var) ^ 0.5
        plt <- array(c(mn, mn+sd, mn-sd, mn+2*sd, mn-2*sd), c(dim(mn),5))
        tf.p <- tframe(obj$simulations)
      }
   on.exit(par(old.par))
   for (i in 1:length(series)) 
    {if (x.sections)
        {for (j in 1:length(periods)) 
          {if (exists("density")) plot(density(data[j,i,]), # -mean[j,i]
                 type="l",xlab=paste("Period",periods[j]),ylab=names[i])
           else if (exists("ksmooth") & !is.R()) plot(ksmooth(data[j,i,], # -mean[j,i],
                              bandwidth=var(data[j,i,])^0.5, kernel="parzen"),
                        type="l",xlab=paste("Period",periods[j]),ylab=names[i])
           else
        stop("Neither ksmooth nor density are available to calculate the plot.")
           if ((i == 1) & (j ==length(periods)%/%2))
              title(main = "kernel estimate of distributions")
           }
        }
     else
        {pl <-plt[,i,]
         tframe(pl) <- tf.p
         tfplot(pl, type="l", lty= c(1,3,3,2,2), ylab=names[i]) #tsplot
         if (i == 1) title(main = "Simulation mean, 1 & 2 S.D. estimates")
        }
    }
  invisible()
  }

############################################################################
#
#       methods for EstEval.  <<<<<<<<<<
#
############################################################################

#e.bb.ar.100 <- EstEval( mod2, replications=100, 
#               estimation.args=list(estimation="estVARXar", verbose=F))

#e.bb.ls.over <- EstEval( simple.mod, replications=100, 
#   estimation.args=list(estimation="estVARXls", max.lag=6, verbose=F), 
#   criterion="coef")


EstEval <- function( model, replications=100, rng=NULL, quiet=FALSE, 
                       simulation.args=NULL,
                       estimation=NULL, estimation.args=NULL, 
                       criterion ="coef", criterion.args =NULL) 
#		       Spawn=if (exists(".SPAWN")) .SPAWN else FALSE)
{
 if(is.null(estimation)) stop("estimation method must be specified.")
 if (is.EstEval(model) | is.MonteCarloSimulations(model))
   {rng  <- getRNG(model)
    model<- TSmodel(model)
   }
  
 truth <- do.call(criterion, append(list(model), criterion.args))

 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
# if (!Spawn) {
    if(!quiet) cat("Calculating ", replications, " estimates.\n")
    result <- vector("list",replications)
    for (i in 1:replications)
       {data <- do.call("simulate", append(list(model), simulation.args))
	m   <-  do.call(estimation, append(list(data),  estimation.args))
        result[[i]]<-do.call(criterion, append(list(m), criterion.args))
       }
#   }
# else {
#     if(!quiet)
#	 cat("Spawning processes to calculate ", replications, " estimates.\n")
#     est.forloop <- function(estimation, estimation.args, model, 
#			      simulation.args, criterion, criterion.args)
#	{data <- do.call("simulate", append(list(model), simulation.args))
#	 m   <-  do.call(estimation, append(list(data),  estimation.args))
#	 do.call(criterion, append(list(m), criterion.args))
#	}
#     assign("est.forloop", est.forloop, where = 1)
#     assign("est.forloop.n", replications, where = 1)
#     assign("est.forloop.result", list(NULL), where = 1)
#     assign("est.forloop.estimation", estimation, where = 1)
#     assign("est.forloop.model", model, where = 1)
#     assign("est.forloop.simulation.args", simulation.args, where = 1)
#     assign("est.forloop.criterion", criterion, where = 1)
#     assign("est.forloop.estimation.args", estimation.args, where = 1)
#     assign("est.forloop.criterion.args", criterion.args, where = 1)
#
#     on.exit(remove(c("est.forloop", "est.forloop.i", "est.forloop.n",
#	  "est.forloop.result",  "est.forloop.estimation","est.forloop.model",
#	  "est.forloop.simulation.args", "est.forloop.criterion", 
#	  "est.forloop.estimation.args","est.forloop.criterion.args"),where = 1))
#
#     For(est.forloop.i = 1:est.forloop.n, est.forloop.result[[est.forloop.i ]]<-
#	est.forloop(est.forloop.estimation, est.forloop.estimation.args, est.forloop.model, 
#	est.forloop.simulation.args, est.forloop.criterion, est.forloop.criterion.args),
#	first=options(warn=-1), sync = TRUE)
#     result<-est.forloop.result
#    }
invisible(classed( # constructor EstEval (EstEvals)
      list(result=result,truth=truth,model=model,
           rng=rng, version=version,
           estimation=estimation, estimation.args=estimation.args,
            criterion=criterion,   criterion.args=criterion.args, 
            simulation.args=simulation.args),
     c(paste(criterion,"EstEval",sep=""), "EstEval")))
}


is.EstEval <- function(obj){inherits(obj,"EstEval")}

testEqual.EstEval <- function(obj1, obj2, fuzz = 0)
 {all(as.character(obj1) == as.character(obj2))}

print.EstEval <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
 cat("Estimation evaluation with model:\n")
 print(x$model, digits=digits)
 cat("Evaluation criterion: ",x$criterion, "\n")
 invisible(x)
}

summary.EstEval <- function(object, ...)
 {#  (... further arguments, currently disregarded)
  classed(list( # constructor summary.EstEval
     class=class(object),
     estimation=object$estimation,
     estimation.args= if(!is.list((object$estimation.args)[[1]]))
                           object$estimation.args    else    NULL,
     criterion=object$criterion,
     labels=object$labels,
     criterion.args=object$criterion.args,
     replications=length(object$result), 
     true.model=object$model,
     rng=getRNG(object)), 
  "summary.EstEval")
 }


print.summary.EstEval <- function(x, digits=options()$digits, ...)
{ #  (... further arguments, currently disregarded)
  cat("Object of class: ", x$class, "\n")
  cat("Evaluation of `",x$estimation,"'")
  if(!is.list((x$estimation.args)[[1]]))
    {cat( " estimation with argument ") 
     cat(labels(x$estimation.args),"= `",x$estimation.args,"'")
    }
  cat("\n")
  cat("using criterion `", x$criterion, "' with argument ")
  cat(x$labels," = `", x$criterion.args, "'\n")
  cat(x$replications, " replications, RNG = ", x$rng, "\n")
  cat("true model:\n")
  print(x$model)
  invisible(x)
 }

tfplot.EstEval <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf),
        truth= if(is.TSdata(x$truth)) outputData(x$truth) else x$truth,
        series = seq(length=nseries(truth)),
	Title="Estimated (and true) results",
        ylab = seriesNames(truth), remove.mean=FALSE,
	graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
 {#  (... further arguments, currently disregarded)
  N<-length(x$result)
  old.par <- par(mfcol = c(min(length(series), graphs.per.page), 1),
           no.readonly=TRUE)
  on.exit(par(old.par))
  tf <- tframe(truth)
  truth <- unclass(truth)
  if (remove.mean) 
    {truth <- sweep(truth,2, colMeans(truth))
     ylab <- paste(ylab, "(mean removed)")
    }
  r <- matrix(NA,Tobs(truth), N)
  for (j in series) {
  	for (i in 1:N)  r[,i] <- x$result[[i]][,j]
	if (remove.mean) r <- sweep(r,2, colMeans(r))
  	tfOnePlot(tframed(tbind(truth[,j], r), tf), ylab= ylab[j],
	              lty=c(1,rep(2, N)), col=c("black",rep("red",N)),
		      start=start,end=end)
        if(!is.null(Title) && (j == series[1]) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = Title)
        }
  invisible()
 }



############################################################################
#
#       methods for rootsEstEval  (EstEval)  <<<<<<<<<<
#
############################################################################

summary.rootsEstEval <- function(object, verbose=TRUE, ...)
{#  (... further arguments, currently disregarded)
  nxt <- if (verbose) NextMethod("summary") else NULL
  if (! verbose) conv <- NULL
  else 
    {if (!is.null(object$result.conv)) conv <- NULL
     else                              conv <- object$result.conv
    }
  N <- length(object$result)
  p <- 0
  for (i in 1:N) p <- max(p, length((object$result)[[i]]))
  r <- matrix(NA, N, p)
  for (i in 1:N) r[i,1:length((object$result)[[i]])] <- (object$result)[[i]]
  m <- colSums(r)/N
  cov <- r- t(matrix(object$truth, p, N))
  cov <- (t(Conj(cov)) %*% cov)/(N-1)
  ecov <- r- t(matrix(m, p, N))
  ecov <- (t(Conj(ecov)) %*% ecov)/(N-1)
  classed(list( # constructor summary.summary.rootsEstEval
     nxt=nxt,
     conv=conv,
     true.criterion=object$truth,
     mean=m,
     cov=cov,
     ecov=ecov),
  "summary.summary.rootsEstEval")
 }

print.summary.rootsEstEval <- function(x, digits=options()$digits, ...)
 {#  (... further arguments, currently disregarded)
  if (!is.null(x$nxt)) print(x$nxt)
  if (!is.null(x$conv))
        {if(all(x$conv)) cat("All estimates converged.\n")
         else cat(sum(!x$conv)," estimates did not converge!\n")
        }
  cat("\nTrue model criterion mean: ",x$true.criterion,"\n")
  cat("Sampling estimate of mean: ",x$mean,"\n")
  cat("Estimate of sampling covariance [e*Conj(t(e))] using true model:\n")
  print(x$cov)
  cat("\nEstimate of sampling covariance (using sample mean and not the true model):\n")
  print(x$ecov)
  invisible(x)
 }


tfplot.rootsEstEval <- function(x, ...)UseMethod("plot.rootsEstEval")

plot.rootsEstEval <- function(x, complex.plane=TRUE, cumulate=TRUE, norm=FALSE,
     bounds=TRUE, transform=NULL, invert=FALSE, Sort=TRUE, ...)
{#  (... further arguments, currently disregarded)
   N<-length(x$result)
   n <- 0
   for (i in 1:N) n <- max(n, length((x$result)[[i]]))
   r <- matrix(0,N, n) 
   for (i in 1:N) r[i,1:length((x$result)[[i]])] <- (x$result)[[i]]
   true.lines <- c(x$truth, rep(0,n-length(x$truth)))
   if (invert)
         {true.lines <- 1/true.lines
          r <- 1/r
         }
   if(!is.null(transform)) 
         {r <- do.call(transform,list(r))
          true.lines <-do.call(transform,list(true.lines))
         }
   if (complex.plane)
    {#plot.roots(x$truth, pch="o")
     plot(x$truth, pch="o")
     for (i in 1:N) addPlotRoots(r[i,], pch="*") 
     addPlotRoots(0, pch="+") # addPlotRoots(0+0i, pch="+")
    }
  else
     {if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      if (cumulate) r <- apply(r,2,cumsum)/matrix(1:N,N,ncol(r))
      else     r <- r 
      if(is.complex(r))
         { r <- cbind(Re(r), Im(r))
          true.lines <-c(Re(true.lines),Im(true.lines))
         }
      r[is.infinite(r)] <- 0
      true.lines <-t(matrix(true.lines, length(true.lines),N))
      matplot(x=seq(nrow(r)), y=cbind(0,true.lines, r), type="l",
             lty=c(1,rep(3,dim(true.lines)[2]), rep(2,dim(r)[2])) )
     }
  invisible(r)
}

roots.rootsEstEval <- function(obj, ...)   {obj}

distribution.rootsEstEval <- function(obj, ..., mod=TRUE, invert=FALSE, Sort=FALSE, 
    bandwidth=0.2, select=NULL)
{#  (... further objectss, currently disregarded)
 # if mod is true the modulus is used, otherwise real and imaginary are separated.
 # if invert is true the reciprical is used.
 # if Sort is true then sort is applied (before cumulate). This is of particular interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
 # If select is not NULL then only the indicated roots are plotted. 
 #     ie - select=c(1,2)  will plot only the two largest roots
      N<-length(obj$result)
      n <- 0
      for (i in 1:N) n <- max(n, length((obj$result)[[i]]))
      r <- matrix(0,N,n)
      for (i in 1:N) r[i,] <- c((obj$result)[[i]], 
                                rep(0,n-length((obj$result)[[i]])))
      true.lines <- c(obj$truth, rep(0,n-length(obj$truth)))
      if (invert)
         {true.lines <- 1/true.lines
          r <- 1/r
         }
      if(mod) 
         {r <- Mod(r) 
          true.lines <-Mod(true.lines)
          xlab <-"Mod root "
         }
      else
         { r <- cbind(Re(r), Im(r))
          true.lines <-c(Re(true.lines),Im(true.lines))
          xlab <-"Real part root "
         }
      r[is.infinite(r)] <- 0
      if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      if(!is.null(select)) r <- r[,select, drop=FALSE]
      par(mfcol=c(dim(r)[2],1))
      for ( i in 1:dim(r)[2])
         {if      (exists("ksmooth") & !is.R()) rd <- ksmooth(r[,i], bandwidth=bandwidth) 
          else if (exists("density")) rd <- density(r[,i], bw= bandwidth)
          else
        stop("Neither ksmooth nor density are available to calculate the plot.")
          if (i > n) xlab <-"Imaginary part root "
          plot(rd, type="l", ylab="density", ylim=c(0, max(rd$y)),
               xlab=paste(xlab, n-(-i%%n)) )
          lines(rep(true.lines[i],2),c(1,0))
         }
      invisible()
}


############################################################################
#
#       methods for coefEstEval (EstEval)  <<<<<<<<<<
#
############################################################################

summary.coefEstEval <- function(object, verbose=TRUE, ...)
  {classed(summary.rootsEstEval(object, verbose=verbose), "summary.coefEstEval")} # constructor
#  (... further arguments, currently disregarded)


print.summary.coefEstEval <- function(x, digits=options()$digits, ...)
  UseMethod("print.summary.rootsEstEval")
#  (... further arguments, currently disregarded)


tfplot.coefEstEval <- function(x, cumulate=TRUE, norm=FALSE, bounds=TRUE,
        invert=FALSE, Sort=FALSE, graphs.per.page = 5, ...){
      #  (... further arguments, currently disregarded)
      N<-length(x$result)
      n <- 0
      for (i in 1:N) n <- max(n, length((x$result)[[i]]))
      r <- matrix(0,N,n)
      plottrue <- TRUE
      for (i in 1:N) {
	ni <- length((x$result)[[i]])
	r[i,1:ni] <- (x$result)[[i]]
	if (ni != n) plottrue <- FALSE
	}
      if (invert) r <- 1/r
      if(norm)    r <- matrix((rowSums(r^2))^.5, N,1)
      r[is.infinite(r)] <- 0
      if (Sort) r <- t(apply(r,1,sort))
      if (n != length(x$truth)) plottrue <- FALSE
      if (plottrue) {
	true.lines <- c(x$truth)
	if (invert) true.lines <- 1/true.lines
	if(norm) true.lines <-sum(true.lines^2)^.5
	if (Sort) true.lines <- sort(true.lines)
	true.lines <-t(matrix(true.lines, length(true.lines),N))
	if (bounds) {
		z  <- r-true.lines
		#Om <- t(z) %*% z/(nrow(z)-1)
		Om <- crossprod(z, z)/(nrow(z)-1)
		Om <- diag(Om)^.5
		Om <- t(matrix(Om, length(Om), N))
		Om <- Om/matrix((1:N)^.5 , N, ncol(Om))
        }
      }
      if (cumulate) r<- apply(r,2,cumsum)/matrix(1:N,N,ncol(r))
      seriesNames(r) <- paste("parm", seq(ncol(r)))
#      matplot(x=matrix(seq(nrow(r)),nrow(r),1), y=cbind(0,true.lines,r, Om), 
#              type="l", lty=c(1,rep(3,dim(true.lines)[2]), rep(4,dim(r)[2]), 
#                     rep(2,2*dim(r)[2])))
      if (plottrue & bounds)
         tfplot(r, true.lines,true.lines + Om, true.lines - Om,
            graphs.per.page = graphs.per.page)
      else if (plottrue) tfplot(r, true.lines, graphs.per.page=graphs.per.page)
      else               tfplot(r, graphs.per.page = graphs.per.page)
      invisible(r)
}

roots.coefEstEval <- function(obj, criterion.args=NULL, ...)
{# extract roots criterion 
  model <- obj$model
  truth <-do.call("roots", append(list(model), criterion.args))
  r <- NULL
  for (m in obj$result)
    {coef(model) <- m
     model <- setArrays(model)
     r <- append(r, 
           list(do.call("roots", append(list(model), criterion.args))))
    }
  ok <- TRUE
  for (m in obj$result)
     ok <- ok & (length(coef(model)) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"roots"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("rootsEstEval","EstEval")))# constructor
}

distribution.coefEstEval <- function(obj, ...,  Sort=FALSE, bandwidth=0.2,
	graphs.per.page=5)
{#  (... further objects, currently disregarded)
 # if Sort is true then sort is applied (before ave). This is of particular interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
      N<-length(obj$result)
      n <- length(obj$truth)
      for (i in 1:N) n <- max(n, length((obj$result)[[i]]))
      if (n == length(obj$truth)) plottrue <- TRUE
      else {
         warning("Number of true parameters does not match number estimated.")
	 plottrue <- FALSE
      }
      r <- matrix(0,N,n)
      for (i in 1:N) r[i,1:length((obj$result)[[i]])] <- (obj$result)[[i]]
      true.lines <- c(obj$truth, rep(0,n-length(obj$truth)))
      if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      xlab <-"parameter "
      old.par <- par(mfcol=c(min(graphs.per.page, ncol(r)),1),
                    mar=c(5.1, 4.1, 4.1, 2.1), no.readonly=TRUE)
      on.exit(par(old.par))
      for ( i in 1:ncol(r))
         {if (is.R())     rd <- density(r[,i], bw=bandwidth)
          else            rd <- ksmooth(r[,i], bandwidth=bandwidth) # Splus
          plot(rd, type="l", ylim=c(0, max(rd$y)),
               ylab="density",  xlab=paste(xlab, i) , main="")
          if (plottrue) lines(rep(true.lines[i],2),c(1,0))
         }
      invisible()
}

TSmodel.coefEstEval <- function(obj, ...)
{# rebuild model from coef
  model <- obj$model
  truth <-TSmodel(model)
  r <- NULL
  for (m in obj$result)
    {coef(model) <- m
     model <- setArrays(model)
     r <- append(r, list(model))
    }
  ok <- TRUE
  for (m in obj$result)
     ok <- ok & (length(coef(model)) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"TSmodel"
  obj$criterion.args <-NULL
  invisible(classed(obj, c("TSmodelEstEval","EstEval")))
}

TSestModel.coefEstEval <- function(obj)
{# rebuild ... 
  model <- obj$model
  truth <-l(TSmodel(model), data)   # need to regenerate data
  r <- NULL
  for (m in obj$result)
    {coef(model) <- m
     model <- l( setArrays(model), data)
     r <- append(r, list(model))
    }
  ok <- TRUE
  for (m in obj$result)
     ok <- ok & (length(coef(model)) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"TSestModel"
  obj$criterion.args <-NULL
  invisible(classed(obj, c("TSestModelEstEval","EstEval")))
}

############################################################################
#
#       methods for TSmodelEstEval (EstEval)  <<<<<<<<<<
#
############################################################################

summary.TSmodelEstEval <- function(object, ...)
 {#  (... further arguments, currently disregarded)
  if (is.null((object$result)[[1]]$converged)) conv <- NULL
  else
    {conv <- rep(NA,length(object$result))
     for (i in 1:length(conv)) conv[i] <- (object$result)[[i]]$converged
    }
 # summary(coef(object))
 #summary(roots(object))  these are slow

  classed(list( # constructor summary.TSmodelEstEval
     class=class(object),
     conv=conv,
     default=summary.default(object)),
  "summary.TSmodelEstEval")
 }

print.summary.TSmodelEstEval <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
  cat("Object of class: ",x$class, "\n")
  if (!is.null(x$conv))
        {if(all(x$conv)) cat("All estimates converged.\n")
         else cat(sum(!x$conv)," estimates did not converge!\n")
        }
  print(x$default)
  invisible(x)
}


coef.TSmodelEstEval <- function(object, criterion.args=NULL, ...)
{#  (... further arguments, currently disregarded)
 # extract parameters from models in the list and 
 #   return a list of class coefEstEval EstEval 
 # criterion.args is not used. It is provided only so calls from 
 #   summary.TSmodelEstEval can provide this argument.
  truth <-coef(object$truth)
  r <- NULL
  for (m in object$result) 
     r <- append(r,list(coef(m)))
  if (! is.null((object$result)[[1]]$converged))
    {rc <- rep(NA,length(object$result))
     for (i in 1:length(rc)) rc[i] <- (object$result)[[i]]$converged
     object$result.conv<-rc
    }
  object$result<-r
  object$truth <-truth
  object$criterion<-"coef"
  object$criterion.args <-criterion.args
  invisible(classed(object, c("coefEstEval","EstEval")))
}

roots.TSmodelEstEval <- function(obj, criterion.args=list( randomize=TRUE), ...)
{# extract roots criterion 
  truth <-do.call("roots", append(list(obj$truth), criterion.args))
  r <- NULL
  for (m in obj$result)
     r <- append(r, 
           list(do.call("roots", append(list(m), criterion.args))))
  if (! is.null((obj$result)[[1]]$converged))
    {rc <- rep(NA,length(obj$result))
     for (i in 1:length(rc)) rc[i] <- (obj$result)[[i]]$converged
     obj$result.conv<-rc
    }
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"roots"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("rootsEstEval","EstEval")))
}


tfplot.TSmodelEstEval <- function(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL, ...){
  #  (... further arguments, currently disregarded)
  # extract criterion and pass to another method with graph.args
  r <- do.call(paste(criterion,".TSmodelEstEval", sep=""), 
               append(list(x), list(criterion.args=criterion.args)))
  do.call("tfplot", append(list(r), graph.args))
  invisible(r)
  }

############################################################################
#
#       methods for TSestModelEstEval (EstEval)   <<<<<<<<<<
#
############################################################################

summary.TSestModelEstEval <- function(object, ...)
  {classed(summary.TSmodelEstEval(object), "summary.TSestModelEstEval") }  # constructor 
#  (... further arguments, currently disregarded)

print.summary.TSestModelEstEval <- function(x, digits=options()$digits, ...)
   UseMethod("print.summary.TSmodelEstEval")
#  (... further arguments, currently disregarded)


coef.TSestModelEstEval <- function(object, criterion.args=NULL, ...)
{#  (... further arguments, currently disregarded)
 # extract parameters from models in the list and convergence info.
 #   return a list of class coefEstEval EstEval
 # criterion.args is not used. It is provided only so calls from 
 #   summary.TSmodelEstEval can provide this argument.
  truth <-coef(object$truth)
  r <- NULL
  for (m in object$result) r <- append(r,list(coef(m)))
  rc <- rep(NA,length(object$result))
  for (i in 1:length(rc)) rc[i] <- (object$result)[[i]]$converged
  object$result<-r
  object$result.conv<-rc
  object$truth <-truth
  object$criterion<-"coef"
  object$criterion.args <-criterion.args
  invisible(classed(object, c("coefEstEval","EstEval")))
}

roots.TSestModelEstEval <- function(obj, criterion.args=NULL, ...)
{# extract roots criterion 
  truth <-do.call("roots", append(list(obj$truth), criterion.args))
  r <- NULL
  for (m in obj$result)
     r <- append(r, 
             list(do.call("roots", append(list(m), criterion.args))))
  rc <- rep(NA,length(obj$result))
  for (i in 1:length(rc)) rc[i] <- (obj$result)[[i]]$converged
  obj$result<-r
  obj$result.conv<-rc
  obj$truth <-truth
  obj$criterion<-"roots"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("rootsEstEval","EstEval")))
}


tfplot.TSestModelEstEval <- function(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL, ...){
  #  (... further arguments, currently disregarded)
  # extract criterion and pass to another method with graph.args
  r <- do.call(paste(criterion,".TSestModelEstEval", sep=""), 
               append(list(x), list(criterion.args=criterion.args)))
  do.call("tfplot", append(list(r), graph.args))
  invisible(r)
  }


############################################################################
#
#  methods for generating test data
#
############################################################################

genMineData <- function(umodel, ymodel, uinput=NULL, sampleT=100, 
	unoise=NULL, usd=1,ynoise=NULL, ysd=1, rng=NULL)
{if (is.TSestModel(umodel)) umodel <- TSmodel(umodel)
 if (is.TSestModel(ymodel)) ymodel <- TSmodel(ymodel)
 if(!is.TSmodel(umodel)) stop("genMineData expecting a TSmodel.")
 if(!is.TSmodel(ymodel)) stop("genMineData expecting a TSmodel.")

 if (nseriesInput(ymodel) != nseriesOutput(umodel))
   stop("umodel output dimension must equal ymodel input dimension.")
 
 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
 input <- outputData(simulate(umodel, input=uinput, sampleT=sampleT, 
                   noise=unoise, sd=usd))
 r <- TSdata(input  = input,
             output = outputData(simulate(ymodel, input=input,
	                             sampleT=sampleT, noise=ynoise, sd=ysd)) )
 r$umodel  <- umodel
 r$ymodel  <- ymodel
 r$uinput  <- uinput
 r$sampleT <- sampleT 
 r$unoise  <- unoise
 r$usd     <- usd
 r$ynoise  <- ynoise
 r$ysd     <- ysd
 r$rng     <- rng
 r
}

build.input.models <- function(data, max.lag=NULL)
{# make a list of univariate models, one for each series in inputData(data)
 #   for use by build.diagonal.model. 
 n <- nseriesInput(data)
 multi.models <- vector("list", n)
 for (i in seq(n))
   {d <-trimNA(TSdata(output= inputData(data, series=i)))
    multi.models[[i]] <- TSmodel(estVARXls(d, max.lag=max.lag))
   }
 multi.models
}

build.diagonal.model <- function(multi.models)
{# build one diagonal model from a list of models as returned  by 
 # build.input.models. Uses the AR part only. This can be used by genMineData.
 n <- length(multi.models)
 lag <- 0
 for (i in seq(n)) lag <- max(lag, dim(multi.models[[i]]$A)[1])
 p <- 0
 for (i in seq(n))  p  <- p + dim(multi.models[[i]]$A)[3]
 model <- array(0, c(lag, p,p))
 p <- 0
 for (i in seq(n))
   {pi <- dim(multi.models[[i]]$A)[3]
    li <- dim(multi.models[[i]]$A)[1]
    model[ 1:li, (p+1):(p+pi), (p+1):(p+pi)] <-  multi.models[[i]]$A
    p <- p + pi
   }
 ARMA(A= model, B=array(diag(1, p), c(1,p,p)))
}
