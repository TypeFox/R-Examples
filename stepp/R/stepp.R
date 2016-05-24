#################################################################
#
# stepp.R
#
#############################
# old stepp interface       #
#############################
# These calls are maintained by backward compatibility reasons:
#
stepp <- function (trttype, coltrt, coltime, colcens=0, coltype=0, colvar, trts, 
			 patspop, minpatspop, timest, noperm)
{
  win.stepp    <- new("stwin", type="sliding", r1=minpatspop, r2=patspop)
  subp.stepp   <- new("stsubpop")
  subp.stepp   <- generate(subp.stepp, win=win.stepp, covariate=colvar)
  if (trttype == "KM")
    mod.stepp  <- new("stmodelKM", coltrt=coltrt, survTime=coltime, censor=colcens, 
				trts=trts, timePoint=timest)
  else
  if (trttype == "CI")
    mod.stepp  <- new("stmodelCI", coltrt=coltrt, coltime=coltime, coltype=coltype, 
				trts=trts, timePoint=timest)
  else
    stop("Treatment type can only be KM or CI !")

  result.stepp <- new("steppes")
  result.stepp <- estimate(result.stepp, subp.stepp, mod.stepp)
  result.stepp <- test(result.stepp, noperm)
  return(result.stepp)
}

stepp_summary <- function(x)
{
  summary(x)
}

stepp_print <- function(x, estimate=TRUE, cov=TRUE, test=TRUE)
{
  print(x@model, x, estimate, cov, test)
}

stepp_plot <- function(x, legendy = 30, pline = -2.5, color = c("red", "black"),
	ylabel= "Specify Timepoint & Endpoint", xlabel="Subpopulations by Median Covariate",
	ncex = 0.7, tlegend=c("Specify 1st Treatment", "Specify 2nd Treatment"), 
	nlas = 0, alpha = 0.05, pointwise = FALSE, diff = TRUE, ci = TRUE, pv = TRUE, 
	showss = TRUE, ylimit=c(0,100,-100,100,0,3), dev="", together=FALSE, noyscale=FALSE, at=NA)
{
  plot(x, legendy=legendy, pline=pline, color=color,
		ylabel=ylabel, xlabel=xlabel, ncex=ncex, tlegend=tlegend, nlas=nlas, alpha=alpha, 
		pointwise=pointwise, diff=diff, ci=ci, pv=pv, showss=showss, ylimit=ylimit, 
		dev=dev, together=together, noyscale=noyscale, at=at)
}

analyze.KM.stepp <- function ( coltrt, coltime, colcens, colvar, trts, patspop, minpatspop, 
				timest, noperm=2500,
			 	ncex = 0.70, legendy = 30, pline = -2.5, color = c("red", "black"),
			 	xlabel="Subpopulations by Median Covariate",
			 	ylabel = "?-year Disease-Free Survival", 
			 	tlegend = c("1st Treatment", "2nd Treatment"),
			 	nlas = 3, pointwise=FALSE)
{
  stepp.KM <- stepp("KM", coltrt=coltrt, coltime=coltime, colcens=colcens, colvar=colvar,
			  trts=trts, patspop=patspop, minpatspop=minpatspop, timest=timest, noperm=noperm)
  stepp_summary(stepp.KM)
  stepp_print(stepp.KM)
  stepp_plot(stepp.KM, ncex=ncex,legendy=legendy, pline=pline, color=color,
		xlabel=xlabel, ylabel=ylabel, tlegend=tlegend, nlas=nlas,
		pointwise=pointwise)
  return(stepp.KM)
}

analyze.CumInc.stepp <- function(coltrt, coltime, coltype, colvar, trts, patspop, minpatspop, 
    				timest, noperm=2500,
				ncex = 0.7, legendy = 30, pline = -2.5, color = c("red", "black"),
    				xlabel = "Subpopulations by Median Covariate",
				ylabel = "?-year Disease-Free Survival", 
    				tlegend = c("1st Treatment", "2nd Treatment"), 
   				nlas = 3, pointwise = FALSE)
{
  stepp.CI <- stepp("CI", coltrt=coltrt, coltime=coltime, coltype=coltype, colvar=colvar,
			  trts=trts, patspop=patspop, minpatspop=minpatspop, timest=timest, noperm=noperm)
  stepp_summary(stepp.CI)
  stepp_print(stepp.CI)
  stepp_plot(stepp.CI, ncex=ncex, legendy=legendy, pline=pline, color=color,
     		xlabel=xlabel, ylabel=ylabel, tlegend=tlegend, nlas=nlas,
		pointwise=pointwise)
  return(stepp.CI)
}


