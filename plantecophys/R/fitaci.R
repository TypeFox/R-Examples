#' Fit the Farquhar-Berry-von Caemmerer model of leaf photosynthesis
#' @description Fits the Farquhar-Berry-von Caemmerer model of photosynthesis to measurements of photosynthesis and intercellular \eqn{CO_2}{CO2} concentration (Ci). Estimates Jmax, Vcmax, Rd and their standard errors. A simple plotting method is also included, as well as the function \code{\link{fitacis}} which quickly fits multiple A-Ci curves. Temperature dependencies are taken into account, see \code{\link{Photosyn}}.
#' @param data Dataframe with Ci, Photo, Tleaf, PPFD (the last two are optional). For \code{fitacis}, also requires a grouping variable.
#' @param varnames List of names of variables in the dataset (see Details).
#' @param Tcorrect If TRUE, Vcmax and Jmax are corrected to 25C. Otherwise, Vcmax and Jmax are estimated at measurement temperature.
#' @param Patm Atmospheric pressure (kPa)
#' @param citransition If provided, fits the Vcmax and Jmax limited regions separately (see Details).
#' @param quiet If TRUE, no messages are written to the screen.
#' @param startValgrid If TRUE (the default), uses a fine grid of starting values to increase the chance of finding a solution.
#' @param algorithm Passed to \code{\link{nls}}, sets the algorithm for finding parameter values.
#' @param useRd If Rd provided in data, and useRd=TRUE (default is FALSE), uses measured Rd in fit. Otherwise it is estimatied from the fit to the A-Ci curve.
#' @param group For batch analysis using \code{fitacis}, the name of the grouping variable in the dataframe.
#' @param object For coef.acifit, and print.acifit, the object returned by \code{fitaci}
#' @param progressbar For \code{fitacis}, whether to display a progress bar (default is TRUE).
#' @param \dots Further arguments passed to \code{\link{Photosyn}}
#' @details Uses non-linear regression to fit an A-Ci curve. No assumptions are made on which part of the curve is Vcmax or Jmax limited. Three parameters are estimated, Jmax, Vcmax and Rd. When \code{Tcorrect=TRUE} (the defualt), Jmax and Vcmax are re-scaled to 25C, using the temperature response parameters provided (but Rd is always at measurement temperature). When \code{Tcorrect=FALSE}, estimates of all parameters are at measurement temperature.
#' 
#' When \code{citransition} is set, it splits the data into a Vcmax-limited (where Ci < citransition), and Jmax-limited region (Ci > citransition). Both parameters are then estimated separately for each region (Rd is estimated only for the Vcmax-limited region). \bold{Note} that the actual transition point as shown in the standard plot of the fitted A-Ci curve may be quite different from that provided, since the fitting method simply decides which part of the dataset to use for which limitation, it does not constrain the actual estimated transition point directly. See the example below.
#' 
#' When plotting the fit, the A-Ci curve is simulated using the \code{\link{Aci}} function, with leaf temperature (Tleaf) and PPFD set to the mean value for the dataset. If PPFD is not provided in the dataset, it is assumed to equal 1800 mu mol m-2 s-1.
#' 
#' Note that atmospheric pressure (Patm) is taken into account, assuming the original data are in molar units (Ci in mu mol mol-1, or ppm). During the fit, Ci is converted to mu bar, and Km and Gammastar are recalculated accounting for the effects of Patm on the partial pressure of oxygen. When plotting the fit, though, molar units are shown on the X-axis. Thus, you should get (nearly) the same fitted curve when Patm was set to a value lower than 100kPa, but the fitted Vcmax and Jmax will be higher. This is because at low Patm, photosynthetic capacity has to be higher to achieve the same measured photosynthesis rate.
#' 
#' Because fitaci returns the fitted nls object (see next section), more details on statistics of the fit can be extracted with standard tools. The Examples below shows the use of the \pkg{nlstools} to extract many details at once.
#' @troubleshooting From time to time, the \code{fitaci} function returns an error, indicating that the A-Ci curve could not be fit. These are usually errors returned by \code{nls}, for example "step size was reduced below ...", or errors about "singular gradients". In the vast majority of cases, this indicates a bad curve that probably could not be fit in any case. Inspect the raw data to check if the curve does not include severe outliers, or large deviations from the expected, typical A-Ci curve. Another solution is to set the \code{ciTransition} parameter directly, which will force \code{fitaci} to fit the Jmax- and Vcmax-limited regions separately, often resulting in a successful fit when \code{fitaci} was not successful.
#' @return A list of class 'acifit', with five components:
#' \describe{
#' \item{df}{A dataframe with the original data, the fitted photosynthetic rate (Amodel), Jmax and Vcmax-limited gross rates (Aj, Ac)}
#' \item{pars}{Contains the parameter estimates and their approximate standard errors}
#' \item{nlsfit}{The object returned by \code{\link{nls}}, and contains more detail on the quality of the fit}
#' \item{Photosyn}{A copy of the \code{\link{Photosyn}} function with the arguments adjusted for the current fit. That is, Vcmax, Jmax and Rd are set to those estimated in the fit, and Tleaf and PPFD are set to the mean value in the dataset.}
#' \item{Ci_transition}{The Ci at which photosynthesis transitions from Vcmax to Jmax limited photosynthesis.}
#' }
#' @examples
#' # Fit an A-Ci curve on a dataframe that contains Ci, Photo and optionally Tleaf and PPFD. 
#' # Here, we use the built-in example dataset 'acidata1'.
#' f <- fitaci(acidata1)
#' 
#' # Note that the default behaviour is to correct Vcmax and Jmax for temperature, 
#' # so the estimated values are at 25C. To turn this off:
#' f2 <- fitaci(acidata1, Tcorrect=FALSE)
#' 
#' # To use different T response parameters (see ?Photosyn),
#' f3 <- fitaci(acidata1, Tcorrect=TRUE, EaV=25000)
#' 
#' # Make a standard plot
#' plot(f)
#' 
#' # Look at a summary of the fit
#' summary(f)
#' 
#' # Extract coefficients only
#' coef(f)
#' 
#' # The object 'f' also contains the original data with predictions.
#' # Here, Amodel are the modelled (fitted) values, Ameas are the measured values.
#' with(f$df, plot(Amodel, Ameas))
#' abline(0,1)
#' 
#' # The fitted values can also be extracted with the fitted() function:
#' fitted(f)
#' 
#' # The non-linear regression (nls) fit is stored as well,
#' summary(f$nlsfit)
#' 
#' # Many more details can be extracted with the nlstools package
#' library(nlstools)
#' overview(f$nlsfit)
#'  
#' # The curve generator is stored as f$Photosyn:
#' # Calculate photosynthesis at some value for Ci, using estimated parameters and mean Tleaf, 
#' # PPFD for the dataset.
#' f$Photosyn(Ci=820)
#' 
#' # Photosynthetic rate at the transition point:
#' f$Photosyn(Ci=f$Ci_transition)$ALEAF
#' 
#' # Set the transition point; this will fit Vcmax and Jmax separately. Note that the *actual* 
#' # transition is quite different from that provided, this is perfectly fine : 
#' # in this case Jmax is estimated from the latter 3 points only (Ci>800), but the actual 
#' # transition point is at ca. 400ppm.
#' g <- fitaci(acidata1, citransition=800)
#' plot(g)
#' g$Ci_transition
#' 
#' # Use measured Rd instead of estimating it from the A-Ci curve. 
#' # The Rd measurement must be added to the dataset used in fitting, 
#' # and you must set useRd=TRUE.
#' acidata1$Rd <- 2
#' f2 <- fitaci(acidata1, useRd=TRUE)
#' f2
#' @export
#' @rdname fitaci
fitaci <- function(data, 
                   varnames=list(ALEAF="Photo", Tleaf="Tleaf", Ci="Ci", PPFD="PARi", Rd="Rd"),
                   Tcorrect=TRUE, 
                   Patm=100,
                   citransition=NULL,
                   quiet=FALSE, startValgrid=TRUE, 
                   algorithm="default", useRd=FALSE,  ...){
  
  # Make sure data is a dataframe; stuff returned by dplyr is no good
  data <- as.data.frame(data)
  
  # Set extra parameters if provided
  m <- as.list(match.call())
  a <- as.list(formals(fitaci))
  f <- names(formals(Photosyn))
  
  extrapars <- setdiff(names(m), c(names(a),""))
  for(i in seq_along(extrapars)){
    if(extrapars[i] %in% f){
      val <- eval(m[[extrapars[i]]])
      formals(Photosyn)[extrapars[i]] <- val
    } else {
      warning("Parameter ", extrapars[i]," not recognized.")
    }
  }
  photpars <- formals(Photosyn)
  removevars <- c("whichA")
  photpars <- photpars[-which(names(photpars) %in% removevars)]
  
  # Check if PAR is provided
  if(!varnames$PPFD %in% names(data)){
    data$PPFD <- 1800
    if(!quiet)warning("PARi not in dataset; assumed PARi = 1800.")
  } else data$PPFD <- data[,varnames$PPFD]
  
  # Check if Tleaf is provided
  if(!varnames$Tleaf %in% names(data)){
    data$Tleaf <- 25
    if(!quiet)warning("Tleaf not in dataset; assumed Tleaf = 25.")
  } else {
    data$Tleaf <- data[,varnames$Tleaf]
  }
  
  # Check if Rd is provided and whether we want to set it as known.
  haveRd <- FALSE
  if(!is.null(varnames$Rd)){ # to avoid breakage with older versions when varnames provided.
    if(varnames$Rd %in% names(data) && useRd){
      
      # Has to be a single unique value for this dataset
      Rd_meas <- data[,varnames$Rd]
      Rd_meas <- unique(Rd_meas)
      if(length(Rd_meas) > 1)
        stop("If Rd provided as measured, it must be a single unique value for an A-Ci curve.")
      
      if(Rd_meas < 0)Rd_meas <- -Rd_meas
      haveRd <- TRUE
      
      if(!is.null(citransition))stop("At the moment cannot provide citransition as well as measured Rd.")
    }
    if(varnames$Rd %in% names(data) && !useRd){
      warning("Rd found in dataset but useRd set to FALSE. Set to TRUE to use measured Rd.")
    }
  }
  
  # Extract Ci and apply pressure correction
  data$Ci_original <- data[,varnames$Ci]
  data$Ci <- data[,varnames$Ci] * Patm/100
  
  data$ALEAF <- data[,varnames$ALEAF]
  
  # Needed to avoid apparent recursion below.
  TcorrectVJ <- Tcorrect
  
  # Wrapper around Photosyn; this wrapper will be sent to nls. 
  acifun_wrap <- function(Ci,..., returnwhat="ALEAF"){
    r <- Photosyn(Ci=Ci,Tcorrect=TcorrectVJ,...)
    if(returnwhat == "ALEAF")return(r$ALEAF)
    if(returnwhat == "Ac")return(r$Ac - r$Rd)
    if(returnwhat == "Aj")return(r$Aj - r$Rd)
  }
  
  # Guess Jmax from max A, T-corrected gammastar
  if(haveRd){
    Rd_guess <- Rd_meas
  } else {
    Rd_guess <- 1.5
  }
  
  maxCi <- max(data$Ci)
  mi <- which.max(data$Ci)
  maxPhoto <- data$ALEAF[mi]
  Tl <- data$Tleaf[mi]
  gammastar <- TGammaStar(Tl,Patm)
  VJ <- (maxPhoto+Rd_guess) / ((maxCi - gammastar) / (maxCi + 2*gammastar))
  Jmax_guess <- VJ*4
  if(Tcorrect){
    Teffect <- TJmax(Tl,  EaJ=39676.89, delsJ=641.3615, EdVJ=200000)
    Jmax_guess <- Jmax_guess / Teffect
  }
  
  # Guess Vcmax, from section of curve that is definitely Vcmax-limited
  dato <- data[data$Ci < 150 & data$Ci > 60 & data$ALEAF > 0,]
  if(nrow(dato) > 0){
    Km <- TKm(dato$Tleaf,Patm)
    gammastar <- TGammaStar(dato$Tleaf,Patm)
    vcmax <- with(dato, (ALEAF + Rd_guess) / ((Ci - gammastar)/(Ci + Km)))
    Vcmax_guess <- median(vcmax)
  } else {
    Vcmax_guess <- Jmax_guess/1.8 
  }
  if(Tcorrect){
    Teffect <- TVcmax(Tl, EaV=82620.87, delsC=645.1013, EdVC=0)
    Vcmax_guess <- Vcmax_guess / Teffect
  }
  
  # Fine-tune starting values; try grid of values around initial estimates.
  if(startValgrid){
    if(!haveRd){
      aciSS <- function(Vcmax, Jmax, Rd){
        Photo_mod <- acifun_wrap(data$Ci, PPFD=data$PPFD, 
                                 Vcmax=Vcmax, Jmax=Jmax, 
                                 Rd=Rd, Tleaf=data$Tleaf, Patm=Patm)
        
        SS <- sum((data$ALEAF - Photo_mod)^2)
        return(SS)
      }
      d <- 0.3
      n <- 20
      gg <- expand.grid(Vcmax=seq(Vcmax_guess*(1-d),Vcmax_guess*(1+d),length=n),
                        Rd=seq(Rd_guess*(1-d),Rd_guess*(1+d),length=n))
    
      m <- with(gg, mapply(aciSS, Vcmax=Vcmax, Jmax=Jmax_guess, Rd=Rd))
      ii <- which.min(m)
      Vcmax_guess <- gg$Vcmax[ii]
      Rd_guess <- gg$Rd[ii]
    } else {
      
      aciSS <- function(Vcmax, Jmax, Rd){
        Photo_mod <- acifun_wrap(data$Ci, PPFD=data$PPFD, 
                                 Vcmax=Vcmax, Jmax=Jmax, 
                                 Rd=Rd, Tleaf=data$Tleaf, Patm=Patm)
        SS <- sum((data$ALEAF - Photo_mod)^2)
        return(SS)
      }
      d <- 0.3
      n <- 20
      gg <- data.frame(Vcmax=seq(Vcmax_guess*(1-d),Vcmax_guess*(1+d),length=n))
      
      m <- with(gg, mapply(aciSS, Vcmax=Vcmax, Jmax=Jmax_guess, Rd=Rd_meas))
      ii <- which.min(m)
      Vcmax_guess <- gg$Vcmax[ii]
      # Rd_guess already set to Rd_meas above
      
    }
    
  }
  
  # Fit curve.
  if(is.null(citransition)){
    
    if(!haveRd){
      nlsfit <- nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, 
                                        Jmax=Jmax, Rd=Rd, Tleaf=Tleaf, Patm=Patm),
                      algorithm=algorithm,
                      data=data, control=nls.control(maxiter=500, minFactor=1/10000),
                      start=list(Vcmax=Vcmax_guess, Jmax=Jmax_guess, Rd=Rd_guess))
      p <- coef(nlsfit)
      pars <- summary(nlsfit)$coefficients[,1:2]
    } else {
      
      nlsfit <- nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, 
                                        Jmax=Jmax, Rd=Rd_meas, Tleaf=Tleaf, Patm=Patm),
                    algorithm=algorithm,
                    data=data, control=nls.control(maxiter=500, minFactor=1/10000),
                    start=list(Vcmax=Vcmax_guess, Jmax=Jmax_guess))
      p <- coef(nlsfit)
      p[[3]] <- Rd_meas
      names(p)[3] <- "Rd"
      
      pars <- summary(nlsfit)$coefficients[,1:2]
      pars <- rbind(pars, c(Rd_meas, NA))
      rownames(pars)[3] <- "Rd"
    }
    
    
    
  } else {
    
    # If citransition provided, fit twice.
    dat_vcmax <- data[data$Ci < citransition,]
    dat_jmax <- data[data$Ci >= citransition,]
    
    if(nrow(dat_vcmax) > 0){
      nlsfit_vcmax <- nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, 
                                        Jmax=10000, Rd=Rd, Tleaf=Tleaf, Patm=Patm, returnwhat="Ac"),
                    algorithm=algorithm,
                    data=dat_vcmax, control=nls.control(maxiter=500, minFactor=1/10000),
                    start=list(Vcmax=Vcmax_guess, Rd=Rd_guess))
      p1 <- coef(nlsfit_vcmax)
    } else {
      nlsfit_vcmax <- NULL
      p1 <- c(Vcmax=NA, Rd=NA)
    }
    
    # Don't re-estimate Rd; it is better estimated from the Vcmax-limited region.
    Rd_vcmaxguess <- if(!is.null(nlsfit_vcmax))coef(nlsfit_vcmax)[["Rd"]] else 1.5
    
    if(nrow(dat_jmax) > 0){
      nlsfit_jmax <- nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=10000, 
                                        Jmax=Jmax, Rd=Rd_vcmaxguess, 
                                        Tleaf=Tleaf, Patm=Patm, returnwhat="Aj"),
                    algorithm=algorithm,
                    data=dat_jmax, control=nls.control(maxiter=500, minFactor=1/10000),
                    start=list(Jmax=Jmax_guess))
      p2 <- coef(nlsfit_jmax)
    } else { 
      nlsfit_jmax <- NULL
      p2 <- c(Jmax=NA)
    }
    
    p <- c(p1[1],p2,p1[2])
    pars1 <- if(!is.null(nlsfit_vcmax))summary(nlsfit_vcmax)$coefficients[,1:2] else matrix(rep(NA,4),ncol=2)
    pars2 <- if(!is.null(nlsfit_jmax))summary(nlsfit_jmax)$coefficients[,1:2] else matrix(rep(NA,2),ncol=2)
    pars <- rbind(pars1[1,],pars2,pars1[2,])
    rownames(pars) <- c("Vcmax","Jmax","Rd")
    nlsfit <- list(nlsfit_vcmax=nlsfit_vcmax,nlsfit_jmax=nlsfit_jmax)
  }
  
  
  # Using fitted coefficients, get predictions from model.
  acirun <- Photosyn(Ci=data$Ci, 
                     Vcmax=p[[1]], Jmax=p[[2]], Rd=p[[3]], 
                     PPFD=data$PPFD, 
                     Tleaf=data$Tleaf,
                     Patm=Patm,
                     Tcorrect=Tcorrect)
  
  acirun$Ameas <- data$ALEAF
  acirun$ELEAF <- NULL
  acirun$GS <- NULL
  acirun$Ca <- NULL
  acirun$Ci_original <- data$Ci_original
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  
  # shuffle
  avars <- match(c("Ci","Ameas","Amodel"),names(acirun))
  acirun <- acirun[,c(avars, setdiff(1:ncol(acirun), avars))]
  
  # Organize output
  l <- list()  
  l$df <- acirun[order(acirun$Ci),]
  l$pars <- pars
  l$nlsfit <- nlsfit
  l$Tcorrect <- Tcorrect
  
  # Save function itself, the formals contain the parameters used to fit the A-Ci curve.
  # First save Tleaf, PPFD in the formals (as the mean of the dataset)
  formals(Photosyn)$Tleaf <- mean(data$Tleaf)
  formals(Photosyn)$Patm <- Patm
  formals(Photosyn)$PPFD <- mean(data$PPFD)
  formals(Photosyn)$Vcmax <- l$pars[1]
  formals(Photosyn)$Jmax <- l$pars[2]
  formals(Photosyn)$Rd <- l$pars[3]
  formals(Photosyn)$Tcorrect <- Tcorrect
  l$Photosyn <- Photosyn
  
  # Store Ci at which photosynthesis transitions from Jmax to Vcmax limitation
  l$Ci_transition <- findCiTransition(l$Photosyn)
  
  l$Vcmax_guess <- Vcmax_guess
  l$Jmax_guess <- Jmax_guess
  l$Rd_guess <- Rd_guess
  l$Rd_measured <- haveRd
  
  # Store GammaStar and Km
  if("GammaStar" %in% extrapars) {
    l$GammaStar <- formals(Photosyn)$GammaStar
    l$gstarinput <- TRUE
  } else {
    l$GammaStar <- TGammaStar(mean(data$Tleaf),Patm)
    l$gstarinput <- FALSE
  }
  if("Km" %in% extrapars){
    l$Km <- formals(Photosyn)$Km
    l$kminput <- TRUE
  } else {
    l$Km <- TKm(mean(data$Tleaf),Patm)
    l$kminput <- FALSE
  }
  
  class(l) <- "acifit"
  
return(l)
}


#' @export print.acifit
#' @S3method print acifit
#' @rdname fitaci
print.acifit <- function(x,...){
  
  cat("Result of fitaci.\n\n")
  
  cat("Data and predictions:\n")
  print(x$df)
  
  cat("\nEstimated parameters:\n")
  
  print(x$pars)
  if(x$Tcorrect)
    cat("Note: Vcmax, Jmax are at 25C, Rd is at measurement T.\n")
  else
    cat("Note: Vcmax, Jmax, Rd are at measurement T.\n")
  
  if(x$Rd_measured)
    cat("Note: measured Rd was provided, only Vcmax and Jmax were fit.\n")
  
  cat("\n")
  
  cat("Parameter settings:\n")
  fm <- formals(x$Photosyn)
  pars <- c("Patm","alpha","theta","EaV","EdVC","delsC","EaJ","EdVJ","delsJ")
  fm <- unlist(fm[pars])
  cat(paste0(pars," = ", fm,"\n"))
  
  if(!x$gstarinput | !x$kminput){
    cat("\nEstimated from Tleaf (shown at mean Tleaf):\n")
    if(!x$gstarinput)cat("GammaStar = ",x$GammaStar,"\n")
    if(!x$kminput)cat("Km = ",x$Km,"\n")
  }
  
  if(x$gstarinput | x$kminput){
    cat("\nSet by user:\n")
    if(x$gstarinput)cat("GammaStar = ",x$GammaStar,"\n")
    if(x$kminput)cat("Km = ",x$Km,"\n")
  }
  
}

#' @export summary.acifit
#' @S3method summary acifit
#' @rdname fitaci
summary.acifit <- function(object,...){
  
  print.acifit(object, ...)
  
}


#' @export coef.acifit
#' @S3method coef acifit
#' @rdname fitaci
coef.acifit <- function(object, ...){
 v <- unname(object$pars[,1])
 names(v) <- rownames(object$pars)
return(v)
}

#' @export fitted.acifit
#' @S3method fitted acifit
#' @rdname fitaci
fitted.acifit <- function(object,...){
  
  object$df$Amodel
  
}


#' @export plot.acifit
#' @S3method plot acifit
#' @param x For plot.acifit, an object returned by \code{fitaci}
#' @param xlim Limits for the X axis, if left blank estimated from data
#' @param ylim Limits for the Y axis, if left blank estimated from data
#' @param whichA By default all three photosynthetic rates are plotted (Aj=Jmax-limited (blue), Ac=Vcmax-limited (red), Hyperbolic minimum (black)). Or, specify one or two of them. 
#' @param what The default is to plot both the data and the model fit, or specify 'data' or 'model' to plot one of them, or 'none' for neither (only the plot region is set up)
#' @param add If TRUE, adds to the current plot
#' @param pch The plotting symbol for the data
#' @param addzeroline If TRUE, the default, adds a dashed line at y=0
#' @param addlegend If TRUE, adds a legend (by default does not add a legend if add=TRUE)
#' @param legendbty Box type for the legend, passed to argument bty in \code{\link{legend}}.
#' @param transitionpoint For plot.acifit, whether to plot a symbol at the transition point.
#' @param linecols Vector of three colours for the lines (limiting rate, Ac, Aj), if one value provided it is used for all three.
#' @param lwd Line widths, can be a vector of length 2 (first element for Ac and Aj, second one for the limiting rate).
#' @rdname fitaci
#' @importFrom graphics points
#' @importFrom graphics abline
#' @importFrom graphics legend
plot.acifit <- function(x, what=c("data","model","none"), xlim=NULL, ylim=NULL, 
                        whichA=c("Ac","Aj","Amin"), add=FALSE, pch=19, 
                        addzeroline=TRUE, addlegend=!add, legendbty='o',
                        transitionpoint=TRUE, linecols=c("black","blue","red"),
                        lwd=c(1,2),
                        ...){
  
  # Note that Ci on the X-axis is in molar units!
  if(is.null(ylim))ylim <- with(x$df, c(min(Ameas), 1.1*max(Ameas)))
  if(is.null(xlim))xlim <- with(x$df,c(0, max(Ci_original)))
  if(length(lwd)==1)lwd <- c(lwd,lwd)
  if(length(linecols)==1)linecols <- rep(linecols,3)
  
  # Vector of Ci values at which to evaluate fitted ACi curve.
  Ci <- with(x$df, seq(min(Ci_original), max(Ci_original), length=101))
  
  # Exact model used to fit the A-Ci curve was saved in the object.
  # (parameter settings etc. are preserved)
  pcor <- mean(x$df$Patm)/100
  pred <- x$Photosyn(Ci=Ci * pcor)
  pred$Ci_original <- pred$Ci / pcor
  
  if(!add){
    with(x$df, plot(Ci_original, Ameas, type='n',
                    ylim=ylim,
                    xlim=xlim,
                    xlab=expression(italic(C)[i]~~(ppm)),
                    ylab=expression(italic(A)[n]~~(mu*mol~m^-2~s^-1)),
                    ...
    ))
  }
  if("data" %in% what)with(x$df, points(Ci_original, Ameas, pch=pch,...))
  
  if("model" %in% what){
    if("Aj" %in% whichA)with(pred, points(Ci_original, Aj-Rd, type='l', col=linecols[2],lwd=lwd[1]))
    if("Ac" %in% whichA)with(pred, points(Ci_original, Ac-Rd, type='l', col=linecols[3],lwd=lwd[1]))
    if("Amin" %in% whichA)with(pred, points(Ci_original, ALEAF, type='l', col=linecols[1], lwd=lwd[2]))
  }
  
  if(transitionpoint && "model" %in% what)
    points(x$Ci_transition / pcor, x$Photosyn(Ci=x$Ci_transition)$ALEAF, pch=21, bg="lightgrey", cex=0.8)
  
  if(addzeroline)
    abline(h=0, lty=3)
  
  if(addlegend){
    legend("bottomright", c(expression(italic(A)[c]),
                            expression(italic(A)[j]),
                            "Limiting rate"), lty=1, lwd=c(lwd[1],lwd[1],lwd[2]), 
           col=linecols[3:1], bty=legendbty)
  }
  
}
