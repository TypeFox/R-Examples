###############################################
#### AUTHOR:    Arnost Komarek             ####
####            (2004)                     ####
####                                       ####
#### FILE:      hazard.smoothSurvReg.R     ####
####                                       ####
#### FUNCTIONS: hazard.smoothSurvReg       ####
###############################################

### ===================================================================================
### hazard.smoothSurvReg: Compute hazard functions for objects of class 'smoothSurvReg'
### ===================================================================================
## x ... object of class 'smoothSurvReg' 
## cov
## time0 .... used when the model was log(T - t0) = alpha + beta'x + siga*epsilon
## plot
## cdf
## by
## xlim
## ylim
## xlab
## ylab
## type
## lty
## main
## sub
## legend
## bty
## ... ....... other parameters passed to plot function
hazard <- function(x, ...){
  UseMethod("hazard")
}  

hazard.smoothSurvReg <-
  function(x, cov, logscale.cov, time0 = 0, plot = TRUE,
           by, xlim, ylim, xlab = "t", ylab = "h(t)", 
           type = "l", lty, main, sub, legend, bty = "n", cex.legend = 1, ...)
{
   if (x$fail >= 99){
        cat("No hazard functions, smoothSurvReg failed.\n")
        return(invisible(x))
   }
   is.intercept <- x$estimated["(Intercept)"]
   common.logscale <- x$estimated["common.logscale"]
   est.scale <- x$estimated["Scale"]
   allregrname <- row.names(x$regres)

## INTERCEPT AND SCALE (if it is common)
## =====================================
   mu0 <- ifelse(is.intercept, x$regres["(Intercept)", "Value"], 0)
   if (common.logscale){
     if (est.scale) s0 <- x$regres["Scale", "Value"]
     else           s0 <- x$init.regres["Scale", "Value"]
   }

## COVARIATES FOR REGRESSION
## =========================
   nx <- x$degree.smooth[1, "Mean param."]   
   ncov <- ifelse(is.intercept, nx - 1, nx)

   ## Manipulate with covariate values from the user
   if (missing(cov) && ncov > 0) cov <- matrix(rep(0, ncov), nrow = 1)
   if (ncov == 0)                cov <- NULL                                                ## only intercept in the model
   if (ncov == 1)                cov <- matrix(cov, ncol = 1)
  
   ## Different covariates combinations
   row.cov <- ifelse(is.null(dim(cov)), 1, dim(cov)[1])
   col.cov <- ifelse(is.null(dim(cov)),
                     ifelse(is.null(cov), 0, length(cov)),
                     dim(cov)[2])

 ## COVARIATES FOR LOG-SCALE
 ## ========================
   nz <- x$degree.smooth[1, "Scale param."]   
   if (!common.logscale){
     is.intercept.inscale <- (allregrname[nx+1] == "LScale.(Intercept)")
     ncovz <- ifelse(is.intercept.inscale, nz - 1, nz)

     ## logscale: Manipulate with covariate values from the user
     if (missing(logscale.cov) && ncovz > 0) logscale.cov <- matrix(rep(0, ncovz), nrow = 1)
     if (ncovz == 0)                         logscale.cov <- NULL                              ## only intercept in the model for log-scale
     if (ncovz == 1)                         logscale.cov <- matrix(logscale.cov, ncol = 1)
  
     ## logscale: Different covariates combinations
     logscale.row.cov <- ifelse(is.null(dim(logscale.cov)), 1, dim(logscale.cov)[1])
     logscale.col.cov <- ifelse(is.null(dim(logscale.cov)),
                                ifelse(is.null(logscale.cov), 0, length(logscale.cov)),
                                dim(logscale.cov)[2])
   }
   else{
     ncovz <- 0
     logscale.row.cov <- row.cov
     logscale.col.cov <- 1
   }    


## LINEAR PREDICTOR
## ================
   beta <- x$regres[1:nx, "Value"]
   if (col.cov != ncov) stop("Incorrect cov parameter ")
   if (ncov > 0){
     if (is.intercept) beta <- matrix(beta[2:nx], nrow = ncov, ncol = 1)
     else              beta <- matrix(beta[1:nx], nrow = ncov, ncol = 1)
     cov <- matrix(cov, nrow = row.cov, ncol = col.cov)
     eta <- mu0 + as.numeric(cov %*% beta)
   }
   else{                          ## only intercept in the model
      eta <- rep(mu0, row.cov)
   }
   

## LINEAR PREDICTORS FOR LOG-SCALE, AND COMPUTATION OF A SCALE
## ===========================================================
   if (!common.logscale){
     pars.scale <- x$regres[(nx+1):(nx+nz), "Value"]
     if (logscale.col.cov != ncovz) stop("Incorrect logscale.cov  parameter ")
     if (row.cov != logscale.row.cov) stop("Different number of covariate combinations for regression and log-scale ")

     if (ncovz > 0){
       if (is.intercept.inscale){
         sint <- pars.scale[1]
         pars.scale <- matrix(pars.scale[2:nz], nrow = ncovz, ncol = 1)
       }
       else{
         sint <- 0
         pars.scale <- matrix(pars.scale[1:nz], nrow = ncovz, ncol = 1)
       }
       logscale.cov <- matrix(logscale.cov, nrow = logscale.row.cov, ncol = logscale.col.cov)
       logscale <- sint + as.numeric(logscale.cov %*% pars.scale)
     }
     else{    ## this should never happen if !common.logscale
        sint <- pars.scale[1]
        logscale <- rep(sint, logscale.row.cov)
     }
     s0 <- exp(logscale)
   }
   else{
     s0 <- rep(s0, row.cov)
   }   

## COMPUTE DESIRED QUANTITIES
## ==========================            
   ccoef <- x$spline[["c coef."]]
   knots <- x$spline$Knot
   sigma0 <- x$spline[["SD basis"]][1]
   shift <- x$error.dist$Mean[1]
   scale <- x$error.dist$SD[1]

   ## Survivor function of the fitted error distribution
   ## (survivor function of epsilon)
   sfitted.un <- function(u){
      normals <- pnorm(u, mean = knots, sd = sigma0)
      value <- 1 - (t(ccoef) %*% normals)[1]
      return(value)
   }

   ## Density function of the fitted error distribution
   ## (density function of epsilon)
   dfitted.un <- function(u){
      normals <- dnorm(u, mean = knots, sd = sigma0)
      value <- (t(ccoef) %*% normals)[1]
      return(value)
   }     

   ## Grid
   if (missing(xlim)){
      xmin <- time0
      xmax <- exp(max(x$y[,1])) + time0
      xlim <- c(xmin, xmax)
   }
   if (missing(by)){
      by <- (xlim[2] - xlim[1])/100
   }
   if (xlim[1] < time0) xlim[1] <- time0
   if (xlim[2] < time0) xlim[2] <- xlim[1] + 0.01

   grid <- seq(xlim[1], xlim[2], by) + 0.01

   ## Values
   etas <- matrix(rep(eta, rep(length(grid), row.cov)), ncol = row.cov)
   s0s <- matrix(rep(s0, rep(length(grid), row.cov)), ncol = row.cov)
   grid2 <- matrix(rep(grid, row.cov), ncol = row.cov)
   grid2 <- (log(grid2 - time0) - etas) / s0s
   haz <- list()
   for (i in 1:row.cov){
      grid3 <- matrix(grid2[,i], ncol = 1)
      Sfun <- apply(grid3, 1, "sfitted.un")
      dfun <- apply(grid3, 1, "dfitted.un")
      Sfun[Sfun <= 0] <- NA
      haz[[i]] <- (1/(s0[i]*(grid - time0))) * (dfun/Sfun)    ### corrected by AK on 14/09/2005 (added s0[i]*)
   }

   ## ylim
   if (missing(ylim)){
     ymax <- max(sapply(haz, max, na.rm = TRUE), na.rm = TRUE)
     ylim <- c(0, ymax)
   }     
   
   ## lty
   if (missing(lty)){
      lty <- 1:row.cov
   }

   ## main and sub
   if (missing(main)) main <- "Fitted Hazard"
   if (missing(sub)){
      aic <- round(x$aic, digits = 3)
      df <- round(x$degree.smooth$df, digits = 2)
      nparam <- x$degree.smooth[["Number of parameters"]]
      sub <- paste("AIC = ", aic, ",   df = ", df, ",   nParam = ", nparam, sep="")
   }

   ## Plot it
   if (plot){
      plot(grid, haz[[1]],
           type = type, lty = lty[1], ylim = ylim, xlab = xlab, ylab = ylab, bty = bty, ...)
      title(main = main, sub = sub)
      if (row.cov > 1){
         for (i in 2:row.cov){
            lines(grid, haz[[i]], lty = lty[i])
         }
      }
      leg <- numeric(2)
      leg[1] <- xlim[1]
      leg[2] <- ylim[2]
      legjust <- numeric(2)
      legjust[1] <- 0
      legjust[2] <- 1
      if (missing(legend)) legend <- paste("cov", 1:row.cov, sep = "")
      legend(leg[1], leg[2], legend = legend, lty = lty, bty = "n", xjust = legjust[1], yjust = legjust[2], cex = cex.legend)
   }
   to.return <- data.frame(grid, haz[[1]])
   if (row.cov > 1)
   for (i in 2:row.cov){
      to.return <- cbind(to.return, haz[[i]])
   }
   names(to.return) <- c("x", paste("y", 1:row.cov, sep = ""))

   if (plot) return(invisible(to.return))
   else      return(to.return)
}   



