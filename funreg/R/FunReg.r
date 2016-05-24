#' @title Perform penalized functional regression
#' @description Performs a penalized functional regression as in Goldsmith et al. (2012)
#' on irregularly measured data such as that found in ecological momentary
#' assessment (see, e.g., Shiffman, Stone, & Hufford, 2008). 
#' @note This function mostly follows code by Jeff Goldsmith and co-workers:
#' the sample code from Goldsmith et al (2011), and the "pfr" function in
#' the "refund" R package.  However, this code is adapted here to allow
#' idiosyncratic measurement times and unequal numbers of observations
#' per subject to be handled easily, and also allows the use of a different
#' estimation method.  Also follows some sample code for penalized
#' B-splines from Eilers and Marx (1996) in implementing B-splines.
#' As the pfr function in refund also does, the function calls
#' the gam function in the mgcv package (Wood 2011) to do much of the
#' internal calculations.
#' @param id An integer or string uniquely identifying the
#' subject to which each observation belongs
#' @param response The response, as a vector, one for each subject
#' @param time The time of observation, as a vector,
#' one for each observation (i.e., each
#' assessment on each person)
#' @param x The functional predictor(s), as a matrix,
#' one row for each observation (i.e., for each
#' assessment on each person)
#' @param basis.method An integer, either 1 or 2, describing how
#' the beta function should be internally
#' modeled. A value of 1 indicates that a
#' truncated power spline basis should be used,
#' and a value of 2 indicates that a B-spline
#' basis should be used.
#' @param deg An integer, either 1, 2, or 3, describing how 
#' complicated the behavior of the beta function between knots 
#' may be.  1, 2, or 3 represent linear, quadratic or cubic 
#' function between knots.
#' @param deg.penalty Only relevant for B-splines.  The
#' difference order used to weight the smoothing
#' penalty (see Eilers and Marx, 1996)
#' @param family The response distribution.  For example,
#' this is  \code{family=gaussian} for normal linear
#' models, \code{family=binomial} for logistic regression models,
#' or \code{family=poisson} for count models.  See the \code{gam}
#' documentation in the \code{mgcv} package, or use \code{help(family)}
#' for details on \code{family} objects.
#' @param other.covariates Subject-level (time-invariant) covariates, 
#' if any, as a matrix, one column per covariate. The default, 
#' \code{NULL}, means that no subject-level covariates will be included 
#' in the model.
#' @param num.bins The number of knots used in the spline basis for the
#'  beta function. The default is based on the  Goldsmith et al. (2011) 
#'  sample code.
#' @param preferred.num.eigenfunctions The number of eigenfunctions to use in approximating the covariance function of x (see Goldsmith et al., 2011)
#' @param preferred.num.knots.for.beta number of knots to use in the spline
#' estimation.  The default, is based on the Goldsmith et al (2011) sample code.
#' @param se.method An integer, either 1 or 2, describing how
#' the standard errors should be calculated. A value
#' of 1 means that the uncertainty related to
#' selecting the smoothing parameter is ignored.
#' Option 2 means that a Bayesian approach is used
#' to try to take this uncertainty into account
#' (see the documentation for Wood's \code{mgcv} package).
#' @param smoothing.method  An integer, either 1 or 2, describing
#' how the weight of the smoothing penalty should be
#' determined.  Option 1 means that the smoothing
#' weight should be estimated using an approach
#' similar to restricted maximum likelihood, and
#' Option 2 means an approach similar to generalized
#' cross-validation.  Option 1 is strongly
#' recommended (based both on our experience and on
#' remarks in the documentation for the gam function
#' in the mgcv package).   
#' @param times.for.fit.grid  Points at which to calculate
#' the estimated beta function.  The default, NULL, means that the
#' code will choose these times automatically.
#' @return An object of type \code{funreg}.  This object 
#' can be examined using \code{summary}, \code{print},
#'  or \code{fitted}.
#' @references Crainiceanu, C., Reiss, P., Goldsmith, J., Huang, L., Huo, L.,
#'    Scheipl, F. (2012). refund: Regression with Functional Data
#'    (version 0.1-6). R package Available online at cran.r-project.org.
#' 
#'  Eilers, P. H. C., and Marx, B. D. (1996). Flexible smoothing with
#'    B-splines and penalties. Statistical Science, 11, 89-121.
#'
#'  Goldsmith, J., Bobb, J., Crainiceanu, C. M., Caffo, B., and Reich, D.
#'    (2011). Penalized functional regression. Journal of Computational
#'    and Graphical Statistics, 20(4), 830-851. The sample code can be
#'    found at www.jeffgoldsmith.com/Downloads/PFR_Web_Appendix.zip;
#'    in writing parts of this function I especially followed "PFR_Example.R",
#'    written on Jan. 15 2010, by Jeff Goldsmith.
#' 
#'  Ruppert, D., Wand, M., and Carroll, R. (2003). Semiparametric regression.
#'    Cambridge, UK: Cambridge University Press.
#' 
#'  Shiffman, S., Stone, A. A., and Hufford, M. R. (2008). Ecological
#'    momentary assessment. Annual Review of Clinical Psychology, 4, 1-32.
#' 
#'  Wood, S.N. (2006) Generalized Additive Models: An Introduction with
#'     R. Chapman and Hall/CRC.
#' 
#'  Wood, S.N. (2011) Fast stable restricted maximum likelihood and
#'    marginal likelihood estimation of semiparametric generalized linear
#'    models. Journal of the Royal Statistical Society (B) 73(1):3-36.
#' @seealso \code{\link{fitted.funreg}}, \code{link{plot.funreg}},
#'          \code{\link{print.funreg}}, \code{link{summary.funreg}}
#' @importFrom mgcv gam
#' @importFrom mgcv s
#' @examples 
#' simple.model <- funreg(id=SampleFunregData$id,
#'                        response=SampleFunregData$y,
#'                        time=SampleFunregData$time,
#'                        x=SampleFunregData$x1,
#'                        family=binomial);
#' print(simple.model);                 
#' par(mfrow=c(2,2));                          
#' plot(x=simple.model$model.for.x[[1]]$bin.midpoints,
#'      y=simple.model$model.for.x[[1]]$mu.x.by.bin,
#'      xlab="Time t",ylab="X(t)",main="Smoothed mean x values"); 
#' # The smoothed average value of the predictor function x(t) at different times t. 
#' # The ``[[1]]'' after model.for.x is there because model.for.x is a list with one entry.
#' # This is because more than one functional covariate is allowed.
#' plot(simple.model,type="correlations");
#' # The marginal correlation of x(t) with y at different times t.
#' # It appears that earlier time points are more strongly related to y.
#' plot(simple.model,type="coefficients");
#' # The functional regression coefficient of y on x(t).  
#' # It also appears that earlier time points are more strongly related to y.
#' plot(simple.model$subject.info$response,
#'      simple.model$subject.info$fitted,
#'      main="Predictive Performance",
#'      xlab="True Y",
#'      ylab="Fitted Y");  
#'@note In the example below, to fit a more complicated model, replace
#'  \code{x=SampleFunregData$x1} with \code{x=cbind(SampleFunregData$x1,
#'   SampleFunregData$x2),other.covariates=cbind(SampleFunregData$s1,
#'   SampleFunregData$s2, SampleFunregData$s3, SampleFunregData$s4)}.  This 
#'   model will take longer to run, perhaps 10 or 20 seconds.  Then try
#'   \code{plot(complex.model)}.      
#'@export
funreg <- function( id,
                    response,  
                    time,     
                    x,      
                    basis.method=1,
                    deg=2,       
                    deg.penalty=2, 
                    family=gaussian, 
                    other.covariates=NULL,
                    num.bins=35,
                    preferred.num.eigenfunctions=30,
                    preferred.num.knots.for.beta=35,
                    se.method=1,
                    smoothing.method=1,  
                    times.for.fit.grid=NULL
                    ) {
    call.info <- match.call();
    ## Pre-process the functional predictors;
        x <- as.matrix(x,drop=FALSE);
        if (ncol(x)>6) {
           print("Please specify no more than six functional predictor variables.");
        }
    ## Handle possible missing data (NA's) within assessments;
        stopifnot(length(time)==nrow(x));
        stopifnot(length(time)==length(id));
        ## Handle possible missing data (NA's) within ID;
        if (is.null(other.covariates)) {
           observations.to.include <- which((!is.na(time))&
                                            (!is.na(id))&
                                            (!is.na(response))&
                                            (apply(is.na(x),1,sum)==0) );
        } else {
           if (length(which(is.na(other.covariates)))>0) {
              warning(paste("At least one of the subject-level covariates had",
                           "missing data.  Therefore, the corresponding rows are",
                           "excluded from the analysis."));
           }
           other.covariates <- as.matrix(other.covariates,drop=FALSE);
           observations.to.include <- which((!is.na(time))&
                                             (!is.na(id))&
                                             (!is.na(response))&
                                             (apply(is.na(x),1,sum)==0)&
                                             (apply(is.na(other.covariates),1,sum)==0) );          
        }
        subjects.to.include <- unique(id[observations.to.include]);
        if (length(subjects.to.include) < length(unique(id))) {
           warning("At least one subject had no complete ",
            "assessments, due to missingness.");
        }
        if (length(observations.to.include)<length(time)) {
            warning(paste("One or more observations have been deleted because of",
                        "missing data."));
        }
        id <- id[observations.to.include];
        time <- time[observations.to.include];
        response <- response[observations.to.include];
        x <- x[observations.to.include,1:ncol(x),drop=FALSE];
        if (!is.null(other.covariates)) {
            other.covariates <- other.covariates[observations.to.include,,drop=FALSE]
        }
    ## Do the eigenfunction decomposition in order to obtain approximated
    ## values for the x functions at a fixed grid of points shared by all
    ## subjects.  As a side effect of this, obtain a list of the unique
    ## values of the id variable.
        eigenfunctions.answers  <- list();
        CJ  <- list();
        M  <- list();
        betafn.estimate.by.bin <- NULL;
        fitted.betafn.by.grid <- NULL;
        se.betafn.estimate.by.bin <- NULL;
        se.fitted.betafn.by.grid <- NULL;
        for (this.xcol in 1:ncol(x)) {
            eigenfunctions.answers[[this.xcol]]  <-
                         funeigen(x=x[,this.xcol,drop=TRUE],
                                           time=time,
                                           id=id,
                                           num.bins=num.bins,
                                           preferred.num.eigenfunctions=
                                                 preferred.num.eigenfunctions);
            id.compressed <- eigenfunctions.answers[[this.xcol]]$id;
            stopifnot(length(id.compressed)==length(unique(id)));
            stopifnot(min(id.compressed==unique(id))==TRUE);
                  # should be the same each time; this is just done for
                  # convenience to get a list of unique ID's included;
        }
    ## Process some other inputs representing user options;
        if (!((deg==1)|(deg==2)|(deg==3))) {
            stop("Please choose a degree of 1, 2, or 3.");
        }
        if (smoothing.method==1) {
            smoothing.type <- "REML";
        } else {
            if (smoothing.method==2) {
                smoothing.type <- "GCV.Cp";
            } else {
            stop("Unrecognized standard errors option");
        }}
        if (basis.method==1) {
           basis.type <- "TruncatedPower";
        } else {
            if (basis.method==2) {
                 basis.type <- "BSpline";
            } else {
            stop("Unrecognized standard errors option");
        }}
    ## Input the non-functional covariates, called
    ## "other.covariates," if any;
        if(!is.null(other.covariates)) {
            other.covariates <- as.matrix(other.covariates,drop=FALSE);
            num.other.covariates <- ncol(other.covariates);
            stopifnot(length(id)==nrow(other.covariates));
            other.covariates.compressed <- matrix(NA,
                                                  nrow=length(id.compressed),
                                                  ncol=ncol(other.covariates));
            if (is.null(colnames(other.covariates))) {
                colnames(other.covariates) <- paste("Covariate",
                                                    1:num.other.covariates);
            }
            rownames(other.covariates.compressed) <- id.compressed;
            colnames(other.covariates.compressed) <-
                colnames(other.covariates);
            for (this.id.num in 1:length(id.compressed)) {
                this.id <- id.compressed[this.id.num];
                these <- which(id==this.id);
                stopifnot(length(these)>0);
                if (min(apply(other.covariates[these,,drop=FALSE],2,var)>1e-20)) {
                    stop(paste("Error: One of the supposedly subject-level",
                               "covariates is not constant within subjects."));
                }
                other.covariates.compressed[this.id.num,] <-
                         other.covariates[min(these),,drop=FALSE];
            }
            stopifnot(nrow(other.covariates.compressed)==length(id.compressed));
            if(min(apply(other.covariates.compressed,2,var)<1e-20)) {
                stop(paste("Error: One of the time-invariant covariates is",
                           "essentially constant between subjects."));
                     # This would be an error because it would make the
                     # model non-identifiable.  There is already
                     # understood to be an intercept column
                     # without specifying it in "other.covariates."
            }
        } else {
            num.other.covariates <- 0;
            other.covariates.compressed <- NULL;
        }
    ## Input the response variable
        if (is.matrix(response)) {stopifnot(ncol(response)==1);}
        response <- as.vector(response);
        response.compressed <- rep(NA,length(id.compressed));
        stopifnot(length(response)==length(id));
        for (this.id.num in 1:length(id.compressed)) {
            this.id <- id.compressed[this.id.num];
            these <- which(id==this.id);
            stopifnot(length(these)>0);
            if (min(var(response[these]))>1e-20) {
                stop(paste("Error: The response must be constant",
                           "within subjects."));
            }
            response.compressed[this.id.num] <- response[min(these)];
        }
        if (identical(family,binomial)) {
            n0 <- sum(response.compressed==0);
            n1 <- sum(response.compressed==1);
            if (n0+n1<length(response.compressed)) {
                warning(paste("The family is given as binomial",
                              "but the responses are not all 0's and 1's."));
            } else {
                preferred.min.n <- max(preferred.num.eigenfunctions,
                                       preferred.num.knots.for.beta);
                if (n0 < preferred.min.n) {
                    warning(paste("There are only",n0,
                                  "subjects with response 0. ",
                                  "Function estimation may be poor."));
                }
                if (n1 < preferred.min.n) {
                    warning(paste("There are only",n1,
                                  "subjects with response 1. ",
                                  "Function estimation may be poor."));
                }
            }
        }
        stopifnot(length(response.compressed)==length(id.compressed));
        if(!is.null(other.covariates)) {
            stopifnot(length(response.compressed)==nrow(other.covariates.compressed));
        }
    ## Read in information from the eigenfunction decomposition;
        bin.midpoints <- eigenfunctions.answers[[1]]$bin.midpoints;
        num.eigenfunctions.to.use <- length(eigenfunctions.answers[[1]]$lambda);
        num.knots.for.betafn <- min(preferred.num.knots.for.beta,
                                  num.eigenfunctions.to.use);
        # num.knots.for.betafn is "kb" in Goldsmith et al.,
        #  "sparse_simulation.R";
        stopifnot(num.knots.for.betafn>5);
    ## Calculate the spline basis
        temp <- make.funreg.basis ( basis.type=basis.type,
                                        deg=deg,
                                        num.knots=num.knots.for.betafn,
                                        times=bin.midpoints );
        interior.knot.locations <- temp$interior.knot.locations;
        basis.for.betafn.by.bin <- temp$basis.for.betafn;
        if (is.null(times.for.fit.grid)) {
            times.for.fit.grid <- seq(min(time,na.rm=TRUE),
                                      max(time,na.rm=TRUE),
                                      length=1000);
        } else {
            times.for.fit.grid <- sort(as.vector(times.for.fit.grid));
        }
    ## Calculate the design matrix and penalty weighting matrix for the intermediate
    ## regression (to some extent this follows the Goldsmith et al. sample code);
        if (is.null(other.covariates.compressed)) {
            X0 <- rep(1,length(id.compressed));
        } else {
            X0 <- cbind(other.covariates.compressed,1);
        }
        for (this.xcol in 1:ncol(x)) {
            J <- t(eigenfunctions.answers[[this.xcol]]$psi) %*%
                      basis.for.betafn.by.bin *
                (bin.midpoints[2]-bin.midpoints[1]);
            CJ[[this.xcol]] <- eigenfunctions.answers[[this.xcol]]$C %*% J;
        num.unpenalized.in.CJ <- deg+1;
           # based on "Sparse_Simulation.r" by Goldsmith;
            stopifnot(ncol(CJ[[this.xcol]])>=num.unpenalized.in.CJ+1);
    ## Calculate the penalty matrices;
        if (basis.type=="TruncatedPower") {
                M[[this.xcol]] <- diag(as.vector(c(rep(0, num.unpenalized.in.CJ),
                                       rep(1,ncol(CJ[[this.xcol]])-num.unpenalized.in.CJ))));
        }
        if (basis.type=="BSpline") {
                temp <- diag(as.vector(rep(1,ncol(CJ[[this.xcol]]))));
                for (d in 1:deg.penalty) { temp <- diff(temp); }
                M[[this.xcol]] <-
                          crossprod(temp);
            }
         }
    ## Do the intermediate regression.  Calls the gam function in mgcv,
    ## in almost the same way as the pfr function in the refund
    ## package does (Ciprian Crainiceanu, Philip Reiss, Jeff Goldsmith,
    ## Lei Huang, Lan Huo, Fabian Scheipl, Sonja Greven, Jaroslaw
    ## Harezlak, Madan Gopal Kundu, Yihong Zhao, Matt McLean, 2013)
    if (ncol(x)==1) {
        X1 <- CJ[[1]]; M1 <- M[[1]];
        fit.model <- gam(response.compressed~X0+X1+0,
                         paraPen=list(X1=list(M1)),
                         method=smoothing.type,
                         family=family);
    }
    if (ncol(x)==2) {
        X1 <- CJ[[1]]; M1 <- M[[1]];
        X2 <- CJ[[2]]; M2 <- M[[2]];
        fit.model <- gam(response.compressed~X0+X1+X2+0,
                         paraPen=list(X1=list(M1),
                                      X2=list(M2)),
                         method=smoothing.type,
                         family=family);
    }
    if (ncol(x)==3) {
        X1 <- CJ[[1]]; M1 <- M[[1]];
        X2 <- CJ[[2]]; M2 <- M[[2]];
        X3 <- CJ[[3]]; M3 <- M[[3]];
        fit.model <- gam(response.compressed~X0+X1+X2+X3+0,
                         paraPen=list(X1=list(M1),
                                      X2=list(M2),
                                      X3=list(M3)),
                         method=smoothing.type,
                         family=family);
    }
    if (ncol(x)==4) {
        X1 <- CJ[[1]]; M1 <- M[[1]];
        X2 <- CJ[[2]]; M2 <- M[[2]];
        X3 <- CJ[[3]]; M3 <- M[[3]];
        X4 <- CJ[[4]]; M4 <- M[[4]];
        fit.model <- gam(response.compressed~X0+X1+X2+X3+X4+0,
                         paraPen=list(X1=list(M1),
                                      X2=list(M2),
                                      X3=list(M3),
                                      X4=list(M4)),
                         method=smoothing.type,
                         family=family);
    }
    if (ncol(x)==5) {
        X1 <- CJ[[1]]; M1 <- M[[1]];
        X2 <- CJ[[2]]; M2 <- M[[2]];
        X3 <- CJ[[3]]; M3 <- M[[3]];
        X4 <- CJ[[4]]; M4 <- M[[4]];
        X5 <- CJ[[5]]; M5 <- M[[5]];
        fit.model <- gam(response.compressed~X0+X1+X2+X3+X4+X5+0,
                         paraPen=list(X1=list(M1),
                                      X2=list(M2),
                                      X3=list(M3),
                                      X4=list(M4),
                                      X5=list(M5)),
                         method=smoothing.type,
                         family=family);
    }
    if (ncol(x)==6) {
        X1 <- CJ[[1]]; M1 <- M[[1]];
        X2 <- CJ[[2]]; M2 <- M[[2]];
        X3 <- CJ[[3]]; M3 <- M[[3]];
        X4 <- CJ[[4]]; M4 <- M[[4]];
        X5 <- CJ[[5]]; M5 <- M[[5]];
        X6 <- CJ[[6]]; M6 <- M[[6]];
        fit.model <- gam(response.compressed~X0+X1+X2+X3+X4+X5+X6+0,
                         paraPen=list(X1=list(M1),
                                      X2=list(M2),
                                      X3=list(M3),
                                      X4=list(M4),
                                      X5=list(M5),
                                      X6=list(M6)),
                         method=smoothing.type,
                         family=family);
        stopifnot(length(fit.model$fitted.values)==length(response.compressed));
    }
    ## Calculate the spline basis for all points at which to fit
    ## the beta function
        temp <- make.funreg.basis( basis.type=basis.type,
                                    deg=deg,
                                    num.knots=num.knots.for.betafn,
                                    times=times.for.fit.grid );
        interior.knot.locations <- temp$interior.knot.locations;
        basis.for.betafn.by.grid <- temp$basis.for.betafn;
    ## Back-transform the coefficients from the intermediate regression
    ## into the functional betafn.
        for (this.xcol in 1:ncol(x)) {
            number.before <- num.other.covariates+1;
            if (this.xcol > 1) {
                for (prev.xcol in 1:(this.xcol-1)) {
                    number.before <- number.before + ncol(CJ[[prev.xcol]]);
        }
            }
            indices.for.this.spline <- (number.before+1):
                                          (number.before+ncol(CJ[[this.xcol]]));
            b.for.this.spline <- fit.model$coefficients[indices.for.this.spline,
                                                   drop=FALSE];
            stopifnot(length(b.for.this.spline)==
                              ncol(basis.for.betafn.by.bin));
            this.fitted.betafn <- basis.for.betafn.by.bin%*%b.for.this.spline;
            betafn.estimate.by.bin <- cbind(betafn.estimate.by.bin,this.fitted.betafn);
            this.fitted.betafn <- basis.for.betafn.by.grid%*%b.for.this.spline;
            fitted.betafn.by.grid <- cbind(fitted.betafn.by.grid,this.fitted.betafn);
        # betafn.estimate.by.bin contains beta function estimates
        # evaluated at the midpoint of each bin used to model the
        # covariance of the x functions.
        # fitted.betafn.by.grid contains beta function estimates
        # evaluated at each point at which the user requested they
        # be calculated in times.for.fit.grid, if this was specified.
        # If the user didn't specify times.for.fit.grid, then
        # these times are chosen by the function.
        }
        rownames(betafn.estimate.by.bin) <- paste("Beta at t=",
                                                 round(bin.midpoints,3));
        rownames(fitted.betafn.by.grid) <- paste("Beta at t=",
                                                 round(times.for.fit.grid,3));
        colnames(betafn.estimate.by.bin) <- paste("Beta",1:ncol(x),sep="");
        colnames(fitted.betafn.by.grid) <- colnames(betafn.estimate.by.bin);
    ## Calculate standard errors (following "PFR_Example.R" with
    ## some modifications);
        if (se.method==1) {
            cov.matrix <- fit.model$Ve;
        } else {
            if (se.method==2) {
                cov.matrix <- fit.model$Vp;
            } else {
            stop("Unrecognized standard errors option");
        }}
        for (this.xcol in 1:ncol(x)) {
            number.before <- num.other.covariates+1;
            if (this.xcol > 1) {
                for (prev.xcol in 1:(this.xcol-1)) {
                    number.before <- number.before + ncol(CJ[[prev.xcol]]);
                }
            }
            indices.for.this.spline <- (number.before+1):
                                          (number.before+ncol(CJ[[this.xcol]]));
            var.betafn.estimate.by.bin <- basis.for.betafn.by.bin %*%
                                      cov.matrix[indices.for.this.spline,
                                                 indices.for.this.spline]%*%
                                      t(basis.for.betafn.by.bin);
            var.fitted.betafn.by.grid <- basis.for.betafn.by.grid %*%
                                         cov.matrix[ indices.for.this.spline,
                                                     indices.for.this.spline]%*%
                                         t(basis.for.betafn.by.grid);
        # "var.fitted.betafn" is "varBetaHat" in Goldsmith et al.,
        # "PFR_Example.R";
            this.se.fitted <- sqrt(diag(var.betafn.estimate.by.bin));
            se.betafn.estimate.by.bin <- cbind(se.betafn.estimate.by.bin,
                                             this.se.fitted);
            this.se.fitted <- sqrt(diag(var.fitted.betafn.by.grid));
            se.fitted.betafn.by.grid <- cbind(se.fitted.betafn.by.grid,
                                             this.se.fitted);
      }
    marg.corr.dataset <- NULL;
    ## Calculate the marginal correlation of x with the linear
    ## predictor at different times;
    nbins <- nrow(eigenfunctions.answers[[1]]$smoothed.cov.mat.between.bins);
    ## Return the results;
    if (num.other.covariates>0) {
        other.covariates.estimate <-
                        fit.model$coefficients[1:num.other.covariates,
                                               drop=FALSE];
        other.covariates.se <- sqrt(diag(cov.matrix[1:num.other.covariates,
                                                    1:num.other.covariates,
                                                    drop=FALSE]));
        if (!is.null(colnames(other.covariates))) {
            names(other.covariates.estimate) <- colnames(other.covariates);
            names(other.covariates.se) <- colnames(other.covariates);
        } else {
            names(other.covariates.estimate) <- paste("Covariate",
                                                 1:num.other.covariates,sep="");
            names(other.covariates.se) <- paste("Covariate",
                                                1:num.other.covariates,sep="");
        }
    } else {
        other.covariates.estimate <- NULL;
        other.covariates.se <- NULL;
    }
    a0c <- fit.model$coefficients[num.other.covariates+1];
    intercept.se <- sqrt(diag(cov.matrix[num.other.covariates+1,
                                         num.other.covariates+1,
                                         drop=FALSE]));
    intercept.back.correction <- rep(NA,ncol(x));
    for (this.coef in 1:ncol(x)) {
        bin.width <- eigenfunctions.answers[[this.coef]]$bin.midpoints[2]-
                  eigenfunctions.answers[[this.coef]]$bin.midpoints[1];
        intercept.back.correction[this.coef] <- bin.width*
           drop(eigenfunctions.answers[[this.coef]]$mu.x.by.bin%*%
           as.matrix(betafn.estimate.by.bin)[,this.coef]);
    }
    a0 <- a0c - sum(intercept.back.correction);
    original.data <- list(id=id,
                          response=response,
                          time=time,
                          x=x,
                          other.covariates=other.covariates);
    settings <- list(basis.method = basis.method,
                     deg = deg,
                     deg.penalty = deg.penalty,
                     family = family,
                     num.bins = num.bins,
                     preferred.num.eigenfunctions = preferred.num.eigenfunctions,
                     preferred.num.knots.for.beta = preferred.num.knots.for.beta,
                     se.method = se.method,
                     smoothing.method = smoothing.method);
    intermediate.work.results <- list(
                basis.for.beta.by.grid=basis.for.betafn.by.grid,
                    # The spline basis constructed for approximating the
                    # beta function.
                CJ=CJ,
                    # The product of the subjects' loadings on the eigenvalues
                    # for x (i.e., C) and the integral of the products of the
                    # eigenfunctions for x times the spline basis function
                    # which was constructed for beta (i.e., J).  This is used
                    # as the design matrix (regression matrix) in the
                    # intermediate calculation which finds the spline
                    # coefficients to estimate beta.
                fit.model=fit.model,
                    # The object returned from the function which was
                    # used to fit the intermediate regression.
                interior.knot.locations=interior.knot.locations
                    # The locations of the interior knots for the spline
                    # approximation of the beta function.
                );
    output <- list(call.info=call.info,
                   data=original.data,
                   intermediate.work.results=intermediate.work.results,
                        # Detailed technical information about the
                        # penalized regression which was used to obtain
                        # the spline coefficients.
                   model.for.x.functions=eigenfunctions.answers,
                        # Detailed technical information about the
                        # estimated mean and covariance functions of the
                        # predictor functions, based on their eigenfunction
                        # decomposition.
                   model.cov.matrix=cov.matrix,
                        # Estimated covariance matrix for the intermediate
                        # fitted model.;
                   betafn.estimate.by.grid=drop(fitted.betafn.by.grid),
                        # The estimated values of the functional regression
                        # coefficient beta at each of a sequence of time points
                        # given by times.for.fit.grid
                   betafn.se.by.grid=drop(se.fitted.betafn.by.grid),
                        # The estimated (pointwise) standard errors for the fitted
                        # beta estimates.  These are underestimates of the true
                        # uncertainty because they ignore bias caused by
                        # smoothing (as Goldsmith et al., 2011, pointed out),
                        # and this bias may be very large for sparse
                        # irregular data, and especially near the edges of the
                        # time interval.
                   betafn.estimate.by.bin=drop(betafn.estimate.by.bin),
                   intercept.estimate.centered=a0c,
                        # Estimate of intercept given centered x;
                   intercept.se=intercept.se,
                        # Estimate of standard error of intercept given centered x;
                   intercept.estimate.uncentered=a0,
                   subject.info=data.frame(cbind(id=id.compressed,
                                      response=response.compressed,
                                      fitted=fit.model$fitted.values,
                                      other.covariates.compressed)),
                   other.covariates.estimate=other.covariates.estimate,
                        # The regression coefficients for the subject-level
                        # time-invariant covariates, if any were present
                   other.covariates.se=other.covariates.se,
                        # The estimated standard errors for the subject-level
                        # time-invariant covariates
                   settings=settings,
                        # The original settings specified by the user;
                   times.for.fit.grid=times.for.fit.grid,
                        # The time points at which the beta function estimate
                        # and its standard errors are evaluated;
                   xnames=colnames(x)
                );
    class(output) <- "funreg";
    invisible(output);
}

