## ->  logLik1: Calculates contribution of one subject to the log-likelihood
#' Calculates contribution of one subject to the log-likelihood
#'
#' This function is generic; method functions can be written to handle specific classes of objects.
#'
#'
#' @param modfit an object representing model fitted to data using ML estimation.
#' @param dt1 a data frame with data for one subject, for whom the log-likelihood function is to be evaluated
#' @param dtInit an optional auxiliary data frame.
#' @return numeric scalar value representing contribution of a given subject to the overall log-likelihood returned by \code{logLik()} function.
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @examples 
#'   require(nlme)
#'   logLik(fm1 <- lme(distance ~ age, data = Orthodont)) # random is ~ age
#'   dt1 <- subset(Orthodont, Subject == "M01")
#'   logLik1(fm1, dt1)
#' @references 
#'   Galecki, A., and Burzykowski, T. (2013) ``Linear Mixed-Effects Models: Step-by-Step Approach'', Springer     
#' @export
logLik1 <-  function(modfit, dt1, dtInit) UseMethod("logLik1")


## logLik1.lme method
#' Calculates contribution of one subject to the log-likelihood for \code{lme} object
#'
#' This is method for \code{logLik1()} generic function.
#'
#' Calculates profile likelihood (with beta profiled out) for *one* subject.
#' Data with *one* level of grouping only.
#' correlation component in modelStruct not implemented.
#'
#' @method logLik1 lme
#' @param modfit an \code{lme} object representing model fitted using maximum likelihood.
#' @param dt1 a data frame with data for one subject, for whom the log-likelihood function is to be evaluated
#' @param dtInit an optional auxiliary data frame.
#' @return numeric scalar value representing contribution of a given subject to
#'   the overall log-likelihood returned by \code{logLik()} function applied to
#'   \code{lme} object defined by \code{modfit} argument.
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @S3method logLik1 lme
#' @examples
#' 
#' require(nlme)
#' lm3.form <- visual ~ visual0 + time + treat.f 
#' (fm16.5ml <-                                # M16.5
#'    lme(lm3.form,                         
#'    random = list(subject = pdDiag(~time)), 
#'    weights = varPower(form = ~time),
#'    data = armd, method = "ML"))
#'  df1 <- subset(armd, subject == "1")       # Panel R20.7
#'  logLik1(fm16.5ml, df1)
#' 
logLik1.lme <- function(modfit, dt1, dtInit){

   .functionName <- "logLik1.lme"                     # Function name
   .traceR <- if (is.null(options()$traceR)) function(...){} else options()$.traceR       

   .traceR(1, lbl = "-> logLik1.lme STARTS")
   if (missing(dtInit)) dtInit <- NULL
 
    m            <- modfit$modelStruct                 # Model structure
    sigma        <- modfit$sigma                       # sigma
 
    D            <- as.matrix(m$reStruct[[1]])         # "subject"
    D            <- D  * sigma^2                       # Matrix D 
    
    clw  <- modfit$call$weights
    vecR <- rep(sigma, nrow(dt1))
    
    
    .traceR(10, lbl = "if len.clw before")
    if (length(clw)){
    .traceR(101, lbl = "if.clw", store = FALSE)
    
    .traceR(102,  lbl = "clw before if.inherits")
    if (inherits(eval(clw), "formula")) clw <- call("varFixed", clw) 
    .traceR(103, lbl = "clw after if.inherits")
    clwl  <-  as.list(clw) 
    .traceR(104, lbl = "clwl as list")
     
    varFun <- as.character(clwl[[1]])                  # VarFun="varPower"
    varSt <- m$varStruct
    vf.coef      <- coef(varSt,                        # Variance function ...
                         unconstrained=FALSE)          # coef. extracted
 
    names(vf.coef) <- NULL
    .traceR(105, lbl = "vf.coef")
 
    ## varFun <- "varPower"
    #args   <- list(form = ~time)  ### <----
    args    <- as.list(clwl[-1]) 
    args    <- c(args, value = vf.coef) # Replace value? In some cases?
    vf      <- do.call(varFun, as.list(args))
    
    .traceR(106, lbl = "if !is.null(dtInit): before")
    if (!is.null(dtInit)){
       .traceR(1061, lbl = "if !is.null(dtInit): executed", store = FALSE)
        dfAug    <- rbind(dtInit,dt1)
        vf.x     <- Initialize(vf, data=dfAug)}      # ... initialized
    .traceR(107, lbl = "if is.null(dtInit): before")     
    if (is.null(dtInit)) vf.x    <- Initialize(vf, data=dt1)


    vecR     <- sigma/varWeights(vf.x)              # AugDiagonal of R_i
    zz      <- sum(log(varWeights(vf.x)))           # Works fine. Work here if dtInit ne NULL
    .traceR(108, lbl = "zz")     
    }
    .traceR(10, lbl = "if len.clw after")                      
   
    
   
    if (!is.null(dtInit)){
        .traceR(222, "IF dtInit", store = FALSE)
        idxInit      <- c(1:nrow(dtInit))           # Indices for dtInit
        vecR         <- vecR[-idxInit]
    }                                               # Diagonal of R_i matrix
    

    ## return(zz)  # Continue to work here
    vecR2        <- vecR^2
    R            <- diag(vecR2, nrow=length(vecR))        # R_i matrix     
    Z            <- model.matrix(m$reStruc, data=dt1)     # Z_i matrix
    V            <- Z %*% D %*% t(Z) + R                  # V_i matrix
    predY        <- predict(modfit, dt1, level=0)         # Predict fixed
    
    ## formula(modfit)[[2]]
    dvName       <- as.character(formula(modfit)[[2]])

    r            <- dt1[[dvName]] - predY                 # Residuals
    n            <- nrow(dt1)                 # No. of obs for subject
    lLik         <- n*log(2*pi) + log(det(V)) + 
                    t(r) %*% solve(V) %*% r
    .traceR(1, lbl = "logLik1.lme ENDS <-")               
    return(-0.5 * as.numeric(lLik))
}
