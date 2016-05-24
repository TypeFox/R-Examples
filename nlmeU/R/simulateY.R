## ->  simulateY  generic function
#' Simulates values of the dependent variable based on a model fit
#'
#' This function is generic; method functions can be written to handle specific classes of objects.
#'
#' @param object an object with a model fit for which dependent variable is to be simulated.
#' @param nsim number of simulations. nsim = 1 by default.
#' @param seed integer scalar used to initiate random numbers generator.
#' @param \dots some methods for this generic function may require additional arguments.
#' @param verbose logical. If TRUE basic information about arguments is provided. By default set to FALSE.
#' @param sigma numeric scalar. Allows to perform simulations employing alternative value of the scale parameter.
#' @return numeric matrix. Number of columns determined by nsim argument.
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @examples
#'
#'  ## simulateY (fm1)
#'
#' @export
simulateY <- function(object, nsim = 1, seed = NULL, ..., verbose = FALSE, sigma) UseMethod("simulateY")



#' @S3method simulateY lme
simulateY.lme <- function (object, nsim = 1, seed = as.integer(runif(1, 0, .Machine$integer.max)), ...,
  verbose = FALSE, sigma) 
{
# Data with one level of grouping only.
   .functionName <- "simulateY.lme"
   .traceR <- if (is.null(options()$traceR)) function(...){} else options()$.traceR      

  .traceR(1, lbl = "-> simulateY.lme STARTS")

if (verbose) print(paste("nsim = ", nsim ,", seed = ", seed, sep = ""))

if (!inherits(object, "lme"))  stop("Object must inherit from class \"lme\" ")

.traceR(110, lbl = "missing(sigma)")

if (!missing(sigma))  object$sigma <- sigma
.traceR(115, lbl = "after !missing(sigma) ")

   groups <- object$groups[[1]]
.traceR(120, lbl = "groups")

        ugroups <- unique(groups)
            individuals <- as.matrix(ugroups)
        if (is.numeric(individuals)) 
            individuals <- ugroups[individuals]
.traceR(130, lbl = "individuals")
Vlist <- nlme::getVarCov(object, individuals, type="marginal")
fitd0 <- fitted(object, level = 0)
chVlist <- lapply(Vlist, chol)

nx <- nsim * length(fitd0)

set.seed(seed)
noise.mtx <- matrix(rnorm(nx), nrow = length(fitd0), ncol = nsim)

.traceR(150, lbl = "lapply STARTS here ***")
dimN   <-   sapply(chVlist, ncol)  # Number of records per subject
cdimN1 <- cumsum(dimN)
cdimN0 <- cumsum(dimN) - dimN + 1
cdimN  <- cbind(cdimN0, cdimN1)
tList <- vector(length(dimN), mode = "list")
tList[] <- 1:length(dimN)
auxF1 <- 
  function(el){ # 1,2, ...
     .traceR(1001, lbl = paste("Local auxF1() STARTS. el =", el) , store = FALSE)
 
     chV <- chVlist[[el]]
     ix <- cdimN[el,]
     i0 <- ix[1]
     i1 <- ix[2]
     noisex <- noise.mtx[i0:i1, ]
     tx <-t(chV) %*% noisex   # Check transpose here
     .traceR(1002, lbl = paste("el=", el))
     .traceR(1001,lbl = paste("Local auxF1() ENDS. el =", el) , store = FALSE)
     tx
}
res <- lapply(tList, FUN= auxF1)
.traceR(150, lbl = "lapply ENDS here ***")

.traceR(160, lbl = "resAll STARTS***", store = FALSE)
resAll <- NULL 
for (i in 1:length(dimN)) resAll <- rbind(resAll, res[[i]])
.traceR(160, lbl = "resAll ENDS***")
.traceR(1, lbl = "simulateY.lme ENDS->")
return(resAll + fitd0)
}

sigmaTolme <- function(object, value){ 
 ### Use this function only with Pwr(), because it  corrupts lme.object
  funNm <- "sigmaTolme"
 
  .functionName <- "simgmaTolme"                                # Function name
  .traceR <- if (is.null(options()$traceR)) function(...){} else options()$.traceR 

  .traceR(1, lbl = "-> sigmaTolme STARTS")
  sigma0 <- object$sigma 
  val <- value * value
  sc  <- sqrt(val)/sigma0  
  object$sigma <- sqrt(val)
  #resids <- object$residuals
  #resids <- resids * sc  
  #std <- attr(resids,"std")*sc
  ## attr(object$residuals,"std") <- std
  
   attr(object$fixDF, "varFixFact") <- 
      sc*attr(object$fixDF, "varFixFact") # Rescaled for anova
  .traceR(5, vcov(object)*sc*sc, funNm)

   object$varFix  <- object$varFix*sc*sc  # vcov rescaled  
   .traceR(1, lbl = "sigmaTolme ENDS <-")
   object
}
