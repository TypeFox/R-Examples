abc <-
function(
x, 
y, 
indices,
family = gaussian,
tuning = c("AIC", "BIC"), 
weights,
offset, 
start, 
control = cat_control(),
plot=FALSE,
...
)

{
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  n <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0

# checks
  if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
      family <- family()
  if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
  }

# offset, start
  if (missing(offset))
      offset <- 0

  if (missing(start))
      start <- rep(0, nvars)
  names(start) <- xnames    
      
  if (missing(control))
      control <- cat_control(...)
    
  if (is.null(weights))
      weights <- rep.int(1, n)
      if (length(weights)!=NROW(x) || !is.vector(weights) || !is.numeric(weights))
      stop("Error in input weights. ") 
      
  if (!is.logical(plot) || !is.matrix(x) || !is.numeric(x) || !is.vector(y) || !is.numeric(y) || nrow(x)!=length(y) || !is.matrix(indices))
       stop ("Error in input arguments. \n")

  if (!is.list(control$K) && control$K > n)
      ("K must be a single integer < dim(data)[1]. \n")

# oml
  try(oml.model <- glm.fit(x, y, weights = weights,
      family = family, intercept = FALSE), silent = TRUE)
  if(!exists("oml.model")) {
      coefs <- rep(NA, nvars)
      names(coefs) <- xnames
      oml.model <- list(coefficients = coefs, residuals = NA, fitted.values = NA,
          rank = NA, family = family, linear.predictors = NA, deviance = NA, aic = NA,
          null.deviance = NA, 
          iter = 0L, weights = NA, prior.weights = weights,
          df.residual = NA, df.null = NA, y = NA, converged = FALSE,
          boundary = TRUE)
      warning("Ordinary maximum likelihood estimate does not exsist. \n")
      control$adapted.weights <- FALSE
      }    
  oml <- round(oml.model$coefficients, control$digits)
  if(any(is.na(oml))) {
      control$adapted.weights <- FALSE
      warning("Ordinary maximum likelihood estimate contains NAs. \n")
      }

# tuning/method
  tuning <- match.arg(tuning)
  if (!(tuning %in% c("AIC", "BIC")))
       stop ("method/tuning is incorrect. \n")

# acoefs / n.sp
  pp <- colSums(as.matrix(indices[-c(1,3),])!=0) # covariates with penalized coefficients
  n.sp <- sum(indices[1,which(pp==1)]) - sum(abs(indices["index2",])*(1-indices["index2b",])) 
    

# model selection
if (tuning == "AIC") {criterion <- function(beta){0} } else
                     {criterion <- function(beta){(log(n)-2)*length(beta)} } 

if (n.sp>0){
    
    # model
    suppressWarnings(try(opt <- abcfit(x, y, weights, family, 
                               tuning, control, indices, start, offset)))
    if(exists("opt")==FALSE) {stop("Variable selection via AIC/BIC failed. \n")}
    
    # plot
    if (plot==TRUE) {plot <- opt$A.models} else {plot <- NA}
    
    }

# nothing to be penalized...
if (n.sp==0){ 
     opt <- oml.model
     plot <- list(NA,NA)
    }

# tuning
tuning <- round(opt$abc.model, digits=2)

# prepare output
reduction <- reduce(opt$coefficients, indices, control$assured.intercept)
x.reduced <- as.matrix(x %*% reduction$C)
x.reduction <- reduction$C
beta.reduced <- as.matrix(reduction$beta)
try(beta.refitted <- suppressWarnings(glm.fit(x.reduced, y, weights, 
     # offset = offset, start = as.vector(t(A.aspirants[[i]])%*%initials), 
     family = family, intercept = FALSE)$coefficients), silent=TRUE) 
if(!exists("beta.refitted")) beta.refitted <- rep(NA, times=length(beta.reduced))
beta.refitted <- round(beta.refitted, digits=control$accuracy)
names(beta.refitted) <- names(beta.reduced)

# output
output <- list(
    coefficients = opt$coefficients,
    coefficients.reduced = beta.reduced,
    coefficients.refitted = beta.refitted,
    coefficients.oml = oml,

    residuals = opt$residuals,
    fitted.values = opt$fitted.values,
    effects = opt$effects, 
    R = opt$R, 
    rank = opt$rank,
    qr = opt$qr,
    family = family,
    linear.predictors = opt$linear.predictors,
    deviance = opt$deviance,
    aic = opt$aic,
    null.deviance = opt$null.deviance,
    iter = opt$iter,
    weights = opt$weights, prior.weights = opt$prior.weights,
    df.residual = opt$df.residual,
    df.null = opt$df.null,
    converged = opt$converged,
    boundary = opt$boundary,
    offset = offset,
    control = control,
    contrasts = options("contrasts"),
    na.action  = "na.omit",
    plot = plot,
    tuning = tuning,
    indices = indices,
    number.selectable.parameters = n.sp,
    number.removed.parameters = nvars-ncol(x.reduced),
    x.reduction = x.reduction,
    beta.reduction = reduction$A
    )

return(output)
   
}

