pest <-    
function(
x,
y,
indices,
family = gaussian,
tuning = list(lambda=TRUE, specific=FALSE, phi=0.5, grouped.fused=0.5, elastic=0.5, vs=0.5, spl=0.5),
weights,
offset,
start = NULL,
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
  if (family$family=="Gamma") family$link <- "log"

# offset, start
  if (missing(offset))
      offset <- 0

  if (missing(start))
      start <- NULL
#  names(start) <- xnames    
      
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
  if(any(abs(oml)<.0001, na.rm = TRUE) && control$adapted.weights) {
      control$adapted.weights <- FALSE
      warning("control$adapted.weights set to FALSE as at least one ML-estimate is too close to zero. \n")
      }
  control$oml <- oml    
      
# tuning parameters
  lambda <- TRUE
  phi <- grouped.fused <- elastic <- vs <- spl <- 0.5
  specific <- FALSE
  tuningnames <- c("lambda", "specific", "phi", "grouped.fused", "elastic", "vs", "spl")
  if (is.null(names(tuning))) names(tuning) <- tuningnames[1:length(tuning)]
  for (i in 1:length(tuning)) assign(names(tuning)[i], tuning[[i]])
  if (elastic <=0 || elastic >= 1 || vs <=0 || vs >= 1 || grouped.fused <=0 || grouped.fused >= 1) 
      stop("error in tuning argument.") 
  control$elastic <- elastic

# a.coefs + number.selectable parameters
  # estimate with small ridge penalty for more stable adaptive weights
  if (control$adapted.weights.ridge==FALSE) {
      est.adapt <- oml
      } else {
      est.adapt <- gvcmcatfitridge(x, y, weights = weights, family = family, control=control, offset = offset)$coefficients
      control$oml <- est.adapt
      }
  acoefs <- a.coefs(indices, control, beta=est.adapt, x)

  pp <- colSums(as.matrix(indices[-c(1,3),])!=0) # covariates with penalized coefficients
  n.sp <- sum(indices[1,which(pp>=1)]) - sum(abs(indices["index2",])*(1-indices["index2b",]))  
     if (any(indices["index6",]!=0)) {n.sp <- n.sp - length(which(colnames(indices)[which(indices["index6",]!=0)]=="L2"))}

# penalized model
if (n.sp>0) { 

    # definitions
    weight.const <- acoefs$w.adaptive * acoefs$w.cases * acoefs$w.categories * acoefs$w.pairwise
    weight.const[which(acoefs$continuous==1)] <- weight.const[which(acoefs$continuous==1)] * control$nu

    cross <- FALSE
    if (!is.logical(specific)) {
                   if (length(rle(acoefs$which.covariate)[[2]]) != length(specific) && length(weight.const)!= length(specific)) 
                               stop("Error in tuning argument specific. ") 
                   if (length(rle(acoefs$which.covariate)[[2]]) == length(specific)) specific <- rep(specific, times=rle(acoefs$which.covariate)[[1]])
                   weight.const <- weight.const * specific
                  } 
    if ( (lambda && is.logical(lambda)) || (length(lambda)>1) ) cross <- TRUE

    phis <- weight.function(phi=phi, grouped.fused, vs, spl, acoefs$phis.v, acoefs$phis.gf, acoefs$phis.vs, acoefs$phis.sp) # vorher phi=.5

    # cross-validation
    if (cross==TRUE) {
        if (control$tuning.criterion %in% c("deviance", "1SE")) {
          if (!is.list(control$K)){
          T.index <- split(sample(1:n), rep(1:control$K, length = n))
          } else {
          T.index <- control$K
          control$K <- length(T.index)
          }
          L.index <- lapply(T.index, function(i) setdiff(1:n, i))
        } else {T.index=NULL; L.index=NULL}
        
        cross <- cv.lambda(x, y, weights, family, control, acoefs, lambda,  
            phis, weight.const, start, offset, L.index, T.index, indices)

        lambda <- cross$lambda
        lambdas <- cross$lambdas
        score <- cross$score
        score.sd <- cross$score.sd
        path <- cross$coefs
        
        if (is.na(lambda)==TRUE) stop("Error in cross validation. \n")
    } else {score <- path <- score.sd <- NA}
    
    # model
    opt    <- gvcmcatfit(x, y, weights, family, control, acoefs$A, lambda, phis, 
                         weight.const, acoefs$which.a, start, offset)

    # plot
    path <- if (plot==TRUE){
      path.matrix(x, y, weights, family, control, acoefs$A, lambda, phis, weight.const, 
             acoefs$which.a, start, offset, opt$coefficients, path, oml)
      } else {NA}
    plot <- list(path=path, score=score, score.sd=score.sd)  

}

# nothing to be penalized...
if (n.sp==0){ 
     opt <- oml.model
     plot <- list(NA,NA)
    }

tuning <- list(lambda=lambda, specific=specific, phi=phi, grouped.fused=grouped.fused, elastic=elastic, vs=vs, spl=spl)
control$elastic <- NULL
       
# prepare output
reduction <- reduce(opt$coefficients, indices, control$assured.intercept)
x.reduced <- as.matrix(x %*% reduction$C)
x.reduction <- reduction$C
beta.reduced <- as.matrix(reduction$beta)
try(beta.refitted <- suppressWarnings(glm.fit(x.reduced, y, weights, family=family,
    intercept = FALSE)$coefficients), silent=TRUE)
if(!exists("beta.refitted")) beta.refitted <- rep(NA, times=length(beta.reduced))
beta.refitted <- round(beta.refitted, digits=control$digits)
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
    , test = est.adapt
    )

return(output)

}

