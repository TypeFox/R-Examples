check.simulation <-
function (n, covariates, correlation = NULL, formula, 
coefficients, sd, seed)
{
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

if ( length(n)!=1 ||n<0 || !is.wholenumber(n))
     stop ("n must be a sinlge integer > 0. \n")

# covariates
b <- 0
for (i in 1:length(covariates)) {
if ( covariates[[i]][[1]]== "multinom" && length(covariates[[i]])!=3)
     stop ("one categorial variable is not well defined. \n")
if ( covariates[[i]][[1]]== "multinom" && length(covariates[[i]])==3 && 
     !(covariates[[i]][[3]] %in% c("nominal","ordinal")))
     stop ("one categorial variable is not well defined. \n")
if ( covariates[[i]][[1]]=="norm" ) { b <- b+1 }
}
if ( length(covariates) == 0 )      
     stop ("There are no covariates defined. \n")
         
# correlation
if (is.matrix(correlation)) {
     if (dim(correlation)[1]!=dim(correlation)[2]) 
         stop ("Error in argument 'correlation'. \n")
     if ((b > 0) && (b != dim(correlation)[1])) 
         stop ("Error in argument 'correlation'. \n")
     if ( length(which(diag(correlation)!=1)) > 0 )    
         stop ("Error in argument 'correlation'. \n")
     d <- correlation
     diag(d) <- 0
     if (max(d)>= 1 || min(correlation)<= -1 || min(correlation)<= -1 ||
         isSymmetric(correlation)==FALSE)    
         stop ("Error in argument 'correlation'. \n")
}
if (!is.null(correlation) && !is.matrix(correlation)) 
      stop ("'correlation' must be either NULL or a quadratic matrix. \n")

# rest
if (("formula" %in% is(formula))==FALSE)
     stop ("formula must be a formula. \n")

coefficients <- as.vector(coefficients)
if (!is.numeric(coefficients) || sum(as.integer(is.na(coefficients))) > 0) 
     stop ("Error in argument 'coefficients'. \n")

if ( !is.numeric(sd) || sd <= 0 ) 
     stop ("'sd' must be numeric and positive. \n")
     
if (!is.numeric(seed)) 
     stop ("seed must be numeric. \n")
           
}

