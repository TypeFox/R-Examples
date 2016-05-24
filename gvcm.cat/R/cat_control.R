cat_control <-
function (
center = FALSE, 
standardize = FALSE, 
accuracy = 2, digits = 4,
g = 0.5, epsilon = 10^(-5), maxi = 250, c = 10^(-5), gama = 20, steps = 25, nu = 1, 
tuning.criterion = "GCV", K = 5,
cv.refit = FALSE, 
lambda.upper=50, lambda.lower=0, lambda.accuracy=.01, scaled.lik=FALSE,
adapted.weights=FALSE, adapted.weights.adj = FALSE, adapted.weights.ridge=FALSE,
assured.intercept=TRUE, 
level.control = FALSE, 
case.control = FALSE,
pairwise = TRUE, # sollen pairwise differences pairwise differences sein? falls nicht: werden durch weighted adjacent differences ersetzt.
grouped.cat.diffs = FALSE,  # bei cat covariate: grouped auf differenzen.
bootstrap = 0, 
start.ml = FALSE,
L0.log = TRUE, # approximation with logistic function!! 
subjspec.gr = FALSE, # grouped fused penalty on effect modifier in model with subject specific coefficients?! betablocker
high = NULL, # unabled if NULL, sonst positive integer, gibt order differenzen bei nominalen variablen in penalty an. wird auf Null gesetzt falls pairwise=FALSE
...)

{

### center = FALSE, standardize = FALSE, accuracy = 2, digits = 4, start = NULL,

if (!is.logical(center))
     stop ("center must be logical. \n")

if (!is.logical(standardize))
     stop ("standardize must be logical. \n")

if (!is.numeric(accuracy) || accuracy < 0)
    stop("'accuracy' must be >= 0")
accuracy <- as.integer(accuracy)

if (!is.numeric(digits) || digits < 0)
    stop("'digits' must be >= 0")
digits <- as.integer(digits)

### g = 0.5, epsilon = 10^(-5), maxi = 250, steps = 25, c = 10^(-5), gama = 30, nu = 1, 

if (!is.numeric(g) || length(as.vector(g))!=1 || g>1 || g<0)
     stop ("g must be a single, numeric value out of ]0;1[. \n")

if (!is.numeric(epsilon) || epsilon <= 0 || epsilon>1 || length(epsilon)!=1)
   stop ("epsilon is ment to be a single, small, positive and numeric value. \n")

if (!is.numeric(maxi) || length(maxi)!=1 || maxi < 1)
     stop ("maxi must be a sinlge integer > 0 \n")
maxi <- as.integer(maxi)

if (!is.numeric(steps) || length(steps)!=1 || steps < 1)
     stop ("steps must be a sinlge integer > 0 \n")
steps <- as.integer(steps)

if (!is.numeric(c) || c<0 || length(c)!=1)
    stop ("c is ment to be a single, small, positive and numeric value. \n")

if (!is.numeric(gama) || length(as.vector(gama))!=1 || gama<0)
     stop ("gama must be a singl, positive and numeric value. \n")

###  tuning.criterion = "deviance", K = 5, cv.refit = FALSE, 

if (!(tuning.criterion %in% c("GCV", "UBRE", "deviance", "1SE")))
     stop ("tuning.criterion must be one out of 'GCV', 'UBRE', 'deviance'. \n")

if (!is.vector(K))
         stop ("K must be a single integer > 1. \n")
if (!is.list(K)) {
    if (!is.numeric(K) || K<2)
         stop ("K must be a single integer > 1. \n")
    }

if (!is.logical(cv.refit))
     stop ("cv.refit must be logical. \n")
     
###

if (!is.numeric(lambda.upper) || length(lambda.upper)!=1)
     stop ("lambda.upper must be numeric and positive \n")

if (!is.numeric(lambda.lower) || length(lambda.lower)!=1)
     stop ("lambda.lower must be numeric and positive \n")
     
if (!is.numeric(lambda.accuracy) || lambda.accuracy <= 0 || length(lambda.accuracy)!=1)
   stop ("lambda.accuracy is ment to be a single, positive and numeric value. \n")

if (!is.logical(scaled.lik))
     stop ("scaled.lik must be logical. \n")

if (!is.logical(adapted.weights.ridge))
     stop ("adapted.weights.ridge must be logical. \n")

if (!is.logical(adapted.weights.adj))
     stop ("adapted.weights.adj must be logical. \n")

if (!is.logical(adapted.weights))
     stop ("adapted.weights must be logical. \n")

if (!is.logical(assured.intercept))
     stop ("assured.intercept must be logical. \n")

if (!is.logical(pairwise))
     stop ("pairwise must be logical. \n")

if (!is.logical(grouped.cat.diffs))
     stop ("grouped.cat.diffs must be logical. \n")

### bootstrap = 0, level.control = FALSE, case.control = FALSE,
    
if (!is.numeric(bootstrap) || bootstrap<0 || length(bootstrap)!=1)
     stop ("bootstrap must be positive and numeric; value zero or >10. \n")
if (bootstrap!=0 && bootstrap<10)
    bootstrap <- 10
    
if (!is.logical(level.control))
     stop ("level.control must be logical. \n")

if (!is.logical(case.control))
     stop ("case.control must be logical. \n")

if (!is.logical(start.ml))
     stop ("start.ml must be logical. \n")

if (!is.logical(L0.log))
     stop ("L0.log must be logical. \n")

if (!is.logical(subjspec.gr))
     stop ("subjspec.gr must be logical. \n")
     
if (pairwise == FALSE) {high <- NULL}     
if (!is.null(high) && !is.numeric(high) )  
     stop ("high must be NULL or a positive integer. \n")
if (is.numeric(high) && high < 1)
     stop ("high must be NULL or a positive integer. \n")
if (is.numeric(high)) high <- as.integer(high)     
    

###  output

list(
center = center, 
standardize = standardize, 
accuracy = accuracy, digits = digits,
g = g, epsilon = epsilon, maxi = maxi, c = c, gama = gama, steps = steps, nu = nu, 
tuning.criterion = tuning.criterion, K = K,
cv.refit = cv.refit, 
lambda.upper = lambda.upper, lambda.lower = lambda.lower,
lambda.accuracy = lambda.accuracy, scaled.lik = scaled.lik,
adapted.weights = adapted.weights, adapted.weights.adj = adapted.weights.adj, adapted.weights.ridge=adapted.weights.ridge,
assured.intercept = assured.intercept, 
level.control = level.control, 
case.control = case.control,
pairwise = pairwise, 
grouped.cat.diffs = grouped.cat.diffs, 
bootstrap = bootstrap,
start.ml = start.ml,
L0.log = L0.log,
subjspec.gr = subjspec.gr, 
high = high
)

}

