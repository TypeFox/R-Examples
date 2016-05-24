wfct <- function(expr)
{
expr <- deparse(substitute(expr))
 
## create new environment
newEnv <- new.env()
 
## get call
mc <- sys.calls()[[1]]
mcL <- as.list(mc)
 
## get data and write to newEnv
DATA <- mcL[["data"]]
DATA <- eval(DATA)
DATA <- as.list(DATA)
NAMES <- names(DATA)
for (i in 1:length(DATA)) assign(NAMES[i], DATA[[i]], envir = newEnv)
 
## get parameter, response and predictor names
formula <- as.formula(mcL[[2]])
VARS <- all.vars(formula)
RESP <- VARS[1]
RHS <- VARS[-1]
PRED <- match(RHS, names(DATA))
PRED <- names(DATA)[na.omit(PRED)]
 
## calculate variances for response values if "error" is in expression
## and write to newEnv
if (length(grep("error", expr)) > 0) {
y <- DATA[[RESP]]
x <- DATA[[PRED]]
## test for replication
if (!any(duplicated(x))) stop("No replicates available to calculate error from!")
## calculate error
error <- tapply(y, x, function(e) var(e, na.rm = TRUE))
error <- as.numeric(sqrt(error))
## convert to original repititions
error <- rep(error, as.numeric(table(x)))
assign("error", error, envir = newEnv)
}
 
## calculate fitted or residual values if "fitted"/"resid" is in expression
## and write to newEnv
if (length(grep("fitted", expr)) > 0 || length(grep("resid", expr)) > 0) {
mc2 <- mc
mc2$weights <- NULL
MODEL <- eval(mc2)
fitted <- fitted(MODEL)
resid <- residuals(MODEL)
assign("fitted", fitted, newEnv)
assign("resid", resid, newEnv)
}
 
## return evaluation in newEnv: vector of weights
OUT <- eval(parse(text = expr), envir = newEnv)
return(OUT)
}
