####  Test that wrong combinations of  varlist / doOne   give *good* error messages
####
require(simsalapar)
source(system.file("xtraR/assertErr-etc.R", package="simsalapar", mustWork=TRUE))
## Want English error/warning messages {works at least on Linux}:
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Must be fast, rather than "interesting":

# Implement the variable list using the funciton varlist
varSmini <- varlist(
  n.sim = list(type="N", expr=quote(N[sim]), value = 3),# 3: "minimal"
  n = list(type="grid", value=c(100,150)),
  # numbers of predictors
  p = list(value = 64),
  s0 = list(type="grid", value=c(3,15)),
  # number of replications (to compute single datapoint)
  N = list(type="frozen", value=50),
  # test level alpha MAYBE ALPHA WITH TYPE INNER?
  alpha = list(type="grid", value=c(0.05, 0.01)), # goal: inner
  # fraction of data used for variable screening in Multi Sample Split Method
  fraction = list(type="grid", expr=quote(mu), value=seq(4, 7, by = 0.5)/10),
  # number of sample splits in the MSS algorithm
  B = list(type="frozen", value=100),
  # the model selector in the MSS algorithm
  model.selector = list(type="frozen", value="lasso.cv")
)

varSmin.1 <- set.n.sim(varSmini, NULL) # No 'n.sim'

## apart from varlist, also has 'DBG' and '...'
do1.Dd <- function(n, p, s0, N, alpha, B, fraction,
                   model.selector, DBG=FALSE, ...)
{
    stopifnot(s0 <= p)
    ## construct the design matrix x and the coeff. vector beta
    ds <- list(x = matrix(rt(n*p, df=3), n,p),
               beta = c(rlnorm(s0), rep(0, p-s0)))
    x <- ds$x
    beta <- ds$beta
    if(DBG) { cat("generate..(): --> (x,beta) \n"); str(ds) }
                                        # check if beta is of the right form
    stopifnot(any(beta[(s0+1):p] == 0))

    if(DBG) {
        cat(sprintf(" .. N=%d, alpha=%g, B=%d, fraction=%g, model.sel.='%s'\n",
                N,alpha,B, fraction, model.selector))
        dots <- list(...)
        if(length(list(...)))
            cat(sprintf(" dots [len. %d]:\n", length(list(...)), deparse(substitute(...))))
        else cat(" '...' is empty\n")
    }
    FWER  <- rbinom( 1, size= 1, prob = 0.1)
    power <- rbinom(s0, size= N, prob = 0.9)/N
    ## initiate the coverage matrix
    coverage.matrix <- matrix(NA, nrow=ncol(x), ncol=N)
    for (i in 1:ncol(x)){
        for (j in 1:N){
            coverage.matrix[i,j] <- (rbinom(1, 10, pr = 0.2) <= 1 &
                                     rbinom(1, 10, pr = 0.2) <= 3)
        }
    }
                                        # calculate coverage rate for each predictor variable 1,..,p
    coverage <- apply(coverage.matrix, 1, mean)
                                        # return FWER, power, power.average and coverage
    c("FWER" = FWER, "power" = power, "coverage" = coverage)
} ## {do1.Dd}

do1.D <- do1.Dd; formals(do1.D) <- head(formals(do1.Dd), -1)
stopifnot(identical(body(do1.D), body(do1.Dd)),
          is.character(print(all.equal(do1.D, do1.Dd))))
##-> do1.D () has 'DBG' but *no* '...


## MM: before doLapply()  use   doCheck()  for a preliminary check:
set.seed(11); cc <- doCheck(do1.Dd, vList = varSmini , n = 1)
set.seed(11); c1 <- doCheck(do1.Dd, vList = varSmin.1, n = 1)
stopifnot(is.numeric(cc), !anyNA(cc),
          identical(c1, cc))
set.seed(12); cc. <- doCheck(do1.D, vList = varSmini , n = 1)
set.seed(12); c1. <- doCheck(do1.D, vList = varSmin.1, n = 1)
stopifnot(is.numeric(cc.), !anyNA(cc.),
          identical(c1., cc.))


options(warn = 2)#--> warnings give errors
system.time(
res1 <- doLapply(varSmin.1, DBG=TRUE, sfile=NULL,#"res_MSS_mini_5.rds",
                 check = FALSE, # <- no check now, did above ..
                 doOne = do1.Dd)
)
system.time(
resN <- doClusterApply(varSmini, cluster = parallel::makeCluster(2, "PSOCK"),
                       sfile=NULL, check=FALSE, doOne = do1.Dd)
)

## do1.D() instead of do1.Dd():
system.time(
res1. <- doLapply(varSmin.1, DBG=TRUE, sfile=NULL,
                 check = FALSE, doOne = do1.D)
)
system.time(
resN. <- doClusterApply(varSmini, cluster = parallel::makeCluster(2, "PSOCK"),
                        sfile=NULL, check=FALSE, doOne = do1.D)
)

stopifnot(
    identical(dim(resN),
              c(n = 2L, s0 = 2L, alpha = 2L, fraction = 7L, n.sim = 3L))
   ,
    identical(dim(res1), dim(resN)[1:4])
   ,
    identical(dimnames(res1), dimnames(resN)[1:4])
   ,
    identical(dimnames(resN),
              list(n = c("100", "150"), s0 = c("3", "15"),
                   alpha = c("0.05", "0.01"),
                   fraction = c("0.40", "0.45", "0.50", "0.55", "0.60", "0.65", "0.70"),
                   n.sim = NULL))
   ,
    identical(dim(resN), dim(resN.)),
    identical(dim(res1), dim(res1.)),
    identical(dimnames(resN), dimnames(resN.)),
    identical(dimnames(res1), dimnames(res1.))
    )


## getArray() does not work [directly], because 's0' influences size of result:
emsg <- tryCatch( getArray(res1), error=function(e) e$message)
stopifnot(identical(emsg, "\"value\" elements of 'x' differ in length"))
v.s03 <- getArray(resN[,"s0" =  "3",,,])
v.s15 <- getArray(resN[,"s0" = "15",,,])
tN <- getArray(resN, "time")

stopifnot(dim(v.s03) == c(68, 2, 2, 7, 3),
          dim(v.s15) == c(80, 2, 2, 7, 3),
          ## more relevantly,  res1 *is* identical to the first round of resN[] -- "seed" worked:
          all.equal(getArray(res1[, "3",,]), v.s03[,,,,1], tol = 1e-14)
          ,
          all.equal(getArray(res1[,"15",,]), v.s15[,,,,1], tol = 1e-14)
          ,
          ## no errors, no warnings:
          !getArray(resN, "error"), !getArray(resN, "warning"),
          ##
          tN >= 0, identical(dimnames(tN)[-2], dimnames(v.s03)[-1])
          )

varSm.no.s0 <- varSmini; varSm.no.s0$s0 <- NULL


## Now, with a typo in the call, we'd want to get a decent error message, but do not (YET ?) here:
system.time(
resE <- doClusterApply(varSmini, cluster = parallel::makeCluster(2, "PSOCK"),
                       sfile=NULL, doOne = do1.Dd, DGB = TRUE)
)
## but actually, there *is* no error, nor warning: the wrongly spelled argument
## is simply "dropped to the floor" ...
eM2 <- tryCatch(getArray(resE), error = function(e) e$message)
stopifnot(identical(eM2, "\"value\" elements of 'x' differ in length"))
## but these work fine
r03 <- getArray(resE[,s0 = "3",,,])
r15 <- getArray(resE[,s0 ="15",,,])
stopifnot(!anyNA(r03), !anyNA(r15),
          all.equal(r03, v.s03, tol=1e-14),
          all.equal(r15, v.s15, tol=1e-14))

##  with  do1.D (no '...'), this now *gives* a relevant warning [warn = 2 above ==> error]
rrr <- tryCatch.W.E(
    doClusterApply(varSmini, cluster = parallel::makeCluster(2, "PSOCK"),
                   sfile=NULL, doOne = do1.D, DGB = TRUE)
)
stopifnot(
    identical(names(rrr), c("value", "warning")),
    grepl("unused argument (DGB = TRUE)", rrr$warning$message, fixed=TRUE),
    ## all NULL values and 'error's are TRUE:
    vapply(lapply(rrr$value, `[[`, "value"), is.null, NA),
    all(getArray(rrr$value, "error"))
)

if(FALSE)# fails
mayplot(r03, vList=varSm.no.s0)

## Does error happen with doLapply()  -- no, not either with 'do1.Dd' (hmm ... FIXME?)
system.time(
reLE <- doLapply(varSmini, sfile=NULL, doOne = do1.Dd, DGB = TRUE)
)
## but it *does* print the 'messages'
## " Have argument names from doOne() not present ........ names:‘DBG’, ‘...’ "
##
## and  with  do1.D()  {no "..." argument}:
options(warn=1)
reL.. <- tryCatch.W.E(
    doLapply(varSmini, sfile=NULL, doOne = do1.D, DGB = TRUE)
)
stopifnot(
    identical(names(reL..), c("value", "warning")),
    grepl("unused argument (DGB = TRUE)", reL..$warning$message, fixed=TRUE),
    ## all NULL values and 'error's are TRUE:
    vapply(lapply(reL..$value, `[[`, "value"), is.null, NA),
    all(getArray(reL..$value, "error"))
)
