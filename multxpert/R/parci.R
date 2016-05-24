# The ParCI function computes one-sided multiplicity-adjusted confidence
# intervals (simultaneous confidence intervals) for the single-step and
# step-down Dunnett procedures in one-sided testing problems with
# a balanced one-way layout and equally weighted null hypotheses
parci<-function(stat,n,est,stderror,covprob=0.975,proc)
# STAT, Vector of test statistics
# N, Common sample size in each treatment group
# EST, Vector of point estimates
# STDERROR, Vector of standard errors associated with the point estimates
# COVPROB, Simultaneous coverage probability
# PROC, Procedure name
    {
    # Number of null hypotheses
    m<-length(stat)

    if (m==0) stop("No test statistics are specified")

    if (m!=length(est)) stop("STAT and EST vectors have different lengths")
    if (m!=length(stderror)) stop("STAT and STDERROR vectors have different lengths")

    if (covprob>=1) stop("Simultaneous coverage probability must be <1")
    if (covprob<=0) stop("Simultaneous coverage probability must be >0")

    if(!all(proc %in% c("Single-step Dunnett", "Step-down Dunnett"))) stop("Procedure name is not recognized. ParCI function supports only the single-step Dunnett and step-down Dunnett procedures")

    if (n<=0) stop("Sample size must be positive")

    # number of procedures specified
    nproc <- length(proc)

    # set up matrix to contain confidence limits
    cimat <- matrix(0,m,nproc)
    dimnames(cimat) <- list(NULL, paste(proc, ".conf.limit", sep=""))

    # set up matrix to contain adjusted p-values
    adjpmat <- matrix(0,m,nproc)
    dimnames(adjpmat) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

    # Degrees of freedon
    nu<-(m+1)*(n-1)

    # Compute adjusted p-values
    result <- paradjp(stat,n,proc)

    #adjpmat <- result[, grep(".adj.pvalue", names(result), value=TRUE)]

    # One-sided familywise error rate
    alpha<-1-covprob

    # Rejection/acceptance of null hypotheses
    #reject<-(adjp<=alpha)

    # Vectors of confidence limits
    ci<-rep(0,m)

    zero<-rep(0,m)

	if (is.element("Single-step Dunnett", proc)) {
            adjpmat[, "Single-step Dunnett.adj.pvalue"] <- round(result[, "Single.step.Dunnett.adj.pvalue"], 4)

            reject <- (result[, "Single.step.Dunnett.adj.pvalue"] <= alpha)
           # Critical value

           c<-qdunnett(1-alpha,nu,m)
           cimat[, "Single-step Dunnett.conf.limit"] <-round(est-c*stderror, 4)
        }

	if (is.element("Step-down Dunnett", proc)) {
            adjpmat[, "Step-down Dunnett.adj.pvalue"] <- round(result[, "Step.down.Dunnett.adj.pvalue"], 4)

            reject <- (result[, "Step.down.Dunnett.adj.pvalue"] <= alpha)

	    # All null hypotheses are rejected
  	    if (sum(reject)==m) {
               # Critical value
               c<-qt(1-alpha,nu)
               cimat[, "Step-down Dunnett.conf.limit"] <- round(pmax(zero,est-c*stderror), 4)
            }

            # Some null hypotheses are accepted
  	    if (sum(reject)<m) {
                for (i in 1:m) {
                   if (reject[i]==1) cimat[i, "Step-down Dunnett.conf.limit"]<-0
                   if (reject[i]==0) {
                      # Critical value
                      c<-qdunnett(1-alpha,nu,m-sum(reject))
                      cimat[i, "Step-down Dunnett.conf.limit"] <- round(est[i]-c*stderror[i],4)
                   }
                }
            }
        }


    # Data frame returned by the function
    result<-data.frame(stat, est, stderror, adjpmat, cimat)
    names(result)[1]<-"Test.statistic"
    names(result)[2]<-"Estimate"
    names(result)[3]<-"Std.error"
    #names(result)[4]<-"Adj.pvalue"
    #names(result)[5]<-"Conf.limit"

    return(result=result)

    }
# End of parci
