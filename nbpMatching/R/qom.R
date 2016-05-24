#'Quality of Match
#'
#'Quality of matches show how well matched pairs differ.  For each variable the
#'average distance is generated.  Each item in a pair is assigned a group and
#'after several iterations the quantile of these average distances is returned.
#'
#'This fuction is useful for determining the effectiveness of your weights
#'(when generating a distance matrix).  Weighting a variable more will lower
#'the average distance, but it could penalize the distance of the other
#'variables. Calculating the standard error requires calling
#'\code{\link{hdquantile}} from \pkg{Hmisc}.  The quantiles may be slighly
#'different when using \code{\link{hdquantile}}.
#'
#'@aliases qom qom,data.frame,data.frame-method qom,data.frame,nonbimatch-method
#'@param covariate A data.frame object.
#'@param matches A data.frame or nonbimatch object.  Contains information on
#'how to match the covariate data set.
#'@param iterations An integer.  Number of iterations to run, defaults to
#'10,000.
#'@param probs A numeric vector.  Probabilities to pass to the quantile
#'function.
#'@param use.se A logical value.  Determines if the standard error should be
#'computed.  Default value of FALSE.
#'@param all.vals A logical value.  Determines if false matches should be
#'included in comparison.  Default value of FALSE.
#'@param seed Seed provided for random-number generation.  Default value of
#'101.
#'@param \dots Additional arguments, not used at the moment.
#'@return a list object containing elements with quality of match information
#'
#'  \item{q}{data.frame with quantiles for each covariate}
#'
#'  \item{se}{data.frame with standard error for each covariate}
#'
#'  \item{sd}{vector with standard deviate for each covariate}
#'@exportMethod qom
#'@author Cole Beck
#'@examples
#'
#'df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
#'df.dist <- gendistance(df, idcol=1)
#'df.mdm <- distancematrix(df.dist)
#'df.match <- nonbimatch(df.mdm)
#'qom(df.dist$cov, df.match)
#'qom(df.dist$cov, df.match$matches)
#'

setGeneric("qom", function(covariate, matches, iterations=10000, probs=NA,
           use.se=FALSE, all.vals=FALSE, seed=101, ...) standardGeneric("qom"))
setMethod("qom", signature(covariate="data.frame", matches="data.frame"), function(covariate, matches, iterations=10000,
          probs=NA, use.se=FALSE, all.vals=FALSE, seed=101, ...) {
    if(exists(".Random.seed", envir = .GlobalEnv)) {
        save.seed <- get(".Random.seed", envir= .GlobalEnv)
        on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    n <- nrow(matches)
    if(n%%2 == 1) {
        stop("There must be an even number of elements")
    }
    if(!all(sapply(matches[,c(2,4)], is.numeric))) {
        stop("matches must contain numeric values in columns two and four")
    }
    # digits can be specified for rounding
    xarg <- list(...)
    if('digits' %in% names(xarg)) {
      digits <- xarg$digits
    } else {
      digits <- 4
    }
    pairs <- matches[matches[,2] < matches[,4], c(2,4)]
    if(all.vals != TRUE) {
        # obtain actual pairs (meaning not matched to phantoms)
        pairs <- pairs[pairs[,2] <= nrow(covariate),]
    }
    npairs <- nrow(pairs)
    ignorecols <- which(sapply(covariate, FUN=function(x) { length(setdiff(x, suppressWarnings(as.numeric(x)))) > 0 }))
    if(length(ignorecols) > 0) covariate <- covariate[, -ignorecols]
    if(is.na(probs) || !is.numeric(probs)) probs <- c(0,25,50,75,90,95,100)/100
    probs <- sort(probs[0 <= probs & probs <= 1])
    if(length(probs) == 0) return(numeric(0))

    worst <- NULL
    pair.sd <- pair.se <- NULL
    # check for ~100% percentile, and remove from probs
    if(1 - probs[length(probs)] < 0.000001) {
        probs <- probs[-length(probs)]
        # 100th percentile is special case, force worst group by placing min value in group 1
        worst <- matrix(colMeans(abs(covariate[pairs[,1],]-covariate[pairs[,2],]), na.rm=TRUE))
        if(all.vals) {
            good.indeces <- which(pairs[,2] <= nrow(covariate))
            # find values that matched phantoms
            others <- covariate[pairs[-good.indeces,1],]
            md <- function(x,y) abs(mean(x, na.rm=TRUE)-mean(y, na.rm=TRUE))
            group.means <- matrix(rowMeans(sapply(apply(pairs[good.indeces,], MARGIN=1, FUN=function(i) covariate[unlist(i),]), FUN=function(j) apply(j, MARGIN=2, sort))), nrow=2)
            for(i in seq(ncol(others))) {
                g1 <- c(rep(group.means[1,i], length(good.indeces)), sort(others[,i]))
                g2 <- c(rep(group.means[2,i], length(good.indeces)), rep(NA, nrow(others)))
                mydiff <- md(g1, g2)
                for(j in seq(length(g1), by=-1, length.out=nrow(others))) {
                    g2[j] <- g1[j]
                    g1[j] <- NA
                    newdiff <- md(g1, g2)
                    if(mydiff >= newdiff) {
                        break
                    } else {
                        mydiff <- newdiff
                    }
                }
                worst[i,1] <- mydiff
            }
        }
        dimnames(worst) <- list(names(covariate), "100%")
    }

    if(length(probs)) {
        if(is.na(iterations) || !is.numeric(iterations) || iterations < 2) iterations <- 10000
        if(!is.numeric(seed)) seed <- 101
        set.seed(seed)
        choices <- matrix(sample(c(-1,1), npairs*iterations, replace=TRUE), ncol=iterations)

        if(all.vals == TRUE) {
            # use alternate, slower method when all values (including NAs) are present
            sp <- seq(n)
            pairdiff.sums.mat <- t(apply(choices, MARGIN=2, FUN=function(x) {
                g1 <- ifelse(x < 0, pairs[,1], pairs[,2])
                g2 <- setdiff(sp, g1)
                abs(colMeans(covariate[g1,], na.rm=TRUE) - colMeans(covariate[g2,], na.rm=TRUE))
            }))
        } else {
            group.one <- sapply(covariate[pairs[,1],], FUN=function(x) {
              colMeans(choices * x, na.rm=TRUE)
            })
            group.two <- sapply(covariate[pairs[,2],], FUN=function(x) {
              colMeans(choices * x, na.rm=TRUE)
            })
            pairdiff.sums.mat <- abs(group.one - group.two)
        }
        if(use.se == TRUE) {
            se.parts <- t(apply(pairdiff.sums.mat, MARGIN=2, FUN=function(x) { tmp<-hdquantile(x, probs=probs, se=TRUE); c(tmp, attr(tmp, "se")) }))
            pairsumm <- se.parts[,seq_along(probs)]
            pair.se <- round(se.parts[,seq(length(probs)+1, length.out=length(probs))], digits)
            # change the format of the colnames to be like "quantile"
            cnames <- sprintf("%.0f%%", as.numeric(colnames(pairsumm))*100)
            colnames(pairsumm) <- cnames
            colnames(pair.se) <- cnames
        } else {
            # it's cheap to add probs=0; do so and remove first row to ensure matrix result
            pairsumm <- t(apply(pairdiff.sums.mat, MARGIN=2, quantile, probs=c(0,probs)))[,-1, drop=FALSE]
        }
        if(!is.null(worst)) {
            pairsumm <- cbind(pairsumm, worst)
        }
        pair.sd <- round(apply(pairdiff.sums.mat, MARGIN=2, sd, na.rm=TRUE), digits)
    } else {
        # only return the 100th percentile
        pairsumm <- worst
    }
    pairsumm <- round(pairsumm, digits)
    list(q=pairsumm, se=pair.se, sd=pair.sd)
})

setMethod("qom", signature(covariate="data.frame", matches="nonbimatch"), function(covariate, matches, iterations=10000,
          probs=NA, use.se=FALSE, all.vals=FALSE, seed=101, ...) {
    qom(covariate, matches$matches, iterations, probs, use.se, all.vals, seed, ...)
})
