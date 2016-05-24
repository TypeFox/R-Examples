#########################################################################
# Categorical Network Class Methods
# Distance between networks

setMethod("initialize", "catNetworkEvaluate", 
          function(.Object, nnodes, numsamples, nnets) {
            .Object@numnodes <- nnodes
            .Object@numsamples <- numsamples
            .Object@nets <- vector("list", nnets)
            .Object@complx <- rep(NA, nnets)
            .Object@loglik <- rep(NA, nnets)
            .Object@KLdist <- rep(NA, nnets)
            .Object@hamm <- rep(NA, nnets)
            .Object@hammexp <- rep(NA, nnets)
            .Object@tp <- rep(NA, nnets)
            .Object@fp <- rep(NA, nnets)
            .Object@fn <- rep(NA, nnets)
            .Object@sp <- rep(NA, nnets)
            .Object@sn <- rep(NA, nnets)
            .Object@fscore <- rep(NA, nnets)
            .Object@skel.tp <- rep(NA, nnets)
            .Object@skel.fp <- rep(NA, nnets)
            .Object@skel.fn <- rep(NA, nnets)
            .Object@order.fp <- rep(NA, nnets)
            .Object@order.fn <- rep(NA, nnets)
            .Object@markov.fp <- rep(NA, nnets)
            .Object@markov.fn <- rep(NA, nnets)
            .Object@time <- 0
            return(.Object)
            })

setMethod("show", "catNetworkEvaluate",
          function(object) {
            if(is(object, "catNetworkEvaluate"))
              str <- sprintf(
" Number of nodes    = %d, 
 Sample size        = %d,
 Number of networks = %d
 Processing time    = %.3f\n",                             
                  object@numnodes, 
                  object@numsamples,
                  length(object@nets),
                  object@time)
            cat(str, "\n")
            return(str)
            })

setMethod("show", "dagEvaluate",
          function(object) {
            if(is(object, "dagEvaluate"))
              str <- sprintf(
" Number of nodes     = %d, 
 Sample size         = %d,
 Number of DAGs      = %d
 Processing time     = %.3f\n",                             
                  object@numnodes, 
                  object@numsamples,
                  object@numDags,
                  object@time)
            cat(str, "\n")
            return(str)
            })

setMethod("cnPlot", "catNetworkEvaluate",
          function(object, file) {

            if(length(object@loglik) > 0 && length(object@complx) > 0 &&
               (length(object@tp) == 0 || is.na(object@tp[1]))) {
              par(mfrow=c(1,1))
              plot(object@complx, object@loglik, xlab="complx", ylab="log(loglik)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))
            }            
            else if(length(object@loglik) > 0 && length(object@complx) > 0 &&
                    length(object@tp) > 0 && !is.na(object@tp[1]) && 
                    length(object@hamm) > 0 && !is.na(object@hamm[1]) &&
                    length(object@markov.fp) > 0 && !is.na(object@markov.fp[1]) && 
                    length(object@hammexp) > 0 && !is.na(object@hammexp[1])) {
              par(mfrow=c(3,2))
              plot(object@complx, object@loglik, xlab="complx", ylab="log(loglik)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))
              ##plot(object@complx, object@KLdist,
              ##     xlab="complx", ylab="KL-dist", lty=1, 
              ##     main="Probability Distance")
              plot(object@complx, object@tp,
                   xlab="complx", ylab="TP", lty=1,
                   main="True Positives Directed Edges")
              plot(object@complx, object@hamm,
                   xlab="complx", ylab="Hamming", lty=1, 
                   main="Parent Matrix Distance")              
              plot(object@complx, object@markov.fp + object@markov.fn,
                   xlab="complx", ylab="Markov", lty=1, 
                   main="Markov Neighbor Distance")
              plot(object@complx, object@fp,
                   xlab="complx", ylab="FP", lty=1, 
                   main="False Positive Directed Edges")
              plot(object@complx, object@fn,
                   xlab="complx", ylab="FN", lty=1,
                   main="False Nagative Directed Edges")
            }
            else if(length(object@loglik) > 0 && length(object@complx) > 0 &&
                    length(object@tp) > 0 && !is.na(object@tp[1]) &&
                    length(object@fp) > 0 && !is.na(object@fp[1]) && 
                    length(object@hamm) > 0 && !is.na(object@hamm[1]) ) {
              par(mfrow=c(2,2))
              xx <- object@complx
              plot(xx, object@loglik, xlab="complx", ylab="log(loglik)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))
              plot(xx, object@hamm,
                   xlab="complx", ylab="hamm", lty=1, 
                   main="Hamming Distance")
              plot(xx, object@tp,
                   xlab="complx", ylab="TP", lty=1,
                   main="True Positive Directed Edges")
              plot(xx, object@fp,
                   xlab="complx", ylab="FP", lty=1,
                   main="False Positive Directed Edges")
              ##plot(xx, object@skel.tp,
              ##     xlab="complx", ylab="skeleton TP", lty=1,
              ##     main="True Positive Edges")
            }
            else if(length(object@nets) > 0) {
              loglik <- sapply(object@nets, function(net) net@loglik)
              complx <- sapply(object@nets, function(net) net@complx)
              par(mfrow=c(1,1))
              plot(complx, loglik, xlab="complx", ylab="log(loglik)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))              
            }
            
            })

setMethod("cnProcTime", "catNetworkEvaluate",
          function(object) {
          return(object@time)
        })

setMethod("cnParHist", "catNetworkEvaluate",
          function(object) {

            numnodes <- object@numnodes

            bfirst <- TRUE
            mhist <- matrix(rep(0, numnodes*numnodes), nrow=numnodes)
            for(bnet in object@nets) {
              if(!is(bnet, "catNetwork"))
                next
              if(bnet@numnodes != numnodes)
                next
              if(bfirst) {
                rownames(mhist) <- bnet@nodes
                colnames(mhist) <- bnet@nodes
                bfirst <- FALSE
              }
              mhist <- mhist + cnMatParents(bnet)
            }
            
            return(mhist)
        })

setMethod("cnParHist", "list",
          function(object) {

            numnodes <- 0
            
            for(bnet in object) {
              if(!is(bnet, "catNetwork"))
                next
              if(numnodes < 1) {
                numnodes <- bnet@numnodes
                mhist <- matrix(rep(0, numnodes*numnodes), nrow=numnodes)
                rownames(mhist) <- bnet@nodes
                colnames(mhist) <- bnet@nodes
              }
              if(bnet@numnodes != numnodes)
                next
              mhist <- mhist + cnMatParents(bnet)
            }
            
            return(mhist)
        })
