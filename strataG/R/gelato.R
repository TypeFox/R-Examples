#' @name gelato
#' @title GELATo - Group ExcLusion and Assignment Test
#' @description Run a GELATo test to evaluate assignment likelihoods of 
#'   groups of samples.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param unknown.strata a character vector listing to assign. Strata must 
#'   occur in \code{g}.
#' @param nrep number of permutation replicates for Fst distribution.
#' @param min.sample.size minimum number of samples to use to characterize 
#'   knowns. If the known sample size would be smaller than this after drawing 
#'   an equivalent number of unknowns for self-assignment, then the comparison 
#'   is not done.
#' @param num.cores number of CPU cores to use. Value is passed to 
#'   \code{\link[parallel]{mclapply}}.
#' @param gelato.result the result of a call to \code{gelato}.
#' @param unknown the name of an unknown stratum in the \code{x$likelihoods} 
#'   element.
#' @param main main label for top of plot.
#' 
#' @return A list with the following elements:
#' \tabular{ll}{
#'   \code{assign.prob} \tab a data.frame of assignment probabilities.\cr
#'   \code{likelihoods} \tab a list of likelihoods.\cr
#' }
#' 
#' @references O'Corry-Crowe, G., W. Lucey, F.I. Archer, and B. Mahoney. 2015. 
#'   The genetic ecology and population origins of the beluga whales of 
#'   Yakutat Bay. Marine Fisheries Review. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' stratify(msats.g, "fine")
#' 
#' gelato.fine <- gelato(msats.g, unk = "Offshore.South", nrep = 10)
#' 
#' gelatoPlot(gelato.fine, "Offshore.South")
#' 
#' @importFrom parallel mclapply
#' @importFrom stats sd dnorm median
#' @export
#' 
gelato <- function(g, unknown.strata, nrep = 1000, min.sample.size = 5, 
                   num.cores = 1) { 
  # Check unknown strata
  all.strata <- strata(g)
  unknown.strata <- unique(as.character(unknown.strata))
  if(!all(unknown.strata %in% all.strata)) {
    stop("Some 'unknown.strata' could not be found in 'g'")
  }
  
  knowns <- sort(setdiff(all.strata, unknown.strata))
  # loop through every unknown strata
  result <- sapply(unknown.strata, function(unknown) {
    unknown.gtypes <- g[, , unknown]
    unknown.mat <- as.matrix(unknown.gtypes)
    unknown.mat <- cbind(strata = unknown, unknown.mat)
    unknown.n <- nInd(unknown.gtypes)
    
    # loop through each known population and calculate distribution
    #   of Fst and log-likelihood of membership
    unknown.result <- sapply(knowns, function(known) {
      known.gtypes <- g[, , known]
      if((nInd(known.gtypes) - unknown.n) >= min.sample.size) {
        fst.dist <- do.call(rbind, mclapply(1:nrep, function(i) {
          # select samples to self assign
          ran.sample <- sample(indNames(known.gtypes), unknown.n)
          
          # extract gtypes of base known strata
          known.to.keep <- setdiff(indNames(known.gtypes), ran.sample)
          known.sample <- known.gtypes[known.to.keep, , ]
          
          # gtypes for observed Fst
          known.mat <- as.matrix(known.sample)
          known.mat <- cbind(strata = known, known.mat)
          obs.gtypes <- df2gtypes(rbind(known.mat, unknown.mat), ploidy = 2,
                                  id.col = NULL, strata.col = 1, loc.col = 2,
                                  sequences = sequences(g))

          # gtypes for null Fst
          st <- as.character(strata(known.gtypes))
          names(st) <- names(strata(known.gtypes))
          st[ran.sample] <- "<gelato.unknown>"
          null.gtypes <- stratify(known.gtypes, st)
          
          c(obs = unname(statFst(obs.gtypes)$result["estimate"]), 
            null = unname(statFst(null.gtypes)$result["estimate"])
          )
        }, mc.cores = num.cores))
        fst.dist <- fst.dist[apply(fst.dist, 1, function(x) all(!is.na(x))), ]

        if(nrow(fst.dist) < 2) {
          NULL
        } else {
          # summarize Fst distribution
          null.mean <- mean(fst.dist[, "null"])
          null.sd <- sd(fst.dist[, "null"])
          null.lik <- dnorm(fst.dist[, "obs"], null.mean, null.sd)
          log.Lik <- sum(log(null.lik), na.rm = T)     
          obs.median <- median(fst.dist[, "obs"], na.rm = T)
          obs.mean <- mean(fst.dist[, "obs"], na.rm = T)
          list(
            fst.dist = fst.dist, 
            log.Lik.smry = c(
              log.Lik = log.Lik, 
              mean.nreps = log.Lik / length(log.Lik),
              median = log(dnorm(obs.median, null.mean, null.sd)), 
              mean = log(dnorm(obs.mean, null.mean, null.sd))
            ),
            norm.coefs = c(mean = null.mean, sd = null.sd)
          )
        }
      } else NULL
    }, simplify = F)
    
    # calculate median logLikehood of assignment to each known
    log.Lik <- sapply(unknown.result, function(x) {
      if(is.null(x)) NA else x$log.Lik.smry["median"]
    })
    
    lik <- exp(log.Lik - max(log.Lik, na.rm = T))
    assign.prob <- lik / sum(lik, na.rm = T) 
    names(assign.prob) <- knowns
    
    list(assign.prob = assign.prob, likelihoods = unknown.result)
  }, simplify = F)
  
  assign.prob <- as.data.frame(t(sapply(result, function(x) x$assign.prob)))
  assign.prob$assignment <- apply(assign.prob, 1, function(x) {
    colnames(assign.prob)[which.max(x)]
  })
  result <- lapply(result, function(x) x$likelihoods)
  
  list(assign.prob = assign.prob, likelihoods = result)
}

#' @rdname gelato
#' @importFrom graphics par hist curve axis mtext
#' @importFrom stats dnorm
#' @export
#' 
gelatoPlot <- function(gelato.result, unknown, main = NULL) { 
  lik <- gelato.result$likelihoods[[unknown]]
  lik <- lik[!sapply(lik, is.null)]
  if(length(lik) == 0) {
    stop(paste("No likelihood distributions available for '", 
               unknown, "'", sep = ""))
  }
  xticks <- pretty(unlist(sapply(lik, function(x) x$fst.dist)))
  xlim <- range(xticks)
  op <- par(mar = c(3, 3, 3, 2) + 0.1, 
            oma = c(2, 2, 0.1, 0.1), 
            mfrow = c(length(lik), 1)
  )
  high.prob <- gelato.result$assign.prob[unknown, "assignment"]
  for(known in names(lik)) {
    known.lik <- lik[[known]]
    null.max <- max(hist(known.lik$fst.dist[, "null"], plot = F)$density)
    obs.max <- max(hist(known.lik$fst.dist[, "obs"], plot = F)$density)
    lik.mean <- known.lik$norm.coefs["mean"]
    lik.sd <- known.lik$norm.coefs["sd"]
    norm.lik <- dnorm(lik.mean, lik.mean, lik.sd)
    ylim <- range(pretty(c(0, null.max, obs.max, norm.lik)))
    hist(known.lik$fst.dist[, "null"], breaks = 10, freq = FALSE, 
         xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "", 
         col = "red", xaxt = "n")
    x <- NULL # To avoid R CMD CHECK warning about no global binding for 'x'
    curve(dnorm(x, lik.mean, lik.sd), from = xlim[1], to = xlim[2],
          add = TRUE, col = "black", lwd = 3, ylim = ylim)
    par(new = TRUE)
    hist(known.lik$fst.dist[, "obs"], breaks = 10, freq = FALSE, 
         xlim = xlim, ylim = ylim, xlab = "", ylab = "", col = "darkgreen", 
         main = "", xaxt = "n", yaxt = "n") 
    axis(1, pretty(xlim))
    ll.median <- known.lik$log.Lik.smry["median"]
    log.lik <- if(!is.infinite(ll.median)) {
      format(ll.median, digits = 4) 
    } else "Inf"
    p.val <- format(gelato.result$assign.prob[unknown, known], digits = 2)
    pop <- paste(known, " (lnL = ", log.lik, ", p = ", p.val, ")", sep = "")
    mtext(pop, side = 3, line = 1, adj = 1, 
          font = ifelse(known == high.prob, 2, 1))
  }
  mtext("Fst", side = 1, outer = T, cex = 1.2)
  mtext("Density", side = 2, outer = T, cex = 1.2)
  par(op)
  if(!is.null(main)) mtext(main, side = 3, line = 3, adj = 0, font = 3)
}