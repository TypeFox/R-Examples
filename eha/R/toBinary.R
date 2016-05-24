toBinary <- function(dat,
                     surv = c("enter", "exit", "event"),
                     strats,
                     max.survs = NROW(dat)) 
{
    if (!is.data.frame(dat))
      stop("dat must be a data frame")
    if (length(surv) != 3)
      stop("surv must have length 3")
    fixed.names <- names(dat)
    surv.indices <- match(surv, fixed.names)
    if (length(which(is.na(surv.indices)))) {
        x <- which(is.na(surv.indices))
        stop(paste(surv[x], " is not a name in the data frame."))
    }
    enter <- dat[, surv.indices[1]]
    exit <- dat[, surv.indices[2]]
    event <- dat[, surv.indices[3]]

    covars <- dat[, -surv.indices, drop = FALSE]

    nn <- NROW(dat)

    if (missing(strats) || is.null(strats)) strats <- rep(1, nn)
    rs <- risksets(Surv(enter, exit, event), strata = strats, max.survs)

    ## Remove this to keep risksets with no survivors:
    ## (Include this to remove risk sets with no survivors):
    ##weg <- (abs(rs$size - rs$n.events) > 0.01)
    ##rs$riskset <- rs$riskset[rep(weg, rs$size)]
    ##rs$eventset <- rs$eventset[rep(weg, rs$n.events)]
    ##rs$n.events <- rs$n.events[weg]
    ##rs$size <- rs$size[weg]
    ###################
    
    n.rs <- length(rs$size)
    ev <- numeric(sum(rs$size))
    start <- 1
    for (i in 1:n.rs) {
        ev[start:(start + rs$n.events[i] - 1)] <- 1
        start <- start + rs$size[i]
    }
    rs$ev <- ev

    out <- data.frame(event = rs$ev,
                      riskset = factor(rep(1:length(rs$size), rs$size)),
                      ##risktime = rep(rs$risktimes[weg], rs$size)
                      risktime = rep(rs$risktimes, rs$size)
                      )
                      
    out <- cbind(out, covars[rs$riskset, , drop = FALSE])
    out$orig.row <- (1:nn)[rs$riskset]
    out
}
