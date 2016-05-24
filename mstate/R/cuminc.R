`Cuminc` <- function(time, status, data, group, failcodes, na.status=c("remove","extra"), variance=TRUE)
{
    ## time
    if (!is.vector(time)) stop("argument \"time\" not of correct type")
    if (is.character(time)) {
        if (length(time) != 1)
            stop("single character string required for \"time\" argument")
        if (missing(data)) stop("argument \"data\" missing")
        time <- data[[time]]
    }
    ## status
    if (!is.vector(status)) stop("argument \"status\" not of correct type")
    if (is.character(status)) {
        if (length(status) == 1) {
            if (missing(data)) stop("argument \"data\" missing")
            status <- data[[status]]
        }
    }

    ## check for compatibility of time and status
    if (length(status) != length(time)) stop("lengths of time and status do not match")
    n <- length(time)

  if (missing(group)) {
    tmp <- data.frame(time=time, status=status)
    # Just call survfit with type="mstate"
    sf <- survfit(Surv(time, status, type="mstate") ~ 1, data=tmp)
    tt <- sf$time
    CIs <- sf$prev
    SSurv <- 1 - apply(CIs, 1, sum)
    ses <- sf$std.err
    sf1 <- survfit(Surv(time, status>0) ~ 1, data=tmp)
    se0 <- sf1$std.err
    res <- cbind(tt, SSurv, CIs, se0, ses)
    res <- as.data.frame(res)
    names(res) <- c("time", "Surv", paste("CI", sf$states, sep="."),
                    "seSurv", paste("seCI", sf$states, sep="."))
    res <- res[!duplicated(res$Surv), ]
    rownames(res) <- 1:nrow(res)
    class(res) <- c("Cuminc", "data.frame")
    attr(res, "survfit") <- sf
  } else {
    if (!(is.matrix(group)|(is.data.frame(group)))) { # then group should be a vector
      if (is.character(group)) { # character vector
        if (missing(data)) stop("argument \"data\" is missing, with no default")
        ngroup <- length(group)
        if (ngroup!=1) stop("only single grouping variable possible")
        kcols <- match(group,names(data))
        if (any(is.na(kcols))) stop("\"group\" column not in data")
        groupname <- group
        group <- data[,kcols]
      }
      else { # other vector, should be of length n
        ngroup <- 1
        groupname <- names(group)
        if (length(group) != n) stop("argument \"group\" has incorrect dimension")
      }
    }
    else {
      ngroup <- ncol(group)
      groupname <- names(group)
      if (nrow(group) != n) stop("argument \"group\" has incorrect dimension")
      if (ngroup!=1) stop("only single grouping variable possible")
      group <- group[,1] # coerce to vector
    }
    group <- as.factor(group)
    if (is.null(groupname)) groupname <- "group"
    tmp <- data.frame(time=time, status=status, group=group)
    # Call survfit with type="mstate"
    sf <- survfit(Surv(time, status, type="mstate") ~ group, data=tmp)
    tt <- sf$time
    CIs <- sf$prev
    SSurv <- 1 - apply(CIs, 1, sum)
    ses <- sf$std.err
    sf1 <- survfit(Surv(time, status>0) ~ group, data=tmp)
    se0 <- sf1$std.err
    res <- cbind(tt, SSurv, CIs, se0, ses)
    res <- as.data.frame(res)
    names(res) <- c("time", "Surv", paste("CI", sf$states, sep="."),
                    "seSurv", paste("seCI", sf$states, sep="."))
    group <- rep(levels(as.factor(tmp$group)), sf$strata)
    res <- cbind(data.frame(group=group), res)
    res <- res[!duplicated(res$Surv), ]
    rownames(res) <- 1:nrow(res)
    class(res) <- c("Cuminc", "data.frame")
    attr(res, "survfit") <- sf
  }
  return(res)
}
