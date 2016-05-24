pairwise <- function(treat,
                     event, n, mean, sd, TE, seTE, time,
                     data=NULL, studlab, ...){
  
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch studlab, treat, event, n, mean, sd, time from data:
  ##
  mf <- match.call()
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  treat <- eval(mf[[match("treat", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  event <- eval(mf[[match("event", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  mean <- eval(mf[[match("mean", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  sd <- eval(mf[[match("sd", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  seTE <- eval(mf[[match("seTE", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  time <- eval(mf[[match("time", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  
  
  if (is.null(treat))
    stop("Argument 'treat' mandatory.")
  ##
  if (is.list(treat))
    chklist(treat)
  ##
  if (!is.null(event))
    if (is.list(event))
      chklist(event)
    else
      meta:::chknumeric(event)
  ##
  if (!is.null(n))
    if (is.list(n))
      chklist(n)
    else
      meta:::chknumeric(n)
  ##
  if (!is.null(mean))
    if (is.list(mean))
      chklist(mean)
    else
      meta:::chknumeric(mean)
  ##
  if (!is.null(sd))
    if (is.list(sd))
      chklist(sd)
    else
      meta:::chknumeric(sd)
  ##
  if (!is.null(TE))
    if (is.list(TE))
      chklist(TE)
    else
      meta:::chknumeric(TE)
  ##
  if (!is.null(seTE))
    if (is.list(seTE))
      chklist(seTE)
    else
      meta:::chknumeric(seTE)
  ##
  if (!is.null(time))
    if (is.list(time))
      chklist(time)
    else
      meta:::chknumeric(time)
  
  
  if (!is.null(event) & !is.null(n) &
      is.null(mean) & is.null(sd) &
      is.null(TE) & is.null(seTE) &
      is.null(time))
    type <- "binary"
  else if (is.null(event) & !is.null(n) &
           !is.null(mean) & !is.null(sd) &
           is.null(TE) & is.null(seTE) &
           is.null(time))
    type <- "continuous"
  else if (!is.null(event) & is.null(n) &
           is.null(mean) & is.null(sd) &
           is.null(TE) & is.null(seTE) &
           !is.null(time))
    type <- "count"
  else if (is.null(event) & is.null(n) &
           is.null(mean) & is.null(sd) &
           !is.null(TE) & !is.null(seTE) &
           is.null(time))
    type <- "generic"
  else
    stop("Type of outcome unclear. Please provide the necessary information:\n  - event, n (binary outcome)\n  - n, mean, sd (continuous outcome)\n  - TE, seTE (generic outcome)\n  - event, time (incidence rates).")
  
  
  
  
  
  ##
  ## Transform long format to list format
  ##
  treat.list <- list()
  event.list <- list()
  n.list     <- list()
  mean.list  <- list()
  sd.list    <- list()
  TE.list    <- list()
  seTE.list  <- list()
  time.list  <- list()
  ##
  if (type == "binary") {
    listformat <- is.list(event) & is.list(n)
    if (!listformat){
      if (is.null(studlab))
        stop("Argument 'studlab' mandatory if argument 'event' is a vector.")
      ##
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      ##
      tdat <- data.frame(studlab, treat, event, n, stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        event.list[[i]] <- tres.i$event
        n.list[[i]]     <- tres.i$n
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      event <- event.list
      n     <- n.list
    }
  }
  else if (type == "continuous") {
    listformat <- is.list(n) & is.list(mean) & is.list(sd)
    if (!listformat){
      if (is.null(studlab))
        stop("Argument 'studlab' mandatory if argument 'mean' is a vector.")
      ##
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      ##
      tdat <- data.frame(studlab, treat, n, mean, sd, stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        n.list[[i]]     <- tres.i$n
        mean.list[[i]]  <- tres.i$mean
        sd.list[[i]]    <- tres.i$sd
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      n     <- n.list
      mean  <- mean.list
      sd    <- sd.list
    }
  }
  else if (type == "count") {
    listformat <- is.list(event) & is.list(time)
    if (!listformat){
      if (is.null(studlab))
        stop("Argument 'studlab' mandatory if argument 'event' is a vector.")
      ##
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      ##
      tdat <- data.frame(studlab, treat, event, time, stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        event.list[[i]] <- tres.i$event
        time.list[[i]]  <- tres.i$time
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      event <- event.list
      time  <- time.list
    }
  }
  else if (type == "generic") {
    listformat <- is.list(TE) & is.list(seTE)
    if (!listformat){
      if (is.null(studlab))
        stop("Argument 'studlab' mandatory if argument 'TE' is a vector.")
      ##
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      ##
      tdat <- data.frame(studlab, treat, TE, seTE, stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        TE.list[[i]]    <- tres.i$TE
        seTE.list[[i]]  <- tres.i$seTE
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      TE    <- TE.list
      seTE  <- seTE.list
    }
  }
  
  
  
  
  
  ##
  ## Check and set study labels
  ##
  if (is.null(studlab))
    studlab <- seq(along=treat[[1]])
  ##
  if (length(studlab) != length(unique(studlab)))
    stop("Study labels must all be distinct.")
  ##
  levs <- studlab
  
  
  narms <- length(treat)
  
  
  if (type=="binary"){
    if (length(event) != narms)
      stop("Different length of lists 'treat' and 'event'.")
    if (length(n) != narms)
      stop("Different length of lists 'treat' and 'n'.")
    ##
    for (i in 1:(narms-1)){
      ##
      if (i==1 & (length(treat[[i]]) != length(event[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'event'.")
      if (i==1 & (length(event[[i]]) != length(n[[i]])))
        stop("Different length of element ", i, " of lists 'event' and 'n'.")
      ##
      for (j in (i+1):narms){
        ##
        if (length(treat[[j]]) != length(event[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'event'.")
        if (length(event[[j]]) != length(n[[j]]))
          stop("Different length of element ", j, " of lists 'event' and 'n'.")
        ##
        dat <- data.frame(TE=NA, seTE=NA,
                          studlab=studlab,
                          treat1=treat[[i]],
                          treat2=treat[[j]],
                          event1=event[[i]], n1=n[[i]],
                          event2=event[[j]], n2=n[[j]])
        ##
        dat <- dat[!(is.na(dat$event1) & is.na(dat$n1)),]
        dat <- dat[!(is.na(dat$event2) & is.na(dat$n2)),]
        ##
        if (nrow(dat) > 0){
          m1 <- metabin(dat$event1, dat$n1,
                        dat$event2, dat$n2, ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          if (i==1 & j==2)
            res <- dat
          else
            res <- rbind(res, dat)
        }
        else
          if (i==1 & j==2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  
  
  if (type=="continuous"){
    if (length(n) != narms)
      stop("Different length of lists 'treat' and 'n'.")
    if (length(mean) != narms)
      stop("Different length of lists 'treat' and 'mean'.")
    if (length(sd) != narms)
      stop("Different length of lists 'treat' and 'sd'.")
    ##
    for (i in 1:(narms-1)){
      ##
      if (i==1 & (length(treat[[i]]) != length(n[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'n'.")
      if (i==1 & (length(treat[[i]]) != length(mean[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'mean'.")
      if (i==1 & (length(treat[[i]]) != length(sd[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'sd'.")
      ##
      for (j in (i+1):narms){
        ##
        if (length(treat[[j]]) != length(n[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'n'.")
        if (length(treat[[j]]) != length(mean[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'mean'.")
        if (length(treat[[j]]) != length(sd[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'sd'.")
        ##
        dat <- data.frame(TE=NA, seTE=NA,
                          studlab=studlab,
                          treat1=treat[[i]],
                          treat2=treat[[j]],
                          n1=n[[i]], mean1=mean[[i]], sd1=sd[[i]],
                          n2=n[[j]], mean2=mean[[j]], sd2=sd[[j]])
        dat <- dat[!(is.na(dat$n1) & is.na(dat$mean1) & is.na(dat$sd1)),]
        dat <- dat[!(is.na(dat$n2) & is.na(dat$mean2) & is.na(dat$sd2)),]
        ##
        if (nrow(dat) > 0){
          m1 <- metacont(dat$n1, dat$mean1, dat$sd1,
                         dat$n2, dat$mean2, dat$sd2,
                         ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          if (i==1 & j==2)
            res <- dat
          else
            res <- rbind(res, dat)
        }
        else
          if (i==1 & j==2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  
  
  if (type=="generic"){
    if (length(TE) != narms)
      stop("Different length of lists 'treat' and 'TE'.")
    if (length(seTE) != narms)
      stop("Different length of lists 'treat' and 'seTE'.")
    ##
    for (i in 1:(narms-1)){
      ##
      if (i==1 & (length(treat[[i]]) != length(TE[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'TE'.")
      if (i==1 & (length(treat[[i]]) != length(seTE[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'seTE'.")
      ##
      for (j in (i+1):narms){
        ##
        if (length(treat[[j]]) != length(TE[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'TE'.")
        if (length(treat[[j]]) != length(seTE[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'seTE'.")
        ##
        dat <- data.frame(TE=NA, seTE=NA,
                          studlab=studlab,
                          treat1=treat[[i]],
                          treat2=treat[[j]],
                          TE1=TE[[i]], seTE1=seTE[[i]],
                          TE2=TE[[j]], seTE2=seTE[[j]])
        dat <- dat[!(is.na(dat$TE1) & is.na(dat$seTE1)),]
        dat <- dat[!(is.na(dat$TE2) & is.na(dat$seTE2)),]
        ##
        if (nrow(dat) > 0){
          m1 <- metagen(dat$TE1 - dat$TE2,
                        sqrt(dat$seTE1^2 + dat$seTE2^2), ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          if (i==1 & j==2)
            res <- dat
          else
            res <- rbind(res, dat)
        }
        else
          if (i==1 & j==2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  
  
  if (type=="count"){
    if (length(event) != narms)
      stop("Different length of lists 'treat' and 'event'.")
    if (length(time) != narms)
      stop("Different length of lists 'treat' and 'time'.")
    ##
    for (i in 1:(narms-1)){
      ##
      if (i==1 & (length(treat[[i]]) != length(event[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'event'.")
      if (i==1 & (length(treat[[i]]) != length(time[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'time'.")
      ##
      for (j in (i+1):narms){
        ##
        if (length(treat[[j]]) != length(event[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'event'.")
        if (length(treat[[j]]) != length(time[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'time'.")
        ##
        dat <- data.frame(TE=NA, seTE=NA,
                          studlab=studlab,
                          treat1=treat[[i]],
                          treat2=treat[[j]],
                          event1=event[[i]], time1=time[[i]],
                          event2=event[[j]], time2=time[[j]])
        dat <- dat[!(is.na(dat$event1) & is.na(dat$time1)),]
        dat <- dat[!(is.na(dat$event2) & is.na(dat$time2)),]
        ##
        if (nrow(dat) > 0){
          m1 <- metainc(dat$event1, dat$time1,
                        dat$event2, dat$time2, ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          if (i==1 & j==2)
            res <- dat
          else
            res <- rbind(res, dat)
        }
        else
          if (i==1 & j==2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }


  ##
  ## Additional checks
  ##
  ##
  ## a) Duplicate treatments ?
  ##
  sel.treat <- as.character(res$treat1) == as.character(res$treat2)
  ##
  if (any(sel.treat)) {
    stop(paste("Identical treatments for the following studies:\n  ",
                paste(paste("'", studlab[sel.treat], "'", sep = ""),
                      collapse = " - "), sep = ""))
  }
  ##
  ## b) Studies missing ?
  ##
  sel.study <- !(studlab %in% unique(as.character(res$studlab)))
  ##
  if (any(sel.study))
    warning(paste("The following studies are not considered in the analysis\n  ",
                  "(due to single study arm or missing values):\n  ",
                  paste(paste("'", studlab[sel.study], "'", sep = ""),
                        collapse = " - "), sep = ""))
  
  
  
  
  
  attr(res, "sm") <- m1$sm
  attr(res, "method") <- m1$method
  attr(res, "version") <- packageDescription("netmeta")$Version
  
  
  res <- res[order(factor(res$studlab, levels=levs), res$treat1, res$treat2),]
  ##
  rownames(res) <- 1:nrow(res)
  res
}
