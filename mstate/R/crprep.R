crprep <- function(Tstop, ...) UseMethod("crprep")

crprep.default <-
function(Tstop, status, data, trans=1, cens=0, Tstart=0, id, strata, keep, shorten=TRUE, rm.na=TRUE, origin=0,
         prec.factor=1000, ...) {

  ## Extract Tstop data if given by column name
  if (!(is.numeric(Tstop))) {
    if (!is.character(Tstop) | length(Tstop)!=1)
      stop("argument \"Tstop\" should be a numeric vector or a character string")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"Tstop\" argument is a character string")
    tcol <- match(Tstop, names(data))
    if (is.na(tcol))
      stop("\"Tstop\" not found in data")
    Tstop <- data[, tcol]
  } else {
    if (!is.vector(Tstop))
      stop("argument should be a numeric vector or a character string")
  }
  nn <- length(Tstop)

  ## Extract Tstart data if given by column name
  if (is.numeric(Tstart)&length(Tstart) == 1) {
    Tstart <- rep(Tstart, nn)
  } else {
    if (!(is.numeric(Tstart))) {
      if (!is.character(Tstart) | length(Tstart)!=1)
        stop("argument \"Tstart\" should be a numeric vector or a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"Tstart\" argument is a character string")
      tcol <- match(Tstart, names(data))
      if (is.na(tcol))
        stop("\"Tstart\" not found in data")
      Tstart <- data[, tcol]
      } else {
      if (!is.vector(Tstart))
        stop("argument should be a numeric vector or a character string")
      }
  }
  if (length(Tstart) != nn)
    stop("Tstop and Tstart have different lengths")
  ## Check whether Tstart is needed
  calc.trunc <- any(Tstart[!is.na(Tstart)] != 0)

  ## Select rows without missing time value
  sel <- !is.na(Tstart) & !is.na(Tstop)
  if (any(Tstart[sel] >= Tstop[sel]))
    stop("Tstop must be greater than Tstart")

  ## Extract status data if given by column name
  if (length(status) == 1) {
    if (!is.character(status))
      stop("argument \"status\" should be a vector or a character string")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"status\" argument is a character string")
    tcol <- match(status, names(data))
    if (is.na(tcol))
      stop("\"status\" not found in data")
    status <- data[ ,tcol]
  }
  if (length(status) != nn)
    stop("Tstop and status have different lengths")

  ## Extract strata data; value 1 if not specified
  if (missing(strata)) {
    strata.val <- rep(1,nn)
  } else {
    if (is.matrix(strata) | is.data.frame(strata))
        stop("only one variable is allowed in \"strata\"")
    if (!(is.vector(as.numeric(factor(strata))) & length(strata) > 1)) {
      if (!is.character(strata))
        stop("argument \"strata\" should be a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"strata\" argument is a character string")
      tcol <- match(strata, names(data))
      if (is.na(tcol))
        stop("\"strata\" not found in data")
      strata.name <- strata
      strata.val <- data[ ,tcol]
    } else {
      if (length(strata) != nn)
        stop("Tstop and strata have different lengths")
      strata.name <- names(strata)
      strata.val <- strata
    }
  }
  strata.num <- as.numeric(factor(strata.val))

  ## Extract id data; values 1:nn if not specified
  if (missing(id)) {
    id.name <- "id"
    id  <- num.id <- 1:nn
  } else {
    if (is.matrix(id) | is.data.frame(id))
      stop("only one variable is allowed in \"id\"")
    if (!(is.vector(id) & length(id) > 1)) { # by name
      if (!is.character(id))
        stop("argument \"id\" should be a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"id\" argument is a character string")
      tcol <- match(id, names(data))
      if (is.na(tcol))
        stop("\"id\" not found in data")
      id.name <- id
      num.id <- 1:nn
      id <- data[, tcol]
    } else {                                 # by value
      if (length(id) != nn)
        stop("Tstop and id have different lengths")
      id.name <- names(id)
      num.id <- 1:nn
    }
  }

  ## Eliminate records with missings in status if rm.na=TRUE
  if(rm.na) sel <- sel & !is.na(status)
  Tstart <- Tstart[sel]
  Tstop <- Tstop[sel]
  status <- status[sel]
  strata.val <- strata.val[sel]
  strata.num <- strata.num[sel]
  id <- id[sel]
  num.id <- num.id[sel]
  n <- length(Tstop)

  ## Extract covariate data
  if(!missing(keep)) {
    if (!(is.matrix(keep) | is.data.frame(keep))) {
      if (is.character(keep)) {  # if given by column name
        if (missing(data))
          stop("missing \"data\" argument not allowed when \"keep\" argument is a character vector")
        nkeep <- length(keep)
        kcols <- match(keep, names(data))
        if (any(is.na(kcols)))
          stop("at least one element of \"keep\" not found in data")
        keep.name <- keep
        keep <- data[, kcols]
      } else {                     # if one column, given by value
        nkeep <- 1
##        keep.name <- names(as.data.frame(keep))
        keep.name <- names(keep)
##        if(is.null(keep.name)) keep.name <- "V1"
        if (length(keep) != nn)
          stop("Tstop and keep have different lengths")
      }
    } else {                         # if a matrix/data.frame
      nkeep <- ncol(keep)
      if(is.data.frame(keep)){
        keep.name <- names(keep)
      } else {
        keep.name <- colnames(keep)
        if(is.null(keep.name)) keep.name <- paste("V",1:nkeep,sep="")
      }
      if (nrow(keep) != nn)
        stop("length Tstop and number of rows in keep are differents")
      if (nkeep == 1)
        keep <- keep[, 1]
    }
  }

  Tstart <- Tstart - origin
  Tstop <- Tstop - origin


  ## Start calculations
  prec <- .Machine$double.eps*prec.factor

  ## Calculate product-limit time-to-censoring distribution, "event" not included in case of ties
  surv.cens <- survival::survfit(Surv(Tstart,Tstop+ifelse(status==cens,prec,0),status==cens)~strata.num)

  ## Calculate time to entry (left truncation) distribution at t-, use 2*prec in order to exclude censorings at same time
  if(calc.trunc) surv.trunc <- survival::survfit(Surv(-Tstop,-(Tstart+2*prec),rep(1,n))~strata.num)
  ## trunc.dist <- summary(surv.trunc)
  ## trunc.dist$time <- rev(-trunc.dist$time)-prec
  ## trunc.dist$surv <- c(rev(trunc.dist$surv)[-1],1)

  ## Create weighted data set for each event type as specified in trans
  data.out <- vector("list",length(trans))
  i.list <- 1
  strat <- sort(unique(strata.num),na.last=TRUE)
  len.strat <- length(strat)
  ## Start weight calculation per event type
  for(failcode in trans) {
    if(len.strat==1){ # no strata
      data.weight <- create.wData.omega(Tstart, Tstop, status, num.id, 1, failcode, cens)
      tmp.time <- data.weight$Tstop
      data.weight$weight.cens[order(tmp.time)] <- summary(surv.cens, times=tmp.time-prec)$surv
      if(calc.trunc) data.weight$weight.trunc[order(-tmp.time)] <- summary(surv.trunc, times=-tmp.time)$surv
    } else {
      data.weight <- vector("list",len.strat)
      if(is.na(strat[len.strat])) {
        tmp.sel <- is.na(strata.num)
        data.weight[[len.strat]] <- data.frame(id=num.id[tmp.sel], Tstart=Tstart[tmp.sel], Tstop=Tstop[tmp.sel], status=status[tmp.sel], strata=NA,  weight.cens=NA)
        if(calc.trunc) data.weight[[len.strat]]$weight.trunc <- NA
      }
      for(tmp.strat in 1:(len.strat-is.na(strat[len.strat]))){
        tmp.sel <- !is.na(strata.num) & strata.num==tmp.strat
        data.weight[[tmp.strat]] <- create.wData.omega(Tstart[tmp.sel], Tstop[tmp.sel], status[tmp.sel], num.id[tmp.sel], tmp.strat, failcode, cens)
        tmp.time <- data.weight[[tmp.strat]]$Tstop
        data.weight[[tmp.strat]]$weight.cens[order(tmp.time)] <- summary(surv.cens[tmp.strat], times=tmp.time-prec)$surv
        if(calc.trunc) data.weight[[tmp.strat]]$weight.trunc[order(-tmp.time)] <- summary(surv.trunc[tmp.strat], times=-tmp.time)$surv
      }
      data.weight <- do.call("rbind", data.weight)
    }
    ## Calculate omega-censoring weights
    data.weight <- data.weight[order(data.weight$id,data.weight$Tstop), ]
    data.weight$weight.cens <- unlist(tapply(data.weight$weight.cens, data.weight$id, FUN=function(x) if(length(x)==1&!is.na(x[1])) 1 else x/x[1]))

    ## Calculate omega-truncation weights
    if(calc.trunc) {
      data.weight$weight.trunc <- unlist(tapply(data.weight$weight.trunc, data.weight$id, FUN=function(x) if(length(x)==1&!is.na(x[1])) 1 else x/x[1]))
    }

    tbl <- table(data.weight$id)

    ## Add covariates
    if(!missing(keep)) {
      ## Extract covariate name from function
      if (is.null(keep.name)) {
        m <- match.call(expand.dots = FALSE)
        m <- m[match("keep", names(m))]
        if(!is.null(m)) {
          keep.name <- as.character(m[1])
          keep.name.split <- strsplit(keep.name, '')[[1]]
          tag <- which(keep.name.split == '$')
          if(length(tag) != 0) {
            keep.name <- substring(keep.name, tag[length(tag)]+1)
          } else {
            tag <- which(keep.name.split == '"')
            if(length(tag) != 0) {
              keep.name <- substring(keep.name, tag[1]+1, tag[2]-1)
            }
          }
        }
      }

      ## Add covariates to the resultset
      if (nkeep > 0) {
        if (nkeep == 1) {
          keep <- keep[sel]
          ddcovs <- rep(keep, tbl)
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- as.character(keep.name)
        } else {
          keep <- keep[sel, ]
          ddcovs <- lapply(1:nkeep, function(i) rep(keep[, i], tbl))
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- keep.name
        }
        data.weight <- cbind(data.weight, ddcovs)
      }
    }

    ## Shorten data set by combining rows with event types without censoring or truncation time in between
    if (shorten) {
      if(calc.trunc) {
        keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0 | diff(weight.trunc)!=0, TRUE))
      } else {
        keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0, TRUE))
      }
      ## First record always included as separate row in order to allow for CSH analysis
      keep.rows[unlist(mapply(seq,1,tbl))==1] <- TRUE
      keep.start <- data.weight$Tstart[unlist(tapply(keep.rows, data.weight$id, FUN=function(x) if(length(x)==1) x else c(TRUE,x[-length(x)])))]
      data.weight <- data.weight[keep.rows,]
      data.weight$Tstart <- keep.start
    }

    ## Recalculate tbl after shorten
    tbl <- table(data.weight$id)
    ## Add count
    data.weight$count <- unlist(mapply(seq,1,tbl))
    data.weight$failcode <- failcode
    ## Return to original id
    data.weight$id <- rep(id,tbl)

    data.out[[i.list]] <- data.weight
    i.list <- i.list+1
  }

  out <- do.call("rbind", data.out)

  if(!missing(strata)) {
  ## Extract strata name if given by value
    if (is.null(strata.name)) {
      m <- match.call(expand.dots = FALSE)
      m <- m[match("strata", names(m))]
      if(!is.null(m)) {
        strata.name <- as.character(m[1])
        strata.name.split <- strsplit(strata.name, '')[[1]]
        tag <- which(strata.name.split == '$')
        if(length(tag) != 0) {
          strata.name <- substring(strata.name, tag[length(tag)]+1)
        } else {
          tag <- which(strata.name.split == '"')
          if(length(tag) != 0) {
            strata.name <- substring(strata.name, tag[1]+1, tag[2]-1)
          }
        }
      }
    }
    ## Use original stratum values
    if(is.factor(strata.val)) {
      out$strata <- factor(out$strata, labels=levels(strata.val))
    } else {
      out$strata <- levels(factor(strata.val))[out$strata]
      if(is.numeric(strata.val)) out$strata <- as.numeric(out$strata)
    }
    ## Use original name of column
    tmp.sel <- match("strata", names(out))
    names(out)[tmp.sel] <- strata.name
  } else {
    out$strata <- NULL
  }

  if (is.null(id.name)) {
    m <- match.call(expand.dots = FALSE)
    m <- m[match("id", names(m))]
    if(!is.null(m)) {
      id.name <- as.character(m[1])
      id.name.split <- strsplit(id.name, '')[[1]]
      tag <- which(id.name.split == '$')
      if(length(tag) != 0) {
        id.name <- substring(id.name, tag[length(tag)]+1)
      } else {
        tag <- which(id.name.split == '"')
        if(length(tag) != 0) {
          id.name <- substring(id.name, tag[1]+1, tag[2]-1)
        }
      }
    }
  }

  row.names(out) <- as.character(1:nrow(out))
  names(out)[1] <- id.name
  class(out) <- c("crprep","data.frame")
  return(out)
}

