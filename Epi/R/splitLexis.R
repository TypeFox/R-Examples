split.lexis.1D <-
function(lex, breaks, time.scale, tol)
{
    time.scale <- check.time.scale(lex, time.scale)

    ## Entry and exit times on the time scale that we are splitting
    time1 <- lex[,time.scale, drop=FALSE]
    time2 <- time1 + lex$lex.dur

    ## Augment break points with +/- infinity
    breaks <- sort( unique( breaks ) )
    I1 <- c(-Inf, breaks)
    I2 <- c(breaks,Inf)

    ## Arrays containing data on each interval (rows) for each subject (cols)
    en <- apply(time1, 1, pmax, I1)     # Entry time
    ex <- apply(time2, 1, pmin, I2)     # Exit time
    NR <- nrow(en)
    NC <- ncol(en)

    ## Does subject contribute follow-up time to this interval?
    ## (intervals shorter than tol are ignored)
    valid <- en < ex - tol
    dur <- ex - en; dur[!valid] <- 0    # Time spent in interval

    ## Cumulative time since entry at the start of each interval
    time.since.entry <- rbind(0, apply(dur,2,cumsum)[-NR,,drop=FALSE])

    cal.new.entry <- function(entry.time) {
        sweep(time.since.entry, 2, entry.time, "+")[valid]
    }

    old.entry <- lex[, timeScales(lex), drop=FALSE]
    new.entry <- lapply(old.entry, cal.new.entry)

    ## Status calculation
    aug.valid <- rbind(valid, rep(FALSE, NC))
    last.valid <- valid & !aug.valid[-1,]
    any.valid <- apply(valid,2,any)

    new.Xst <- matrix( lex$lex.Cst, NR, NC, byrow=TRUE)
    new.Xst[last.valid] <- lex$lex.Xst[any.valid]

    n.interval <- apply(valid, 2, sum)
    new.lex <- Lexis("entry" = new.entry,
                     "duration" = dur[valid],
                     "id" = rep(lex$lex.id, n.interval),
                     "entry.status" = rep(lex$lex.Cst, n.interval),
                     "exit.status" = new.Xst[valid])

    ## Update breaks attribute and tranfer time.since attribute
    breaks.attr <- attr(lex, "breaks")
    breaks.attr[[time.scale]] <- sort(c(breaks.attr[[time.scale]], breaks))
    attr(new.lex, "breaks") <- breaks.attr
    attr(new.lex, "time.since") <- attr(lex, "time.since")
    return(new.lex)
}


splitLexis <- function(lex, breaks, time.scale=1, tol= .Machine$double.eps^0.5)
{
  ## Advise the uninformed user...
  if( inherits(lex,"stacked.Lexis") )
    stop( "It makes no sense to time-split after stacking ---\n",
    "split your original Lexis object and stack that to get what you want.\n")


  ## Set temporary, unique, id variable
  lex$lex.tempid <- lex$lex.id
  lex$lex.id <- 1:nrow(lex)

  ## Save auxiliary data
  aux.data.names <- setdiff(names(lex), timeScales(lex))
  aux.data.names <- aux.data.names[substr(aux.data.names,1,4) != "lex."]
  aux.data <- lex[, c("lex.id","lex.tempid", aux.data.names), drop=FALSE]

  ## Check for NAs in the timescale
  ts <- check.time.scale(lex, time.scale)
  ts.miss <- any(is.na(lex[,ts]))
  if( ts.miss )
    {
    na.lex <- lex[ is.na(lex[,ts]),]
       lex <- lex[!is.na(lex[,ts]),]
    cat( "Note: NAs in the time-scale \"", ts, "\", you split on\n")
    }

  ## If states are factors convert to numeric while splitting
  factor.states <- is.factor( lex$lex.Cst )
  if( factor.states )
    {
    state.levels <- levels( lex$lex.Cst )
    nstates     <- nlevels( lex$lex.Cst )
    lex$lex.Cst <- as.integer( lex$lex.Cst )
    lex$lex.Xst <- as.integer( lex$lex.Xst )
    }

  ## Split the data
  lex <- split.lexis.1D(lex, breaks, time.scale, tol)

  ## Reinstitute the factor levels
  if( factor.states )
    {
    lex$lex.Cst <- factor( lex$lex.Cst, levels=1:nstates, labels=state.levels )
    lex$lex.Xst <- factor( lex$lex.Xst, levels=1:nstates, labels=state.levels )
    }

  ## Put the NA-rows back
  if( ts.miss ) lex <- rbind( lex, na.lex[,colnames(lex)] )

  ## Save attributes
  lex.attr <- attributes(lex)
  ## Merge
  lex <- merge.data.frame(lex, aux.data, by="lex.id")
  ## Restore attributes
  attr(lex,"breaks") <- lex.attr$breaks
  attr(lex,"time.scales") <- lex.attr$time.scales
  attr(lex,"time.since") <- lex.attr$time.since
  class(lex) <- c("Lexis", "data.frame")
  ## Restore id variable
  lex$lex.id <- lex$lex.tempid
  lex$lex.tempid <- NULL

  return(lex)
}
