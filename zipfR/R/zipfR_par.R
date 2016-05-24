zipfR.par <- function (..., bw.mode=FALSE)
{
  args <- list(...)
  known.pars <- ls(.PAR) 

  if (! length(args)) {                 # --> zipfR.par()
    return(as.list(.PAR))
  }
  else if (is.null(names(args))) {      

    if (length(args) == 1 && is.list(args[[1]])) { # --> zipfR.par(par.save)
      args <- args[[1]]                 # then fall through to default case below

    } else {                            # --> zipfR.par("lty", "lwd")
      pars <- as.character(unlist(args))
      if (bw.mode) {                    # bw.mode=TRUE -> return corresponding B/W mode parameters if possible
        pars.bw <- paste(pars, "bw", sep=".")
        idx <- pars.bw %in% known.pars
        if (any(idx)) pars[idx] <- pars.bw[idx]
      }
      return( if (length(pars) == 1) .PAR[[ pars[1] ]] else as.list(.PAR)[pars] )
    }
  }               

  invalid <- !(names(args) %in% known.pars) # --> zipfR.par(lwd.bw=2, bw=TRUE, ...)
  if (any(invalid)) {
    warning("invalid zipfR graphics parameter(s) ", paste(names(args)[invalid], collapse=", "), " ignored")
    args <- args[! invalid]
  }
  old <- as.list(.PAR)[names(args)] # make sure we get a list (can't call zipfR.arg(names))
  for (key in names(args)) { .PAR[[ key ]] <- check.par(key, args[[ key ]]) }

  return(invisible(old))
}

## private environment for graphics parameter data
## (I wonder whether this is the way it's supposed to work, but at least it works ...)
.PAR <- new.env()

## styles for colour plots

.PAR$lty <- rep( c("solid", "solid", "33", "solid", "1232"), 2)
.PAR$lwd <- rep( c(3,3,3,3,3), 2)
.PAR$col <- rep( c("black", "#FF0000", "#00AA00", "#4444FF", "black"), 2)
.PAR$pch <- rep( c(1, 3, 15, 2, 20), 2)
.PAR$barcol <- rep( c("black", "red", "blue", "yellow", "green"), 2 )

## styles for b/w plots 
.PAR$lty.bw <- rep( c("solid", "dashed", "12", "solid", "dashed"), 2)
.PAR$lwd.bw <- rep( c(2,2,3,1,2), 2)
.PAR$col.bw <- rep( c("grey30", "black", "black", "black", "grey30"), 2)
.PAR$pch.bw <- rep( c(1, 3, 15, 2, 20), 2)
.PAR$barcol.bw <- rep( c("black", "grey50", "white", "grey70", "grey20"), 2 )

## whether to produce b/w graphics by default
.PAR$bw <- FALSE

## for the zipfR plotutils functions
.PAR$device <- if (capabilities()["aqua"]) "quartz" else "x11"
.PAR$init.par <- list()
.PAR$width <- 6
.PAR$height <- 6
.PAR$bg <- "white"
.PAR$pointsize <- 12

## interal helper function to check validity of graphics parameters
check.par <- function (name, value) {
  if (name %in% c("bw")) {              # Boolean parameters
    value <- as.logical(value)
    if (is.na(value) || length(value) != 1) stop("parameter '",name,"' must be a Boolean value (logical)")
  }
  else if (name == "device") {
    supported <- c("x11", "png", "eps", "pdf")
    if (.Platform$GUI == "AQUA") supported <- c("quartz", supported)
    value <- match.arg(value, supported)
  }
  else if (name == "init.par") {
    if (is.null(value)) value <- list() # init.par=NULL translates to empty list
    if (! is.list(value)) stop("parameter 'init.par' must be a list of name=value pairs")
  }
  else if (name %in% c("lwd", "lwd.bw")) {
    if (! is.numeric(value)) stop("parameter '",name,"' must be a numeric vector")
    if (length(value) > 10) warning("extra style options for parameter '",name,"' ignored")
    value <- rep(value, length.out=10)
  }
  else if (name %in% c("lty", "col", "lty.bw", "col.bw")) {
    if (length(value) > 10) warning("extra style options for parameter '",name,"' ignored")
    value <- rep(value, length.out=10)
  }

  value
}
