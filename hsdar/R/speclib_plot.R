setMethod("plot", signature(x = "Speclib"), 
          definition = function(x,
                                FUN = NULL,
                                new = TRUE,
#                                 byattributes = NULL,
#                                 legend = NULL,
                                ...)
{
  
  if (dim(x)[1] == 1) FUN <- 1
  
  if (!is.null(attr(x, "setmask")))
  {
    if (attr(x, "setmask"))
    {
      dropped <- attr(x, "dropped")
      if (nspectra(x) == 1)
      {
        duplicate_first <- TRUE
        spectra(x) <- rbind(spectra(x), spectra(x))
      } else {
        duplicate_first <- FALSE
      }
      for (i in 1:nrow(dropped))
      {
        wav_back <- wavelength(x)
        wavelength(x) <- c(wavelength(x)[wav_back<dropped[i,1]], 
                           (dropped[i,1]+dropped[i,2])/2,
                           wavelength(x)[wav_back>dropped[i,2]])
        spectra(x) <- cbind(spectra(x)[,wav_back<dropped[i,1]],
                            matrix(data = NaN, nrow = dim(x)[1], ncol = 1),
                            spectra(x)[,wav_back>dropped[i,2]])        
      }
      if (duplicate_first)
        spectra(x) <- matrix(spectra(x)[1,], nrow = 1)
    }
  }
  
  ## Defaults
  lty <- par()$lty
  col <- par()$col
  lwd <- par()$lwd
  type <- "l"
  lend <- par()$lend
  ljoin <- par()$ljoin 
  lmitre <- par()$lmitre 
  pch <- par()$pch
  
  if (is.null(FUN))
  {
    mean_spec <- apply(x, FUN = mean, na.rm = TRUE)
    sd_spec   <- apply(x, FUN = sd, na.rm = TRUE)
    spectra2plot <- rbind(spectra(mean_spec) + spectra(sd_spec), 
                          spectra(mean_spec),                           
                          spectra(mean_spec) - spectra(sd_spec))    
    lty <- c("dashed", "solid", "dashed")
  } else {
    if (is.numeric(FUN))
    {
      spectra2plot <- matrix(spectra(x)[FUN,], ncol = length(wavelength(x)))
    } else {
      spectra2plot <- spectra(apply(x, FUN = FUN))
    }
  }
  
  call_fu <- list(...)
  
  if (any(names(call_fu) == "lty")) 
    lty <- call_fu$lty  
  if (any(names(call_fu) == "col")) 
    col <- call_fu$col
  if (any(names(call_fu) == "lwd")) 
    lwd <- call_fu$lwd
  if (any(names(call_fu) == "type")) 
    type <- call_fu$type  
  if (any(names(call_fu) == "pch")) 
    pch <- call_fu$pch
  if (any(names(call_fu) == "lend")) 
    lend <- call_fu$lend  
  if (any(names(call_fu) == "ljoin")) 
    ljoin <- call_fu$ljoin  
  if (any(names(call_fu) == "lmitre")) 
    lmitre <- call_fu$lmitre
    
  nlines <- nrow(spectra2plot) 
  
  lty    <- .adopt_dotsvar(lty, nlines)
  col    <- .adopt_dotsvar(col, nlines)
  lwd    <- .adopt_dotsvar(lwd, nlines)
  type   <- .adopt_dotsvar(type, nlines)
  pch    <- .adopt_dotsvar(pch, nlines)
  lend   <- .adopt_dotsvar(lend, nlines)
  ljoin  <- .adopt_dotsvar(ljoin, nlines)
  lmitre <- .adopt_dotsvar(lmitre, nlines)
  
  if (new)
  {      
    if (any(names(call_fu)=="xlim")) xlim <- call_fu$xlim
    else xlim <- range(wavelength(x), na.rm = TRUE)
   
    if (any(names(call_fu)=="ylim")) ylim <- call_fu$ylim
    else ylim <- range(spectra2plot, na.rm = TRUE)

    if (any(names(call_fu)=="xlab")) xlab <- call_fu$xlab
    else xlab <- paste(x@xlabel," (", x@wlunit,")", sep = "")
    
    if (any(names(call_fu)=="ylab")) ylab <- call_fu$ylab
    else ylab <- x@ylabel 
      
    if (any(names(call_fu)=="main")) main <- call_fu$main
    else main <- ""
      
    if (any(names(call_fu)=="xaxt")) xaxt <- call_fu$xaxt
    else xaxt <- "s" 
    
    if (any(names(call_fu)=="yaxt")) yaxt <- call_fu$yaxt
    else yaxt <- "s" 
    
      
    plot(xlim, ylim, type = "n", xlab = "", ylab = "",
         xaxt = xaxt, yaxt = yaxt)
    title(xlab = xlab, ylab = ylab, main = main)
  }
  
  status <- sapply(c(1:nrow(spectra2plot)), .plot_spec_curves,
                   wavelength(x), spectra2plot, lty, col, lwd, 
                   type, pch, lend, ljoin, lmitre)
}
)

.plot_spec_curves <- function(i, x, y, lty, col, lwd, type, pch, 
                             lend, ljoin, lmitre)
{
  if (type[i] == "l")
  {
    lines(x, y[i,], lty = lty[i], col = col[i], lwd = lwd[i],
          lend = lend[i], ljoin = ljoin[i], lmitre = lmitre[i])  
  } else {
    if (type[i] == "b")
      lines(x, y[i,], lty = lty[i], col = col[i], lwd = lwd[i], 
            type = "b", pch = pch[i], lend = lend[i], 
            ljoin = ljoin[i], lmitre = lmitre[i])  
  }
}

.adopt_dotsvar <- function(x, n)
{
  if (length(x) != n)
    x <- rep.int(unlist(x)[1], n)
  return(x)
}



legendSpeclib <- function(x, speclib, ...)
{
  legend.args <- c("y", "legend", "fill", "border", "lty", "lwd", "pch", "angle", "density", 
    "bty", "bg", "box.lwd", "box.lty", "box.col", "pt.bg", "cex", "pt.cex", "pt.lwd", "xjust", 
    "yjust", "x.intersp", "y.intersp", "adj", "text.width", "text.col", "merge", "trace", 
    "plot", "ncol", "horiz", "title", "inset", "xpd", "title.col", "title.adj","seg.len")
  call_fu <- match.call()
  if (!any(names(x)=="x")) stop("Position of legend must be specified via variable 'x'")
  if (!any(names(x)=="legend")) 
  {
    if (is.speclib(speclib))
    {
      x$legend <- idSpeclib(speclib)
    }
  }
  if (!any(names(x)=="col")) 
    x$col <- if (any(names(call_fu)=="col")) eval(parse(text = call_fu[which(names(call_fu)=="col")])) else "black"
  if (!any(names(x)=="lty")) 
    x$lty <- "solid"
  if (length(names(call_fu))>2)
  {
    for (i in 3:length(names(call_fu)))
    {
      if (any(legend.args==names(call_fu)[i]))
      {
        x[[i-1]] <- eval(parse(text = call_fu[i]))
        names(x)[i-1] <- names(call_fu)[i]
      }
    }
  }
  if (is.null(x$legend)) stop("argument 'legend' is missing")
  do.call("legend",x)
}