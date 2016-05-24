
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



regression <- function(x, y, main, scr,
                       value=function(regr) regr$coeff[2], variable="var",
                       mode=c("nographics", "plot", "interactive"),
                       regr.select=c("coefficients"),
                       averaging=FALSE, cex.main=0.85,
                       col.passive = "grey",
                       col.active = "green",
                       col.main.passive = "black",
                       col.main.active = "red",
                       col.points.chosen = "lightblue",
                       col.abline = "yellow",
                       col.abline.chosen = "darkblue",
                       col.smooth = "red",
                       ...) {
  
  bg <- par()$bg
  par(bg="white")
  mode <- match.arg(mode)
  if (is.array(y)) {
    stopifnot(length(x)==dim(y)[1])
    y <- as.vector(y)
  } else stopifnot(length(y) %% length(x) == 0)
  if (length(y)>length(x)) {
    if (averaging) y <- rowMeans(matrix(y, nrow=length(x)))
    else x <- rep(x, len=length(y))
  }
  regr <- lm(y ~ x)[regr.select]
  val <- value(regr)
  names(val) <- NULL
  x.u <- y.u <- regr.u <- val.u <- NULL
  ## curmain must be set since within repeat loop curmain is not set if
  ## the user immediately breaks the loop
  curmain <- main <- paste(main,": ",variable,"=", format(val, dig=3), sep="")
  
  sm <- ksmooth(x, y, n.points=100, kernel="box", bandwidth=0.5)[c("x", "y")]
  if (mode!="nographics") { 
    screen(scr)
    par(mar=c(4.1, 4.1, 3, 0.1))
    plot(x, y,
         col=if (mode=="plot") col.passive else col.active, main=main,
         col.main=if (mode=="plot")  col.main.passive else col.main.active,
         cex.main=cex.main, ...)
    abline(regr, col=col.abline)
    lines(sm$x, sm$y, col=col.smooth)
    if (mode=="interactive") {
      repeat {
        loc <- try(locator(2), silent=TRUE)
        if (is(loc, "try-error")) loc <- NULL
        if (length(loc)==0) break;
        r <- range(loc$x)
        x.u <- x[idx <- x >= r[1] & x <= r[2]]
        y.u <- y[idx]
        if (diff(r)==0 || length(x.u) <= 1) next;
        regr.u <- lm(y.u ~ x.u)[regr.select]
        val.u <- value(regr.u)
        names(val.u) <- NULL
        curmain <- paste(main, ", ",variable,".u=", format(val.u, dig=3), sep="")
        screen(scr)
        plot(x, y, col=col.active, main=curmain, col.main=col.main.active,
             cex.main=cex.main, ...)
        points(x.u, y.u, col=col.points.chosen, ...)
        abline(regr, col=col.abline)
        lines(sm$x, sm$y, col=col.smooth)
        abline(regr.u, col=col.abline.chosen)
      }
      screen(scr)
      plot(x, y, col=col.passive, main=curmain, col.main=col.main.passive,
           cex.main=cex.main, ...)
      if (!is.null(x.u)) points(x.u, y.u, col=col.points.chosen, ...)
      abline(regr, col=col.abline)
      lines(sm$x, sm$y, col=col.smooth)
      if (!is.null(x.u)) abline(regr.u, col=col.abline.chosen)
    }
  }
  par(bg=bg)
  return(list(regr=regr, val=val,
         regr.u=regr.u, val.u=val.u, x.u=x.u, y.u=y.u,
         sm=sm))
}

  
RFhurst <- function(x, y = NULL, z = NULL, data, sort=TRUE,
                   block.sequ =
                    unique(round(exp(seq(log(min(3000, dimen[1] / 5)),
                     log(dimen[1]), len=min(100, dimen[1]))))),
                   fft.m = c(1, min(1000, (fft.len - 1) / 10)),
                   fft.max.length = Inf, ## longer ts are cut down
                   method=c("dfa", "fft", "var"),
                   mode = if (interactive()) c("plot", "interactive")
                   else "nographics",
                   pch=16, cex=0.2, cex.main=0.85,
                   printlevel=RFoptions()$general$printlevel,
                   height=3.5,
                   ...
                   ) {
  l.method <- eval(formals()$method)
  pch <- rep(pch, len=length(l.method))
  cex <- rep(cex, len=length(l.method))
  T <- NULL # if (length(T)!=0) stop("time not programmed yet")
  grid <- TRUE # if (!grid) stop("only grids are possible")
  
  method <- l.method[pmatch(method, l.method)]
  do.dfa <- any(method=="dfa")
  do.fft <- any(method=="fft")
  do.var <- any(method=="var")

  l.block.sequ <- l.dfa.var <- dfa <-
    l.varmeth.sequ <- l.varmeth.var <- varmeth <- 
      l.lambda <- l.I.lambda <- fft <- NULL
  
  modes <- c("nographics", "plot", "interactive")
  if (any(is.na(mode <- modes[pmatch(mode, modes)])))
    stop("unknown values of `mode'")
  

  ct <- CheckXT(x=x, y=y, z=z, T=T, grid=TRUE)
  stopifnot(ct$grid)
  dimen <- cbind(ct$x, ct$T)[3, ] 

  if (block.sequ[1] < 2500)
    warning("results may show high variation due to short sequence(s).")

  if (printlevel>=PL_STRUCTURE) cat("(formatting) ")
  if (ncol(ct$x)>1) {
    if (!is.array(data)) data <- array(as.double(data), dim=dimen)
    if (sort) { ## so that the first dimension gives the side with
      ##           the greatest length (in points)
      ord <- order(dimen, decreasing=TRUE)
      if (any(diff(ord)<0)) {
        data <- aperm(data, ord)
        ct$x <- ct$x[, ord, drop=FALSE]
      }
    }
    gc()
  } else {
    data <- as.matrix(data)
    dimen <- c(dimen, 1)
  }
  repet <- prod(dimen[-1])

  ## periodogram, code taken from jean-francois coeurjolly
  if (do.fft) {


     
    if (printlevel>=PL_STRUCTURE) cat("periodogramm")
    fft.len <- min(dimen[1], fft.max.length)
    fft.m <- round(fft.m)
    stopifnot(diff(fft.m)>0, all(fft.m>0), all(fft.m<=fft.len))
#    Print("periodogramm")
    l.I.lambda <-
      .Call("periodogram",
            data,
            as.integer(dimen[1]), 
            as.integer(repet),# Produkt der anderen Dimensionen
            as.integer(fft.m),# Ausschnitt aus Fourier-Trafo aus Stueck
            ##                  nachf. Laenge
            as.integer(fft.len),# Reihe zerhackt in Stuecke dieser Laenge 
            as.integer(fft.len / 2), ## WOSO(?)-Sch\"aetzer
            PACKAGE="RandomFields")
    l.lambda <-  log((2 * pi * (fft.m[1]:fft.m[2])) / fft.len)
  } 


  ## detrended fluctuation analysis
  if (do.dfa || do.var) {
    if (printlevel>=PL_STRUCTURE) {
      if (do.dfa) cat("detrended fluctuation; ")
      if (do.var) cat("aggregated variation; ")
    }
    
    stopifnot(all(diff(block.sequ)>0))
    l.block.sequ <- log(block.sequ)
    dfa.len <- length(block.sequ)
    ##    l.dfa.var <- double(dfa.len * repet)
    ##    l.varmeth.var <- double(dfa.len * repet)
    ## wird data zerstoert ?!
    ## log already returned!
#Print("detrendedfluc")
    l.var <- ## rbind(l.varmeth.var, l.dfa.var)
      .Call("detrendedfluc", as.double(data), as.integer(dimen[1]),
            as.integer(repet),
            as.integer(block.sequ), as.integer(dfa.len),
            PACKAGE="RandomFields")
    ## 1:dfa.len since data could be a matrix; l.block.sequ has length dfa.len
    ## and is then shorter than l.varmeth.var!
    l.dfa.var <- l.var[2, ]
    varmeth.idx <- is.finite(l.var[1, 1:dfa.len])
    l.varmeth.var <- l.var[1, varmeth.idx]
    l.varmeth.sequ <- l.block.sequ[varmeth.idx]
  }

  gc() 
  if (printlevel>=PL_STRUCTURE) cat("\n")

  if (any(mode=="plot" | mode=="interactive"))
  {
    cat("\nuse left mouse for locator and right mouse to exit\n")
    plots <- do.dfa + do.fft + do.var
    ScreenDevice(height=height, width=height * plots)
    par(bg="white")
    screens <- seq(0, 1, len=plots+1)
    screens <- split.screen(figs=cbind(screens[-plots-1], screens[-1], 0, 1))
    on.exit(close.screen(screens))  ## nervig sonst bei Fehlern
  }

  ## calculation and plot are separated into blocks since calculation time
  ## is sometimes huge. So here, the user might wait for any result quite
  ## a long time, but is than bothered only once.
  for (m in mode) {
    scr <- 0
    if (do.dfa) {
      scr <- scr + 1
      dfa <- regression(l.block.sequ, l.dfa.var, variable="H", pch=pch[1],
                        main="detrended fluctuation analysis", scr=scr,
                        cex=cex[1], value=function(regr) regr$coeff[2]/2, mode=m,
                        cex.main=cex.main, ...)
    }
    if (do.fft) {
      scr <- scr + 1
      fft <- regression(l.lambda, l.I.lambda, main="periodogram", pch=pch[2],
                        scr=scr,  value=function(regr) 0.5 * (1-regr$coeff[2]),
                        averaging = repet > 1,
                        mode=m, cex=cex[2], variable="H",
                        cex.main=cex.main, ...)
    }
    if (do.var) {
      scr <- scr + 1
      varmeth <- regression(l.varmeth.sequ, l.varmeth.var, variable="H",
                            pch=pch[3],
                            main="aggregated variation", scr=scr, cex=cex[3],
                            value=function(regr) 1 + 0.5 * (regr$coeff[2]-2),
                            mode=m, cex.main=cex.main,...)    
    }
  }
  
  if (printlevel>PL_SUBIMPORTANT ) {
    cat("#################### Hurst #################### \n")
    cat(c(dfa.H=dfa$val, varmeth.H=varmeth$val, fft.H=fft$val))
    cat("---- by interactively defined regression interval:\n")
    cat(c(dfa.H=dfa$val.u, varmeth.H=varmeth$val.u, fft.H=fft$val.u))
    #cat("---- beta (Cauchy model) ----\n")
    #Print(2 * (1 - c(dfa.H=dfa$val, fft.H=fft$val, varmeth.H=varmeth$val)))
    #Print(2 * (1 - c(dfa.H=dfa$val.u, fft.H=fft$val.u,varmeth.H=varmeth$val.u)))
    cat("############################################### \n")
  }
  #if (any(mode=="plot" | mode=="interactive")) close.screen(screens)
  if (any(mode=="plot" | mode=="interactive")) {
    close.screen(screens)
    dev.off() # OK
  }
   
  return(list(dfa=list(x=l.block.sequ, y=l.dfa.var, regr=dfa$regr,
                sm=dfa$sm,
                x.u=dfa$x.u, y.u=dfa$y.u,  regr.u=dfa$regr.u,
                H=dfa$val, H.u=dfa$val.u),
              varmeth=list(x=l.varmeth.sequ, y=l.varmeth.var, regr=varmeth$regr,
                sm=varmeth$sm,
                x.u=varmeth$x.u, y.u=varmeth$y.u, regr.u=varmeth$regr.u,
                H=varmeth$val, H.u=varmeth$val.u),
              fft=list(x=l.lambda, y=l.I.lambda, regr=fft$regr,
                sm=fft$sm,
                x.u=fft$x.u, y.u=fft$y.u, regr.u=fft$regr.u,
                H = fft$val, H.u = fft$val.u)
              )
         )
}


RFfractaldim <-
  function(x, y = NULL, z = NULL, data, grid, 
           bin=NULL,
           vario.n=5,
           sort=TRUE,
           #box.sequ=unique(round(exp(seq(log(1),
           #  log(min(dimen - 1, 50)), len=100)))),
           #box.enlarge.y=1,
           #range.sequ=unique(round(exp(seq(log(1),
           #  log(min(dimen - 1, 50)), len=100)))),
           fft.m = c(65, 86), ## in % of range of l.lambda
           fft.max.length=Inf,
           fft.max.regr=150000,
           fft.shift = 50, # in %; 50:WOSA; 100: no overlapping
           method=c("variogram", "fft"),#"box","range", not correctly implement.
           mode = if (interactive()) c("plot", "interactive") else "nographics",
           pch=16, cex=0.2, cex.main=0.85,
           printlevel = RFoptions()$general$printlevel,
           height=3.5,
           ...) {
  l.method <- eval(formals()$method)
  pch <- rep(pch, len=length(l.method))
  cex <- rep(cex, len=length(l.method))
  T <- NULL # if (length(T)!=0) stop("time not programmed yet")
  
  method <- l.method[pmatch(method, l.method)]
  do.vario <- any(method=="variogram")
  do.box <- any(method=="box")     # physicists box counting method
  do.range <- any(method=="range") # Dubuc et al.
  do.fft <- any(method=="fft")

  ev <- l.midbins <- l.binvario <- vario <-
    Ml.box.sequ <- l.box.count <- box <-
      Ml.range.sequ <- l.range.count <- rnge <-
        l.lambda <- l.I.lambda <- fft <- NULL

  range.sequ <- box.enlarge.y <- box.sequ <- NULL ## not used but otherwise
  ##  check will complain

  modes <- c("nographics", "plot", "interactive")
  if (any(is.na(mode <- modes[pmatch(mode, modes)])))
    stop("unknown values of `mode'")

  if (isSpObj(data)) data <- sp2RF(data)
  if (is(data, "RFsp")) {
    if (!(missing(x) && length(y)==0 && length(z)==0 && length(T)==0))
      stop("x, y, z, T may not be given if 'data' is of class 'RFsp'")
    gridtmp <- isGridded(data)
    compareGridBooleans(grid, gridtmp)
    grid <- gridtmp
    tmp <- RFspDataFrame2conventional(data)
    x <- tmp$x
    y <- NULL
    z <- NULL
    T <- tmp$T
    data <- tmp$data
    rm(tmp)
  }

  ct <- CheckXT(x=x, y=y, z=z, T=T, grid=grid, length.data=length(data))

  ## variogram method (grid or arbitray locations)
  if (do.vario) {
    stopifnot(vario.n > 0)
    if (printlevel>PL_STRUCTURE) cat("variogram; ")
    if (bin.not.given <- length(bin) == 0) {
      if (ct$grid) {
        step <- min(ct$x[2, ])
        end <- ct$x[2, ] * ct$x[3, ] / 4
        bin <- seq(step / 2, end, step)
      } else {        
        step <- apply(ct$x, 2, function(x) list(unique(diff(sort(unique(x))))))
        edge.lengths <- apply(ct$x, 2, function(x) diff(range(x)))
        end <- sqrt(sum(edge.lengths^2)) / 4
        #Print(step)
        if (max(sapply(step, length)) <= 20) {
          if (printlevel>=PL_IMPORTANT) cat("locations on a grid.\n")
          step <- min(unlist(step))
        } else {          
          if (printlevel>PL_IMPORTANT) cat("locations not on a grid.\n")
          dim <- ct$Zeit + ct$spatialdim
          inv.lambda <- prod(edge.lengths) / ct$restotal
          step <- inv.lambda^(1 / dim) # Ann: glm Gitter
          ## Ann: Poisson Punktprozess, Leerwk != p (heuristisch =0.5)
          ##      nach radius R aufgeloest (heuristisch mit 2 multipliziert)
          ## step <- 2 * (invlambda * log(2) * gamma(dim / 2 + 1))^(1 / dim) / sqrt(pi)
        }
      }
      bin <- seq(step / 2, end, step)
    }

 #   Print(bin, step, end, edge.lengths)
     
    ev <- RFempiricalvariogram(x=x, y=y, z=z, T=T, data=data, grid=ct$grid,
                             bin=bin, spConform=FALSE)

    idx <-  which(is.finite(l.binvario <- log(ev$emp)))
    if (length(idx) < vario.n) {
      note <- if (length(idx) <= 1) stop else warning
      if (bin.not.given) 
        note("A good choice for 'bin' was not found. Please try out the arguments manually. Current value is 'bin = ", paste(bin, collapse=", "), "' so that only the bin(s) with center value ", paste(l.binvario[idx], collapse=", "), " possess values of the variogram cloud.")
      else
        note("Current value of 'bin' leads to ", length(idx), " values for the empirical variogram. Consider trying out different arguments.")
    }

    idx <- idx[1:vario.n]
    l.midbins <- log(ev$centers[idx])
    l.binvario <- l.binvario[idx]

  }

  if (ct$grid) {
    dimen <- cbind(ct$x, ct$T)[3, ] 
    
    if (printlevel>=PL_STRUCTURE) cat("(formatting) ")
    if (ncol(ct$x)>1) {
      if (!is.array(data)) data <- array(as.double(data), dim=dimen)
      if (sort) { ## see fractal
        ord <- order(dimen, decreasing=TRUE)
        if (any(diff(ord)<0)) {
          data <- aperm(data, ord)
          ct$x <- ct$x[, ord, drop=FALSE]
        }
      }
    } else {
      data <- as.matrix(data)
      dimen <- c(dimen, 1)
    }
    repet <- prod(dimen[-1])
    gc()
    
    if (do.range) {
      if (printlevel>=PL_STRUCTURE) cat("ranges")
      lrs <- length(range.sequ) ## range.sequ is bound right here!
      #   l.range.count <- double(lrs * repet)
      ## logarithm is already taken within minmax
      storage.mode(data) <- "double"
      l.range.count <-
        .Call("minmax", data, as.integer(dimen[1]),
              as.integer(repet), as.integer(range.sequ), as.integer(lrs),
              PACKAGE="RandomFields")
      box.length.correction <- 0 ## might be set differently
      ##                            for testing or development
      Ml.range.sequ <- -log(range.sequ + box.length.correction)
    }
    
    if (do.box) {
      stop("this point cannot be reached")
      if (printlevel>=PL_STRUCTURE) cat("box counting; ")
      x <- ct$x[,1] / ct$x[3,1] ## alles auf integer
      factor <- diff(x[c(1,2)]) / diff(range(data)) * box.enlarge.y
      ## note: data changes its shape! -- should be irrelevant to any procedure!
      data <- matrix(data, nrow=dimen[1])
     # l.box.count <- double(length(box.sequ) * repet)
      l.box.count <-
        .Call("boxcounting",
              as.double(rbind(data[1,], data, data[nrow(data),])),
              as.integer(dimen[1]),
              as.integer(repet), as.double(factor),
              as.integer(box.sequ),
              PACKAGE="RandomFields")
      gc()
      
      box.length.correction <- 0 ## might be set differently
      ##                            for testing or development
      Ml.box.sequ <- -log(box.sequ + box.length.correction)
      l.box.count <- log(l.box.count)
    }

    if (do.fft) {
      ## periodogram, code taken from jean-francois coeurjolly and WOSA
      if (printlevel>=PL_STRUCTURE) cat("periodogramm; ")
      fft.len <- min(dimen[1],  fft.max.length) ## the virtual length
     
      ## fft.m is given in percent [in log-scale]; it is now rewriten
      ## in absolute indices for the fft vector
      stopifnot(all(fft.m>=0), all(fft.m<=100))
      ## 0.5 * fft.len not fft.len since only first half should be considered
      ## in any case! (Since the fourier transform returns a symmetric vector.)
      spann <- log(2 * pi / c(0.5 * fft.len, 1))
      fft.m <- round(2 * pi / exp(spann[1] + diff(spann) * (1 - fft.m / 100)))
     
      l.lambda <- log((2 * pi * (fft.m[1]:fft.m[2])) / fft.len)
      # l.I.lambda <- double(repet * (fft.m[2] - fft.m[1] + 1));
      l.I.lambda <- .Call("periodogram",
         data,
         as.integer(dimen[1]),
         as.integer(repet),## Produkt der anderen Dimensionen
         as.integer(fft.m),## Ausschnitt aus Fourier-Tr aus Stueck nachf. Laenge
         as.integer(fft.len),## Reihe zerhackt in Stuecke dieser Laenge 
         as.integer(fft.len / 100 * fft.shift), ## shift (WOSA-Sch\"aetzer)
         PACKAGE="RandomFields")
    }
  } else { # not grid
    if (do.box || do.range || do.fft) {
      if (printlevel>=PL_IMPORTANT) {
        cat("\nThe following methods are available only for grids:\n")
        if (do.box)   cat("  * box counting\n")
        if (do.range) cat("  * max-min method\n")
        if (do.fft)   cat("  * periodogram/WOSA\n")
      }
      do.box <- do.range <- do.fft <- FALSE
    }
  }
  if (printlevel>PL_STRUCTURE) cat("\n")


  if (any(mode=="plot" | mode=="interactive")) {
    plots <- do.vario + do.box + do.range + do.fft
    ScreenDevice(height=height, width = height * min(3.4, plots))
    par(bg="white")
    screens <- seq(0, 1, len=plots+1)
    screens <- split.screen(figs=cbind(screens[-plots-1], screens[-1], 0, 1))
    on.exit(close.screen(screens))
  }

  for (m in mode) {
    scr <- 0
    if (do.vario) {
      scr <- scr + 1
      vario <-
        regression(l.midbins, l.binvario, variable="D", pch=pch[1],
                   main="variogram method", scr=scr, cex=cex[1],
                   value =function(regr) ncol(ct$x) + 1 - regr$coeff[2]/2,
                   mode=m, cex.main=cex.main, ...)
    }
    if (do.box) {
      scr <- scr + 1
      box <- regression(Ml.box.sequ, l.box.count,
                        variable="D", main="box counting", scr=scr, pch=pch[2],
                        value =function(regr) ncol(ct$x) - 1 + regr$coeff[2],
                        mode=m, cex=cex[2], cex.main=cex.main, ...)
    }
    if (do.range) {
      scr <- scr + 1
      rnge <- regression(Ml.range.sequ, l.range.count, variable="D",
                         pch=pch[3],
                         main="max-min method", scr=scr, mode=m, cex=cex[3],
                         value = function(regr) ncol(ct$x) - 1 + regr$coeff[2],
                         cex.main=cex.main,...)
    }
    if (do.fft) {
      scr <- scr + 1
      if (length(l.lambda)>fft.max.regr && printlevel>=PL_IMPORTANT)
        cat("too large data set; no fft regression fit")
      else
        fft <- regression(l.lambda, l.I.lambda, main="periodogram", pch=pch[4],
                          scr=scr,
                          value=function(regr) 2.5 + 0.5 * regr$coeff[2],
                          averaging=fft.len!=dimen[1],
                          mode=m, variable="D", cex=cex[4], cex.main=cex.main,
                          ...)
    }
  }

  
  if (printlevel>=PL_SUBIMPORTANT) {
    cat("##################### fractal D ############### \n")
    cat(c(D.vario=vario$val,# D.box=box$val, D.range=rnge$val,
          D.fft=fft$val))
    cat("---- by interactively defined regression interval:\n")
    cat(c(D.vario=vario$val.u, # D.box=box$val.u, D.range=rnge$val.u,
            D.fft=fft$val.u))
    #cat("----alpha (Cauchy)----\n")
    #DIM <- 1ev
    #Print(2 * (DIM + 1 - c(D.vario=vario$val, D.box=box$val, D.range=rnge$val)))
    #Print(2 * (DIM + 1 - c(D.vario=vario$val.u, D.box=box$val.u,
    #                       D.range=rnge$val.u))) 
    cat("############################################### \n")
  }
  if (any(mode=="plot" | mode=="interactive")) close.screen(screens)
 
  return(list(vario=list(ev=ev, vario.n=vario.n,
                x=l.midbins, y=l.binvario, regr=vario$regr,
                sm=vario$sm,
                x.u=vario$x.u, y.u=vario$y.u, regr.u=vario$regr.u,
                D = vario$val, D.u = vario$val.u),
#              box = list(x=Ml.box.sequ, y=l.box.count, regr=box$regr,
#                sm=box$sm,
#                x.u=box$x.u, y.u=box$y.u, regr.u=box$regr.u,
#                D = box$val, D.u = box$val.u),
#              range=list(x=Ml.range.sequ, y=l.range.count, regr=rnge$regr,
#                sm=rnge$sm,
#                x.u=rnge$x.u, y.u=rnge$y.u, regr.u=rnge$regr.u,
#                D = rnge$val, D.u = rnge$val.u),
              fft=list(x=l.lambda, y=l.I.lambda, regr=fft$regr,
                sm=fft$sm,
                x.u=fft$x.u, y.u=fft$y.u, regr.u=fft$regr.u,
                D = fft$val, D.u = fft$val.u)
              )
         )
}








