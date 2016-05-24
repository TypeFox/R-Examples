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



## Methods for classes 'RFempVario' and 'RFfit'  #######################

summary_RFempVariog <- function(OP, object, ...) {
  variogram <- do.call(OP, list(object, "emp.vario"))
  if (is.array(variogram)) {
    dims <- dim(variogram)
    dims <- dims[dims > 1]
    dim(variogram) <- dims
  }
  l <- list(centers=do.call(OP, list(object, "centers")),
            variogram=variogram)
  obj <- do.call(OP, list(object, "sd"))
  if (length(obj) > 0) {
    l$sd <- obj
    if (is.array(variogram)) dim(l$sd) <- dims
  }

  obj <- do.call(OP, list(object, "phi.centers"))
  if (length(obj) > 1) l$phi <- obj
  obj <- do.call(OP, list(object, "theta.centers"))
  if (length(obj) > 1) l$theta <- obj
  obj <- do.call(OP, list(object, "T"))
  if (length(obj)>0) l$T <- obj
  obj <- do.call(OP, list(object, "vdim"))
  if (obj > 1) l$vdim <- obj
  class(l) <- "summary.RFempVariog"
  l
}

summary.RFempVariog <- function(object, ...) summary_RFempVariog("@",object,...)
summary.RF_empVariog <- function(object,...) summary_RFempVariog("$",object,...)



print.summary.RFempVariog <- function(x, ...) {
  cat("Object of class 'RFempVariog'\n")
  str(x, no.list=TRUE) #
  invisible(x)
}

print.RFempVariog <- function(x, ...) {
  print.summary.RFempVariog(summary.RFempVariog(x, ...))
}
print.RF_empVariog <- function(x, ...) {
  print.summary.RFempVariog(summary.RF_empVariog(x, ...))
}


setMethod(f="show", signature="RFempVariog",
          definition=function(object) print.RFempVariog(object)) #

## coercion methods
setAs("RFempVariog", "list",
      def=function(from) {
        li <- list()
        for (name in slotNames(from)) {
          li[[name]] <- eval(parse(text=paste("from@", name, sep="")))
        }
        return(li)
      })

setAs("RFfit", "RFempVariog", def=function(from) from@ev)

list2RFempVariog <- function(li) {
  RFopt.coords <- RFoptions()$coords
  return(new("RFempVariog",
             centers=li$centers,
             emp.vario=li$emp.vario,
             var=var,
             sd=li$sd,
             n.bin=li$n.bin,
             phi.centers=li$phi.centers,
             theta.centers=li$theta.centers,
             T=li$Tbins,
             vdim = li$vdim,
             coordunits = RFopt.coords$coordunits,
             varunits = RFopt.coords$varunits,
             call=li$call ))
}


## plot method


plotRFempVariogUnbinned <- function(x, coordunits, varunits, varnames,
                                    legend, ..., plotmethod="plot") {
  stop("currently not used ") # !!
  graphics <- RFoptions()$graphics
  dots = mergeWithGlobal(list(...))
  dotnames <- names(dots)
  coords <- GridTopology2gridVectors(cbind(x@centers$x, x@centers$T))  
  if (length(coords)>1) {    
    coords[[1]] <- sort(unique((coords[[1]] - min(coords[[1]])) *
                               rep(c(-1, 1), each=length(coords[[1]]))))
    lab.names <- dimnames(x@centers$x)[[2]]
    
    if (!(x@centers$spatialdim==1 && x@centers$Zeit) &&
        !(x@centers$x[3,2]==dim(x@emp.vario)[2]))
      coords[[2]] <- sort(unique((coords[[2]] - min(coords[[2]])) *
                                 rep(c(-1, 1), each=length(coords[[2]]))))
    else
      lab.names[2] <- "T"
  } else {
    x@centers <- coords[[1]]-min(coords[[1]])
    do.call(graphics::plot, args=c(dots, list(x=x, plot.nbin=FALSE)))
    return()
  }
  
  if (!("main" %in% dotnames)) {
    main <- "Variogram image plot"
    if (length(varnames)>0) main <- paste(main, "for", varnames)
    dots$main <- main
  }
  lab.names <- paste(lab.names, "-distance", sep="")
  
  idx <- lab.names != "T-distance"
  if (any(idx) && all(coordunits[idx] != ""))
      lab.names[idx] <-
        paste(lab.names[idx], " [", coordunits[idx], "]", sep="")
  if (!all(idx) && all(coordunits[!idx] != ""))
    lab.names[!idx] <-
      paste(lab.names[!idx], " [", coordunits[!idx], "]", sep="")
  
  if (!("xlab" %in% dotnames)) dots$xlab <- lab.names[1]
  if (!("ylab" %in% dotnames)) dots$ylab <- lab.names[2]
  if (!("xlim" %in% dotnames)) dots$xlim <- range(coords[[1]])
    if (!("ylim" %in% dotnames)) dots$ylim <- range(coords[[2]])
  
  idx1 <- coords[[1]] >= dots$xlim[1] & coords[[1]] <= dots$xlim[2]
  idx2 <- coords[[2]] >= dots$ylim[1] & coords[[2]] <= dots$ylim[2]
  coords[[1]] <- coords[[1]][idx1]
  coords[[2]] <- coords[[2]][idx2]
  
  dims <- dim(x@emp.vario)
  ev.plot <- matrix(x@emp.vario, nrow=dims[1], ncol=dims[2])
  ev.plot <- ev.plot[idx1, idx2]
  
  range.ev <- range(ev.plot)#x@emp.vario)
  col <- if ("col" %in% dotnames) dots$col else
  default.image.par(NULL, NULL)$default.col

  if (graphics$split.screen && legend) {
    scr <- split.screen( rbind(c(0,.85,0,1), c(.85,1,0,1)))
    screen(scr[1])
    par(mar=c(4,4,3,1))
  } else {
    scr <- NULL
    par(mar=c(4,4,3,1))
  }

  dots$type <- NULL
  do.call(plotmethod,
          args=c(dots, list(x=coords[[1]], y=coords[[2]], z=ev.plot)))
                                        #      xlab=lab.names[1],
                                        #      ylab=lab.names[2],
                                        #      main="Variogram image plot")
  if (legend) {
    stop("das ist doch quark?!")
   
    
    if (graphics$split.screen) {
      screen(scr[2])
      par(mar=c(4,0,3,3))
      z.legend <- seq(range.ev[1], range.ev[2], length=length(col))
      image(y=z.legend, x=c(0,1), z=rbind(z.legend, z.legend), axes=F, col=col,
            xlab="")
      axis(4)
      box()
    } else {
      my.legend(min(coords[[1]]), max(coords[[2]]), range(ev.plot),
                col=col, bg="white")
    }
  }

  if (graphics$close_screen) {
    close.screen(scr)
    scr <- NULL
  }
  return(invisible(scr))
}


RFplotEmpVariogram <- function(x, model = NULL, nmax.phi = NA, nmax.theta = NA,
                               nmax.T = NA,
                               plot.nbin = TRUE, plot.sd=FALSE, method = "ml",
                               variogram=TRUE,
                               boundaries = TRUE,
                               ...) {
  if (!variogram)
    stop("plot of estimated covariance functions not programmed yet.")

  graphics <- RFoptions()$graphics
  newx <- list()
  methodnames <- double()
  if (is(x, "RFfit") || is(x, "RF_fit")) {
    OP <- c("$", "@")[1 + is(x, "RFfit")]
    if(length(do.call(OP, list(x, "ev")))==0) stop("The fit does not contain an empirical variogram.")
    newx$autostart <- do.call(OP, list(x, "autostart"))
    newx$self <- do.call(OP, list(x, "self"))
    newx$plain <- do.call(OP, list(x, "plain"))
    newx$sqrt.nr <- do.call(OP, list(x, "sqrt.nr"))
    newx$sd.inv <- do.call(OP, list(x, "sd.inv"))
    newx$internal1 <- do.call(OP, list(x, "internal1"))
    newx$internal2 <- do.call(OP, list(x, "internal2"))
    newx$internal3 <- do.call(OP, list(x, "internal3"))
    newx$ml <- do.call(OP, list(x, "ml"))

    x <- do.call(OP, list(x, "ev"))
    if( is.null(method) )
      method <- if (length(newx$ml@name) > 0) "ml" else "plain"

    methodidx<-names(newx[unlist(lapply(newx, function(x) is(x, ZF_MODELEXT)))])
    methodidx <- method %in% methodidx
    methodnames <- method[methodidx]
    nomethodnames <- method[!methodidx]

    if( !all(methodidx) )
      warning( paste("The following method does not exist: ", nomethodnames) )
  } else {
    if (!is(x, "RFempVariog") && !is(x, "RF_empVariog"))
      stop("method only for objects of class 'RFempVariog' or 'RFfit'")
    OP <- c("$", "@")[1 + is(x, "RFempVariog")]
  }
 
  ev <- do.call(OP, list(x, "emp.vario"))
  sd <- do.call(OP, list(x, "sd"))
  centers <- do.call(OP, list(x, "centers"))
  phi.centers <- do.call(OP, list(x, "phi.centers"))
  theta.centers <- do.call(OP, list(x, "theta.centers"))
  Time <- do.call(OP, list(x, "T"))
  n.bin <- do.call(OP, list(x, "n.bin"))
  coordunits <- do.call(OP, list(x, "coordunits"))

  varnames <- if (is.matrix(ev)) dimnames(ev)[[2]][1] else names(ev)[1]
 
  dots = list(...)
  dotnames <- names(dots)
   if (!("type" %in% dotnames)) dots$type <- "b" 
  cex <- if ("cex" %in% dotnames) dots$cex else .8
  par(cex=cex, xaxs="i")
  dots$cex <- NULL
  if (!("pch" %in% dotnames)) dots$pch <- 19
             
  if(!("xlim" %in% dotnames)) dots$xlim <- range(centers)
  ylim.not.in.dotnames <- !("ylim" %in% dotnames)
  xlab <- if ("xlab" %in% dotnames) dots$xlab else "distance"
  dots$xlab <- NULL
  ylab.ev <- if ("ylab" %in% dotnames) dots$ylab else "semivariance"
  dots$ylab <- NULL
  
  main0 <- if ("main" %in% dotnames) dots$main else "Variogram plot"  
  dots$main <- NULL
  if (!is.null(main0)) oma.top <- 2 else oma.top <- 0
  

  has.sd <- !is.null(sd)
  
 if(!is.null(model)) {
   if (!is.list(model)) model <- list(model)
    if (!all(unlist(lapply(model, FUN=function(x) is(x, ZF_MODEL)))))
      stop("model must be (a list of elements) of class 'ZF_MODEL'")
    modelnames <-
      if(length(names(model)) > 0) names(model)
      else paste("model", 1:length(model))
    methodnames <- c(methodnames, modelnames)
    names(model) <- modelnames
    newx <- c(newx, model)
 }
  n.methods <- length(methodnames)


  n.phi <- min(nmax.phi, l.phi <- max(1,length(phi.centers)), na.rm=TRUE)
  n.theta <- min(nmax.theta, l.theta <- max(1,length(theta.centers)),
                 na.rm=TRUE)
  n.T <- min(nmax.T, l.T <- max(1,length(Time)), na.rm=TRUE)
  vdim <- dim(ev)
  vdim <- vdim[length(vdim)]

  if (n.phi > 6 || n.theta > 3 || n.T > 3)
    message("'If you feel the picture is overloaded, set the parameters 'nmax.phi', 'nmax.theta' and 'nmax.T'")
  
  halfphi.sector <- pi/(2*l.phi)
  halftheta.sector <- pi/(2*l.theta)
  phi.angles <- c(halfphi.sector, 0, -halfphi.sector) 
  theta.angles <- seq(-halftheta.sector, halftheta.sector, len=5) # len=odd!
  if (n.phi>1 && boundaries) {
    phi.angles <- phi.angles * 0.96
    theta.angles <- theta.angles * 0.96
  } 
  
  TandV <- (n.T > 1 && vdim > 1) && graphics$split_screen
  if (vdim>1 && length(varnames)==0)
    varnames <- paste("v", 1:vdim, sep="")

  range.nbin <- range(c(0, n.bin), na.rm=TRUE)
  ylim.nbin <- range.nbin * c(1,1.2)
  
  col.v <- col <- 
    if ("col" %in% dotnames) rep(dots$col, len=n.phi) else 1:max(n.phi)
  dots$col <- NULL

  if (n.methods > 0){
    dotsRFfit <- dots
    dotsRFfit$type <- "l"
    dotsRFfit$lwd <- 2
    ltyRFfit.v <- 1:n.methods
    dotsRFfit$lty <- NULL
  }


  dir2vario <- function(dir.vec, x.eval, x.time, method.model, v1, v2){
    x.space <- as.matrix(x.eval) %*% t(dir.vec)            
    for (j in (1:4)) { ## orig 20
      vario.vals <- try(RFvariogram(x = cbind(x.space, x.time),
                                    model = method.model,
                                    grid = FALSE,
                                    internal.examples_reduced=FALSE),
                        silent = FALSE)
       if(!is(vario.vals, "try-error")) {
        if (is.array(vario.vals)) {
          return(vario.vals[, v1, v2])
        } else {
          return(vario.vals)
        }
      }      
      x.space <- cbind(x.space, 0)    
    }

     
    return(NULL)
  }

  oma.left <- 6
  Screens <- if (TandV) c(n.T, n.theta) else c(n.T * vdim * vdim, n.theta)
  ArrangeDevice(graphics, Screens)

  all.scr <- scr <- (if (!graphics$split_screen) NULL else
              split.screen(if (TandV) c(vdim, vdim) else Screens))
    
  for (v1 in 1:vdim) {
    for (v2 in 1:vdim) {
      if (TandV) {
        scr <- split.screen(Screens, all.scr[v1 + (v2 - 1) * vdim])
        all.scr <- c(all.scr, scr)
      }
      if (vdim == 1) {
        main <-
          if (!is.null(main0) && length(varnames)>0)
            paste(main0, "for", varnames) else main0
      } else {
        main <-
          if (!TandV) main0
          else paste(main0, "for", varnames[v1], "vs.",  varnames[v2])
      }      
      if (ylim.not.in.dotnames)
        dots$ylim <- range(ev[,,,, v1, v2], na.rm=TRUE)
      for (iT in 1:n.T) {
        for (ith in 1:n.theta) {
          ## plot n.bin
          if (plot.nbin) {
            screen(scr[1])
            par(oma=c(4,oma.left,oma.top,0))
            scr2 <- split.screen(rbind(c(0,1,.2,1), c(0,1,0,.2)), screen=scr[1])
            all.scr <- c(all.scr, scr2)
            screen(scr2[2])
            par(mar=c(0,.5,0,.5))
            for (iph in 1:n.phi) {
              if (n.phi > 1) col <- col.v[iph]
              if (iph==1) {
                lab <- xylabs("bin centers", NULL, units=coordunits)
                plot(centers,
                     n.bin[ ,iph, ith, iT, v1, v2],
                     xlim=dots$xlim, ylim=ylim.nbin,
                     type=if (n.phi>1) "p" else "h",
                     col =if (n.phi>1) col else "darkgray", lwd=4,
                     pch=16, axes=FALSE, ylab="", xlab = lab$x
                     )
                box()
                at <- seq(range.nbin[1], range.nbin[2], len=3)
                if (ith==1)
                  axis(2, las=1, at=at, labels=format(at, scient=-1, dig=2),
                       outer=TRUE)
                else  axis(2, las=1, at=at, labels=FALSE)
                
                if (iT==n.T && (n.T > 1 || (v1==vdim && v2==vdim))) axis(1)
                if (ith==1) title(ylab="n.bin", line=5, outer=TRUE, adj=0)
                
              } else {
                points(centers, n.bin[ ,iph, ith, iT, v1, v2],
                       type="p", col=col, pch=16)
              }
            }
            screen(scr2[1])
          } else {
            screen(scr[1])
            par(oma=c(4,oma.left,oma.top,0))
          }
          
          ## plot emp.vario
          ##if (ith==1) par(mar=c(0,6,1,1)) else
          par(mar=c(0,.5,1,.5))
          
          plotted.meth <- NULL  # needed for legend
          
          for (iph in 1:n.phi) {
            if (n.phi>1) col <- col.v[iph]
            if (iph==1) {
              do.call(graphics::plot,
                      c(dots, list(
                          x=centers,
                          y=ev[,iph,ith,iT,v1,v2],
                                        #ylim=ylim.ev, type=type, pch=19,
                          col=col, axes=FALSE, ylab=""))
                      )
              box()
              axis(2, las=1, labels=(ith==1), outer=(ith==1))
              if (!plot.nbin) axis(1)          
              if (l.theta > 1 || l.T > 1 || vdim > 1) {
                L <- character(3)
                if (!TandV && vdim > 1) L[1] <- paste(varnames[c(v1,v2)],
                                                      collapse=":")
                if (l.T>1) L[2] <-
                  paste("T=",signif(Time[iT], 3), " ", sep="")
                if (l.theta>1) L[3] <-
                  paste(sep="", "theta=", signif(theta.centers[ith],3))
                legend("topleft", legend=paste(L[L!=""], collapse=","))
              }
              if (ith == 1) title(ylab=ylab.ev, line=5, outer=TRUE)
              if (has.sd && plot.sd)
                legend("topright", bty="n", #col=NA, lwd=NA,
                       legend="arrows represent +-1 sd intervals")
              
            } else {  ## iph > 1
              do.call(graphics::points,
                      c(dots, list(x=centers,
                                   y=ev[,iph,ith,iT,v1,v2],
                                   col=col))
                      )   #type=type,  pch=19)
            } ## if iph 
            
            k <- 1  
            if (n.methods > 0) {          
              for(i in 1:n.methods) {        
                method.model <- newx[[methodnames[i]]]
              
                if(length(method.model@name) == 0){
                  warning("The method '", methodnames[i], "' was not fitted.")
                  next
                }
                if (!is.null(phi.centers)) {
                  x.radial <- cbind(cos(phi.centers[iph]+phi.angles),
                                    sin(phi.centers[iph]+phi.angles))
                  if (!is.null(phi.centers))
                    x.radial <-
                      cbind(x.radial,
                            cos(rep(theta.angles, each = length(phi.angles))
                                + theta.centers[ith]))
                  
                } else x.radial <- matrix(1, nrow=1, ncol=1)
                
                x.time <- NULL
                if (!is.null(Time)) x.time <- Time[iT]
                x.eval <- seq(from = max(dotsRFfit$xlim[1],1e-3),
                              to = dotsRFfit$xlim[2], len = 150) 
                
                ## sehr genau Abschaetzung, indem mehrere (3)
                ## Winkel angeschaut werden und dann der mittlere
                ## Wert der Variogramme angeschaut wird, da ja auch
                ## das emp. Variogramm ueber ein Winkelintervall gemittelt wird
                dummy.vals <- as.matrix(apply(x.radial, 1, FUN = dir2vario,
                                              x.eval=x.eval, x.time=x.time,
                                              method.model=method.model,
                                              v1=v1, v2=v2))

               
               # Print(i, x.eval, dummy.vals);print.RMmodelFit(method.model)
                
                if(!is.null(dummy.vals)){
                 # Print(n.phi, boundaries)
                  if (n.phi>1 && boundaries) {
                    do.call(graphics::matplot,
                            args=c(dotsRFfit, list(x=x.eval,
                                y= if (ncol(dummy.vals) == 1) dummy.vals
                                else t(apply(dummy.vals, 1, range)),
                                add = TRUE, col=col, lty = 3)))
                  } else {
                     do.call(graphics::points,
                             args=c(dotsRFfit, list( x=x.eval,
                                 y = rowMeans(dummy.vals)  , 
                                 col=col, lty = ltyRFfit.v[k])))
                  }             
                  k <- k+1
                  if(iph == 1) plotted.meth <- c(plotted.meth, methodnames[i])
                }
                
              }
            }  
            
            
            if (has.sd && plot.sd) {
              sdnot0 <-  sd[ ,iph, ith, iT] != 0
              arrows(centers[sdnot0],
                     ev[sdnot0 ,iph, ith,iT] - sd[sdnot0,iph,ith,iT],
                     centers[sdnot0],
                     ev[sdnot0 ,iph, ith,iT] + sd[sdnot0,iph,ith,iT],
                     code=2, angle=90, length=0.05, col=col)
            }
            
          } # nphi
          
          
          pifactor <- signif((phi.centers[1:n.phi]) / pi, 2)
          
          len.mnames <- length(plotted.meth)
          string.emp <- "empirical"
           if(len.mnames > 0) {
            labels <-
              if (n.phi>1) paste('"phi=",', rep(pifactor, each = len.mnames+1),
                            ', pi, ', '", ')
              else '" '
            labels <- paste("c(",
                            paste('expression(paste(',labels,
                                  rep(c(string.emp, plotted.meth), l.phi),
                                  '"', '))', collapse=","),
                            ")")
            labels <- eval(parse(text=labels))
          } else
            labels <- if (n.phi > 1) 
              eval(parse(text=paste(
                           "c(",
                           paste('expression(paste("phi=",', pifactor, ', pi))',
                                 collapse=","),
                           ")")))
          
          
          if (l.phi > 1 || len.mnames > 0)# && iT==1 && ith==1) # auskommentiert auf Sebs wunsch        
            legend("bottomright", col=rep(col.v, each = len.mnames+1),
                   lwd=1, pch=rep(c(19,rep(NA, len.mnames)), l.phi), bty="n",
                   legend=labels, lty =
                   rep(c(1, if(len.mnames==0) NULL else ltyRFfit.v[1:len.mnames]),
                   l.phi))
          
          scr <- scr[-1]
        } # n.theta
      } # T
    } # vdim 1
  } # vdim 2
  dots$type <- NULL
  
  if (!is.null(main))
    do.call(graphics::title, args=c(dots, main=main, outer=TRUE))
  if (!is.null(xlab))
    do.call(graphics::title, args=c(dots, xlab=xlab, outer=TRUE))

  if (graphics$close_screen) {
    close.screen(all.scr)
    all.scr <- NULL
  }
 
  return(invisible(all.scr))
}

setMethod(f="plot", signature(x="RFempVariog", y="missing"),
          function(x, y, ...) RFplotEmpVariogram(x, ...))
