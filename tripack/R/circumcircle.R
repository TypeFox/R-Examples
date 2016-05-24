circumcircle <- function(x,y=NULL,num.touch=2,plot=FALSE,debug=FALSE)
  {

    circumcircle3 <- function(ch,i,j,...){
      # helper function, given a convex hull ch and two indices i,j of
      # its points find a third point k on ch so that the circumcircle
      # of the triangle i,j,k contains all points of ch and has minimum
      # radius.
      kmin <- NULL
      for(k in (1:npch)[-c(i,j)]){
        # cat(paste("<",i,",",j,",",k,">\n"))
        circ <- circum(ch$x[c(i,j,k)],ch$y[c(i,j,k)])
        if(debug){
          plot(tri.mesh(ch$x[c(i,j,k)],ch$y[c(i,j,k)]),
               add=TRUE,lty="dotted")
          circles(circ$x,circ$y,circ$radius,lty="dotted",col="grey")
        }
        d <- as.matrix(dist(rbind(cbind(circ$x,circ$y),
                                  cbind(ch$x[-c(i,j,k)],
                                        ch$y[-c(i,j,k)]))
                            ))[-1,1]
        if(all(d<=circ$radius)){
          if(is.null(rmin)){
            rmin <- circ$radius
          }            
          if(debug){
            plot(tri.mesh(ch$x[c(i,j,k)],ch$y[c(i,j,k)]),
                 add=TRUE,lty="dotted")
            circles(circ$x,circ$y,circ$radius,lty="dotted")
          }
          if(circ$radius<=rmin){
            rmin <- circ$radius
            kmin <- k
            xmin <- circ$x
            ymin <- circ$y
            rmin <- circ$radius
          }
        }
      }
      if(!is.null(kmin)){
        if(debug)
          plot(tri.mesh(ch$x[c(i,j,k)],
                        ch$y[c(i,j,k)]),
               add=TRUE,col="red")
        xtri <- ch$x[c(i,j,kmin)]
        ytri <- ch$y[c(i,j,kmin)]
      } else {
        xtri <- ytri <- NULL
      }
      list(x=xmin,y=ymin,radius=rmin,xtri=xtri,ytri=ytri)
    }
    
    rmin <- NULL
    xmin <- NULL
    ymin <- NULL

    if(debug) plot <- TRUE

    if(is.null(x))
      stop("argument x missing.")
    if(is.null(y)){
      y<-x$y
      x<-x$x
      if (is.null(x) || is.null(y))
        stop("argument y missing and x contains no $x or $y component.")
    }

    
    n <- length(x)
    tri <- tri.mesh(x,y,duplicate="remove")
    ch <- convex.hull(tri)
    npch <- length(ch$x)

    if(num.touch!=2 & num.touch!=3)
      stop("num.touch can only take values 2 or 3!")

    # get the diameter, it's somewhat tricky to extracrt the index of the
    # maximum out of the return value of dist():
    chdist <- dist(cbind(ch$x,ch$y))
    idmax <- which.max(as.matrix(chdist))[1] # take only the first!
    jmax <- (idmax-1)%/%npch+1
    imax <- (idmax-1)%%npch+1
    if(imax==0) imax <- npch
    
    if(plot){
      # taken partially from eqscplot
      ratio <- 1
      tol <- 0.02
      xlim <- range(x[is.finite(x)])
      ylim <- range(y[is.finite(y)])
      midx <- 0.5 * (xlim[2L] + xlim[1L])
      xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2L] - xlim[1L])
      midy <- 0.5 * (ylim[2L] + ylim[1L])
      ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2L] - ylim[1L])
      oldpin <- par("pin")
      xuin <- oxuin <- oldpin[1L]/abs(diff(xlim))
      yuin <- oyuin <- oldpin[2L]/abs(diff(ylim))
      if (yuin > xuin * ratio) 
        yuin <- xuin * ratio
      else xuin <- yuin/ratio
      xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
      ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5

      xrange <- diff(range(x))
      plot(x,y,
           xlim=xlim,
           ylim=ylim)
    }
    # (ch$x[imax],ch$y[imax]) and (ch$x[jmax],ch$y[jmax])
    # are the two points where the minumum circle touches if it
    # touches in only two points.

    # check if the circle with center between (ch$x[imax],ch$y[imax]) and
    # (ch$x[jmax],ch$y[jmax]) and diameter = dist((ch$x[imax],ch$y[imax]),
    # (ch$x[jmax],ch$y[jmax])) already contains all points, ...
    xc <- (ch$x[imax]+ch$x[jmax])/2
    yc <- (ch$y[imax]+ch$y[jmax])/2
    radiusc <- sqrt((ch$x[imax]-ch$x[jmax])^2+(ch$y[imax]-ch$y[jmax])^2)/2
    # this radius gives also a lower bound for searching later:
    rlower <- radiusc
    
    if(all(as.matrix(dist(rbind(cbind(xc,yc),
                                cbind(ch$x[-c(imax,jmax)],
                                      ch$y[-c(imax,jmax)]))
                          ))[-1,1]<=radiusc)){
      if(num.touch==2){
        # then if num.touch==2: this is the solution, ...
        if(plot)
          points(ch$x[c(imax,jmax)],ch$y[c(imax,jmax)],col="red",pch="X")
        ret <- list(x=xc,y=yc,radius=radiusc)
        ret3 <- NULL
      } else {
        #      if num.touch==3: find exactly one third point on the hull
        #                       so that the circumcircle of the resulting
        #                       triangle contains all points and is minimal.
        #                       ( circumcircle3() ) ...
        ret3 <- circumcircle3(ch,imax,jmax,plot,debug)
        ret <- list(x=ret3$x,y=ret3$y,radius=ret3$radius)
      }
    } else {
      # else find exactly one third point on the hull so that the circumcircle
      #      of the resulting triangle contains all points and is minimal.
      #      ( circumcircle3() ) ...
      ret3 <- circumcircle3(ch,imax,jmax,plot,debug)
      ret <- list(x=ret3$x,y=ret3$y,radius=ret3$radius)
    }
    # if circumcircle3() doesnt find a third point to imax and jmax
    # (see e.g. circtest2) search all remaining combinations on the convex hull.

    # enclosing rectangle (r1,r2,r3,r4):
    r1 <- c(min(ch$x),min(ch$y))
    r2 <- c(max(ch$x),min(ch$y))
    r3 <- c(max(ch$x),max(ch$y))
    r4 <- c(min(ch$x),max(ch$y))
    # upper bound for the minimum circle radius taken from the
    # enclosing rectangle, used later:
    cupper <- circum(c(r1[1],r2[1],r3[1]),c(r1[2],r2[2],r3[2]))
    rupper <- cupper$radius
    xmin <- cupper$x
    ymin <- cupper$y
    rmin <- rupper
    if(is.null(ret$x)){
      for(i1 in 1:npch){
        for (i2 in (1:npch)[-(1:i1)]){
          for(i3 in (1:npch)[-c(i1,i2)]){
            # skip combinations which contain imax and jmax:
            if(imax %in% c(i1,i2,i3) & jmax %in% c(i1,i2,i3)){
#              cat(paste("skipping: <",i1,",",i2,",",i3,"> contains ",imax,",",jmax,"\n"))
            } else {
#              cat(paste("<",i1,",",i2,",",i3,">\n"))
              circ <- circum(ch$x[c(i1,i2,i3)],ch$y[c(i1,i2,i3)])
              # use bounds to avoid calling dist() with no need:
              if(circ$radius<=rupper & circ$radius>=rlower){
                d <- as.matrix(dist(rbind(cbind(circ$x,circ$y),
                                          cbind(ch$x[-c(i1,i2,i3)],
                                                ch$y[-c(i1,i2,i3)]))
                                    ))[-1,1]
                if(all(d<=circ$radius)){
                  if(circ$radius<=rmin){
                    if(debug){
                      plot(tri.mesh(ch$x[c(i1,i2,i3)],ch$y[c(i1,i2,i3)]),
                           add=TRUE,lty="dotted")
                      circles(circ$x,circ$y,circ$radius,lty="dotted")
                    }
                    rmin <- circ$radius
                    i1min <- i1
                    i2min <- i2
                    i3min <- i3
                    xmin <- circ$x
                    ymin <- circ$y
                    rmin <- circ$radius
                  }
                }
              }
            }
          }
        }
      }
      if(debug)
        plot(tri.mesh(ch$x[c(i1min,i2min,i3min)],
                      ch$y[c(i1min,i2min,i3min)]),
             add=TRUE,col="blue")
      if(plot)
        points(ch$x[c(i1min,i2min,i3min)],
               ch$y[c(i1min,i2min,i3min)],col="red",pch="X")
      ret <- list(x=xmin,y=ymin,radius=rmin)
    } else {
      if(plot & !is.null(ret3))
        points(ret3$xtri,
               ret3$ytri,col="red",pch="X")
    }
    if(plot)
      circles(ret$x,ret$y,ret$radius,col="red")
    ret
  }
    


    
