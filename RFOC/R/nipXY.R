`nipXY` <-
function(MEC, x=x, y=y, size=c(1,1), fcol=gray(0.90), nipcol="black", cex=.4)
  {
    if(missing(fcol)) { fcol = gray(0.90) }
    if(missing(nipcol)) { nipcol =  fcol  }
    if(is.null(MEC$UP)) MEC$UP = TRUE
    if(missing(cex)) { cex=.4  }

    LP1 = lowplane( MEC$az1, MEC$dip1, col=1, UP=MEC$UP, PLOT=FALSE)

    lines(x+size[1]*LP1$x, y+size[2]*LP1$y, col=fcol, xpd=TRUE )

    A1 = RSEIS::TOCART(MEC$U$az, 90-MEC$U$dip)
    A2 = RSEIS::TOCART(MEC$V$az,  90-MEC$V$dip)
    N = CROSSL(A1, A2)
    Q = qpoint(N$az, N$dip , lab="N", UP=MEC$UP, PLOT=FALSE)
    points(x+size[1]*Q$x, y+size[2]*Q$y, pch=16, col=nipcol, cex=cex, xpd=TRUE)

    invisible(list(Q=Q, N=N) )
  }

