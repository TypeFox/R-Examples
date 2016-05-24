sizelegend<-function(se, am, pch=pch)
  {
    if(missing(pch)) pch= 1

    u = par('usr')
    ex = c(u[1]+ .05*(u[2]-u[1]), u[1]+ .2*(u[2]-u[1]))
    why = u[3]+.95*(u[4]-u[3])
    N = length(se)

    rect(u[1], u[3]+.9*(u[4]-u[3]) , u[1]+ .25*(u[2]-u[1]) , u[4], col="white", border=NA, xpd=TRUE)

    points(seq(from=ex[1], to=ex[2], length=N), rep(why, length=N), pch=pch, cex=se, xpd=TRUE)
    text(seq(from=ex[1], to=ex[2], length=N), rep(why, length=N),labels=am, pos=3, xpd=TRUE)

  }

