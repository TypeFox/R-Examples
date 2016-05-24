`ov.p.grid` <-
function (typ=NA, fun=NA, dis=NA, p1=1:49/50, p2=p1, a1=0, a0=1, plt=FALSE, invisible=FALSE, wire=FALSE, round=FALSE) { 
  p                  <- matrix(nrow=length(p1), ncol=length(p2)) 
  dimnames(p)        <- list(p1, p2) 
  names(dimnames(p)) <- c("p1", "p2") 
  i <- 0 
  for (x in p1) { 
    i <- i+1 
    j <- 0 
    for (y in p2) { 
      p[i, j <- j+1] <- ov.p(typ=typ, fun=fun, dis=dis, p1=x, p2=y, a1=a1, a0=a0, round=FALSE)
      }
    }
  if (plt) {
    fig1      <- expand.grid(p1=p1, p2=p2) 
    fig1$p    <- as.vector(p) 
    maintitle <- "Overall p-values" 
    subtitle  <- paste("typ =",
                       typ,
                       ", fun =",
                       fun,
                       ", dis =",
                       dis,
                       ", a1 =",
                       round(a1, digits=round),
                       ", a0 =",
                       round(a0, digits=round)
                       )
    if (!wire) {
      surface <- cloud(p ~ p1 * p2, 
                       fig1, 
                       xlim   = c(0, 1), 
                       ylim   = c(0, 1), 
                       zlim   = c(0, 1), 
                       scales = list(arrows = FALSE), 
                       pch    = ".", 
                       main   = maintitle, 
                       sub    = subtitle, 
                       zlab   = "p"
                       )
      }
    else {
      surface <- wireframe(p ~ p1 * p2, 
                           fig1, 
                           xlim   = c(0, 1), 
                           ylim   = c(0, 1), 
                           zlim   = c(0, 1), 
                           scales = list(arrows = FALSE), 
                           main   = maintitle, 
                           sub    = subtitle, 
                           zlab   = "p",
                           drape=TRUE,
                           colorkey=FALSE
                           )
      }
    print(surface)
    }
  if (round) {
    p <- round(p, round)
    }
  if (invisible) {
    invisible(p)
    }
  else {
    p
    }
  }

