### Arrows.R

Arrows <- function (x1, y1, x2, y2,
                    size= 1,
                    width= 1.2/4/cin,
                    open=TRUE,
                    sh.adj=0.1,
                    sh.lwd=1,
                    sh.col=if(is.R()) par("fg") else 1,
                    sh.lty=1,
                    h.col=sh.col,
                    h.col.bo=sh.col,
                    h.lwd=sh.lwd,
                    h.lty=sh.lty,
                    verbose=FALSE)
  ## Author: Andreas Ruckstuhl, refined by Rene Locher
  ## Version: 2007-05-28
{
  arr <- ((x2-x1)^2+(y2-y1)^2) > (10*.Machine$double.eps)
  if (any(!arr)) {
    warning("Some of the arrows have undefined directions and are therefore not drawn.\n")
    x1 <- x1[arr]
    x2 <- x2[arr]
    y1 <- y1[arr]
    y2 <- y2[arr]
  }

  cin <- size * par("cin")[2]
  uin <- if (is.R())
    1/xyinch()
  else par("uin")
  x <- sqrt(seq(0, cin^2, length.out = floor(35 * cin) + 2))
  delta <-  sqrt(h.lwd)*par("cin")[2]*0.005
  x.arr <- c(-rev(x), -x)
  wx2 <- width * x^2
  y.arr <- c(-rev(wx2 + delta), wx2 + delta)
  deg.arr <- c(atan2(y.arr, x.arr), NA)
  r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)
  theta <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
  lx <- length(x1)
  Rep <- rep(length(deg.arr), lx)
  p.x2 <- rep(x2, Rep)
  p.y2 <- rep(y2, Rep)
  ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
  r.arr <- rep(r.arr, lx)

  if(open) lines((p.x2 + r.arr * cos(ttheta)/uin[1]),
                 (p.y2 + r.arr*sin(ttheta)/uin[2]),
                 lwd=h.lwd, col = h.col.bo, lty=h.lty) else
  polygon(p.x2 + r.arr * cos(ttheta)/uin[1], p.y2 + r.arr*sin(ttheta)/uin[2],
          col = h.col, lwd=h.lwd,
          border=h.col.bo, lty=h.lty)
  r.seg <- rep(cin*sh.adj, lx)
  th.seg <- theta + rep(atan2(0, -cin), lx)
  segments(x1, y1, x2+r.seg*cos(th.seg)/uin[1], y2+r.seg*sin(th.seg)/uin[2],
           lwd=sh.lwd, col=sh.col, lty=sh.lty)
  if(verbose) return(list(width=width))
} # Arrows
