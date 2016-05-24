fleches <- function(x=TRUE,y=TRUE){
  extr <- par("usr")
  if (isTRUE(x)) arrows(extr[1],extr[3],extr[2]+0.05*diff(extr[1:2]),extr[3],
                        xpd=T,length=0.1)
  if (isTRUE(y)) arrows(extr[1],extr[3],extr[1],extr[4]+0.05*diff(extr[3:4]),
                        xpd=T,length=0.1)
}

