`pimg` <-
function(xmin,xmax,xn,
                 ymin,ymax,yn) {
  res <- list(xbound=c(xmin,xmax,xn),
       ybound=c(ymin,ymax,yn),
       offset=c(1,1),
       image=NULL)
       
       class(res) <- c("pimg", "list")
  res
}

