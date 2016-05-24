`coords.pimg` <-
function(pimg) {
  list(x=seq(pimg$xbound[1],pimg$xbound[2],length=pimg$xbound[3]),
       y=seq(pimg$ybound[1],pimg$ybound[2],length=pimg$ybound[3]))
}

