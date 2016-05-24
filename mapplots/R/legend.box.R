legend.box <-
function(x,y=NULL,maxradius,mab=1.2,inset=0,double=F){
  auto <- if (is.character(x)) 
      match.arg(x, c("bottomright", "bottom", "bottomleft", 
          "left", "topleft", "top", "topright", "right", "center")) else NA
  asp <- get.asp()
  h <- mab*2*maxradius
  w <- h*asp
  if(double) h <- h*2
  usr <- par("usr")
  inset <- rep(inset, length.out = 2)
  if(!is.na(auto)){
    insetx <- inset[1L] * (usr[2L] - usr[1L])
    left <- switch(auto, bottomright = , topright = , 
        right = usr[2L] - w - insetx, bottomleft = , 
        left = , topleft = usr[1L] + insetx, bottom = , 
        top = , center = (usr[1L] + usr[2L] - w)/2)
    insety <- inset[2L] * (usr[4L] - usr[3L])
    top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] + 
        h + insety, topleft = , top = , topright = usr[4L] - 
        insety, left = , right = , center = (usr[3L] + 
        usr[4L] + h)/2)
  } else {
    left <- x-1.2*asp*maxradius
    top <- y+1.2*maxradius
  }
return(c(left,top,left+w,top-h))
}

