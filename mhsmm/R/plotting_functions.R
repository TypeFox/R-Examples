plot.hsmm <- function(x,...) {
  tmp = x$model$d
  plot(1:nrow(tmp),tmp[,1],type='l',...,ylab="d(u)",xlab="u",ylim=range(tmp))
  for(i in 2:x$J)
    lines(tmp[,i],type='l',col=i)
  legend("topright",legend=1:x$J,col=1:x$J,lty=1)
}


plot.hsmm.data <- function(x,...) {
  	plot(ts(x$x),...)
 	  if(!is.null(x$s)) .add.states(x$s,ht=axTicks(2)[1],time.scale=1)
  	if(length(x$N)>1) abline(v=cumsum(x$N),lty=2) 	  
}

addStates <- function (states,x=NULL,ybot = axTicks(2)[1], ytop=ybot + (axTicks(2)[2] - axTicks(2)[1])/5,dy  = ytop - ybot,greyscale = FALSE, leg = NA, 
                J = length(unique(states)), time.scale = 1, shiftx = 0) 
{

  draw.it <- function(hats, ybot, ytop, cols, greyscale){
    ##cat("ybot", ybot, "ytop", ytop, "\n")
    for (ii in 1:length(hats$state)){
      if (greyscale) { 
        rect(xleft   = hats$intervals[ii],
             ybottom = ybot,
             xright  = hats$intervals[ii + 1], 
             ytop    = ytop,
             col = cols[hats$state[ii]], border = 1)
      } else {
        rect(xleft   = hats$intervals[ii],
             ybottom = ybot,
             xright  = hats$intervals[ii + 1], 
             ytop    = ytop,
             col = cols[hats$state[ii]], border = cols[hats$state[ii]])
      }
    }
  }

  
  if (is.null(states)){
    states <- x
    if (!is.list(states))
      states <- list(states)
    x <- seq_along(states[[1]])
  } else {
    if (!is.list(states))
      states <- list(states)
    if(is.null(x)) x <- seq_along(states[[1]])        
  }
  
##   cat("states:\n");
##   print(states)
##   cat("x:\n");
##   print(x)

  x <- as.numeric(x)
  rr  <- range(x)
  
  J = length(unique(states))
  if (greyscale) { 
    cols <- c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252", "#252525")
  } else {
    cols <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  }


  st.list <- lapply(states,
                    function(st){
                      runs = rle(st)
                      cs  <- cumsum(c(0, runs$lengths))
                      hats <- list(intervals=rr[1]+ diff(rr)*cs/max(cs), states=runs$values)
                      hats
                    })



  ##cat("dy:", dy, "\n")
  for (ii in seq_along(st.list)){
    draw.it (st.list[[ii]], ybot, ytop, cols, greyscale)
    ybot <- ytop + .2*dy
    ytop <- ybot + dy
  }
  
  if (any(!is.na(leg))) 
    legend("topleft", legend = leg, fill = cols, bg = "white")
}
