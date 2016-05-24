`PLTpicks` <-
function(picks, labs=NA, cols=NA)
{
##X##  plot a list of picks on a seismogram
##X## 
##X##   picks = vector of times relative to the start of the plot (seismogram)
##X##   labs = labels to plot next to picks
  ##X##  cols = vector of colors to plot line and label
  
  if(missing(labs)) { labs = rep(NA, length(picks)) }
  if(missing(cols)) { cols = rep(NA, length(picks)) }

  u=par("usr")
     pyt=u[4]-0.05*(u[4]-u[3])
  if(length(cols)<length(picks))  {cols = rep( cols[1],length(picks))  }
  if(length(labs)<length(picks))  {labs = rep( labs[1],length(picks))  }

  
     for(i in 1:length(picks))
     {
       if(!is.na(cols[i]))
         {
           col=cols[i]
         }
       else
         {
           col = 1
         }
       abline(v=picks[i], col=col, lty=3)
       if(!is.na(labs[i])) { text( picks[i], pyt, labs[i], adj=0, col=col) }
     }
     

}

