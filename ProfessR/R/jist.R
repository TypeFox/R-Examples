`jist` <-
function(h, Z=1, L=1, col=2)
  {
    if(missing(col)) { col = 1 }
    if(missing(Z))  {  Z = h$grades  }
    if(missing(L))  {  L = h$lett  }
    
    
    ###  h = histogram structure
    ###  Z = scores
    #### L = letter grades or labels

    u = par("usr")
    

    if( length(h$lett)>0)
      {
        Z = h$grades
        L = h$lett
        h = h$hist
      }

  ###  plot(h)
   ## box()
    
  ###  abline(v=seq(from=round(min(Z[Z>0])), to=max(Z), by=2))
     axis(1, at=seq(from=10*round(min(Z[Z>0])/10), to=max(Z), by=10), labels=TRUE)
   
    axis(1, at=seq(from=10*round(min(Z[Z>0])/10), to=max(Z), by=2), labels=FALSE)

    
    
    J = rep(1, length(Z))

    Kounts = c(h$counts)
    
    for(i in 1:length(Kounts))
      {
        if(Kounts[i]==0) next
        x1 = h$breaks[i]
        x2 = h$breaks[i+1]
        w = which(Z>=x1&Z<x2)
        k = Kounts[i]

        v = Kounts[i]*seq(from=1, to=length(w))/length(w)
        J[w] = v
        
      }

    #############  deal with the last group 
    x2 = h$breaks[length(h$breaks)]
    w = which(Z>=x2)
    if(length(w)>0)
      {
        v = length(w)*seq(from=1, to=length(w))/length(w)
        J[w] = v
      }
    ############

    #####  o = order(Z); cbind(Z[o], J[o])
    
    text(Z, J, labels=L, col=col, xpd=TRUE, cex=.8, font=2)

    
    invisible(J)
  }

