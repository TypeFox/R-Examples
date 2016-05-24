
translation = function(pp,Length,Breadth,direction){
  
  if(direction==1){
    doubleSize = 2*Breadth
    # translation at right border
    newxvals = c(length(pp$x),1)
    newyvals = pp$y
    
    for(i in 1:length(pp$x)){
      newxvals[i] =  (Breadth + pp$x[i])
    }
    
    plot(pp$x,pp$y,pch=18,cex=2,col="blue",xlim=c(0,doubleSize),ylim=c(0,Length),
         xlab="x",ylab="y")
    points(newxvals,newyvals,pch=16,cex=2,col="red",xlim=c(0,doubleSize),ylim=c(0,Length))
  }
  if(direction==2){
    doubleSize = 2*Breadth
    # translation at left border
    newxvals = c(length(pp$x),1)
    newyvals = pp$y
    
    for(i in 1:length(pp$x)){
      newxvals[i] =  ( pp$x[i] - Breadth)
    }
    
    plot(pp$x,pp$y,pch=18,cex=2,col="blue",xlim=c(-Breadth,Breadth),ylim=c(0,Length),
         xlab="x",ylab="y")
    points(newxvals,newyvals,pch=16,cex=2,col="red",xlim=c(-Breadth,Breadth),ylim=c(0,Length))
    
  }
  if(direction==3){
    # translation at top border
    doubleSize = 2*Length
    newxvals = pp$x
    newyvals = c(length(pp$x),1)
    
    for(i in 1:length(pp$x)){
      newyvals[i] =  ( pp$y[i] + Length)
    }
    
    plot(pp$x,pp$y,pch=18,cex=2,col="blue",xlim=c(0,Breadth),ylim=c(0,doubleSize),
         xlab="x",ylab="y")
    points(newxvals,newyvals,pch=16,cex=2,col="red",xlim=c(0,Breadth),ylim=c(0,doubleSize))
    
  }
  
  if(direction==4){
    # translation at bottom border
    newxvals = pp$x
    newyvals = c(length(pp$x),1)
    
    for(i in 1:length(pp$x)){
      newyvals[i] =  -(Length - pp$y[i])
    }
    
    plot(pp$x,pp$y,pch=18,cex=2,col="blue",xlim=c(0,Breadth),ylim=c(-Length,Length),
         xlab="x",ylab="y")
    points(newxvals,newyvals,pch=16,cex=2,col="red",xlim=c(0,Breadth),ylim=c(-Length,Length))
  }
  list(newxvals,newyvals)
}