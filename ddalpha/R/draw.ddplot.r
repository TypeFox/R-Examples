draw.ddplot <- function(ddalpha, depth.space, cardinalities,
                        main = "DD plot", xlab = "C1", ylab = "C2", classes = c(1,2), colors = c("red", "blue", "green")){
  
  if(!missing(ddalpha)){
    
    col = c()
    points = NULL
    for (c in classes){
    points = rbind(points, ddalpha$patterns[[c]]$depths[,classes])
    col = c(col, rep(colors[c], ddalpha$patterns[[c]]$cardinality))
    }
    
    if (xlab == "C1" && ylab == "C2"){
      xlab = ddalpha$patterns[[1]]$name
      ylab = ddalpha$patterns[[2]]$name
    }
    
    plot(points, col = col, main = main, xlab = xlab, ylab = ylab)
    
  } else if(!missing(depth.space)){
    
    col = c()
    for (c in classes){
      col = c(col, rep(colors[c], cardinalities[c]))
    }
    
    plot(depth.space[,classes], col = col, main = main, xlab = xlab, ylab = ylab)
    
  } else stop("Both 'ddalpha' and 'depth.space' area missing")
  
}