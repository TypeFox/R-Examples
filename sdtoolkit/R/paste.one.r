`paste.one` <-
function(x, y, x.init, y.init, box, paste.alpha,
                      mass.min, threshold, d, n)
{
  box.new <- box
  mass <- length(y)/n
  y.mean <- mean(y)
  n.box <- length(y)

  box.init <- apply(x.init, 2, range)
  
  if (is.vector(x)) x <- as.matrix(t(x))
  
  y.mean.paste <- matrix(0, nrow=2, ncol=d)
  mass.paste <- matrix(0, nrow=2, ncol=d)
  box.paste <- matrix(0, nrow=2, ncol=d)
  x.paste1.list <- list()
  x.paste2.list <- list()
  y.paste1.list <- list()
  y.paste2.list <- list()

  box.paste1 <- box
  box.paste2 <- box
  
  for (j in 1:d){   
    ## candidates for pasting
    box.diff <- (box.init[2,] - box.init[1,])[j]
		
    box.paste1[1,j] <- box[1,j] - box.diff*paste.alpha
    box.paste2[2,j] <- box[2,j] + box.diff*paste.alpha
    
    x.paste1.ind <- in.box(x=x.init, box=box.paste1, d=d, boolean=TRUE)
    x.paste1 <- x.init[x.paste1.ind,]
    y.paste1 <- y.init[x.paste1.ind]
    
    x.paste2.ind <- in.box(x=x.init, box=box.paste2, d=d, boolean=TRUE)
    x.paste2 <- x.init[x.paste2.ind,]
    y.paste2 <- y.init[x.paste2.ind]
    
		if(box.diff > 0){
#		if(TRUE){
			while (length(y.paste1) <= length(y) & box.paste1[1,j] >= box.init[1,j])
			{
				box.paste1[1,j] <- box.paste1[1,j] - box.diff*paste.alpha
				x.paste1.ind <- in.box(x=x.init, box=box.paste1, d=d, boolean=TRUE)
				x.paste1 <- x.init[x.paste1.ind,]
				y.paste1 <- y.init[x.paste1.ind]
				#print(c(length(y.paste1), length(y), box.paste[1,j], box.init[1,j]))
			}
			
			while (length(y.paste2) <= length(y) & box.paste2[2,j] <= box.init[2,j])
			{
				box.paste2[2,j] <- box.paste2[2,j] + box.diff*paste.alpha
				x.paste2.ind <- in.box(x=x.init, box=box.paste2, d=d, boolean=TRUE)
				x.paste2 <- x.init[x.paste2.ind,]
				y.paste2 <- y.init[x.paste2.ind]
			}
		}
   
    ## y means of pasted boxes
    y.mean.paste[1,j] <- mean(y.paste1)
    y.mean.paste[2,j] <- mean(y.paste2)

    ## mass of pasted boxes
    mass.paste[1,j] <- length(y.paste1)/n
    mass.paste[2,j] <- length(y.paste2)/n
    
    x.paste1.list[[j]] <- x.paste1
    y.paste1.list[[j]] <- y.paste1
    x.paste2.list[[j]] <- x.paste2
    y.paste2.list[[j]] <- y.paste2
    box.paste[1,j] <- box.paste1[1,j]
    box.paste[2,j] <- box.paste2[2,j]
  }

  ## break ties by choosing box with largest mass
  
  y.mean.paste.max <- which(y.mean.paste==max(y.mean.paste, na.rm=TRUE), arr.ind=TRUE)
  
  if (nrow(y.mean.paste.max)>1) {
     y.mean.paste.max <- cbind(y.mean.paste.max, mass.paste[y.mean.paste.max])
     y.mean.paste.max.ind <- y.mean.paste.max[order(y.mean.paste.max[,3], decreasing=TRUE),][1,1:2]
  } else {
    y.mean.paste.max.ind <- as.vector(y.mean.paste.max)       
  }
  ## paste along dimension j.max
  j.max <- y.mean.paste.max.ind[2]
	
  ## paste lower 
  if (y.mean.paste.max.ind[1]==1){
     x.new <- x.paste1.list[[j.max]] 
     y.new <- y.paste1.list[[j.max]]
     box.new[1,j.max] <- box.paste[1,j.max]
  }  else if (y.mean.paste.max.ind[1]==2){ ## paste upper
     x.new <- x.paste2.list[[j.max]] 
     y.new <- y.paste2.list[[j.max]]
     box.new[2,j.max] <- box.paste[2,j.max]
  } 
  
  mass.new <- length(y.new)/n
  y.mean.new <- mean(y.new)

  if ((y.mean.new > threshold) & (mass.new >= mass.min) & (y.mean.new >= y.mean)
      & (mass.new > mass))
    return(list(x=x.new, y=y.new, y.mean=y.mean.new, box=box.new, mass=mass.new))
 
}

