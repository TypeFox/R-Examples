
#peel.one - 

#the main change I made was adding a new quantile function pquantile, that handles
#things more appropriately for PRIMs purposes (see that function).  

`peel.one` <-
function(x, y, box, peel.alpha, mass.min, threshold, d, n,
					peel_crit)
{
  
  box.new <- box
  mass <- length(y)/n
  
  peel.alpha <- max(peel.alpha, 1/nrow(x)) #Added 2009-10-25

  if (is.vector(x)) return(NULL)
  
  y.mean <- mean(y)
	
	#Matrix to store the means first row is for boxes restricted from below,
	#second row is for boxes restricted from above:
  y.mean.peel <- matrix(0, nrow=2, ncol=d)
	#Equivalent for boxes that get peeled off, only for use in alternate peeling
	#criteria:
	y.mean.peeled <- matrix(0, nrow=2, ncol=d)
	
	#Matrix to store volume of resulting boxes - same indexing as y.mean.peel
  box.vol.peel <- matrix(NA, nrow=2, ncol=d) 
  
	#Matrix to store volume of the "little b" boxes that get peeled off, only for 
	#use in alternate peeling criteria:
	box.vol.peeled  <- matrix(0, nrow=2, ncol=d) 
#	box.supp.peeled <- matrix(NA, nrow=2, ncol=d)  
	box.supp.peeled <- matrix(0, nrow=2, ncol=d)  
	
	ranges	<- apply(x,2,range)
	spans 	<- ranges[2,]-ranges[1,]
#  if(length(which(spans!=0))==0){
#		#browser()
#		stop("Out of variables with any variation.") #would be better to make it 
#  }                                              #gracefully exit and end the 
																								 #peeling process...
	
	#for (j in which(spans!=0)){
	#print("peeled")
	for (j in 1:d){
   # print(j)
		if(TRUE){
#		if(j%in%which(spans!=0)){
			box.min.new <- pquantile(x[,j], peel.alpha, ptype="lowend")
			box.max.new <- pquantile(x[,j], 1-peel.alpha, ptype="highend")
		#Possibly switch to this as faster version in the future, though seems to 
		#cause problems for some reason, even with type=2)
		#	box.min.new <- quantile(x[,j], peel.alpha, type=2)
		#	box.max.new <- quantile(x[,j], 1-peel.alpha, type=2)
		} else {
			box.min.new <- min(x[,j])
			box.max.new <- max(x[,j])
			if(box.min.new!=box.max.new) stop("what?")
		}
		
    y.mean.peel[1,j] <- mean(y[x[,j] >= box.min.new])
    y.mean.peel[2,j] <- mean(y[x[,j] <= box.max.new])
		
		y.mean.peeled[1,j] <- mean(y[x[,j] < box.min.new])
    y.mean.peeled[2,j] <- mean(y[x[,j] > box.max.new])
		
		box.supp.peeled[1,j] <- sum(x[,j] < box.min.new)
		box.supp.peeled[2,j] <- sum(x[,j] > box.min.new)
		
		#Initialize boxes that will be subset to equal the current full box:
    
		#These are the boxes will hold the actual new box of interest:;
		#(1 and 2) indicate restriction from below and above, respectively)
		box.temp1 <- box
    box.temp2 <- box
		
		#These are the boxes that hold the boxes that get peeled away:
		box.temp1.peeled <- box
		box.temp2.peeled <- box
		
    if (!is.na(box.min.new)){
      box.temp1[1,j] <- box.min.new  #new box
			box.temp1.peeled[2,j] <- box.min.new #box to get peeled away
    }
    if (!is.na(box.max.new)){
      box.temp2[2,j] <- box.max.new #new box
			box.temp2.peeled[1,j] <- box.max.new #box to get peeled away
    }   
		
		
		
#    box.vol.peel[1,j] <- vol.box(box.temp1)
#    box.vol.peel[2,j] <- vol.box(box.temp2)   
		
#		box.vol.peeled[1,j] <- vol.box(box.temp1.peeled)
#		box.vol.peeled[2,j] <- vol.box(box.temp2.peeled)
		
		box.vol.peel[1,j] <- vol.box(box.temp1, x)
    box.vol.peel[2,j] <- vol.box(box.temp2, x)   
		
		box.vol.peeled[1,j] <- vol.box(box.temp1.peeled, x)
		box.vol.peeled[2,j] <- vol.box(box.temp2.peeled, x)
		
  }
  
	#Some code for starting to implement alternative peeling criteria -- 
	#specifically Eq 14.5 of F&F 1999 paper.
	toprint <- y.mean.peel
	#colnames(toprint) <- colnames(x)
	#print(toprint)
	if(peel_crit==1){
		y.mean.peel.max.ind <- which(y.mean.peel==max(y.mean.peel, na.rm=TRUE), arr.ind=TRUE)
		#print(y.mean.peel.max.ind)
	} else if(peel_crit==2){
		#Final FF form: evaled.peel.crit2 <- y.mean.peel - y.mean.peeled
		#How it is conceived, and original form:
		evaled.peel.crit2 <- (y.mean.peel - y.mean)/box.supp.peeled
#		if(all(is.na(evaled.peel.crit2))) evaled.peel.crit2[] <- 1
		y.mean.peel.max.ind <- which(evaled.peel.crit2==max(evaled.peel.crit2, na.rm=TRUE), arr.ind=TRUE)
	} else if(peel_crit==3){
		evaled.peel.crit3 <- (y.mean.peel - y.mean)*
													(nrow(x) - box.supp.peeled)/box.supp.peeled
#		if(all(is.na(evaled.peel.crit3))) evaled.peel.crit3[] <- 1
		y.mean.peel.max.ind <- which(evaled.peel.crit3==max(evaled.peel.crit3, na.rm=TRUE), arr.ind=TRUE)											
	} else {
		stop("There is no implemented peeling criteria associated with the value 
					that was passed to the peel_crit argument")
	}
	
  ## break ties by choosing box with largest volume
	
  nrr <- nrow(y.mean.peel.max.ind) 
  if (nrr > 1){
    box.vol.peel2 <- rep(0, nrr)
    for (j in 1:nrr)
      box.vol.peel2[j] <- box.vol.peel[y.mean.peel.max.ind[j,1],
                                       y.mean.peel.max.ind[j,2]]
#   print(box.vol.peel2)
    row.ind <- which(max(box.vol.peel2)==box.vol.peel2)[1] #added 5/28/08 the [1] to deal with ties
#	 print(row.ind)
	} else {
    row.ind <- 1
  }
	
#	print("gotabove")
  y.mean.peel.max.ind <- y.mean.peel.max.ind[row.ind,]
#  print("gotbelow")
	
	## peel along dimension j.max
  j.max <- y.mean.peel.max.ind[2]  #why is this [2]? Asked 2011-10-13

  
  if (y.mean.peel.max.ind[1]==1){ 				## peel lower 
    box.new[1,j.max] <- pquantile(x[,j.max], peel.alpha, ptype="lowend")
    x.index <- x[,j.max] >= box.new[1,j.max] 
  } else if (y.mean.peel.max.ind[1]==2) { ## peel upper 
    box.new[2,j.max] <- pquantile(x[,j.max], 1-peel.alpha, ptype="highend")
    x.index <- x[,j.max] <= box.new[2,j.max]
  }
 
  x.new <- x[x.index,]
  y.new <- y[x.index]
  mass.new <- length(y.new)/n
  y.mean.new <- mean(y.new)

  ## if min. y mean and min. mass conditions are still true, update
  ## o/w return NULL  
  # cat("mm:",mass.min, mass.new,mass,y.mean,y.mean.new,"\n")     #diagnostic
  if ((y.mean.new >= threshold) & (mass.new >= mass.min) & (mass.new < mass) & (y.mean < 1))
    return(list(x=x.new, y=y.new, y.mean=y.mean.new, box=box.new,
                mass=mass.new))
}

