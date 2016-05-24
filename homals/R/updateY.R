`updateY` <-
function(dframe,x,y,active,rank,level,sets,verbose=0){
  nobj<-dim(x)[1]
  ndim<-dim(x)[2]
  nset<-length(sets)
  for (l in 1:nset) {                     #number of sets (if no sets specified, sets = nvar)
  	indi <- sets[[l]]                     #set index
    jndi <- indi[which(active[indi])]     #index for active variables in set 
    if (length(jndi) == 0) next()         #no active variables
  	
    ii <- which(!is.na(dframe[,jndi[1]])) #non-NA values in active variables
    if (length(ii) == 0) next()           #only NA
  	
    ss <- sumSet(dframe,nobj,ndim,y,jndi) #scores for variable (set), 0 where NA
    
    for (j in jndi) {                     #runs over active variables in set
		  gg <- dframe[ii,j]
      yy <- y[[j]]
      d <- as.vector(table(gg))           
		  s1 <- sum((x[ii,]-ss[ii,])^2)
		  ss[ii,] <- ss[ii,]-yy[gg,]
		  yc <- computeY(gg,x[ii,]-ss[ii,])   #compute scores
		  yy <- restrictY(d,yc,rank[j],level[j],verbose=verbose)$y
		  ss[ii,] <- ss[ii,]+yy[gg,]
		  s2 <- sum((x[ii,]-ss[ii,])^2)
		  y[[j]]<-yy
		  if (verbose > 1) {
                    cat("\nSet:", l)
                    cat("\nAfter Variable:", j)
                    cat("\nLoss:", c(s1, s2))
         		#at("Set: ",formatC(l,digits=3,width=3),
		  	#" After Variable: ",formatC(j,digits=3,width=3),
		  	#" Loss: ", formatC(c(s1,s2),digits=6,width=9, format="f"),"\n")	
                  } 	
   	}
  }
return(y)
}

