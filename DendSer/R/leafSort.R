# leafSort <-
# function(x,w) {
    # merges <- x$merge
    # n <- nrow(merges)
    # avg <- matrix(0,n,2)
    # size <- matrix(1,n,2)
    # for (i in 1:n) {
	# j <- merges[i,1]
	# k <- merges[i,2]
	# if (j<0){
	  # avg[i,1] <- w[-j]
	# } else {
	# size[i,1] <- size[j,1]+size[j,2]
	# avg[i,1] <- (avg[j,1]*size[j,1]+avg[j,2]*size[j,2])/size[i,1]
	# }
	# if (k<0){
	  # avg[i,2] <- w[-k]
	# } else {
		# size[i,2] <- size[k,1]+size[k,2]
	    # avg[i,2] <- (avg[k,1]*size[k,1]+avg[k,2]*size[k,2])/size[i,2]
	# }
	# }
    
    # clusters <- as.list(1:n)
    # for (i in 1:n) {
	# j <- merges[i,1]
	# k <- merges[i,2]
	
	# cj <- if(j < 0) -j else clusters[[j]]
	# ck <- if(k < 0) -k else clusters[[k]]

	# clusters[[i]]<- if (avg[i,1] <= avg[i,2]) c(cj,ck) else c(ck,cj)
	# } 
		# clusters[[n]]   
    
# }


leafSort <-
function(h,w) {
    merges <- h$merge
    n <- nrow(merges)
    
    if (length(w) != n +1 || !is.numeric(w))
           stop(paste("'w' must be numeric vector of length",n+1))
         
    avg <- matrix(0,n,2)
    size <- matrix(1,n,2)
    oldres <- h$order
    ord<-sort.int(oldres,method="quick",na.last=NA,index.return=TRUE)$ix
    startp <- vector("numeric",n)
    for (i in 1:n) {
	  j <- merges[i,1]
	  k <- merges[i,2]
	if (j<0){
	   startp[i] <- ord[-j]
	  avg[i,1] <- w[-j]
	} else {
	startp[i] <- startp[j]
	size[i,1] <- size[j,1]+size[j,2]
	avg[i,1] <- (avg[j,1]*size[j,1]+avg[j,2]*size[j,2])/size[i,1]
	}
	if (k<0){
	  avg[i,2] <- w[-k]
	} else {
		size[i,2] <- size[k,1]+size[k,2]
	    avg[i,2] <- (avg[k,1]*size[k,1]+avg[k,2]*size[k,2])/size[i,2]
	}
	}
	res <- oldres
   for (i in 1:n) {
   	   if (avg[i,1] > avg[i,2]){
   	   	   x <- startp[i]
   	   	   s1 <- size[i,1]
   	   	   s2 <- size[i,2]
   	   	  r <- x:(x+s1+s2-1)
   	   	  r1 <- c(s1:(s1+s2-1),0:(s1-1))+x
   	   	  res[r]<- res[r1]

    	   }
   	}
   	res
   	}
   	   

leafSort <- cmpfun(leafSort)
# leafSort1 <- cmpfun(leafSort1)
