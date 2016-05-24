


preph <- function(h,o){
	m <- h$merge
	n <- nrow(m)
	se <- matrix(1L,nrow=4,ncol=n)
	mat <- 	sort.int(o,method="quick",na.last=NA,index.return=TRUE)$ix
	#mat <- order(o)
	for (i in 1:n){
		j <- m[i,1]
		if (j < 0 ) {
			a1 <- a2 <- mat[-j]
			}
		else {
			a1 <- se[1,j]
			a2 <- se[4,j]
			}
		j <- m[i,2]
		if (j < 0 ) {
			b1 <- b2<- mat[-j]
			}
		else {
			b1 <- se[1,j]
			b2 <- se[4,j]
			}
		if (a1 < b1)
		  se[,i] <- c(a1,a2,b1,b2)
		else se[,i]	<-	c(b1,b2,a1,a2)	
			
	}
	list(se=se,sort.int=mat,order=o)
}



# preph <- function(h,o){
	# m <- h$merge
	# n <- nrow(m)
	# se <- matrix(1L,nrow=4,ncol=n)
	# mat <- 	sort.int(o,method="quick",na.last=NA,index.return=TRUE)$ix
	# #mat <- order(o)
	# for (i in 1L:n){
		# j <- m[i,1]
		# if (j < 0 ) {
			# a1 <- a2 <- mat[-j]
			# }
		# else {
			# a1 <- se[1,j]
			# a2 <- se[4,j]
			# }
		# j <- m[i,2]
		# if (j < 0 ) {
			# b1 <- b2<- mat[-j]
			# }
		# else {
			# b1 <- se[1,j]
			# b2 <- se[4,j]
			# }
		# if (a1 < b1)
		  # se[,i] <- c(a1,a2,b1,b2)
		# else se[,i]	<-	c(b1,b2,a1,a2)			
	# }	
	 # mode(se) <- "integer"	
	# se
# }





# preph1 <- function(h,o,oldse,kmax){
	# m <- h$merge
	# n <- nrow(m)
	# se <- oldse
	# if (m[kmax,2] > 0 || m[kmax,1] > 0){
	# sek <- se[,kmax]
	# #mat <- 	sort.int(o,method="quick",na.last=NA,index.return=TRUE)$ix
	# #mat <- order(o)
	# for (i in 1L:kmax){
		# if (se[1,i] >= sek[1] & se[4,i] <= sek[4]){
		# j <- m[i,1]
		# if (j < 0 ) {
			# a1 <- a2 <-  match(-j,o)
			# }
		# else {
			# a1 <- se[1,j]
			# a2 <- se[4,j]
			# }
		# j <- m[i,2]
		# if (j < 0 ) {
			# b1 <- b2<-  match(-j,o)
			# }
		# else {
			# b1 <- se[1,j]
			# b2 <- se[4,j]
			# }
		# if (a1 < b1)
		  # se[,i] <- c(a1,a2,b1,b2)
		# else se[,i]	<-	c(b1,b2,a1,a2)			
	# }
	# }	
	 # mode(se) <- "integer"	
	# }
	# se
# }


preph1 <- function(h,o,oldse,kmax){
	m <- h$merge
	n <- nrow(m)
	se <- oldse$se
	mat <- oldse$sort.int
	y <- which(oldse$order !=o)
    mat[o[y]]<-y

    if (m[kmax,2] > 0 || m[kmax,1] > 0){
    sek1 <- se[1,kmax]
    sek4 <- se[4,kmax]
	for (i in 1:kmax){
		if (se[1,i] >= sek1 & se[4,i] <= sek4){
		j <- m[i,1]
		if (j < 0 ) {
			a1 <- a2 <-  mat[-j]
			}
		else {
			a1 <- se[1,j]
			a2 <- se[4,j]
			}
		j <- m[i,2]
		if (j < 0 ) {
			b1 <- b2<-  mat[-j]
			}
		else {
			b1 <- se[1,j]
			b2 <- se[4,j]
			}
		if (a1 < b1)
		  se[,i] <- c(a1,a2,b1,b2)
		else se[,i]	<-	c(b1,b2,a1,a2)			
	}
	}	
	}
	list(se=se,sort.int=mat,order=o)
}



#-----------DendSer works for R0, T0, R1, T1, C0, R01 and T01--------------






defaultOp <- function(costfn,n) {
	if (identical(costfn,costLS))
	  list(t0,"up",1, FALSE)
	 else if (identical(costfn,costED))
	   list(r1,"up",1, FALSE)
	else if (identical (costfn,costARc)){
		if (n >= 300)
		 list(t0,"down",10,FALSE)
		else list(r01,"down",10,TRUE)
	}
	else if (identical (costfn,costBAR)){
		if (n >= 300)
		 list(c0,"down",10,TRUE)
		else list(r01,"down",10,TRUE)
	}
	else { # PL, LPL and others
		if (n >= 1000)
		 list(c0,"down",10,TRUE)
		else list(r01,"down",10,TRUE)
	}
}


defaultOpE <- function(costfn) {
	if (identical(costfn,costLS))
	  list(t0,"up",1, FALSE)
	 else if (identical(costfn,costED))
	   list(r1,"up",1, FALSE)
	else if (identical (costfn,costARc))
	  list(t0,"down",10, FALSE)
      else  list(c0,"down",10, TRUE) #PL,:LPL, BAR 
}



DendSer <- function(h, ser_weight,  cost=costBAR, node_op=NULL, costArg=NULL, maxloops=NULL, saveinfo=FALSE,direction=NULL,GW=NULL,...){
	
	if (identical(cost,costED)) return (gclus::reorder.hclust(h,ser_weight)$order)
	if (identical(cost,costLS)) return (leafSort(h,ser_weight))
	
	if (is.character(node_op))
	   node_op <- getFromNamespace(node_op, "DendSer")
	    
	if (! inherits(h,"hclust"))
	   stop("'h' must inherit from hclust")
	   
	   
    
    n <- length(h$order)
       
	if (isTRUE(all.equal(cost, costLS))){
		if (missing(ser_weight) || length(ser_weight) != n  || !is.numeric(ser_weight))
           stop(paste("'ser_weight' must be numeric vector of length",n))
        names(ser_weight) <- NULL
        }
           
    if (identical(cost, costBAR) || identical(cost, costPL) || identical(cost, costLPL) || 
                    identical(cost, costARc) || identical(cost, costED)) {
     if (missing(ser_weight) ||  !is.numeric(ser_weight))
	    stop(paste("'ser_weight' must be a numeric symmetric ",n ,"x" ,n," matrix or a dist",sep=""))

	 if (inherits(ser_weight,"dist"))
		   ser_weight <- as.matrix(ser_weight)
	 if (!is.matrix(ser_weight) || nrow(ser_weight) != ncol(ser_weight))
		  stop(paste("'ser_weight' must be a numeric symmetric ",n ,"x" ,n," matrix or a dist",sep=""))
		 
	 rownames(ser_weight) <- colnames(ser_weight) <- NULL
	 }
	 
    
    
	
	
  if (is.null(node_op) || is.null(direction) || is.null(maxloops)){
        nop <- defaultOp(cost,n)
     if (is.null(node_op)) node_op <- nop[[1]]  
     if (is.null(direction)) direction <- nop[[2]] 
     if (is.null(maxloops)) maxloops <- nop[[3]] 
     if (is.null(GW)) GW <- nop[[4]]         
     }
     
     
  if (GW && is.matrix(ser_weight)) h<- reorder.hclust(h,ser_weight)
	
  ord <- h$order

  weight <- ser_weight


  if (is.null(costArg)) costArg<- defaultcostArg(cost,weight)
	reprep <- direction=="down" || isTRUE(all.equal(node_op, t01))
  # saveScore <- !isTRUE(all.equal(cost,costED))
  

  stoploop <- FALSE
  denditer <- seq_len(n-1)
  if (direction=="down") denditer <- rev(denditer)
  nloops <- 0
  nperms <- 0
  nswitches <- 0
  stoploop <- FALSE
  score_cur <- cost(weight,ord,costArg)
  # if (saveScore) score_cur <- cost(weight,ord,costArg)
  
   while (!stoploop) {
  	nloops <- nloops+1
  	ord0 <- ord
  	se <- preph(h,ord)
   
    for (i in denditer){
    	p <- node_op(i,ord,se$se,h)
    		
	    # update order
	    if (!is.null(p)){
  	  	  score1 <- score_cur
  	  	  
	    	
	  	  if (!is.matrix(p)) {
             # score2 <- cost(weight,p,costArg,node=i,se=se) 
             score2 <- cost(weight,p,costArg) 
             nperms <- nperms+1
 		      if (score2 < score1){
 		      	nswitches<-nswitches+1
 		      	 ord <- p
 		      	  score_cur <- score2
	             if (reprep) se <- preph1(h,ord,se,i)
                 
              }
  		    }
	  	  else {
            scores <- apply(p,2, function(pj) cost(weight,pj,costArg))
             # scores <- apply(p,2, function(pj) cost(weight,pj,costArg,node=i,se=se))
	        j <- which.min(scores)
	        nperms <- nperms+ncol(p)
	        score2 <- scores[j]
            if (score2 < score1) {
            	nswitches<-nswitches+1
               ord <- p[,j]
 		        score_cur <- score2
 		         if (reprep) se <- preph1(h,ord,se,i)
	   	        
               }  
           }
       }
    }
     
     stoploop <- all(ord0==ord) || nloops == maxloops
 
  }
  if (saveinfo) 
    return(list("order" = ord, "nloops" =nloops, "nperms" = nperms,nswitches=nswitches,score=score_cur))
  else
    return(ord)
}


preph <- cmpfun(preph)
preph1 <- cmpfun(preph1)


DendSer <- cmpfun(DendSer)


