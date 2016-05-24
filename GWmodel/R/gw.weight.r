################## NAME: gw.weight ###################
## AUTHOR: IG (and MC gw_ridge_regression_02.R)
## DESCRIPTION: The function finds the weights for one particular point, given the distances and the desired kernel function. 
## FUNCTION DEPENDING ON: gw.dist
## ARGUMENTS IN:
#### vdist: numeric vector or matrix of distances (from gw.dist.r)
#### bw: scalar, bandwidth or number of nearest neighbours
#### kernel: text vector of function chosen
########## boxcar: wgt=1 if dist < bw, wgt=0 otherwise  
########## gaussian: wgt = exp(-.5*(vdist/bw)^2)
########## bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise   
########## tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise    
########## adaptive: if TRUE calulate the adaptive kernel, and bw correspond to the number of nearest neighbour
####
## ARGUMENTS OUT:
#### if vdist is a vector the output is a vector of weights of the length of vdist
#### if vdist is a matrix the output is a matrix of weights of the dim of vdist, containing in column i the weights for observation i
## REFERENCES: Book pg 56-57
########################

####### Fixed Kernel

# exponential kernel

gw.weight.exponential <- function(vdist,bw){
	 exp(-vdist/bw)
	}


# Boxcar kernel

gw.weight.box <- function(vdist,bw){
	{vdist<= bw}*1
	}

# Gaussian kernel

gw.weight.gau <- function(vdist,bw){
	exp(vdist*vdist/{-2*bw*bw})
	}

# Fixed Bisquare kernel
	
gw.weight.bis<- function(vdist,bw){
	cond<- which(vdist< bw)          # condition for locations within the bandwidth
	
	wgt<-numeric(length(vdist))
	sqwgt<- 1-vdist[cond]^2/{bw*bw}
	wgt[cond]<-sqwgt*sqwgt
	
	if (is.matrix(vdist)) wgt<-matrix(wgt,dim(vdist))
	
	wgt
	}

# Fixed Tricube kernel
	
gw.weight.tri<- function(vdist,bw){
	cond<- which(vdist<= bw)          # condition for locations within the bandwidth
	
	wgt<-numeric(length(vdist))
	sqwgt<- 1-vdist[cond]^3/{bw*bw*bw}
	wgt[cond]<-sqwgt*sqwgt*sqwgt
	
	if (is.matrix(vdist)) wgt<-matrix(wgt,dim(vdist))
	
	wgt
	}

####### Adaptive Kernel

# Adaptive Bisquare kernel 
## bw correspond to the number of nearest neighbours

#####

gw.weight.gau.ad<- function(vdist,bw){
	
	if (is.matrix(vdist)){
		rnk<-apply(vdist,2,rank,ties.method='first')
		bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance	
    if(bw>0)
    	wgt<-t(exp(t(vdist*vdist)/{-2*bw*bw}))
    else
      wgt <- diag(1, dim(vdist)[1], dim(vdist)[2])
	}else{
    rnk<-rank(vdist,ties.method='first') # ranked distances
    cond<- which(rnk <= bw) 
    bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance
    if(bw>0)
      wgt<- exp(vdist*vdist/{-2*bw*bw})
    else
    {
      wgt <- numeric(length(vdist))
      wgt[cond] <- 1
    }  
	}
	wgt
	}

#####

gw.weight.exp.ad<- function(vdist,bw){
  if (is.matrix(vdist)){
		rnk<-apply(vdist,2,rank,ties.method='first')
		bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance	
    if(bw>0)
		  wgt<-exp(-vdist/bw)
    else
      wgt <- diag(1, dim(vdist)[1], dim(vdist)[2])
	}
	else
	{
	  rnk<-rank(vdist,ties.method='first') # ranked distances
    cond<- which(rnk <= bw) 
	  bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance
    if(bw>0)
	    wgt<- exp(-vdist/bw)
    else
    {
      wgt <- numeric(length(vdist))
      wgt[cond] <- 1
    }
	}
	wgt
	}

#####


gw.weight.bis.ad<- function(vdist,bw){
	
	if (is.matrix(vdist))
  {
		rnk<-apply(vdist,2,rank,ties.method='first')
		cond<- rnk <=bw  
		bw<- vdist[rnk == bw]
	  wgt<- matrix(0,nrow(vdist),ncol(vdist))
    if(bw>0)
    {
		  mdist<- matrix(vdist[cond==1],nrow=ncol(vdist),byrow=TRUE)
		  sqwgt<- t(1-mdist^2/{bw*bw})
		  wgt[cond==1]<-sqwgt*sqwgt
    }
    else
      diag(wgt) <- 1
		
  }
  else
  {
    rnk<-rank(vdist,ties.method='first') # ranked distances
    cond<- which(rnk <= bw)               # condition for locations less than bw-th
    bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance
    wgt<-numeric(length(vdist))
    if(bw >0)
    {
      bw<-bw*bw
      sqwgt<- 1-vdist[cond]^2/bw
      wgt[cond]<-sqwgt**2
    }
    else
    {
      wgt[cond]<- 1
    }
	}
	wgt

	}


#####

gw.weight.tri.ad<- function(vdist,bw){
	
	if (is.matrix(vdist)){
		rnk<-apply(vdist,2,rank,ties.method='first')
		cond<- rnk<= bw 
		bw<- vdist[rnk == bw]
		wgt<- matrix(0,nrow(vdist),ncol(vdist))
		mdist<- matrix(vdist[cond==1],nrow=ncol(vdist),byrow=TRUE)
    if(bw>0)
    {
		  sqwgt<- t(1-mdist^3/{bw*bw*bw})
		  wgt[cond==1]<-sqwgt*sqwgt*sqwgt
    }
    else
      diag(wgt) <- 1
		
	}else{
	rnk<-rank(vdist,ties.method='first') # ranked distances
	cond<- which(rnk <= bw)               # condition for locations less than bw-th
	bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance	
	wgt<-numeric(length(vdist))
  if(bw>0)
  {
	  sqwgt<- 1-vdist[cond]^3/{bw*bw*bw}
	  wgt[cond]<-sqwgt*sqwgt*sqwgt
  }
  else
    wgt[cond]<-1
 }
	wgt
	}


#####

gw.weight.box.ad<- function(vdist,bw){
	
	if (is.matrix(vdist)) rnk<-apply(vdist,2,rank,ties.method='first')
	
	else rnk<-rank(vdist,ties.method='first') # ranked distances

	{rnk <= bw}*1
	}


# MAIN FUNCTION

gw.weight<-function(vdist,bw,kernel,adaptive=FALSE){
	
	    if(adaptive==FALSE) switch(kernel,
	    gaussian = gw.weight.gau(vdist,bw),
        bisquare = gw.weight.bis(vdist,bw),
        tricube  = gw.weight.tri(vdist,bw),
        boxcar   = gw.weight.box(vdist,bw),
        exponential = gw.weight.exponential(vdist,bw))
        
        else switch(kernel,
	      gaussian = gw.weight.gau.ad(vdist,bw),
        bisquare = gw.weight.bis.ad(vdist,bw),
        tricube  = gw.weight.tri.ad(vdist,bw),
        boxcar   = gw.weight.box.ad(vdist,bw),
        exponential = gw.weight.exp.ad(vdist,bw))
}


######### End of the Code
