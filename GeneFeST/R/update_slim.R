#source("math.R")

#constants
#e_f 		       <-  0.05
#sd_prior_beta          <-  1 #3
#mean_prior_beta        <- -1  #-5.15
#lambda 	       <-  1


X_update_d_alphaco_i <- function(id){ # id


 	if(id[1]==FALSE){ # Pilot runs
 	  id      <- 1:BBB2$locnum
 	  locnum2 <- BBB2$locnum
 	}else{                # Main Loop
 	  locnum2 <- length(id)
 	}
##

	 sub_freq_locus    <- BBB2$freq_locus[id,]
	 sub_freq_pop      <- lapply(BBB2$freq_pop,function(yy){return(yy[id,])})
         old_alpha         <- BBB2$d_alpha[id]
         sub_sample_size   <- BBB2$sample_size[id,]
         sub_GROUP         <- BBB2$GROUP[id]

##

 	old_theta     <- outer(old_alpha,BBB2$d_beta,"+")
 	old_theta     <- exp(-(old_theta))
 	old1          <- lgamma(old_theta)-lgamma(sub_sample_size + old_theta)
 	
  val <- matrix(0,locnum2,dim(BBB2$freq_locus)[2])
  for(xx in 1:BBB2$popnum){
    xyz <- old_theta[,xx]*sub_freq_locus
    val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
  }
  
  old1 <- rowSums(old1)
  oldL <- old1 + rowSums(val,na.rm=TRUE)

  ## modification ------------------------------------  
    oldL <- tapply(oldL,sub_GROUP,sum)
  ##--------------------------------------------------
  
    # e             <- rnorm(locnum2,0,sqrt(BBB2$var_prop_alpha))
    # new_alpha     <- old_alpha + e
  
  ### modification -----------------------------------------------

   if(!BBB2$PILOT){

    whatgroup       <- unique(sub_GROUP)
    xxx             <- BASIX.match(whatgroup,BBB2$GROUP)
    e_group         <- rnorm(length(whatgroup),0,sqrt(BBB2$var_prop_alpha[xxx]))

  #
    BBB2$new_alpha   <- old_alpha
    BBB2$jjj <- 1
    sapply(whatgroup,function(x){
    what            <- x==sub_GROUP
    BBB2$new_alpha[what] <- old_alpha[what] + e_group[BBB2$jjj]
    BBB2$jjj <- BBB2$jjj + 1
    })

    new_alpha <- BBB2$new_alpha
  #

 }else{

    e_group    <- rnorm(BBB2$N.REGIONS,0,sqrt(BBB2$var_prop_alpha))
    new_alpha  <- old_alpha + e_group[BBB2$GROUP] 
    whatgroup  <- 1:BBB2$N.REGIONS

 }

  #
  ##--------------------------------------------------------------
  
  new_theta <- outer(new_alpha,BBB2$d_beta,"+")
  new_theta <- exp(-(new_theta))
  new1      <- lgamma(new_theta)-lgamma(sub_sample_size+new_theta)
  
  val <- matrix(0,locnum2,dim(BBB2$freq_locus)[2])
 	for(xx in 1:BBB2$popnum){
          xyz <- new_theta[,xx]*sub_freq_locus
 	  val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
 	}
  
  new1 <- rowSums(new1)
  newL <- new1 + rowSums(val,na.rm=TRUE)
 
  ## modification ----------------------------
  newL            <- tapply(newL,sub_GROUP,sum)
   
  old_alpha_group <- tapply(old_alpha,sub_GROUP,function(x){return(x[1])})
  new_alpha_group <- tapply(new_alpha,sub_GROUP,function(x){return(x[1])})
  
 # print(old_alpha)
 #  print(new_alpha)

  ##-----------------------------------------
 

#A        	       <-  -oldL + newL + (old_alpha^2-new_alpha^2)/(2*BBB2$sd_prior_alpha^2)

## mod ------------------------
A         	       <-  -oldL + newL + (old_alpha_group^2-new_alpha_group^2)/(2*BBB2$sd_prior_alpha^2)
##-----------------------------

# r        	      <-  runif(locnum2)
  
## mod ---------------------------  
  r                    <-  runif(length(whatgroup))
## mod ---------------------------

check    	       <-  log(r) < A

 # BBB2$d_alpha[id[check]]      <- new_alpha[check]  
 # BBB2$acc_alpha[id[check]]    <- BBB2$acc_alpha[id[check]] + 1


# modification -------------------------------------
  # update groups
  	
  group_idX   <- which(check)
  group_id    <- whatgroup[group_idX]
  
  
#  make this also in C

  if(length(group_idX)>0){
   
  #res <- .Call("update_group2",group_id,BBB2$GROUP,BBB2$d_alpha,new_alpha,BBB2$acc_alpha)
  #BBB2$d_alpha   <- res[[1]]
  #BBB2$acc_alpha <- res[[2]]

   mm         <- sapply(group_id,function(xx){xx==BBB2$GROUP}) 
   BBB2$iii   <- 1
   apply(mm,2,function(x){
   BBB2$d_alpha[x]   <- new_alpha_group[group_idX[BBB2$iii]] # new_alpha[x]  
   BBB2$acc_alpha[x] <- BBB2$acc_alpha[x] + 1	  
   BBB2$iii <- BBB2$iii + 1

  })
 }

# ---------------------------------------------------	

}# End of Function



update_d_betaco <- function(){

   
  # Init 
  #val <- rep(NA,2)
   val <- BBB2$GLOBAL_INIT2
  # ---------------
 
  old_beta    <- BBB2$d_beta 
  old_theta   <- outer(BBB2$d_alpha,old_beta,"+")
  old_theta   <- exp(-(old_theta))
  
  old1 <- lgamma(old_theta)-lgamma(BBB2$sample_size+old_theta)
  
  for(xx in 1:BBB2$popnum){
      xyz     <- old_theta[,xx]*BBB2$freq_locus
      val[xx] <- sum(lgamma(BBB2$freq_pop[[xx]]+xyz)-lgamma(xyz),na.rm=TRUE)
  }

  old1 <- colSums(old1)
  oldL <- old1 + val
  #--------- calculate new likelihood
  e          <- rnorm(BBB2$popnum,0,BBB2$var_prop_beta)
  BBB2$d_beta <- BBB2$d_beta + e
  
  new_theta <- outer(BBB2$d_alpha,BBB2$d_beta,"+")
  new_theta <- exp(-(new_theta))
  new1      <- lgamma(new_theta)-lgamma(BBB2$sample_size+new_theta)
  for(xx in 1:BBB2$popnum){
      xyz     <- new_theta[,xx]*BBB2$freq_locus
      val[xx] <- sum(lgamma(BBB2$freq_pop[[xx]]+xyz)-lgamma(xyz),na.rm=TRUE)
  }
  new1 <- colSums(new1)
  newL <- new1 + val

A                <-  -oldL + newL #+ (old_beta^2-BBB2$d_beta^2)/(2*BBB2$sd_prior_beta^2)
#+ (old_alpha_group^2-new_alpha_group^2)/(2*BBB2$sd_prior_alpha^2)
#+ (old_beta-BBB2$d_beta)*(old_beta+BBB2$d_beta-2*BBB2$mean_prior_beta)/(2*BBB2$sd_prior_beta^2)
r                <-  runif(BBB2$popnum)
check            <-  log(r)>A
BBB2$d_beta[check]        <- old_beta[check]
BBB2$acc_beta[!check]     <- BBB2$acc_beta[!check] + 1

 
}# End of Function


update_freq_codominant <- function(){


lambda    <- 1
# choose randomly two haplotype of each locus
# samp      <-  sapply(BBB2$hapcount,function(x){return(sample(1:x,size=2)) })

theta     <-  outer(BBB2$d_alpha,BBB2$d_beta,"+")
theta     <-  exp(-theta)
# ok

#BBB2$iii       <- 1
# frequencies of haplotypes of the two random haps
#samp_freq <- apply(samp,2,function(x){
#	      r   <- BBB2$freq_locus[BBB2$iii,x]
#	      BBB2$iii <- BBB2$iii+1
#	      return(r)
#          })
#ok

## mod
samp_freq <- t(BBB2$freq_locus)
## end mod


#print(dim(samp_freq))
#print(dim(BBB2$freq_locus))

# haplotype distribution for each pop of the two random haps
# BBB2$iii           <- 1
samp_pop_freq <- vector("list",BBB2$popnum)
for(yy in 1:BBB2$popnum){

	samp_pop_freq[[yy]] <- t(BBB2$freq_pop[[yy]])

	#samp_pop  <- apply(samp,2,function(x){
	#	      r       <- BBB2$freq_pop[[yy]][BBB2$iii,x]
	#	      BBB2$iii <- BBB2$iii+1
	#	      return(r)
        #	   })
	#samp_pop_freq[[yy]] <- samp_pop
	#BBB2$iii <- 1
}
#ok
#print(dim(samp_pop_freq[[1]]))

old_L <- BBB2$GLOBAL_INIT3
# 1 <- m
# 2 <- n
for(xx in 1:BBB2$popnum){

	old_L <- old_L + lgamma(samp_pop_freq[[xx]][1,] + theta[,xx] * samp_freq[1,]) - 
	         lfactorial (samp_pop_freq[[xx]][1,]) - lgamma(theta[,xx]*samp_freq[1,])+
	         lgamma(samp_pop_freq[[xx]][2,]+theta[,xx]*samp_freq[2,]) - 
                 lfactorial(samp_pop_freq[[xx]][2,]) - lgamma(theta[,xx]*samp_freq[2,])		

}

old_m <- samp_freq[1,]
old_n <- samp_freq[2,]



#### BwPrp
old_m_n 	 <- old_m + old_n
old_m_ancestral  <- old_m + BBB2$e_ancestral
old_m_ancestral2 <- old_m - BBB2$e_ancestral

checkmin <- cbind(old_m_n,old_m_ancestral)
checkmax <- cbind(BBB2$GLOBAL_INIT3,old_m_ancestral2)
minn     <- apply(checkmin,1,min)
maxx     <- apply(checkmax,1,max)
BwPrp    <- minn - maxx
###################


# veränder Frequenzen
u     <- runif(BBB2$locnum,0,BwPrp) + maxx 
# samp_freq_new
samp_freq_new     <- samp_freq
samp_freq_new[1,] <- u
samp_freq_new[2,] <- old_n + old_m - u

### FwPrp
ss1   <- colSums(samp_freq_new)
ss2   <- samp_freq_new[1,] + BBB2$e_ancestral
ss3   <- samp_freq_new[1,] - BBB2$e_ancestral
minn  <- apply(cbind(ss1,ss2),1,min)
maxx  <- apply(cbind(BBB2$GLOBAL_INIT3,ss3),1,max)
FwPrp <- minn - maxx
###################

new_L <- BBB2$GLOBAL_INIT3
for(xx in 1:BBB2$popnum){

	new_L <- new_L + lgamma(samp_pop_freq[[xx]][1,] + theta[,xx] * samp_freq_new[1,]) - 
	         lfactorial (samp_pop_freq[[xx]][1,]) - lgamma(theta[,xx]*samp_freq_new[1,])+
	         lgamma(samp_pop_freq[[xx]][2,]+theta[,xx]*samp_freq_new[2,]) - 
                 lfactorial(samp_pop_freq[[xx]][2,]) - lgamma(theta[,xx]*samp_freq_new[2,])		

}

A     <- -old_L + new_L + log(BwPrp) - log(FwPrp) +
      (lambda-1)*(log(samp_freq_new[1,]) + log(samp_freq_new[2,]) - log(old_m) - log(old_n) )	
r     <- runif(BBB2$locnum)
check <- log(r) < A

#for(xx in which(check)){
# BBB2$freq_locus[xx,samp[1,xx]]    <- samp_freq_new[1,xx]
# BBB2$freq_locus[xx,samp[2,xx]]    <- samp_freq_new[2,xx]
# BBB2$acc_freq_ancestral[xx]       <- BBB2$acc_freq_ancestral[xx] + 1
#}

for(xx in which(check)){
 BBB2$freq_locus[xx,1]             <- samp_freq_new[1,xx]
 BBB2$freq_locus[xx,2]             <- samp_freq_new[2,xx]
 BBB2$acc_freq_ancestral[xx]       <- BBB2$acc_freq_ancestral[xx] + 1
}


}# End of Function



X_jump_model_codominant <- function(){

old_alpha 	   <- BBB2$d_alpha
excluded   	   <- which(!BBB2$alpha_included)

BBB2$d_alpha[BBB2$alpha_included]  <-  0 #rnorm(1,0,1)

#cat(1)

if(length(excluded)!=0){


 ### mod -------------------------------------------
  whatgroups    <- unique(BBB2$GROUP[excluded])
  group_id      <- match(whatgroups,BBB2$GROUP)
  e_group       <- rnorm(length(whatgroups),mean=BBB2$mean_alpha[group_id],sd=sqrt(BBB2$var_alpha[group_id]))
  #e_group      <- rnorm(length(whatgroups),mean=BBB2$d_alpha[group_id],1)
  #e_group      <- rlaplace(length(whatgroups),BBB2$mean_alpha[group_id],BBB2$var_alpha[group_id])

   .Call("update_group_Jump",	
	BBB2$d_alpha,
	as.integer(BBB2$GROUP),
	as.integer(whatgroups),
	e_group
   )  

}

BBB2$alpha_included     <- !BBB2$alpha_included


## old_theta
old_theta          <- outer(old_alpha,BBB2$d_beta,"+")
old_theta          <- exp(-(old_theta))
old1               <- lgamma(old_theta)-lgamma(BBB2$sample_size+old_theta)
 	
val  <- BBB2$GLOBAL_INIT1 

for(xx in 1:BBB2$popnum){
    xyz <- old_theta[,xx]*BBB2$freq_locus
    val <- val + lgamma(BBB2$freq_pop[[xx]]+xyz) - lgamma(xyz)
}

old1 <- rowSums(old1)
oldL <- old1 + rowSums(val,na.rm=TRUE)  
oldL <- tapply(oldL,BBB2$GROUP,sum) 

#+ tapply(oldL,BBB2$GROUP,min)
#- log(1-exp(tapply(oldL,BBB2$GROUP,max)))

  
  

## new_theta
new_theta          <- outer(BBB2$d_alpha,BBB2$d_beta,"+")
new_theta          <- exp(-(new_theta))
new1               <- lgamma(new_theta)-lgamma(BBB2$sample_size+new_theta)
 	
val  <- BBB2$GLOBAL_INIT1
for(xx in 1:BBB2$popnum){
    xyz <- new_theta[,xx]*BBB2$freq_locus
    val <- val + lgamma(BBB2$freq_pop[[xx]]+xyz) - lgamma(xyz)
}


new1 <- rowSums(new1)
newL <- new1 + rowSums(val,na.rm=TRUE) 
newL <- tapply(newL,BBB2$GROUP,sum) 

#+ tapply(newL,BBB2$GROUP,min) 
#- log(1-exp(tapply(newL,BBB2$GROUP,max)))


A_group <- -oldL + newL
A 	<- A_group[BBB2$GROUP]


included       		<-  BBB2$alpha_included #==TRUE
excluded       		<-  !included

A_included <- 0
A_excluded <- 0

# include in the selection model 
if(any(included)){
A_included     		<- - log(BBB2$tempting) 
			   - log(dnorm(BBB2$d_alpha[included],BBB2$mean_alpha[included],sqrt(BBB2$var_alpha[included]))) 
			   + log(dnorm(BBB2$d_alpha[included],0,1))
}

# exclude from the selection model 
if(any(excluded)){
A_excluded     		<-   log(BBB2$tempting) 
			   + log(dnorm(BBB2$d_alpha[excluded],BBB2$mean_alpha[excluded],sqrt(BBB2$var_alpha[excluded]))) 
			   - log(dnorm(BBB2$d_alpha[excluded],0,1))
}

A[included] 	           <-  A[included] + A_included 
A[excluded]                <-  A[excluded] + A_excluded


   A                       <-  tapply(A,BBB2$GROUP,unique)


write.table(A,"Jump_A")


   r                       <-  runif(BBB2$N.REGIONS)


check            	   <-  log(r) > A


   group_id       <- which(check)

   # tempting
   #  BBB2$acc.ratio <- BBB2$acc.ratio + length(which(!check))/length(A)
   #

   if(length(group_id)>0){

        .Call("update_group_Jump2",
	BBB2$d_alpha,
	old_alpha,
	as.integer(BBB2$GROUP),
	as.integer(group_id),
	BBB2$alpha_included
        )
  }

}# End of Function


X_update_d_alphaco_i_new <- function(){

	
##

	 sub_freq_locus    <- BBB2$freq_locus
	 sub_freq_pop      <- BBB2$freq_pop
         old_alpha         <- BBB2$d_alpha
         sub_sample_size   <- BBB2$sample_size
         sub_GROUP         <- BBB2$GROUP
         locnum2           <- length(BBB2$d_alpha)

##

 	old_theta     <- outer(old_alpha,BBB2$d_beta,"+")
 	old_theta     <- exp(-(old_theta))
 	old1          <- lgamma(old_theta)-lgamma(sub_sample_size + old_theta)
 	
  val <- matrix(0,locnum2,dim(BBB2$freq_locus)[2])
  for(xx in 1:BBB2$popnum){
    xyz <- old_theta[,xx]*sub_freq_locus
    val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
  }
  
  old1 <- rowSums(old1)
  oldL <- old1 + rowSums(val,na.rm=TRUE)
  oldL <- tapply(oldL,sub_GROUP,sum) 

  #+ tapply(oldL,sub_GROUP,min)
  #- log(1-exp(tapply(oldL,sub_GROUP,max))) 
  
  
    e_group   <- rnorm(BBB2$N.REGIONS,0,tapply(BBB2$var_prop_alpha,BBB2$GROUP,unique)) # BUG FIXED
   #e_group   <- rnorm(BBB2$N.REGIONS,0,1)

    new_alpha  <- old_alpha + e_group[BBB2$GROUP] 
    whatgroup  <- 1:BBB2$N.REGIONS


  
  new_theta <- outer(new_alpha,BBB2$d_beta,"+")
  new_theta <- exp(-(new_theta))
  new1      <- lgamma(new_theta)-lgamma(sub_sample_size+new_theta)
  
  val <- matrix(0,locnum2,dim(BBB2$freq_locus)[2])
 	for(xx in 1:BBB2$popnum){
          xyz <- new_theta[,xx]*sub_freq_locus
 	  val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
 	}
  
  new1 <- rowSums(new1)
  newL <- new1 + rowSums(val,na.rm=TRUE)
  newL <- tapply(newL,sub_GROUP,sum) 

  #+ tapply(newL,sub_GROUP,min)
  #- log(1-exp(tapply(newL,sub_GROUP,max)))
    
 
  old_alpha_group <- tapply(old_alpha,sub_GROUP,function(x){return(x[1])})
  new_alpha_group <- tapply(new_alpha,sub_GROUP,function(x){return(x[1])})
  

A         	     <-  -oldL + newL + (old_alpha_group^2-new_alpha_group^2)/(2*BBB2$sd_prior_alpha^2)
r                    <-  runif(length(whatgroup))

check    	     <-  log(r) < A

  	
  group_idX   <- which(check)
  group_id    <- whatgroup[group_idX]
  

  if(length(group_idX)>0){
   
   #CCode
   .Call("update_group2",
	as.integer(group_id),
	as.integer(BBB2$GROUP),
	BBB2$d_alpha,
	new_alpha,
	BBB2$acc_alpha
   )	

  }

}# End of function
X_update_d_alphaco_i_new_id <- function(id){

	
         sub_freq_locus    <- BBB2$freq_locus[id,]
	 sub_freq_pop      <- lapply(BBB2$freq_pop,function(yy){return(yy[id,])})
         old_alpha         <- BBB2$d_alpha[id]
         new_alpha         <- old_alpha
         sub_sample_size   <- BBB2$sample_size[id,]
         sub_GROUP         <- BBB2$GROUP[id]
	 locnum2           <- length(id)


 	old_theta     <- outer(old_alpha,BBB2$d_beta,"+")
 	old_theta     <- exp(-(old_theta))
 	old1          <- lgamma(old_theta)-lgamma(sub_sample_size + old_theta)
 	
  val <- matrix(0,locnum2,dim(BBB2$freq_locus)[2])
  for(xx in 1:BBB2$popnum){
    xyz <- old_theta[,xx]*sub_freq_locus
    val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
  }
  
  old1 <- rowSums(old1)
  oldL <- old1 + rowSums(val,na.rm=TRUE)
  oldL <- tapply(oldL,sub_GROUP,sum)  

    # create new alpha 
    group_id   <- unique(sub_GROUP) # what groups have to be updated
    xxx        <- BASIX.match(group_id,BBB2$GROUP)
    e_group    <- rnorm(length(group_id),0,sqrt(BBB2$var_prop_alpha[xxx])) 

    # update the whole groups
    .Call("update_groupJumpX",
	new_alpha,
        as.integer(sub_GROUP),
	as.integer(group_id),
	e_group
    )  

    
  #whatgroup  <- 1:BBB2$N.REGIONS

  new_theta <- outer(new_alpha,BBB2$d_beta,"+")
  new_theta <- exp(-(new_theta))
  new1      <- lgamma(new_theta)-lgamma(sub_sample_size+new_theta)
  
  val <- matrix(0,locnum2,dim(BBB2$freq_locus)[2])
 	for(xx in 1:BBB2$popnum){
          xyz <- new_theta[,xx]*sub_freq_locus
 	  val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
 	}
  
  new1 <- rowSums(new1)
  newL <- new1 + rowSums(val,na.rm=TRUE)
  newL <- tapply(newL,sub_GROUP,sum)  
 
  old_alpha_group <- tapply(old_alpha,sub_GROUP,function(x){return(x[1])})
  new_alpha_group <- tapply(new_alpha,sub_GROUP,function(x){return(x[1])})
  

A         	     <-  -oldL + newL + (old_alpha_group^2-new_alpha_group^2)/(2*BBB2$sd_prior_alpha^2)
r                    <-  runif(length(group_id))

check    	     <-  log(r) < A

  	
  group_idX   <- which(check)
  group_id    <- group_id[group_idX]
  

  if(length(group_idX)>0){
   
   #CCode
   .Call("update_group2",
	as.integer(group_id),
	as.integer(BBB2$GROUP),
	BBB2$d_alpha,
	new_alpha,
	BBB2$acc_alpha
   )	

  }

}# End of function update alpha new id

