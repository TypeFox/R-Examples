#source("math.R")

#constants
#e_f 		       <-  0.05
#sd_prior_beta          <-  1 #3
#mean_prior_beta        <- -1  #-5.15
#lambda 	       <-  1


update_d_alphaco_i <- function(id){ # id


 	if(id[1]==FALSE){ # Pilot runs
 	  id      <- 1:BBB$locnum
 	  locnum2 <- BBB$locnum
 	}else{                # Main Loop
 	  locnum2 <- length(id)
 	}
##

	 sub_freq_locus    <- BBB$freq_locus[id,]
	 sub_freq_pop      <- lapply(BBB$freq_pop,function(yy){return(yy[id,])})
         old_alpha         <- BBB$d_alpha[id]
         sub_sample_size   <- BBB$sample_size[id,]
    # mod #     sub_GROUP         <- GROUP[id]

##

 	old_theta     <- outer(old_alpha,BBB$d_beta,"+")
 	old_theta     <- exp(-(old_theta))
 	old1          <- lgamma(old_theta)-lgamma(sub_sample_size + old_theta)
 	
  val <- matrix(0,locnum2,dim(BBB$freq_locus)[2])
  for(xx in 1:BBB$popnum){
    xyz <- old_theta[,xx]*sub_freq_locus
    val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
  }
  
  old1 <- rowSums(old1)
  oldL <- old1 + rowSums(val,na.rm=TRUE)

  ## modification ------------------------------------  
  # oldL <- tapply(oldL,sub_GROUP,sum)
  ##--------------------------------------------------
  
    e             <- rnorm(locnum2,0,sqrt(BBB$var_prop_alpha))
    new_alpha     <- old_alpha + e
  
  ### modification -----------------------------------------------

 # if(!PILOT){

 #  whatgroup       <- unique(sub_GROUP)
 #  e_group         <- rnorm(length(whatgroup),0,sqrt(var_prop_alpha))

  #
 #  new_alpha <<- old_alpha
 #  jjj <<- 1
 #  sapply(whatgroup,function(x){
 #  what            <- x==sub_GROUP
 #  new_alpha[what] <<- old_alpha[what] + e_group[jjj]
 #  jjj <<- jjj + 1
 #  })
 #  new_alpha <- new_alpha

 #}else{

 #  e_group    <- rnorm(N.REGIONS,0,sqrt(var_prop_alpha))
 #  new_alpha  <- old_alpha + e_group[GROUP] 
 #  whatgroup  <- 1:N.REGIONS

 #}

  #
  ##--------------------------------------------------------------
  
  new_theta <- outer(new_alpha,BBB$d_beta,"+")
  new_theta <- exp(-(new_theta))
  new1      <- lgamma(new_theta)-lgamma(sub_sample_size+new_theta)
  
  val <- matrix(0,locnum2,dim(BBB$freq_locus)[2])
 	for(xx in 1:BBB$popnum){
          xyz <- new_theta[,xx]*sub_freq_locus
 	  val <- val + lgamma(sub_freq_pop[[xx]] + xyz) - lgamma(xyz)
 	}
  
  new1 <- rowSums(new1)
  newL <- new1 + rowSums(val,na.rm=TRUE)
  
  ## modification ----------------------------
  #newL            <- tapply(newL,sub_GROUP,sum)
   
  

  #old_alpha_group <- tapply(old_alpha,sub_GROUP,function(x){return(x[1])})
  #new_alpha_group <- tapply(new_alpha,sub_GROUP,function(x){return(x[1])})
  ##-----------------------------------------
 

 
A        	       <-  -oldL + newL + (old_alpha^2-new_alpha^2)/(2*BBB$sd_prior_alpha^2)

## mod ------------------------
#A         	       <-  -oldL + newL + (old_alpha_group^2-new_alpha_group^2)/(2*sd_prior_alpha^2)
##-----------------------------

 r        	      <-  runif(locnum2)
  
## mod ---------------------------  
#  r                    <-  runif(length(whatgroup))
## mod ---------------------------

check    	       <-  log(r) < A

 BBB$d_alpha[id[check]]      <- new_alpha[check]  
 BBB$acc_alpha[id[check]]    <- BBB$acc_alpha[id[check]] + 1


# modification -------------------------------------
  # update groups
  	
#  group_idX   <- which(check)
#  group_id    <- whatgroup[group_idX]
  
  
#  make this also in C

#  if(length(group_idX)>0){
   
#   mm         <- sapply(group_id,function(xx){xx==GROUP}) 
#   iii <<- 1
#   apply(mm,2,function(x){
#   d_alpha[x]   <<- new_alpha_group[group_idX[iii]]   # d_alpha[x] <<- new_alpha[x] geht auch !
  
#   acc_alpha[x] <<- acc_alpha[x] + 1	  
#   iii <<- iii + 1

#  })
#  }



# ---------------------------------------------------	

}# End of Function




update_d_betaco <- function(){

   
  # Init 
  #val <- rep(NA,2)
   val <- BBB$GLOBAL_INIT2
  # ---------------
 
  old_beta    <- BBB$d_beta 
  old_theta   <- outer(BBB$d_alpha,old_beta,"+")
  old_theta   <- exp(-(old_theta))
  
  old1 <- lgamma(old_theta)-lgamma(BBB$sample_size+old_theta)
  
  for(xx in 1:BBB$popnum){
      xyz     <- old_theta[,xx]*BBB$freq_locus
      val[xx] <- sum(lgamma(BBB$freq_pop[[xx]]+xyz)-lgamma(xyz),na.rm=TRUE)
  }

  old1 <- colSums(old1)
  oldL <- old1 + val
  #--------- calculate new likelihood
  e          <- rnorm(BBB$popnum,0,sqrt(BBB$var_prop_beta))
  BBB$d_beta <- BBB$d_beta + e
  
  new_theta <- outer(BBB$d_alpha,BBB$d_beta,"+")
  new_theta <- exp(-(new_theta))
  new1      <- lgamma(new_theta)-lgamma(BBB$sample_size+new_theta)
  for(xx in 1:BBB$popnum){
      xyz     <- new_theta[,xx]*BBB$freq_locus
      val[xx] <- sum(lgamma(BBB$freq_pop[[xx]]+xyz)-lgamma(xyz),na.rm=TRUE)
  }
  new1 <- colSums(new1)
  newL <- new1 + val

A                <-  -oldL + newL + (old_beta-BBB$d_beta)*(old_beta+BBB$d_beta-2*BBB$mean_prior_beta)/(2*BBB$sd_prior_beta^2)
r                <-  runif(BBB$popnum)
check            <-  log(r)>A
BBB$d_beta[check]        <- old_beta[check]
BBB$acc_beta[!check]     <- BBB$acc_beta[!check] + 1

 
}# End of Function


update_freq_codominant <- function(){


lambda    <- 1
# choose randomly two haplotype of each locus
samp      <-  sapply(BBB$hapcount,function(x){return(sample(1:x,size=2)) })
theta     <-  outer(BBB$d_alpha,BBB$d_beta,"+")
theta     <-  exp(-theta)
# ok

BBB$iii       <- 1
# frequencies of haplotypes of the two random haps
samp_freq <- apply(samp,2,function(x){
	      r   <- BBB$freq_locus[BBB$iii,x]
	      BBB$iii <- BBB$iii+1
	      return(r)
           })
#ok

# haplotype distribution for each pop of the two random haps
BBB$iii           <- 1
samp_pop_freq <- vector("list",BBB$popnum)
for(yy in 1:BBB$popnum){

	samp_pop  <- apply(samp,2,function(x){
		      r       <- BBB$freq_pop[[yy]][BBB$iii,x]
		      BBB$iii <- BBB$iii+1
		      return(r)
        	   })
	samp_pop_freq[[yy]] <- samp_pop
	BBB$iii <- 1
}
#ok



old_L <- BBB$GLOBAL_INIT3
# 1 <- m
# 2 <- n
for(xx in 1:BBB$popnum){

	old_L <- old_L + lgamma(samp_pop_freq[[xx]][1,] + theta[,xx] * samp_freq[1,]) - 
	         lfactorial (samp_pop_freq[[xx]][1,]) - lgamma(theta[,xx]*samp_freq[1,])+
	         lgamma(samp_pop_freq[[xx]][2,]+theta[,xx]*samp_freq[2,]) - 
                 lfactorial(samp_pop_freq[[xx]][2,]) - lgamma(theta[,xx]*samp_freq[2,])		

}

old_m <- samp_freq[1,]
old_n <- samp_freq[2,]


#### BwPrp
old_m_n 	 <- old_m + old_n
old_m_ancestral  <- old_m + BBB$e_ancestral
old_m_ancestral2 <- old_m - BBB$e_ancestral

checkmin <- cbind(old_m_n,old_m_ancestral)
checkmax <- cbind(BBB$GLOBAL_INIT3,old_m_ancestral2)
minn     <- apply(checkmin,1,min)
maxx     <- apply(checkmax,1,max)
BwPrp    <- minn - maxx
###################


# veraender Frequenzen
u     <- runif(BBB$locnum,0,BwPrp) + maxx 
# samp_freq_new
samp_freq_new     <- samp_freq
samp_freq_new[1,] <- u
samp_freq_new[2,] <- old_n + old_m - u

### FwPrp
ss1   <- colSums(samp_freq_new)
ss2   <- samp_freq_new[1,] + BBB$e_ancestral
ss3   <- samp_freq_new[1,] - BBB$e_ancestral
minn  <- apply(cbind(ss1,ss2),1,min)
maxx  <- apply(cbind(BBB$GLOBAL_INIT3,ss3),1,max)
FwPrp <- minn - maxx
###################

new_L <- BBB$GLOBAL_INIT3
for(xx in 1:BBB$popnum){

	new_L <- new_L + lgamma(samp_pop_freq[[xx]][1,] + theta[,xx] * samp_freq_new[1,]) - 
	         lfactorial (samp_pop_freq[[xx]][1,]) - lgamma(theta[,xx]*samp_freq_new[1,])+
	         lgamma(samp_pop_freq[[xx]][2,]+theta[,xx]*samp_freq_new[2,]) - 
                 lfactorial(samp_pop_freq[[xx]][2,]) - lgamma(theta[,xx]*samp_freq_new[2,])		

}

A     <- -old_L + new_L + log(BwPrp) - log(FwPrp) +
      (lambda-1)*(log(samp_freq_new[1,]) + log(samp_freq_new[2,]) - log(old_m) - log(old_n) )	
r     <- runif(BBB$locnum)
check <- log(r) < A

for(xx in which(check)){
 BBB$freq_locus[xx,samp[1,xx]]    <- samp_freq_new[1,xx]
 BBB$freq_locus[xx,samp[2,xx]]    <- samp_freq_new[2,xx]
 BBB$acc_freq_ancestral[xx]       <- BBB$acc_freq_ancestral[xx] + 1
}


}# End of Function



jump_model_codominant <- function(){



old_alpha 	   <- BBB$d_alpha
excluded   	   <- which(!BBB$alpha_included)

BBB$d_alpha[BBB$alpha_included]  <- 0

if(length(excluded)!=0){

 ### mod -------------------------------------------
 # whatgroups    <- unique(GROUP[excluded])
 # group_id      <- match(whatgroups,GROUP)
 # e_group       <- rnorm(length(whatgroups),mean=mean_alpha[group_id],sd=sqrt(var_alpha[group_id]))
 #


 # jjj <<- 1
 # sapply(whatgroups,function(x){
 # what            <- x==GROUP
 # d_alpha[what]  <<- old_alpha[what] + e_group[jjj]
 # jjj <<- jjj + 1
 # })

  #
### ------------------------------------------------


 BBB$d_alpha[excluded] <- rnorm(length(excluded),mean=BBB$mean_alpha[excluded],sd=sqrt(BBB$var_alpha[excluded])) 

}

BBB$alpha_included     <- !BBB$alpha_included


## old_theta
old_theta          <- outer(old_alpha,BBB$d_beta,"+")
old_theta          <- exp(-(old_theta))
old1               <- lgamma(old_theta)-lgamma(BBB$sample_size+old_theta)
 	
#val <- matrix(0,locnum,dim(freq_locus)[2])
val  <- BBB$GLOBAL_INIT1 
for(xx in 1:BBB$popnum){
    xyz <- old_theta[,xx]*BBB$freq_locus
    val <- val + lgamma(BBB$freq_pop[[xx]]+xyz) - lgamma(xyz)
}

old1 <- rowSums(old1)
oldL <- old1 + rowSums(val,na.rm=TRUE)  

#  mod
#   oldL <- tapply(oldL,GROUP,sum)
#--------------------------------------------------------------

## new_theta
new_theta          <- outer(BBB$d_alpha,BBB$d_beta,"+")
new_theta          <- exp(-(new_theta))
new1               <- lgamma(new_theta)-lgamma(BBB$sample_size+new_theta)
 	
# val <- matrix(0,locnum,dim(freq_locus)[2])
val  <- BBB$GLOBAL_INIT1
for(xx in 1:BBB$popnum){
    xyz <- new_theta[,xx]*BBB$freq_locus
    val <- val + lgamma(BBB$freq_pop[[xx]]+xyz) - lgamma(xyz)
}

new1 <- rowSums(new1)
newL <- new1 + rowSums(val,na.rm=TRUE) 

# mod
#  newL <- tapply(newL,GROUP,sum)
#------------------------------------------------------------


 A                     <-  -oldL + newL

### mod --------------------------------
#A_group              		<-  -oldL + newL
### ------------------------------------

# mod -------------------------------
#  A   <<- numeric(locnum)
#  jjj <<- 1
#  sapply(1:N.REGIONS,function(x){
#  what        <- x==GROUP
#  A[what]     <<- A_group[jjj]
#  jjj         <<- jjj + 1
#  })
#A <- A
#------------------------------------

included       		<-  BBB$alpha_included==TRUE
excluded       		<-  !included


A_included <- 0
A_excluded <- 0

if(any(included)){
A_included     		<-  log_prior_alpha(BBB$d_alpha[included]) - ( -0.5*log(2*pi*BBB$var_alpha[included])-((BBB$d_alpha[included]-BBB$mean_alpha[included])^2)/(2*BBB$var_alpha[included])) - log(10)
}
 #log_prior_alpha(alpha[i]) -(-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))-log(prior_odds);

if(any(excluded)){
A_excluded     		<-  (-0.5*log(2*pi*BBB$var_alpha[excluded]) - ((BBB$d_alpha[excluded]-BBB$mean_alpha[excluded])^2)/(2*BBB$var_alpha[excluded])) - log_prior_alpha(BBB$d_alpha[excluded]) + log(10)
}

#(-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))
#            -log_prior_alpha(alpha[i])+log(prior_odds);


A[included] 	          <-  A[included] + A_included 
A[excluded]               <-  A[excluded] + A_excluded

 r                       <-  runif(BBB$locnum)


### mod ---------------------------------------------
#  A                       <-  tapply(A,GROUP,unique)
#  r                       <-  runif(N.REGIONS)
### mod ---------------------------------------------  
  
  
check            	        <-  log(r) > A

# modification -------------------------------------

  # update groups	
#  group_id    <- which(check)

#  if(length(group_id)>0){
   
#   mm         <- sapply(group_id,function(xx){xx==GROUP}) 
#   iii <<- 1
#   apply(mm,2,function(x){
#   d_alpha[x]         <<- old_alpha[x]    
#   alpha_included[x]  <<- !alpha_included[x]
#   iii <<- iii + 1
#  })
#  }
  # nb_alpha_included fehlt noch
  
# ---------------------------------------------------	




BBB$d_alpha[check]         	<- old_alpha[check]
BBB$alpha_included[check]       <- !BBB$alpha_included[check]

#nb_alpha_included[!check &  alpha_included]  <<- nb_alpha_included[!check &  alpha_included] + 1
#nb_alpha_included[!check & !alpha_included]  <<- nb_alpha_included[!check & !alpha_included] - 1



}# End of Function

