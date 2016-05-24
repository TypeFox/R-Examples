X_update_d_alphaco_i_new3 <- function(){
	
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
  oldL <- tapply(oldL,BBB2$GROUP,sum)
  
  e          <- rnorm(locnum2,0,BBB2$var_prop_alpha)
 #e          <- rnorm(locnum2,0,1)	

  new_alpha  <- old_alpha + e
    
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
  newL <- tapply(newL,BBB2$GROUP,sum)
 
newL      <- newL[BBB2$GROUP]
oldL      <- oldL[BBB2$GROUP]
new_alpha <- tapply(new_alpha,BBB2$GROUP,mean)
new_alpha <- new_alpha[BBB2$GROUP]
  
A        	         <-  -oldL + newL + (old_alpha^2-new_alpha^2)/(2*BBB2$sd_prior_alpha^2)
r        	         <-  runif(length(unique(A)))
r                        <-  r[BBB2$GROUP]
  
check    	         <- log(r) < A
BBB2$d_alpha[check]      <- new_alpha[check]  
BBB2$acc_alpha[check]    <- BBB2$acc_alpha[check] + 1

}
# END of Function
X_jump_model_codominant3 <- function(){

old_alpha 	   <- BBB2$d_alpha
excluded   	   <- which(!BBB2$alpha_included)

BBB2$d_alpha[BBB2$alpha_included]  <-  0


if(length(excluded)!=0){
 # choose random SNP of each group excluded
 #groups <- BBB2$GROUP[excluded]
 #ids    <- tapply(excluded,groups,function(x){rand <- sample(1:length(x),1); return(x[rand]) })

 #e_group               <- rnorm(length(ids),mean=BBB2$mean_alpha[ids],sd=sqrt(BBB2$var_alpha[ids]))
 #BBB2$d_alpha[ids]     <- e_group 
 e_group                <- rnorm(length(excluded),mean=BBB2$mean_alpha[excluded],sd=sqrt(BBB2$var_alpha[excluded]))
 BBB2$d_alpha[excluded] <- e_group 

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


A_group              		<- -oldL + newL
A 	                        <- A_group[BBB2$GROUP]


included       		<-  BBB2$alpha_included #==TRUE
excluded       		<-  !included


A_included <- 0
A_excluded <- 0

if(any(included)){
A_included     		<- -log(BBB2$tempting) 
}

if(any(excluded)){
A_excluded     		<- log(BBB2$tempting)  
}

A[included] 	        <-  A[included] + A_included 
A[excluded]             <-  A[excluded] + A_excluded
A                       <-  tapply(A,BBB2$GROUP,unique)


write.table(A,"Jump_A")

   r                       <-  runif(BBB2$N.REGIONS)
   check            	   <-  log(r) > A
  
   # tempting
   # BBB2$acc.ratio <- BBB2$acc.ratio + length(which(!check))/length(A)
   #

   group_id     <- which(check)
   BBB2$d_alpha <- tapply(BBB2$d_alpha,BBB2$GROUP,mean)
   BBB2$d_alpha <- BBB2$d_alpha[BBB2$GROUP]

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

