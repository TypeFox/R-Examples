"HWE2sided" <-
function(geno,qplot=F,title=NULL)
   {

	#### small alteration by AL!
	obs_hets <- geno[1]
	obs_hom1 <- geno[2]
	obs_hom2 <- geno[3]


   if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
      return(-1.0)

   # total number of genotypes
   N <- obs_hom1 + obs_hom2 + obs_hets
   
   # rare homozygotes, common homozygotes
   obs_homr <- min(obs_hom1, obs_hom2)
   obs_homc <- max(obs_hom1, obs_hom2)

   # number of rare allele copies
   rare  <- obs_homr * 2 + obs_hets

   # Initialize probability array
   probs <- rep(0, 1 + rare)

   # Find midpoint of the distribution
   mid <- floor(rare * ( 2 * N - rare) / (2 * N))
   if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

   probs[mid + 1] <- 1.0
   mysum <- 1.0

   # Calculate probablities from midpoint down 
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr

   while ( curr_hets >=  2)
      {
      probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
      mysum <- mysum + probs[curr_hets - 1]

      # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
      curr_hets <- curr_hets - 2
      curr_homr <- curr_homr + 1
      curr_homc <- curr_homc + 1
      }    

   # Calculate probabilities from midpoint up
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr
   
   while ( curr_hets <= rare - 2)
      {
      probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
      mysum <- mysum + probs[curr_hets + 3]
         
      # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
      curr_hets <- curr_hets + 2
      curr_homr <- curr_homr - 1
      curr_homc <- curr_homc - 1
      }    

    # P-value calculation
    target <- probs[obs_hets + 1]

    #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

    #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

    # This assignment is the last statement in the fuction to ensure 
    # that it is used as the return value
    pval.H <- min(1.0, sum(probs[probs <= target])/ mysum)
    #print(pval.H)

####################################### added by AL 

probs <- probs/mysum
#print(probs)

if (rare %% 2 == 0) null_hets <- seq(0,rare,2)
else null_hets <- seq(1,rare,2)

probs <- probs[null_hets+1]
#print(probs)

mean_hets <- rare*(2*N-rare)/(2*N-1)

#print(mean_hets)
#print(null_hets)

weightleft <- sum(probs[null_hets<=mean_hets])
weightright <- sum(probs[null_hets>=mean_hets])

#print(weightleft)
#print(weightright)

if(obs_hets==mean_hets) pval.cond <- 1
else if(obs_hets<mean_hets) pval.cond <- sum(probs[null_hets<=obs_hets])/weightleft
else pval.cond <- sum(probs[null_hets>=obs_hets])/weightright

#print(pval.cond)

################# one-sided p-values for inbreeding

if (rare %% 2 == 0) pval.inbreed <- sum(probs[1:((obs_hets+2)/2)])
else pval.inbreed <- sum(probs[1:((obs_hets+1)/2)])

###########################################

if(qplot){
 plot(null_hets,probs,type="b",xlab="n_Aa",ylab="prob",main=title)
 #abline(v=obs_hets)
 abline(v=mean_hets)
}

return(list(pval.cond=pval.cond,pval.H=pval.H,pval.inbreed=pval.inbreed))

    }

