sharing3secondcousins = function()
{
# Function designed specifically to compute the probability that all
# three second cousins receive a variant and that none of three second 
# cousins receive a variant introduced by two founders of their pedigree

# Value:
# numa : vector of probabilities that all three second cousins receive a variant introduced by each of the 7 different types of pairs of founders of their pedigree
# p0a : vector of probabilities that none of three second cousins receive a variant introduced by each of the 7 different types of pairs of founders of their pedigree
# npairs : the number of each type of founder pairs among the founders of the pedigree
Ds = 9
nfd=3

# 7 different cases to consider (see Solving kinship of founders for family 15157)
numa = p0a = numeric(7)

# case of founder couple 101 and 102
# Expression for probability that all receive variant from equation 18 of Bureau and Ruczinski
numa[1] = 1/2^(Ds-nfd)

# Vector of terms from equation 20 of Bureau and Ruczinski
vt = numeric(10)
# 0 child gets 2, 0 gets 1
vt[1] = 1/4^3
# 0 child gets 2, 1 gets 1
vt[2] = 3 * 3/4 * 1/(4*4*2)
# 0 child gets 2, 2 get 1
vt[3] = 3 * (3/4)^2 * 1/(4*2*2)
# 0 child gets 2, 3 get 1
vt[4] = (3/4)^3 * 1/2^3
# 1 child gets 2, 0 gets 1
vt[5] = 3 * 1/2 * 1/4^3
# 1 child gets 2, 1 gets 1
vt[6] = 6 * 1/2 * 3/4 * 1/(4*4*2)
# 1 child gets 2, 2 get 1
vt[7] = 3 * 1/2 * (3/4)^2 * 1/(4*2*2)
# 2 children get 2, 1 gets 0
vt[8] = 3 * 1/2^2 * 1/4^3
# 2 children get 2, 1 gets 1
vt[9] = 3 * 1/2^2 * 3/4 * 1/(4*4*2)
# 3 children get 2
vt[10] = 1/2^3 * 1/4^3

p0a[1] = sum(vt)

# case of founder 101 or 102 with founders 202, 205 or 207
# Expression for probability that all receive variant from equation 21 of Bureau and Ruczinski
  numa[2] = 1/2^7 *(1/2^2 + 1/2)
# Expression for probability that none receive variant from equation 23 of Bureau and Ruczinski
  p0a[2] = (1/2^3 + 3/4*1/2 + 1/2^2)*(1 - 1/2^3)^2

# case of founder 101 or 102 with founders 302, 307 or 309
# Expression for probability that all receive variant from equation 22 of Bureau and Ruczinski
  numa[3] = 1/2^(Ds+1) *(1 + 2^3)
# Expression for probability that none receive variant from equation 24 of Bureau and Ruczinski
  p0a[3] = 1/2*(1 - 1/2^3)^3

# case of pairs of founders among 202, 205 and 207
# numerator is 0
# Expression for probability that none receive variant from equation ?? of Bureau and Ruczinski
  p0a[4] = (1 - 1/2^2)^2

# case of pairs of founders 202-307, 202-309, 205-302, 205-309, 207-302 and 207-307
# numerator is 0
# Expression for probability that none receive variant from equation ?? of Bureau and Ruczinski
  p0a[5] = (1 - 1/2^2)*1/2

# case of pairs of founders among 202, 205 and 207
# numerator is 0
# Expression for probability that none receive variant from equation ?? of Bureau and Ruczinski
  p0a[6] = 1/2^2

# case of pairs of founders 202-302, 205-307 and 207-309
# numerator is 0
  p0a[7] = 1/4

npairs = c(1,6,6,3,6,3,3)

list(numa=numa,p0a=p0a,npairs=npairs)
}

sample.sharing3secondcousins = function()
  {
# Function designed specifically to sample the event that a variant is shared by all three second cousins
# and the event that a variant is seen in at least one of three second cousins, given that it was
# introduced only once in their pedigree

# Value:
# RVshared : takes value 1 (TRUE) if the variant is shared by all three second cousins, 0 (FALSE) otherwise
# RVseen : takes value 1 (TRUE) if the variant is seen in at least one of three second cousins, 0 (FALSE) otherwise

    # Sampling the degree of relatedness of the founder introducing the RV to final descendant
    Df = sample(rep(1:3,c(3,3,2)),1)
    # Sample whether the RV is shared or not
    # top founders with Df=3 have 1/2^9 prob of transmitting it to all descendents, others have prob 0
    RVshared = ifelse(Df==3,rbinom(1,1,1/2^9),FALSE)
    # If the variant is shared it is necessarily seen
    if (RVshared) RVseen = TRUE
    # Else we compute the probability that the variant was seen given that it was not shared and sample RVseen with that prob
    else
      {
      pseennotshared = ifelse(Df==3,(1-(1-1/2^Df)^Df - 1/2^9)/(1 - 1/2^9),1-(1-1/2^Df))
      RVseen = rbinom(1,1,pseennotshared)
      }
    list(RVshared=RVshared,RVseen=RVseen)
  }

approxsharing = function(nf,ord,phihat,sharingobj,fam.id,fam.dadid,fam.momid,nsim=1000)
{
# Main function to compute analytical and Monte Carlo approximations to the probability 
# a rare variant is shared by all sequenced members of a pedigree using Method 2 of Bureau et al.
# The analytical computation is written to be general, but the Monte Carlo approximation
# is specific to three second cousins

# Arguments:
# nf: number of founders of the pedigree
# ord : order of the polynomial approximation to the distribution of the number 
#       of distinct alleles in the founders (noted d in Bureau et al.) Must be <= 5.
# phihat : vector of the mean estimated kinship between founders
#          The computation is repeated with every value of phi in the vector phihat
# sharingobj: object returned by a function computing the probability that all
#             sequenced subject receive a variant and that none of the sequenced
#             subjects receive a variant introduced by two founders of their pedigree
#             (e.g. the function sharing3secondcousins)
# fam.id : vector of subject IDs
# fam.dadid : vector of father IDs. Founders' parents should be coded to 0
# fam.momid : vector of mother IDs. Founders' parents should be coded to 0
# nsim : number of Monte Carlo draws

# Value:
# papprox: vector of approximations of the sharing probability obtained with the analytical method for every value in phihat
# pMC: vector of approximations of the sharing probability obtained by Monte Carlo simulation for every value in phihat
# PFU.vec: vector of approximations of the probability that a founder is the unique founder to introduce the RV in the pedigree

nrep = length(phihat)
papprox = pMC = PFU.vec = rep(NA,nrep)

# Retrieving elements from sharing obj
numa = sharingobj$numa
p0a = sharingobj$p0a
npairs = sharingobj$npairs

# Vector of phia
phi.vec = compute.phi.vec(nf,2*nf-ord)

# Loop over the estimates of phi in the vector phihat 
for (r in 1:nrep)
  {
  # Solve the polynomial equation to obtain the estimate of the theta parameter of the distribution of the number of distinct alleles in the founders
  theta = max(infer.theta(phihat[r],phi.vec))
  # The approximation can be obtained only if a theta is strictly positive
  if (theta>0)
  {
  # Computation of the approximation of the probability that a founder is the unique founder to introduce the RV in the pedigree
  PFU.vec[r] = PFU = PFU.direct(nf,theta,ord)
  
  # Analytical approximation (for general pedigrees)
  
  numam = weighted.mean(numa,npairs)
  p0am2 = weighted.mean(p0a,npairs)
  
  p.approx2 = RVsharing.approx2(fam.id,fam.dadid,fam.momid,PFU=PFU)
  
  num.approx2 = numam*(1-nf*PFU) + p.approx2$num
  denom.approx2 = 1-p0am2*(1-nf*PFU) - p.approx2$p0

  papprox[r] = num.approx2/denom.approx2
  
  # Monte Carlo approximation (some steps are for general pedigrees, others are specific to three second cousins)
  
  # Computing distribution of the number of duplicated alleles (for general pedigrees)
  distri = c(1,theta,theta^2/2,theta^3/6,theta^4/24,theta^5/120)[1:(ord+1)]
  distri = distri/sum(distri)

  # Initialize counters of the number of Monte Carlo draws where the variant was shared and 
  # where the variant was seen among the three second cousins
  sumshared = sumseen = 0
  # Loop over the Monte Carlo draws
  for (m in 1:nsim)
  {
  # Sample number of duplicated alleles (for general pedigrees)
  nduplic = sample(0:ord,1,prob=distri)
  
  # Sample RV among the alleles (for general pedigrees)
  RV = sample.int(2*nf-nduplic,1)
  
  # The rest of the MC sampling steps are specific to three second cousins
  
  # If RV is among the duplicated alleles (the first nduplic alleles),
  if (RV <= nduplic)
    {
    # Sample the type of founder pair that introduces the RV 
    pairtype = sample(rep(1:length(npairs),npairs),1)

    # Sample whether the RV is shared or not
    RVshared = rbinom(1,1,numa[pairtype])
    # If RV is not shared, sample whether is was seen
    if (RVshared) RVseen = TRUE
    else RVseen = rbinom(1,1,(1 - p0a[pairtype] - numa[pairtype])/(1 - numa[pairtype]))
    }
  # Else the RV is introduced only once
  # We call sample.sharing3secondcousins to perfom the sampling
  else
    {
    obj = sample.sharing3secondcousins()
    RVshared = obj$RVshared
    RVseen = obj$RVseen
    }
  sumshared = sumshared + RVshared
  sumseen = sumseen + RVseen
  }
  # Record the MC estimate for the current value of phihat
  pMC[r] = sumshared/sumseen
  }
  }
list (papprox=papprox,pMC=pMC,PFU.vec=PFU.vec)
}


