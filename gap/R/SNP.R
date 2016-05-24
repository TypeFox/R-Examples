# 9-3-2008 MRC-Epid JHZ

PARn <- function(p,RRlist)
{
   if (abs(sum(p)-1)>.Machine$double.eps) print("The frequencies do not sum up to one")
   1 - 1/(p%*%RRlist)
}

snp.ES <- function(beta,SE,N)
{
   if (SE<=0) stop("Incorrect input of SE(beta)")
   t <- beta/SE
   df <- N - 2
   t^2/(t^2+df)
}

snp.HWE <- function(g)
{
   obs_hom1 <- g[1]
   obs_hets <- g[2]
   obs_hom2 <- g[3]
   if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) return(-1.0)
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
   #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)
   #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)
   # P-value calculation
   target <- probs[obs_hets + 1]
   min(1.0, sum(probs[probs <= target])/ mysum)
}

snp.PAR <- function(RR,MAF,unit=2)
{
   RR2 <- RR^unit
   RRlist <- c(1,RR,RR2)
   q <- MAF
   p <-  1 - MAF
   pq <- c(p^2,2*p*q,q^2)
   PARn(pq,RRlist)
}
