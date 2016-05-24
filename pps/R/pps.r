#########################################################################
# Copyright (C) 2003 Jack G. Gambino
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.	You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA 02111-1307  USA.
#########################################################################

#########################################################################
#
# stratsrs: stratified simple random sampling function
#
# Inputs:
#         - a vector of stratum indicators for the population units
#           e.g., stratum <- c(1,1,1,1,2,2,2,3,3,3,3,3,8,8,8,8)
#         - a vector of sample sizes (one per stratum)
#           e.g., nh <- c(2,1,3,2) (i.e., here there are four strata)
# It is assumed that the stratum vector is sorted (or at least grouped).
#
# Returns: indices of units in sample
#
# Usage:
# 	samplesubframe <- pop[stratsrs(pop$stratum,nh),]
#
# Here, pop is a frame that includes a variable stratum. This example
# takes pop and returns a stratified simple random sample of units (rows)
# from pop. The vector nh is the stratum sample sizes in the same order
# as the strata in pop.
#
#########################################################################

stratsrs <- function(stratum,nh)
{

   H <- length(nh)				# number of strata

   # H should equal length(unique(stratum))

   Nh <- as.vector(table(stratum))		# compute stratum sizes

   Nhcum <- cumsum(Nh)				# cumulate stratum sizes

   strsample <- sample(Nh[1],nh[1])		# sample in stratum 1

   if (H>1)
	   for (h in 2:H)			# samples in strata 2 to H
	   {
	      strsample <- c(strsample, sample(Nh[h],nh[h])+Nhcum[h-1])
	   }
   strsample
}

############################################################
#
# ppswr: select n units out of N with probability 
# proportional to the size of the units; 
# units are selected with replacement;
# the sizes are input as a vector
#
# Usage: sampleindeces <- ppswr(sizes,n)
#
# Returns: the indeces of the selected units
#
############################################################

ppswr <- function(sizes,n)
{

   N <- length(sizes)		# number of units in the population
   cumsizes <- cumsum(sizes)
   totsize <- cumsizes[N]	# the sum of all N sizes
   s <- numeric(n)		# initialize the sample indeces vector

   for (u in 1:n) {		# determine the uth sample unit
      r <- runif(1,0,totsize)
      i <- 1
      while (cumsizes[i] < r) {i <- i+1}
      s[u] <- i
   }

   s

}

############################################################
#
# pps1: select one unit out of N with probability 
# proportional to the size of the units; 
# the sizes are input as a vector
#
# Usage: sampleindex <- pps1(sizes)
#
# Returns: the index of the selected unit
#
############################################################

pps1 <- function(sizes)
{

   N <- length(sizes)		# number of units in the population
   cumsizes <- cumsum(sizes)
   totsize <- cumsizes[N]	# the sum of all N sizes
   r <- runif(1,0,totsize)
   
   i <- 1
   while (cumsizes[i] < r) {i <- i+1}

   i

}

############################################################
#
# ppss: use pps systematic sampling to select a sample of 
# n units out of N with probability proportional to the
# size of each unit; the sizes are input as a vector
#
# Usage: sampleindeces <- ppss(sizes,n)
#
# Returns: a vector corresponding to the indeces of the 
#          units that fall in the sample
#
############################################################

ppss <- function(sizes,n)
{

   N <- length(sizes)		# number of units in the population
   cumsizes <- cumsum(sizes)
   totsize <- cumsizes[N]	# the sum of all N sizes
   int <- totsize/n		# the sampling interval
   r <- runif(1,0,int)		# the random starting point
   s <- numeric(n)		# initialize sample indeces
   
   i <- 1
   for (j in 1:n)		# determine the jth sample unit
   {
      u <- r + (j-1)*int
      while (cumsizes[i] < u) {i <- i+1}
      s[j] <- i
   }
   s

}

#############################################################################
# ppssstrat: stratified pps systematic sampling
# In each stratum, select a sample using pps systematic sampling.
#
# Requires: ppss(), found in ppss.r
#
# Usage: sampleindeces <- ppssstrat(unitsizes,stratumcodes,samplesizes)
# unitsizes[i] is the size of the ith unit (assumed to be sorted by stratum),
# stratumcodes[i] is the stratum unit i belongs to,
# samplesizes[h] is the number of units to be selected from stratum h
#
# Returns: a vector corresponding to the indeces of the units that
#          fall in the sample
#
# Example:
#	sizes <- c(1:5,10:6)*10 (i.e., 10,20,30,40,50,100,90,80,70,60)
#	strat <- c(1,1,1,2,2,3,3,3,3,3) 
#	n <- c(2,1,3) (i.e., select 2 units in stratum 1, 1 in 2, 3 in 3)
#	print(ppssstrat(sizes,strat,n))
#############################################################################

ppssstrat <- function(sizes,stratum,n)
{
	H <- length(n)			# number of strata; should equal
					#    length(unique(stratum))
	
	Nh <- as.vector(table(stratum))	# compute stratum sizes
	
	s <- ppss(sizes[1:Nh[1]], n[1])	# sample indeces for first stratum
	stratstart <- 1+Nh[1]		# start of stratum 2
	if (H>1) 
	   for (h in 2:H) {
	
		stratend <- stratstart + Nh[h]-1   # end of stratum h
		
		# sample indeces for stratum h :
		tmp <- ppss(sizes[stratstart:stratend], n[h]) + stratstart-1	
		
		s <- c(s,tmp)
		
		stratstart <- stratend + 1	# start of stratum h+1
		
	   }
	s
}

#########################################################################
# sizesok: check that unit sizes are not too big 
# Returns: The number of "bad" units (i.e., that are too big)
# In pps systematic sampling such units are selected with 
# certainty. In Sampford''s method, such units have negative
# adjusted sizes/probabilities.
#########################################################################

sizesok <- function(size,n)
{
	N <- length(size)
	totsize <- sum(size)
	limit <- totsize/n
	toobig <- 0
	for (i in 1:N) {
		if (size[i] >= limit) {
			toobig <- toobig + 1
			cat("Unit",i,"is too big.\n")
		}
	}
	toobig	# the number of bad units
}

#########################################################################
# sampford: Sampford''s method for selecting n units out of N with
# probability proportional to size. The method fails if
# a unit has a size greater than or equal to totsize/n
#
# Usage: sampleindeces <- sampford(size,n)
#
# Returns: the indeces of the selected units
#########################################################################

sampford <- function(size,n)
{

	toobig <- sizesok(size,n)		    # check if sizes are ok
	if (toobig == 0) {
		totsize <- sum(size)
		adjsize <- size/(totsize-n*size)    # adjusted sizes
		s <- numeric(n)

		while ( length(unique(s)) < n ) {   # keep selecting samples
						    # until one with n distinct
						    # units appears
	
			s <- ppswr(adjsize,n-1)	    # select n-1 units using the
						    #   adjusted sizes
			s[n] <- pps1(size)	    # select 1 unit using the
						    #   original sizes

		}
	
		s
	}
	else cat("Some units are too big. Aborting.\n")
}

################################################################
# sampfordpi: compute pi_i and pi_ij values of PPS sampling with 
# Sampford''s method
# See sampford.r to actually select a sample
#
# Usage: jointprobs <- sampfordpi(sizes,n)
#   where sizes is a vector of sizes of the population units and
#   n is the desired sampli size
#
# Returns: a matrix with the inclusion probability pi(i) for
#   for each unit i in the population and with the joint
#   inclusion probability pi(i,j) of units i and j in 
#   position (i,j) in the matrix, where i and j are not equal
#
# Test data: from Sampford article (Biometrika 1967)
#   sizes <- c(18,14,13,11,10,10,8,7,5,4)
#   n <- 5
################################################################
#
# In the following, L[1], L[2], etc. correspond to L_0, L_1, etc.
# in Sampford''s paper. Similarly, LL[1,i,j], LL[2,i,j], etc.
# correspond to L_0(i,j), L_1(i,j), etc.
#
################################################################

sampfordpi <-function(sizes,n)
{
   toobig <- sizesok(sizes,n)			# check if sizes are ok
   if (toobig == 0 & n>1) {			# do nothing if n=1
   N <- length(sizes)
   N1 <- N+1; n1 <- n+1
   sizetot <- sum(sizes)
   p <- sizes/sizetot
   lambda <- p/(1-n*p)
   
   L <- numeric(n); R <- numeric (n); LM <- numeric(n); KI <- numeric(n)
   LL <- array(dim=c(n-1,N,N)); f <- numeric(n); pi2 <- array(dim=c(N,N))
   
   L[1] <- 1.
   L[2] <- sum(lambda)
   R[2] <- L[2]

   if (n>2) {			# compute L[3], ..., L[n]
   for (m in 2:(n-1)) {
      R[m+1] <- sum(lambda^m)
      for (r in 1:m) LM[r] <- (-1)^(r-1)*R[r+1]*L[m-r+1]
      L[m+1] <- sum(LM)/m
   } }
   for (t in 1:n) KI[t] <- t*L[n-t+1]/n^t
   K <- 1/sum(KI)
   
   for (i in 1:N) {
      for (j in i:N) {
         LL[1,i,j] <- 1.
         if (n>2) { LL[2,i,j] <- L[2] - lambda[i] - lambda[j] }
         if (n>3) {		# compute LL[3, i,j], ..., LL[n-1,i,j]
	 for (m in 2:(n-2)) {
          LL[m+1,i,j] <- 
          L[m+1]-(lambda[i]+lambda[j])*LL[m,i,j]-lambda[i]*lambda[j]*LL[m-1,i,j]
         } }
         for (t in 2:n) f[t] <- ( t - n*(p[i]+p[j]) ) * LL[n-t+1,i,j]/n^(t-2)
         pi2[i,j] <- K * lambda[i] * lambda[j] * sum(f) # joint inclusion prob.
	 pi2[j,i] <- pi2[i,j]	# symmetric matrix
      }
      pi2[i,i] <- n*p[i]	# put inclusion probabilities p_i on diagonal
   }
   pi2
   }
   else if (n>1) cat("Some units are too big. Aborting.\n")
}

#########################################################################
# stratumsizes: given a vector of sorted stratum indicators, returns the 
# number of units in each stratum.
# Example: stratumsizes(c(1,1,1,1,2,2,2,3,3,3,3,3,8,8,8,8))
#	   returns 4 3 5 4
#########################################################################

stratumsizes <- function(stratum)
{

	Nh <- as.vector(table(stratum))	# compute stratum sizes
	Nh

}

# Alternative method: Nh <- as.data.frame(table(stratum))$Freq

#########################################################################
# permuteinstrata: randomize order of elements *within* each stratum
# Example: there are three strata with 9, 10 and 10 elements, respectively;
#	   thus original indeces are 1 2 ... 9 10 ... 19 20 ... 29;
#	   then permuteinstrata(c(9,10,10)) returns the permuted indeces
#	   within each stratum (e.g., 3 1 ... 7 16 12 ... 18 28 21 ... 25)
#########################################################################

permuteinstrata <- function(stratsizes)
{
	H <- length(stratsizes)		# number of strata
	
	neworder <- sample(1:stratsizes[1])	# reorder first stratum

	stratend <- stratsizes[1]

	if (H>1) 
		for (h in 2:H) {
		   neworder <- c(neworder, sample(1:stratsizes[h])+stratend)
		   stratend <- stratend + stratsizes[h]
		}

	neworder
}

