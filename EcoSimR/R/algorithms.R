

#' Vector Sample Function
#' @description Takes an input binary vector and a weight vector.
#' Reassigns 1s randomly in proportion to vector weights.
#' @details This function takes an input vector of binary presence-absence
#' values and a vector of non-negative probability weights. Both vectors
#' must be of identical length. 
#' @param speciesData binary vector representing species presences and 
#' absences.
#' @param weights a vector of non-negative read numbers representing
#' probabilistic weights for species occurrencs of the same length 
#' as speciesData.
#' @return Returns a re-ordered binary vector in which the occurrences
#' are placed in cells with probabilities proportional to values given 
#' in weights.
#' @references Gotelli, N.J., G.R. Graves, and C. Rahbek. 2010. Macroecological 
#' signals of species interactions in the Danish avifauna. Proceedings of the 
#' National Academy of Sciences, U.S.A. 107: 530-535.
#' @note Several of the randomization algorithms use this function
#' to assign species occurrences with probabilities that reflect
#' species or site differences. It is an effective method for 
#' conditioning the marginal probabilities of a null matrix on
#' independent measurements of site or species characteristics.
#' @seealso \code{\link{sim10}} randomization algorithm.
#' @examples 
#' myColonizer <- vector_sample(speciesData=rbinom(10,1,0.5),weights=runif(10))

#' @export

vector_sample <- function(speciesData,weights) 

{
x <- mat.or.vec(length(speciesData),1)
x[sample(1:length(speciesData),size=sum(speciesData),prob=weights)] <- 1
return(x)
}


#' Sim1 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling 
#' all of its elements equiprobably.
#' @details This algorithm assumes species and sites are equiprobable.
#' It preserves the total matrix fill, but places no other constraints on 
#' row or column totals.
#' @param speciesData binary presence-absence matrix
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and percent fill as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. Ecology 81: 2606-2621. 
#' @note This is the simplest of all randomization algorithms for a presence-
#' absence matrix. However, it assumes that both species and sites are
#' equiprobable, and has poor Type I error frequencies when tested with 
#' purely random matrices. If the input matrix is sparse, it will often
#' generate null matrices with empty rows or columns. Not recommended for 
#' co-occurrence analysis.
#' @examples 
#' randomMatrix <- sim1(speciesData=matrix(rbinom(40,1,0.5),nrow=8))
#' @export

sim1 <- function(speciesData) 

{
matrix(sample(speciesData), ncol=ncol(speciesData))
}



#' Sim2 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling elements 
#' within each row equiprobably.
#' @details This algorithm assumes sites are equiprobable, but preserves 
#' differences among species (= row sums).	
#' @param speciesData binary presence-absence matrix
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and rowsums as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. ecology 81: 2606-2621. 
#' @note This algorithm preserves differences in the commonness and rarity of 
#' species (= rowsums), but assumes that all sites are equiprobable. It would 
#' not be appropriate for islands that vary greatly in area, but it would be 
#' appropriate for quadrat censuses in a relatively homogeneous environment. 
#' sim2 can sometimes generate matrices with empty columns, but this is 
#' unlikely unless the matrix is very sparse. sim2 has good Type I error 
#' frequenceis when tested against random matrices. However, if sites do vary 
#' in their suitability or habitat quality, it will often identify aggregated 
#' patterns of species co-occurrence. sim2 and sim9 have the best overall 
#' performance for species co-occurrence analyses. However, because they differ 
#' in their assumptions about site quality, they often differ in their results, 
#' with sim9 often detecting random or segregated patterns for matrices in 
#' which sim2 detects aggregated patterns.

#' @seealso \code{\link{sim9}} co-occurrence algorithm.
#' @examples 
#' randomMatrix <- sim2(speciesData=matrix(rbinom(40,1,0.5),nrow=8))
#' @export

sim2 <- function(speciesData) 

{
t(apply(speciesData,1,sample))
}


#' Sim3 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling elements 
#' within each column equiprobably. 
#' @details This algorithm assumes species are equiprobable, but preserves 
#' differences among sites (= column sums).	
#' @param speciesData binary presence-absence matrix
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and colsums as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. Ecology 81: 2606-2621. 
#' @note This algorithm preserves differences in species richness among sites 
#' (= colsums), but assumes that all species are equiprobable. This assumption 
#' is usually unrealistic, and sim3 has a high frequency of Type I errors with 
#' random matrices, so it is not recommended for co-occurrence analysis.

#' @examples 
#' randomMatrix <- sim3(speciesData=matrix(rbinom(40,1,0.5),nrow=8))
#' @export

sim3 <- function(speciesData) 

{
apply(speciesData,2,sample)
}



#' Sim4 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling elements 
#' within each row.  Sampling weights for each column are proportional to 
#' column sums. Makes a call to the vector_sample function.
#' @details This algorithm preserves differences among species in occurrence
#' frequencies, but assumes differences among sites in suitability are 
#' proportional to observed species richness (= column sums).
#' @param speciesData binary presence-absence matrix 
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and rowsums as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. Ecology 81: 2606-2621. 
#' @note This algorithm preserves differences in the commonness and rarity of 
#' species (= rowsums), but assumes differences among sites in suitability are 
#' proportional to observed species richness (= column sums). sim4 has a 
#' somewhat high frequency of Type I errors with random matrices, so it is not 
#' recommended for co-occurrence analysis.

#' @examples 
#' randomMatrix <- sim4(speciesData = matrix(rbinom(40,1,0.5),nrow=8))
#' @export

sim4 <- function(speciesData) 

{
t(apply(speciesData,1,vector_sample,weights=colSums(speciesData)))
}



#' Sim5 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling elements 
#' within each column. Sampling weights for each row are proportional to 
#' row sums. Makes a call to the vector_sample function.
#' @details This algorithm preserves differences among sites in species 
#' richness, but assumes differences among species in commonness and rarity 
#' are proportional to observed species occurrences (= row sums).
#' @param speciesData binary presence-absence matrix
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and colsums as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. Ecology 81: 2606-2621. 
#' @note This algorithm preserves differences among sites in species richness
#' (= colsums), but assumes differences among species in commonness and rarity 
#' are proportional to observed species occurrences (= rowsums). sim5 has a 
#' high frequency of Type I errors with random matrices, so it is not 
#' recommended for co-occurrence analysis.

#' @examples 
#' randomMatrix <- sim5(speciesData = matrix(rbinom(40,1,0.5),nrow=8))
#' @export

sim5 <- function(speciesData) 

{
apply(speciesData,2,vector_sample,weights=rowSums(speciesData))
}


#' Sim6 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling all 
#' elements. Rows are equiprobable, and columns are proportional to column sums.
#' Makes a call to the vector_sample function.
#' @details This algorithm assumes that species are equiprobable, but that 
#' differences in suitability among sites are proportional to observed species
#' richness (=colsums).
#' @param speciesData binary presence-absence matrix 
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and fill as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. Ecology 81: 2606-2621. 
#' @note This algorithm assumes that species are equiprobable, and that
#' differences among sites are proportional to observed species richness
#' (=colsums). sim6 has a high frequency of Type I errors with random matrices, 
#' so it is not recommended for co-occurrence analysis.

#' @examples 
#' randomMatrix <- sim6(speciesData = matrix(rbinom(40,1,0.5),nrow=8))
#' @export

sim6 <- function(speciesData) 

{
matrixWeights <- outer(rep(1,nrow(speciesData)),colSums(speciesData))
out <-matrix(vector_sample(speciesData, weights=matrixWeights),ncol=ncol(speciesData))
}


#' Sim7 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling all 
#' elements. Columns are equiprobable, and rows are proportional to row sums.
#' Makes a call to the vector_sample function.
#' @details This algorithm assumes that sites are equiprobable, but that 
#' differences in frequency of occurrence among species are proportional to 
#' observed species richness (=colsums).
#' @param speciesData binary presence-absence matrix
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and fill as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. Ecology 81: 2606-2621. 
#' @note This algorithm assumes that species are equiprobable, and that
#' differences among sites are proportional to observed species richness
#' (=colsums). sim7 has a high frequency of Type I errors with random matrices, 
#' so it is not recommended for co-occurrence analysis.

#' @examples 
#' randomMatrix <- sim7(speciesData = matrix(rbinom(40,1,0.5),nrow=8))
#' @export


sim7 <- function(speciesData) 

{
matrixWeights <- outer(rowSums(speciesData),rep(1,ncol(speciesData)))
matrix(vector_sample(speciesData, weights=matrixWeights),ncol=ncol(speciesData))
}



#' Sim8 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling all 
#' elements. Columns are proportional to column sums, and rows are proportional 
#' to row sums. Makes a call to the vector_sample function.
#' @details This algorithm assumes that the probability that a species occurs
#' in a site is depends on the joint independent probability of randomly 
#' selecting the species and randomly selecting the site, with these
#' probabilities set proportional to row and column sums of the matrix. 

#' @param speciesData binary presence-absence matrix 
#' (rows = species, columns = sites).
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and fill as the input matrix.
#' @references Gotelli, N.J. 2000. Null model analysis of species co-occurrence 
#' patterns. Ecology 81: 2606-2621. 
#'
#' Ulrich, W. and N.J. Gotelli. 2012. A null model algorithm for presence-
#' absence matrices based on proportional resampling. Ecological Modelling 
#' 244:20-27.
#' @note This algorithm is theoretically attractive because it incorporates 
#' heterogeneity in species occurrences and species richness per site in a 
#' probabilistic way that does not fix row and column frquencies. However,
#' in spite of its appeal, sim8 does not generate average row and column sums 
#' that match the original matrix, and it is susceptible to Type I errors when 
#' tested with random matrices. It is not recommended for co-occurrence 
#' analysis. See Ulrich and Gotelli (2012) for a more complicated algorithm 
#' for probabilistic row and column totals that has better statistical behavior. 

#' @examples 
#' randomMatrix <- sim8(speciesData = matrix(rbinom(40,1,0.5),nrow=8))

#' @export


sim8 <- function(speciesData) 

{
matrixWeights <- outer(rowSums(speciesData),colSums(speciesData))
matrix(vector_sample(speciesData,weights=matrixWeights),ncol=ncol(speciesData))
}



#' Sim10 Co-occurrence Randomization Algorithm
#' @description Randomizes a binary matrix speciesData by reshuffling all 
#' elements. Rows and column probabilities are proportional to user-supplied 
#' row and column weights, which define relative suitability probabilities for
#' species and sites. Makes a call to the vector_sample function.

#' @details This function incorporates vectors of weights for species and/or 
#' sites to condition the simulation. These two vectors are used as outer 
#' products to set cell probabilities for the entire matrix. Thus:  

#' \deqn{p(cell_{ij})=p(row_i)p(col_j)}{p(cell_ij)=p(row_i)p(col_j)}
#' Weights must be positive real numbers. The algorithm will scale them so they
#' sum to 1.0, so they can be used in their natural units (e.g. island area, 
#' species abudance), and will be scaled properly. If all species (or sites)
#' are assumed to be equally likely, the weight vector should be set to the
#' same constant for all elements.
#' @param speciesData binary presence-absence matrix 
#' (rows = species, columns = sites).
#' @param rowWeights vector of positive values representing species weights.
#' @param colWeights vector of positive values representing site weights.
#' @return Returns a binary presence-absence matrix with the same 
#' dimensions and fill as the input matrix.
#' @references Jenkins, D.G. 2006. In search of quorum effects in metacommunity
#' structure: species co-occurrence analyses. Ecology 87:1523-1531 
#'
#' Gotelli, N.J., G.R. Graves, and C. Rahbek. 2010. Macroecological
#' signals of species interactions in the Danish avifauna. Proceedings of the 
#' National Academy of Sciences, U.S.A. 107: 530-535.

#' @note sim10 allows users to incorporate independent data on species
#' occurrence probabilities and site suitabilities. This represents an
#' important conceptual advance over standard co-occurrence analyses, which
#' must infer these probabilities from the matrix itself. sim10 may generate
#' empty rows or columns, especially if weights are very small for some species
#' or sites. Also, the results may be sensitive to algebraic transformations 
#' of the weights (x, x^2, log(x), etc.), and these transformations may be hard 
#' to justify biologically. Nevertheless, sim10 is worth exploring for rich 
#' data sets with site and species attributes.
#' @seealso \code{\link{vector_sample}} for weighted vector sampling.
#' @examples 
#' randomMatrix <- sim10(speciesData=matrix(rbinom(40,1,0.5),nrow=8))
#' @export

sim10 <- function(speciesData,rowWeights = runif(dim(speciesData)[1]),colWeights = runif(dim(speciesData)[2])) 

{
matrixWeights <- outer(rowWeights,colWeights)
matrix(vector_sample(speciesData, weights =matrixWeights),ncol=ncol(speciesData))
} 

#' RA1 Niche Overlap Randomization Algorithm 
#' @description Randomizes a numeric utilization matrix speciesData by 
#' replacing all elements with a random uniform [0,1] value.
#' @details The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @param speciesData a resource utilization matrix (rows = species, columns = discrete 
#' resource states) filled with non-negative real numbers. 
#' @return Returns a random utilization matrix with the same dimensions as the 
#' input matrix.
#' @references Kobayashi, S. 1991. Interspecific relations in forest floor 
#' coleopteran assemblages: niche overlap and guild structure. Researches 
#' in Population Ecology 33: 345-360.
#' @note Because all matrix elements, including zeroes, are replaced with a 
#' random uniform distribution, the null expectation is based on an assemblage
#' of generalist species with maximum niche breadth. This algorithm retains
#' neither the niche breadth of the individuals species nor the placement of 0
#' values (= unutilized resource states) in the matrix. These assumptions are
#' unrealistic, and a random matrix with zeroes will generate significantly low
#' niche overlap values with this metric. It is not recommended for niche 
#' overlap analysis.
#' @examples 
#' ranUtil <- ra1(speciesData=matrix(rpois(40,0.5),nrow=8))
#' @export


ra1 <- function(speciesData=matrix(rpois(80,1),nrow=10)){

matrix(runif(prod(dim(speciesData))),ncol=ncol(speciesData))
}

#' RA2 Niche Overlap Randomization Algorithm
#' @description Randomizes a numeric utilization matrix speciesData by 
#' replacing all non-zero elements with a random uniform [0,1] value.
#' @description Randomizes a numeric utilization matrix speciesData by 
#' replacing all elements with a random uniform [0,1].
#' @details The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @param speciesData a resource utilization matrix (rows = species, columns = discrete 
#' resource states) filled with non-negative real numbers. 
#' @return Returns a random utilization matrix with the same dimensions as the 
#' input matrix.
#' @references Kobayashi, S. 1991. Interspecific relations in forest floor 
#' coleopteran assemblages: niche overlap and guild structure. Researches 
#' in Population Ecology 33: 345-360.
#'	
#' Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural
#' assemblages of desert lizards and tropical fishes. Ecological Monographs
#' 60: 27-55.
#' @note This algorithm retains the number and position of zero states in the 
#' original matrix. However, all non-zero values are again replaced by a random 
#' [0,1] value, which tends to inflate niche breadths of the simulated 
#' assemblage. Although the results are not as severe as for RA1, this 
#' algorithm is still prone to Type I errors, and is not recommended for niche 
#' overlap analysis.

#' @examples 
#' ranUtil <- ra2(speciesData=matrix(rpois(40,0.5),nrow=8))
#' @export

ra2 <- function(speciesData=matrix(rpois(80,1),nrow=10)) 

{
z <- which(speciesData > 0)
RM <- speciesData
RM[z] <- runif(length(z))
return(RM)
}

   

#' RA3 Niche Overlap Randomization Algorithm
#' @description Randomizes a numeric utilization matrix speciesData by 
#' reshuffling the elements within each row.
#' @details The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @param speciesData a resource utilization matrix (rows = species, columns = discrete 
#' resource states) filled with non-negative real numbers. 
#' @return Returns a random utilization matrix with the same dimensions as the 
#' input matrix.
#' @references Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural
#' assemblages of desert lizards and tropical fishes. Ecological Monographs
#' 60: 27-55.
#' @note This algorithm retains the niche breadth and zero states for each 
#' species, but randomizes the assignment of each utilization value to a 
#' different niche category. It performs effectively in simulation studies and 
#' is recommended for analysis of niche overlap patterns. 
#' @examples 
#' ranUtil <- ra3(speciesData=matrix(rpois(40,0.5),nrow=8))
#' @export

ra3 <- function(speciesData=matrix(rpois(80,1),nrow=10)) 

{
RM <- apply(speciesData,1,sample)
RM <- t(as.matrix(RM,nrow=nrow(speciesData),ncol=ncol(speciesData)))
return(RM)
}

   

#' RA4 Niche Overlap Randomization Algorithm
#' @description Randomizes a numeric utilization matrix speciesData by 
#' reshuffling the non-zero elements within each row. 
#' @details The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @param speciesData a resource utilization matrix (rows = species, columns = discrete 
#' resource states) filled with non-negative real numbers. 
#' @return Returns a random utilization matrix with the same dimensions as the 
#' input matrix.
#' @references Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural
#' assemblages of desert lizards and tropical fishes. Ecological Monographs
#' 60: 27-55.
#' @note This algorithm is similar to RA3, but adds the additional constraint of
#' retaining the positions of all of the zero elements of the matrix, and 
#' reshuffling only the non-zero elements of the matrix within each row. It is
#' more conservative than RA3, but has a low Type I error rate, and, along with 
#' RA3, is recommended for null model analysis of niche overlap. 
#' @examples 
#' ranUtil <- ra4(speciesData=matrix(rpois(40,0.5),nrow=8))
#' @export

ra4 <- function(speciesData=matrix(rpois(80,1),nrow=10)) 
  
{
    #.....................................
    ## Nonzero row shuffle function for RA4
    NonZeroRowShuffle <- function(vec=runif(10)) {
      nonzero <- which(vec > 0)
      shuffledvec <- vec
      shuffledvec[nonzero] <- vec[sample(nonzero)]
      return(shuffledvec)
    }
    #....................................
    
    RM <- t(apply(speciesData,1,NonZeroRowShuffle))
    return(RM)
    
  }

#' SizeUniform Size Overlap Randomization Algorithm
#' @description Function to randomize body sizes within a uniform distribution
#' with boundaries set by the largest and smallest species in the assemblage.
#' @details If the assemblage contains n species,
#' only the body sizes of the inner n - 2 species are randomized.
#' @param speciesData a vector of positive real values representing the body
#' sizes or trait values for each species.
#' @return Returns a vector of body sizes that have been randomly assigned. The
#' largest and smallest body sizes in the randomized assemblage match those in
#' the empirical data.
#' @references Simberloff, D. and W. Boecklen. 1981. Santa Rosalia
#' reconsidered: size ratios and competition. Evolution 35: 1206-1228.
#' 
#' Tonkyn, D.W. and B.J. Cole. 1986. The statistical analysis of size ratios.
#' American Naturalist 128: 66-81.
#' @note Although the distribution of body sizes may not be truly uniform,
#' it may be approximately uniform within the range of observed values,
#' particularly for small assemblages. 
#' @seealso \code{\link{size_gamma}} size distribution function.
#' @examples 
#' nullSizes <-size_uniform(speciesData=runif(20))
#'@export

size_uniform <- function(speciesData=runif(20)) {
  
  endpoints <- c(min(speciesData),max(speciesData))  # save max and min boundaries
  sim <- runif(n=(length(speciesData)-2),min=endpoints[1],max=endpoints[2])
  randomVec <- c(endpoints,sim)
  return(randomVec)
}

#' SizeUser Size Overlap Randomization Algorithm
#' @description Observed body sizes are randomized with a uniform distribution
#' for which the user has defined the minimum and maximum possible body size.
#' @details Within the user-defined limits, body sizes of all n species are
#' randomized, whereas uniform_size randomizes only n - 2 of the body
#' sizes and uses the extreme values to set the endpoints.
#' @param speciesData a vector of observed body sizes.
#' @param userLow a user-defined lower limit.
#' @param userHigh a user-defined upper limit.
#' @return Returns a vector of randomized body sizes.	
#' @note As the difference between the lower and upper boundaries is increased
#' the test will yield results that are random or aggregated, even though the 
#' same data might yield a segregated pattern when the uniform_size algorithm
#' is used. For this reason, this algorithm is not recommended for size ratio
#' analyses.
#' @seealso \code{\link{size_uniform}} size distribution algorithm.
#' @examples 
#' nullSizes <- size_uniform_user(speciesData=runif(20,min=10,max=20),userLow=8,userHigh=24)
#' @export

size_uniform_user <- function(speciesData=runif(n=20),userLow=0.9*min(speciesData),
                              userHigh=1.1*max(speciesData)){
#  if(!is.null(Param.List$Special)){User.low <- Param.List$Special[1]
#                                    User.high <- Param.List$Special[2]}
  randomVec <- runif(n=length(speciesData),min=userLow,max=userHigh)
  
  return(randomVec)
}

#' SizeSourcePoolDraw Size Overlap Randomization Algorithm
#' @description Function to randomize body sizes by drawing species from a 
#' user-defined source pool. Species are drawn without replacement, 
#' and there is a specified probability vector for the source pool species

#' @param speciesData a vector of observed body sizes.
#' @param sourcePool a vector of body sizes of species in the user-defined
#' pool of potential colonists.
#' @param speciesProbs a vector of relative colonization weights of 
#' length 'sourcePool'.
#' @return Returns a vector of body sizes of an assemblage randomly drawn
#' from a user-defined source pool.
#' @references Strong, D.R. Jr., L.A. Szyska, and D. Simberloff. 1979. Tests of
#' community-wide character displacement against null hypotheses. Evolution 33:
#' 897-913.

#' Schluter, D. and P.R. Grant. 1984. Determinants of morphological patterns in
#' communities of Darwin's finches. American Naturalist 123: 175-196.
#' @note Although delineating a source pool of species and estimating their
#' relative colonization probabilities is difficult, this is the most realistic
#' approach to constructing a null distribution.
#' @examples 
#' obsOverlap <- size_source_pool(dataRodents$Sonoran)
#' @export
size_source_pool <- function(speciesData=21:30,sourcePool=
                               runif(n=2*length(speciesData),min=10,max=50),
                             speciesProbs=rep(1,length(sourcePool))) {

  speciesDraw <- sample(sourcePool,size=length(speciesData),replace=FALSE,
                         prob=speciesProbs)
  
  return(speciesDraw)
}

#' SizeGamma Size Overlap Randomization Algorithm
#' @description Function to generate a random distribution of body sizes by
#' drawing from a gamma distribution. Shape and rate parameters of the gamma
#' are estimated from the empirical data.
#' @param speciesData a vector of body sizes or other trait measurements of
#' species. All values must be positive real numbers.
#' @return Returns a vector of simulated body sizes as the same length as speciesData.
#' @note The shape and rate parameters are estimated from the real data using
#' the maximum likelihood estimators generated from the fitdr function in 
#' the MASS library. The flexible gamma distribution can be fit to a variety of 
#' normal, log-normal, and exponential distributions that are typical for trait 
#' data measured on a continuous non-negative scale. 
#' @seealso \code{\link{fitdistr}} in the MASS library.
#' @examples 
#' obsOverlap <- size_gamma(speciesData=rnorm(50,mean=100,sd=1))
#' @import MASS
#' @export

size_gamma <- function (speciesData=rnorm(50,mean=100,sd=1)) {
 
  mleFit <- fitdistr(speciesData, 'gamma')
  a <- mleFit$estimate["shape"]
  b <- mleFit$estimate["rate"]
  gammaDraw <- rgamma(length(speciesData),shape=a,rate=b)

  return(gammaDraw)
}
