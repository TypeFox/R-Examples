
#' Pianka Niche Overlap Metric
#' @description Takes a resource utilization matrix as input and
#' returns the average pairwise Pianka's niche overlap index.
#' @details Pianka's niche overlap index is averaged 
#' over each unique species pair. The index is symmetric, 
#' with a normalization term in the denominator for the overlap 
#' between species 1 and 2. Values of Pianka's niche overlap index close to 0.0
#' reflect usage of exclusive resource categories, whereas values close to 1.0
#' reflect similar resource utilization spectra.
#' \deqn{O_{jk} = O_{kj} = \frac{\sum_{n}^{i}p_{ij}p_{jk}}{\sqrt{\sum_{n}^{i}p_{ij}^2\sum_{n}^{i}p_{ik}^2}}}{O_jk = O_kj = sum(p_ij*p_jk) / sqrt(sum((p_ij)^2)sum((p_jk)^2))}
#' 
#' @param m a matrix of resource utilization values. 
#' @return Returns the average pairwise niche overlap.
#' @references Pianka, E. 1973. The structure of lizard communities.
#' Annual Review of Ecology and Systematics 4:53-74.
#' 
#' Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural assemblages
#' of desert lizards and tropical fishes. Ecological Monographs 60: 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{czekanowski}} niche overlap index.
#' @examples 
#' obsOverlap <- pianka(m=matrix(rpois(40,0.5),nrow=8))
#' @export
#'

pianka <- function(m=matrix(rpois(80,1),nrow=10)) 

{
m <- m/rowSums(m)
pairwise <- cbind(t(combn(nrow(m),2)),0)	# set up pairwise species list
for (i in 1:nrow(pairwise)) 
	pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
	sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
return(mean(pairwise[,3]))
}


#' Czekanowski Niche Overlap Metric
#' @description Takes a resource utilization matrix as input and
#' returns the average pairwise Czekanowki niche overlap index.
#' @details The Czekanowski niche overlap index is averaged 
#' over each unique species pair. The index measures the area of intersection
#' of the resource utilization histograms of each species pair. 
#' Values of Czekanowski niche overlap index close to 0.0
#' reflect usage of exclusive resource categories, whereas values close to 1.0
#' reflect similar resource utilization spectra.
#' 
#' \deqn{O_{jk} = O_{kj} = 1 - 0.5\sum_{i=1}^n|p_{ij} - p_{ik}|}{O_jk = O_kj = 1 - 0.5*sum|p_ij - p_ik|}' 
#' 
#' @param m a matrix of resource utilization values. 
#' @return Returns the average pairwise niche overlap.
#' @references Feinsinger, P., E.E. Spears, and R. Poole. 1981. A simple measure 
#' of niche breadth. Ecology 62: 27-32. 
#' 
#' Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural assemblages
#' of desert lizards and tropical fishes. Ecological Monographs 60: 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{pianka}} niche overlap index.
#' @examples 
#' obsOverlap <- czekanowski(m=matrix(rpois(40,0.5),nrow=8))
#' @export
#'

czekanowski <- function(m=matrix(rpois(80,1),nrow=10)) 

{
m <- m/rowSums(m)
pairwise <- cbind(t(combn(nrow(m),2)),0)	# set up pairwise species list
for (i in 1:nrow(pairwise)) 
	pairwise[i,3] <- 1 - 0.5*sum(abs((m[pairwise[i,1],] - m[pairwise[i,2],])))

return(mean(pairwise[,3]))
}


#' PiankaVariance Niche Overlap Metric
#' @description Takes a niche utilization matrix as in put and 
#' returns the variance of Pianka's niche overlap index.
#' @details A large value for variance implies that some species pairs show high
#' niche overlap and others show low niche overlap. A low value for variance 
#' implies that niche overlap (high or low) is very similar among all species pairs.
#'
#' @param m a matrix of resource utilization values. 
#' @return Returns the variance of the average pairwise niche overlap.
#' @references Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural 
#' assemblages of desert lizards and tropical fishes. Ecological Monographs 60: 
#' 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{pianka}} niche overlap index.
#' @examples 
#' obsVar <- pianka_var(m=matrix(rpois(40,0.5),nrow=8))
#' @export

pianka_var <- function(m=matrix(rpois(80,1),nrow=10))
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
      sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
  
  return(var(pairwise[,3]))
}

##
##
#' CzekanowskiVariance Niche Overlap Metric
#' @description Takes a niche utilization matrix returns the variance of the 
#' Czekanowski niche overlap index
#' 
#' @details A large value for variance implies that some species pairs show high
#' niche overlap and others show low niche overlap. A low value for variance 
#' implies that niche overlap (high or low) is very similar among all species pairs.
#'
#' @param m a matrix of resource utilization values. 
#' @return Returns the variance of the average pairwise niche overlap.
#' @references Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural 
#' assemblages of desert lizards and tropical fishes. Ecological Monographs 60: 
#' 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{czekanowski}} niche overlap index.
#' 
#' @examples 
#' obsVar <- czekanowski_var(m=matrix(rpois(40,0.5),nrow=8))
#' @export

czekanowski_var <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)	# set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- 1 - 0.5*sum(abs((m[pairwise[i,1],] - m[pairwise[i,2],])))
  
  return(var(pairwise[,3]))
}


#' PiankaSkew Niche Overlap Metric
#' @description Takes a niche utilization matrix returns the skewness of the 
#' Pianka pairwise niche overlap index.
#' 
#' @details A large positive value for skewness implies that there are more species pairs
#' with high than low niche overlap. A large negative value for skewness implies there are more 
#' species pairs with low than high niche overlap. The performance of this algorithm
#' has not been thoroughly tested with real data sets.
#'
#' @param m a matrix of resource utilization values. 
#' @return Returns the skewness of the average pairwise niche overlap.
#' @references Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural assemblages
#' of desert lizards and tropical fishes. Ecological Monographs 60: 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{pianka}} niche overlap index.
#' 
#' @examples 
#' obsSkew<- pianka_skew(m=matrix(rpois(40,0.5),nrow=8)) 
#' @export

pianka_skew <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
    sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
  
    m3 <- mean((pairwise[,3]-mean(pairwise[,3]))^3)
    piankaSkew <- m3/(sd(pairwise[,3])^3)
    
  
  return(piankaSkew)
}


#' CzekanowskiSkew Niche Overlap Metric
#' @description Takes a niche utilization matrix returns the skew of the 
#' Czekanowski pairwise niche overlap index.
#' 
#' @details A large positive value for skewness implies that there are more species pairs
#' with high than low niche overlap. A large negative  value for skewness implies there are more 
#' species pairs with low than high niche overlap. The performance of this algorithm
#' has not been thoroughly tested with real data sets.
#'
#' @param m a matrix of resource utilization values. 
#' @return Returns the skewness of the average pairwise niche overlap.
#' @references Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural assemblages
#' of desert lizards and tropical fishes. Ecological Monographs 60: 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{czekanowski}} niche overlap index.
#' 
#' @examples 
#' obsSkew <- czekanowski_skew(m=matrix(rpois(40,0.5),nrow=8)) 
#' @export

czekanowski_skew <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- 1 - 0.5*sum(abs((m[pairwise[i,1],] - m[pairwise[i,2],])))
  
  m3 <- mean((pairwise[,3]-mean(pairwise[,3]))^3)
  czekanowskiSkew <- m3/(sd(pairwise[,3])^3)
  
  
  return(czekanowskiSkew)
}

#' SpeciesCombo Co-occurrence Metric
#' @description Function to calculate number of unique species combinations in a matrix
#' @details In Diamond's (1975) assembly rules model, species interactions lead to 
#' certain "forbidden combinations" of species. A set of communities structured 
#' this way should contain fewer species combinations than expected by chance.
#'
#' @param m a binary presence-absence matrix in which rows are species and columns
#' are sites. 
#' @return Returns the number of unique species combinations represented by 
#' the different columns (= sites) in the matrix.
#'  
#' @references Diamond, J.M. 1975. Assembly of species communities. p. 342-444 in:
#' Ecology and Evolutoin of Communities. M.L. Cody and J.M. Diamond (eds.). 
#' Harvard University Press, Cambridge. 
#' 
#' Pielou, D.P. and E.C. Pielou. 1968. Association among species of infrequent
#' occurrence: the insect and spider fauna of Polyporus betulinus (Bulliard) Fries.
#' Journal of Theoretical Biology 21: 202-216.
#' 
#' @note This metric is most useful when the number of sites (= columns) is relatively
#' large compared to the number of species (= rows). Empty sites are excluded 
#' from the matrix and are not counted as a unique species combination.
#' 
#' 
#' @examples 
#' obsCombo <- species_combo(m=matrix(rbinom(100,1,0.5),nrow=10)) 
#' @export
species_combo <- function(m=matrix(rbinom(100,1,0.5),nrow=10))
{
  
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  return(ncol(unique(m,MARGIN=2))) # number of columns in submatrix of uniques
  
}

#' Checker Co-occurrence Metric
#' @description Function to calculate number of unique pairs of species 
#' that never co-occur and form a "checkerboard pair".
#' @details In Diamond's (1975) assembly rules model, pairs of species that 
#' never co-occur in any site are interpreted as examples of interspecific competition.  
#' A set of communities structured this way should contain more checkerboard
#' pairs than expected by chance.
#'
#' @param m A binary presence-absence matrix in which rows are species and columns
#' are sites. 
#' @return Returns the number of unique species pairs that never co-occur.
#'  
#' @references Diamond, J.M. 1975. Assembly of species communities. p. 342-444 in:
#' Ecology and Evolution of Communities. M.L. Cody and J.M. Diamond (eds.). 
#' Harvard University Press, Cambridge. 
#' 
#' Connor, E.F. and D. Simberloff. 1979. The assembly of species communities: chance
#' or competition? Ecology 60: 1132-1140.
#' 
#' @examples 
#' obsChecker <- checker(m=matrix(rbinom(100,1,0.5),nrow=10)) 
#' @export
#' 
checker <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  shared <- mat.or.vec(1,nrow(pairwise))
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
  }
  
  
  
  return(sum(shared==0)) # return number of pairs with no shared sites
  
}

#' CScore Co-occurrence Metric
#' @description Takes a binary presence-absence matrix and returns 
#' Stone and Roberts' (1990) C-score.
#' @details For each unique pair of species, the C-score is calculated as
#' 
#' \deqn{C_{ij} = (R_i - S)(R_j - S)}{C_ij = (R_i - S)(R_j - S)}
#' 
#' where R_i and R_j are the row sums for species i and j, and S is the number 
#' of shared sites in which both species i and species j are present. For any 
#' particular species pair, the larger the C-score, the more segregated the 
#' pair, with fewer shared sites. However, the index can be difficult to 
#' interpret when calculated as a matrix-wide average, because a single matrix
#' can contain individual pairs of species that are segregated, random, or aggregated.
#' 
#' Degenerate matrices result from simulations where a row or column sum may be 0. <nick can you fill in the implications as to what this means if they are included or not?>
#'
#' @param m a binary presence-absence matrix in which rows are species and columns
#' are sites. 
#' @return Returns the average C-score calculated across all possible species pairs
#' in the matrix.
#'  
#' @references Stone. L. and A. Roberts. 1990. The checkerboard score and species
#' distributions. Oecologia 85: 74-79.
#' 
#' Gotelli, N.J. and W. Ulrich. 2010. The empirical Bayes approach as a tool to 
#' identify non-random species associations. Oecologia 162:463-477.
#' 
#' @note The matrix-wide C-score is not calculated for missing species, so empty
#' rows in the matrix do not affect the result.
#' @examples 
#' obsCScore <- c_score(m=matrix(rbinom(100,1,0.5),nrow=10)) 
#' @export
#' 

c_score <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  shared = tcrossprod(m)
  sums = rowSums(m)
  
  upper = upper.tri(shared)
  
  scores = (sums[row(shared)[upper]] - shared[upper])*
      (sums[col(shared)[upper]] - shared[upper])
  
  return(mean(scores))
}


#' CScoreVariance Co-occurrence Metric
#' @description Takes a binary presence-absence matrix and returns 
#' the variance of the Stone and Roberts' (1990) C-score.
#' @details A large value of this variance implies that some species pairs 
#' in the matrix are strongly segregated (large C-score) and other species pairs 
#' are random or aggregated.
#'
#' @param m a binary presence-absence matrix in which rows are species and columns
#' are sites. 
#' @return Returns the variance of the C-score calculated across all possible 
#' species pairs in the matrix.
#'  
#' @references Stone, L. and A. Roberts. 1990. The checkerboard score and species
#' distributions. Oecologia 85: 74-79.
#' 
#' Stone, L. and A. Roberts. 1992. Competitive exclusion, or species aggregation?
#' An aid in deciding. Oecologia 91: 419-424.
#' 
#' @note The matrix-wide C-score is not calculated for missing species, so empty
#' rows in the matrix do not affect the result. This index has not been 
#' thoroughly tested with real data sets.
#' 
#' @seealso \code{\link{c_score}} co-occurrence index.
#' 
#' @examples 
#' varCScore <- c_score_var(m=matrix(rbinom(100,1,0.5),nrow=10)) 
#' @export

c_score_var <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  shared = tcrossprod(m)
  sums = rowSums(m)
  
  upper = upper.tri(shared)
  
  scores = (sums[row(shared)[upper]] - shared[upper])*
    (sums[col(shared)[upper]] - shared[upper])
  
  return(var(scores))  
}


#' CScoreSkew Co-occurrence Metric
#' @description Takes a binary presence-absence matrix and returns 
#' the skewness of the Stone and Roberts' (1990) C-score.
#' @details A large positive value of skewness implies a preponderance of species pairs 
#' with large C-score values (segregated), whereas a large negative value of 
#' skewness implies a preponderance of species pairs with small C-score values 
#' (aggregated). 
#'
#' @param m a binary presence-absence matrix in which rows are species and columns
#' are sites. 
#' @return Returns the skewness of the C-score calculated across all possible 
#' species pairs in the matrix.
#'  
#' @references Stone, L. and A. Roberts. 1990. The checkerboard score and species
#' distributions. Oecologia 85: 74-79.
#' 
#' Stone, L. and A. Roberts. 1992. Competitive exclusion, or species aggregation?
#' An aid in deciding. Oecologia 91: 419-424.
#' 
#' @note The matrix-wide C-score is not calculated for missing species, so empty
#' rows in the matrix do not affect the result. This index has not been 
#' thoroughly tested with real data sets.
#' 
#' @seealso \code{\link{c_score}} co-occurrence index.
#' 
#' @examples 
#' skewCScore <- c_score_skew(m=matrix(rbinom(100,1,0.5),nrow=10)) 
#' @export

#' 
c_score_skew <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  shared = tcrossprod(m)
  sums = rowSums(m)
  
  upper = upper.tri(shared)
  
  scores = (sums[row(shared)[upper]] - shared[upper])*
    (sums[col(shared)[upper]] - shared[upper])
  
  
  m3 <- mean((scores-mean(scores))^3)
  cScoreSkew <- m3/(sd(scores)^3)
  
  
  return(cScoreSkew)  # return skewness of pairwise C-score
  
}

#' SchlutersVRatio Co-occurrence Metric
#' @description Takes a binary presence-absence matrix or a matrix of 
#' abundances and returns Schluter's (1984) variance ratio. 
#' @details The variance ratio is the ratio of the variance in species number 
#' among sites to the sum of the variance of the species occurrences. If the average
#' covariation in abundance (or occurrence) of each species pair is close to zero,
#' the expected value for this ratio is approximately 1.0. V-ratios larger than 1.0
#' imply positive average covariation in the abundance of species pairs, whereas V-ratios
#' significantly smaller than 1.0 imply negative average covariation.
#'
#' @param m a binary presence-absence matrix in which rows are species and columns
#' are sites. The entries may be either abundances or occurrences of indivdual species.
#' @return Returns the variance ratio of the matrix.
#'  
#' @references Schluter, D. 1984. A variance test for detecting species associations,
#' with some example applications. Ecology 65: 998-1005.
#' 
#' McCulloch, C.E. 1985. Variance tests for species association. Ecology 66: 1676-1681.
#' 
#' @note This index is determined exclusively by the row and column sums of the 
#' matrix, so it cannot be used with null model algorithms that hold both of those
#' elements fixed. A simple randomization of the rows of the matrix (see sim2)
#' assumes that all sites are equiprobable, so it may generate large values (= 
#' positive covariance) that reflect heterogeneity among sites.
#' 
#'
#' 
#' @examples 
#' varCScore <- v_ratio(m=matrix(rbinom(100,1,0.5),nrow=10)) 
#' @export

v_ratio <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  
  v <- var(colSums(m))/sum(apply(m,1,var))
  return(v)
}


#' MinDiff Size Overlap Metric
#' @description Function to calculate the minimum absolute size difference
#' between species within an assemblage.
#' @details Although this index is typically used to examine body size differences
#' in an animal assemblage, it could be used for any morphological index,	 
#' such as beak size, or for a phenological "trait", such as peak flowering time in a plant
#' assemblage.
#'
#' @param m a vector of non-negative trait measures, one for each species
#' @return Returns the minimum difference between adjacent, ordered values.
#' 
#' @references Simberloff, D. and W.J. Boecklen. 1981. Santa Rosalia reconsidered: size
#' ratios and competition. Evolution 35: 1206-1228.
#' 
#' @examples 
#' MinSizeDif <- min_diff(rgamma(20,shape=3,scale=2)) 
#' @export
#' 
min_diff <- function(m=runif(20)) {
  m <- sort(m)
  md <- min(diff(m))
  return(md)
}

#' MinRatio Size Overlap Ratio Metric
#' @description Function to calculate the minimum size ratio (larger/next larger)
#' between species within an assemblage.
#' @details This index is based on the minimum size ratio (larger/next larger) difference between
#' consecutively ordered species in an assemblage. It is appropriate for 
#' morphological traits, but not phenological ones.
#'
#' @param m a vector of non-negative trait measures, one for each species
#' @return Returns the minimum size ratio difference between adjacent, ordered values.
#' 
#' @references Simberloff, D. and W.J. Boecklen. 1981. Santa Rosalia reconsidered: size
#' ratios and competition. Evolution 35: 1206-1228.
#' 
#' @examples 
#' MinSizeDif <- min_ratio(rgamma(20,shape=3,scale=2)) 
#' @export

min_ratio <- function(m=runif(20)) {
  m <- sort(log(m))
  mr <- exp(min(diff(m)))
  return(mr)
}



#' VarDiff Size Overlap Ratio Metric
#' @description Function to calculate the variance in size differences
#' between adjacent, ordered species. If there is a tendency towards a
#' constant absolute size difference between adjacent species, this 
#' variance will be relatively small. Alternatively, if some adjacent species
#' are close in size, but others are very distant, this variance will be 
#' large. Small variances might be indicative of assemblages in which 
#' there is a competitively-based limit to similarity.
#'
#' @param m a vector of non-negative trait measures, one for each species
#' @return Returns the variance of the absolute difference between adjacent, ordered values.
#' 
#' @references Poole, R.W. and B.J. Rathcke. 1979. Regularity, randomness,
#' and aggregation in flowering phenologies. Science 203:470-471.
#'
#' Simberloff, D. and W.J. Boecklen. 1981. Santa Rosalia reconsidered: size
#' ratios and competition. Evolution 35: 1206-1228.
#' 
#' @examples 
#' SizeDifVar <- var_diff(rgamma(20,shape=3,scale=2)) 
#' @export

var_diff <- function(m=runif(20)){
  m <- sort(m)
  vd <- var(diff(m))
  return(vd)
}

#' VarRatio Size Overlap Ratio Metric
#' @description  Function to calculate the variance in size ratios
#' between adjacent, ordered species. If there is a tendency towards a
#' constant size ratio between adjacent species, this 
#' variance will be relatively small. Alternatively, if some adjacent species
#' are close in size, but others are very distant, this variance will be 
#' large. Small variances might be indicative of assemblages in which 
#' there is a competitively-based limit to similarity.
#'
#' @param m a vector of non-negative trait measures, one for each species
#' @return Returns the variance of the size ratios between adjacent, ordered values.
#' 
#' @references Poole, R.W. and B.J. Rathcke. 1979. Regularity, randomness,
#' and aggregation in flowering phenologies. Science 203:470-471.
#'
#' Simberloff, D. and W.J. Boecklen. 1981. Santa Rosalia reconsidered: size
#' ratios and competition. Evolution 35: 1206-1228.
#' 
#' @examples 
#' SizeRatioVar <- var_ratio(rgamma(20,shape=3,scale=2)) 
#' @export

var_ratio <- function(m=runif(20)){
  m <- sort(log(m))
  vr <- var(exp(diff(m)))
  return(vr)
}
