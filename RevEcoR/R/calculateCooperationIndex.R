##==============================================================================
##==============================================================================
#' Calculating the species interactions
#' 
#' Calculating the metabolic complementarity index and complementarity index of 
#' based on species metabolic network.
#' 
#'@param g1, igraph object, a species-specific metabolic network.
#'
#'@param g2, igraph object, a species-specific metabolic network, the complementary
#' network of g1
#'
#'@param seed.set1, seeds slot of a seed-set object, seeds of the metabolic 
#'network g1, more details see \code{\link{seedset-class}}.
#'
#'@param seed.set2, seeds slot of a seed-set object, seeds of the metabolic 
#'network g2, more details see \code{\link{seedset-class}}.
#'
#'
#'@param threshold, the cutoff of confidence score to be serve as a seed set, 
#'default is 0.
#'
#'@param p, a logical value which determins whether the calculated index is 
#'statistical or biological significant. default is FALSE.
#'
#'@param nperm, the number of permuations of metabolic network node labes, which
#'is used for complementarity index's P value calculating, default is 1000.
#'
#'@return a two length list: complementarity index or competition index: range 
#'from 0 to 1, p value of complementarity index. Or a single value of 
#'complementarity or competition index while p is FALSE.
#'
#'@details Metabolic competition index is defined as the fraction of compounds 
#'in a species seed set of metabolic network that are also included in its partner; 
#'However, metabolic complementarity index is the fraction of compounds in one 
#'species seed set of metabolic network appearing in the metabolic network but 
#'not in the seed set of its partner. However, seed compounds are associated 
#'with a confidence score (1/size of SCC), so this fraction is calculated as a 
#'normalized weighted sum.
#'
#' Based on the metabolic network and seed sets of species, this functions help 
#' us to predict the species interactions of species1 on the presence of species2. 
#'
#'@name Interactions
#'@rdname Interactions
#'
#'@seealso \code{\link{getSeedSets}},
#'\code{\link{calculateCooperationIndex}}
#'
#'@examples 
#'\dontrun{
#'## metabolic network reconstruction and seed set identity of sample data anno.species
#'net <- lapply(anno.species,reconstructGsMN)
#'seed.sets <- lapply(net, getSeedSets) 
#'seed.sets <- lapply(seed.sets, function(x)x@@seeds)
#'
#'## calculate the complementarity index of the first species
#'complementarity.index <- complementarityIndex(net[[1]],net[[2]], 
#'  seed.sets[[1]], seed.sets[[2]])
#'competition.index <- competitionIndex(net[[1]],net[[2]], 
#'  seed.sets[[1]], seed.sets[[2]])
#'}
NULL


#'@export
#'@rdname Interactions

complementarityIndex <- function(g1,g2, seed.set1, seed.set2, threshold = 0, 
  p  = FALSE, nperm = 1000){
  
  ## calculating the normalized complementarity index
  normalize_complementarity <- function(complement.cpd, seed.set1) {
    if (length(complement.cpd)) {
      norm.seed <- lapply(seed.set1, 
        function(x)match(complement.cpd,x, nomatch = 0)) %>% 
        sapply(sum) %>% is_greater_than(0) %>% sum
      complementary.index <- norm.seed / length(seed.set1)
    } else
      complementary.index <- 0
  }
    
  ## random permutation complementarity index, random permutation the nodes 
  ## label of metabolic network, which can preserve the original graph's 
  ## seed set number and size.
  perm_complementary <- function(seed.set1, seed.set2, g1, nonseed2){
    seed1.len <- lengths(seed.set1)
    sample.seed1 <- sample(V(g1)$name,sum(seed1.len)) %>% 
      split(rep(1:length(seed.set1), seed1.len))
    ##nonseed2 <- setdiff(V(g2)$name, unlist(seed.set2))
    complement.cpd <- setdiff(unlist(sample.seed1),unlist(seed.set2)) %>%
      intersect(., nonseed2)
    perm.complementary <- normalize_complementarity(complement.cpd, sample.seed1)
    return(perm.complementary)
  }
  
  if (!is.igraph(g1) || !is.igraph(g2))
    stop("Both g1 and g2 must be igraph object")
  ##seed.set1 <- getSeedSets(g1,threshold)@seeds
  ##seed.set2 <- getSeedSets(g2, threshold)@seeds
  nonseed2 <- setdiff(V(g2)$name, unlist(seed.set2))
  complement.cpd <- setdiff(unlist(seed.set1),unlist(seed.set2)) %>%
    intersect(., nonseed2)
  complementary.index <- normalize_complementarity(complement.cpd, seed.set1)
  if (p){
    perm.complementary <- replicate(nperm,
      perm_complementary(seed.set1, seed.set2, g1, nonseed2))
    p.value <- sum(perm.complementary < complementary.index, 1)/(nperm + 1)
      p.value <- round(p.value, digits = 5)
    return(list(complementary.index = complementary.index,
      p.value = p.value))
  }else
    return(complementary.index)
}


#'@export
#'@rdname Interactions
competitionIndex <- function(g1,g2, seed.set1, seed.set2, threshold=0, 
  p = FALSE, nperm = 1000){
  
  ## calculating the normalzied competition index
  normalize_competition <- function(intersect.seed, seed.set1){ 
    if(length(intersect.seed)){
      norm.seed <- lapply(seed.set1, 
        function(x)match(intersect.seed, x, nomatch = 0)) %>% 
        sapply(sum) %>%  is_greater_than(0) %>% sum
      competition.index <- norm.seed/length(seed.set1)
    }else{
      competition.index <- 0
    }
  }
  
  ## random permutation competition index
  perm_competition <- function(seed.set1, seed.set2, g1){
    seed1.len <- lengths(seed.set1)
    sample.seed1 <- sample(V(g1)$name,sum(seed1.len)) %>% 
      split(rep(1:length(seed.set1), seed1.len))
    intersect.seed <- intersect(unlist(sample.seed1),unlist(seed.set2))
    return(normalize_competition(intersect.seed, sample.seed1))
  }
  
  if (!is.igraph(g1) || !is.igraph(g2))
    stop("Both g1 and g2 must be igraph object")
  ##seed.set1 <- getSeedSets(g1,threshold)@seeds
  ##seed.set2 <- getSeedSets(g2, threshold)@seeds
  intersect.seed <- intersect(unlist(seed.set1),unlist(seed.set2))
  if(length(intersect.seed)){
    norm.seed <- lapply(intersect.seed,function(x)lapply(seed.set1,
      function(y)match(x,y,nomatch=0))) %>%
      lapply(.,function(x)which(x>0)) %>%
      unique %>%
      length
    competition.index <- norm.seed/length(seed.set1)
  }else{
    competition.index <- 0
  }
  if (p){
    perm.competition <- replicate(nperm, 
      perm_competition(seed.set1, seed.set2, g1))
    p.value <- sum(perm.competition > competition.index, 1)/(nperm + 1) 
    p.value <- round(p.value, digits = 5)
    return(list(competition.index = competition.index, p.value = p.value))
  }else
    return(competition.index)
}

##==============================================================================
##==============================================================================
#'Calculating the metabolic competition and complementarity index
#'
#'Calculating the metabolic competition complementarity index among all metabolic 
#'networks
#'
#'@param g, igraph that represents a metabolic network, see \code{\link{reconstructGsMN}}
#'
#'@param ..., a list of metabolic networks or a network append to g
#'
#'@param threshold threshold, the cutoff of confidence score to be serve as a 
#'seed set, default is 0.2
#'
#'@param p, a logical value which determins whether the calculated index is 
#'statistical or biological significant. default is FALSE
#'
#'@param nperm, the number of permuations of metabolic network node labes, which
#'is used for p value calculation, default is 1000.
#'
#'@export
#'
#'@details Metabolic competition index is defined as the fraction of compounds 
#'in a species seed set of metabolic network that are alse included in its 
#'partner; However, metabolic complementarity index is the fraction of 
#'compounds in one species seed set of metabolic network appearing in the 
#'metabolic network but not in the seed set of its partner; The biosynthetic 
#'support score represents the extent to which the metabolic requirements of a 
#'potential parasitic organism can be supported by the biosynthetic capacity of
#'a potential host. It is measured by calculating the fraction of the source 
#'components of a, in which at least one of the compounds can be found in the 
#'network of b. However, seed compounds are associated with a confidence score 
#'(1/size of SCC), so this fraction is calculated as a mormalized weighted sum. 
#'
#'The ith row and jth col elements of the returnd matrix represents the 
#'metabolic competition index or complementarity index of the ith network on the
#'jth metabolic network.
#'
#'@return a cooperation index matrix whose nrow and ncol is equal to the number 
#'of species to be compared, for more see details.
#'
#'@seealso \code{\link{complementarityIndex}},
#'\code{\link{competitionIndex}}
#'
#'@examples
#'\dontrun{
#'## metabolic network reconstruction and seed set identity of sample data anno.species
#'net <- lapply(anno.species,reconstructGsMN)
#'interactions <- calculateCooperationIndex(net)
#'}

calculateCooperationIndex <- function(g, ...,threshold=0, p = FALSE, nperm = 1000){
  if (is.igraph(g)){
    g <- c(list(g), list(...))
  }else{
    g <- c(g, list(...))
  }
  l <- length(g)
  if (l < 2)
    stop("At least two species to compare")
  seed.set <- lapply(g, getSeedSets) %>% lapply(function(x)x@seeds)
  competition.index <- matrix(1,l,l)
  complementarity.index <- matrix(0,l,l)
  index <- permutations(l,2) %>% t %>% as.data.frame
  perm.seed.set <- lapply(index,function(x)seed.set[x])
  perm.g <- lapply(index, function(x)g[x])
  competition.v <- map2(perm.g,perm.seed.set,
    function(x,y)competitionIndex(x[[1]],x[[2]],y[[1]],y[[2]], 
      threshold, p, nperm))
  complementarity.v <- map2(perm.g, perm.seed.set,
    function(x,y)complementarityIndex(x[[1]],x[[2]],y[[1]],y[[2]], 
      threshold, p, nperm))
  for(i in 1:length(index)){
    competition.index[index[1,i],index[2,i]] = unlist(competition.v[i])[1]
    complementarity.index[index[1,i],index[2,i]] = unlist(complementarity.v[i])[1]
  }
  row.names(competition.index)  <- names(g)
  colnames(competition.index)  <- names(g)
  row.names(complementarity.index) = colnames(complementarity.index)  <- names(g)
  if (p){ 
    competition.index.p <- matrix(0,l,l)
    row.names(competition.index.p) = colnames(competition.index.p) <- names(g)
    complementarity.index.p <- matrix(0,l,l)
    row.names(complementarity.index.p) = colnames(complementarity.index.p)  <- names(g)
    for(i in 1:length(index)){
      competition.index.p[index[1,i],index[2,i]] = unlist(competition.v[i])[2]
      complementarity.index.p[index[1,i],index[2,i]] = unlist(complementarity.v[i])[2]
    }
    return(list(competition.index = competition.index, 
      competition.index.p = competition.index.p,
      complementarity.index = complementarity.index, 
      complementarity.index.p = complementarity.index.p)) 
  } else
    return(list(competition.index = competition.index,
      complementarity.index = complementarity.index))
}
##==============================================================================
##==============================================================================