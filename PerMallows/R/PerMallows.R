# -----------------------------------------
# Author: Ekhine Irurozki
#         University of the Basque Country
# -----------------------------------------

###################   GENERAL   ###################

.CAYLEY.DISTANCE  <- 0
.KENDALL.DISTANCE <- 1
.HAMMING.DISTANCE <- 2
.ULAM.DISTANCE    <- 3

.MALLOWS.MODEL <- 0
.GENERALIZED.MALLOWS.MODEL <-1

#' Generate every permutation of perm.length item
#' 
#' This functions returns a matrix in thich each of rows
#' is a different permutation of the specified number of items
#'
#' @param perm.length number of items in the permutation
#' @param alert optional ask for confirmation when the number of permtuations to show is very large
#' @return A collection of every permutation of the specified number of items
#' @export
#' @examples
#' permutations.of(3)
#' permutations.of(10)
permutations.of <- function(perm.length, alert=TRUE){
  if ( alert && perm.length > 9) {
    cat("Do you really want to generate all ",factorial(perm.length)," permutations? (Y/N)")
    answer = readline()
    if ( answer != 'Y' ) return("Process cancelled")
  }
  return (.permutations.all.recursive(c(), c(1:perm.length)))
}

.permutations.all.recursive <- function(vec1, vec2){
  len <- length(vec2)
  perms = matrix(NA, ncol=length(vec2)+length(vec1), nrow=0)
  if(len == 0) return (aperm(as.matrix(vec1))) 
  else
    for(i in 1:len)
      perms <- rbind(perms, .permutations.all.recursive(c(vec1, vec2[i]), vec2[-i]))
  return (perms)
}

.same.length.perms <- function(s1, s2 ){
  return(is.permutation(s1) &&  is.permutation(s2) && length(s1) == length(s2))
}

#' Generate identity the permutation
#'
#' This function generates the identity permutation of a given number of items
#'
#' @param perm.length number of items in the permutation
#' @return The identity permutation of the specified number of items
#' @export
#' @examples
#' identity.permutation(3)
#' identity.permutation(7)
identity.permutation<-function(perm.length){
  return(c(1:perm.length));
}

#' Generate inverse permutation
#'
#' This function generates the inverse of a given permutation. If the input 
#' is a matrix of permutations, invert all the permutations in the input. 
#'
#' @param perm a permutation or matrix of permutations
#' @return The inverse permutation. 
#' If the input is a matrix,the matrix with the inverses
#' @export
#' @examples
#' inverse.perm(c(1,2,3,4))
#' inverse.perm(c(2,3,4,1))
#' data <- matrix(c(1,2,3, 4,1,4,3,2,1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
#' inverse.perm(data)
inverse.perm <- function(perm){
  if (is.permutation(perm))
    if (is.vector( perm ) ) 
      return (order(perm))
    else if (is.matrix(perm))
      return (t(apply(perm, MARGIN=1,FUN=order)))
  stop("Error in the input parameters: the input must be a permutation or a matrix of permutations")
}

#' Compose permutations
#'
#' This function composes two given permtuations. One of the arguments 
#' can be a collection of permutations, but not both at the same time. In this 
#' case, every permutation in the collection is composed with the 
#' other argument
#'
#' @param perm1 a permtuation or a collection of permutations
#' @param perm2 a permtuation or a collection of permutations
#' @return The composition of the permutations
#' @export
#' @examples
#' compose(c(3,1,2,4), c(4,1,3,2))
compose <- function(perm1, perm2){
  if (!is.permutation(perm1) || ! is.permutation(perm2))
    stop("Check parameters")
  if(is.matrix(perm1) && is.matrix(perm2)) 
    stop("Both params can not be samples\n") 
  if(is.vector(perm1) && is.vector(perm2) && .same.length.perms(perm1, perm2)) 
    return (.compose.perms(perm1,perm2))
  if(is.matrix(perm1) && dim(perm1)[2] == length(perm2) )
    return (t(apply(perm1, MARGIN=1,FUN=function(x){.compose.perms(x,perm2)})))
  if(is.matrix(perm2) && dim(perm2)[2] == length(perm1))
    return (t(apply(perm2, MARGIN=1,FUN=function(x){.compose.perms(perm1,x)})))
  stop ("Check input parameters")
}

.compose.perms <- function(s,p){
  return(s[p])
}

#' Check if its argument is a permutation
#'
#' This function tests if the given argument is a permutation of the first 
#' n natural integers (excluding 0)
#'
#' @param perm a vector (or a bidimensional matrix)
#' @return TRUE iff perm is a valid permutation (or a matrix of valid permutations)
#' @export
#' @examples
#' is.permutation(c(3,1,2,4))
#' is.permutation(c(6,1,2,3))
#' is.permutation(matrix(c(1,2,3, 4,1,4,3,2,1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE))
is.permutation <- function(perm){
  if (is.vector(perm)) return (.is.permutation.one(perm))
  else if (is.matrix(perm))
    return(all(apply(perm, MARGIN=1, .is.permutation.one)))
  else stop("Error in the input parameters: the input must be a permutation or a matrix of permutations")
}

.is.permutation.one <- function(perm){
  bid<-logical(length = length(perm)) #bid: logical vector, all false
  lapply(X=perm, function(x){
    bid[x] <<- TRUE
  })
  return(all(bid))
}

#' Read a text file with a collection of permtuations
#'
#' This function reads the text file in the specified path and 
#' checks if each row is a proper permutation
#'
#' @param path string with a path
#' @return A collection of permutations in matrix form
#' @export
#' @examples
#' path = system.file("test.txt", package="PerMallows")
#' sample = read.perms(path)
read.perms <- function(path){
  sample = as.matrix( read.table(path, sep=" "))
  if (! is.permutation(sample) )
    stop("Not valid permtuations in the matrix")
  return(sample)
}

#' Convert rating to permutation
#'
#' This function is given a collection of ratings and converts each row to
#' a permutation
#'
#' @param ratings a matrix in which each row is a vector of ratings of several items
#' @return A matrix in which each row is the corresponding permutation of the items
#' @export
#' @examples
#' order.ratings(c(0.1, 4, 0.5, -4))
order.ratings <- function(ratings){
  if (is.vector(ratings))
    return (order(ratings))
  return (t(apply(ratings, MARGIN=1, FUN=order)))
  stop("The input parameter must be a vector")
}

#' Random permutation
#'
#' Generate a collection of n permutations uniformly at random
#'
#' @param n optional number of permutations to generate
#' @param perm.length length of the permutations generated
#' @return A single permutation or a matrix with n rows, each being a permutation. 
#' Every permutation is drawn uniformly at random and has length perm.length
#' @export
#' @examples
#' runif.permutation(1,5)
runif.permutation <- function(n=1 , perm.length){
  x <- c(1:perm.length)
  if ( n == 1 ) {
    x <- .C("random_permutation", as.integer(perm.length), as.integer(x))[[2]]
    return (x)
  }
  return (t (replicate(n, .C("random_permutation", as.integer(perm.length), as.integer(x))[[2]])))
}

#' Compute the frequency matrix
#'
#' Compute the first order marginal probability. In other words, given at least one 
#' permutation, calculate the proportion of them that have each item in each position
#'
#' @param perm a permutation or a collection of them
#' @return A matrix with n rows and n columns with the proportion of the permutations 
#' in the input that have each item in each position
#' @export
#' @examples
#' freq.matrix(c(1,3,2,4,5))
freq.matrix <- function(perm){
  if( !is.permutation(perm) ) stop ("Check parameters")
  if ( is.vector(perm) ) {
    perm.length <- length(perm)
    #mat <- matrix(0, perm.length, perm.length)
    mat<-diag(perm.length)[,perm]
  }else if (is.matrix(perm)){
    perm.length <- dim(perm)[2]
    mat <- matrix(0, perm.length, perm.length)
    apply(X=perm, MARGIN=1, FUN=function(x){
      #mat <<- mat + diag(perm.length)[,x]   # fast
      for (i in 1:perm.length) mat[i,x[i]] <<- mat[i,x[i]] + 1 #low mem
    })
    mat = mat / dim(perm)[1]
  }else stop("Check arguments")
  return (mat)
}

###################   PROBABILITY DISTRIBUTIONS   ###################
.check.distance.name <- function(dist.name, for.gmm = FALSE){
  if ( for.gmm ) err.msg = "Choose one of these distances: Cayley, Kendall, Hamming"
  else err.msg = "Choose one of these distances: Cayley, Kendall, Hamming, Ulam"
  
  if      (dist.name == "c" || dist.name == "cayley"  || dist.name == "Cayley")  dist_id = .CAYLEY.DISTANCE
  else if (dist.name == "k" || dist.name == "Kendall" || dist.name == "kendall") dist_id = .KENDALL.DISTANCE
  else if (dist.name == "h" || dist.name == "Hamming" || dist.name == "hamming") dist_id = .HAMMING.DISTANCE
  else if ((dist.name == "u" || dist.name == "Ulam"    || dist.name == "ulam") &&  !for.gmm )    dist_id = .ULAM.DISTANCE
  else stop("nnot valid distance: ", dist.name, ". ",err.msg)
  return (dist_id)
}

.check.theta.length.gmm <- function(dist_id, theta.length, perm.length){
  if ( dist_id == 0 && ! ( theta.length == perm.length - 1))      
    stop ("Check Theta length, must be ", perm.length - 1, ".")
  else if ( dist_id == 1 && ! ( theta.length == perm.length - 1)) 
    stop ("Check Theta length, must be ", perm.length - 1, ".")
  else if ( dist_id == 2 && ! ( theta.length == perm.length ))    
    stop ("Check Theta length, must be ", perm.length, ".")
  return (theta.length == 1)
}

.get.theta.length.gmm <- function(dist.id, perm.length){
  if ( dist.id == .HAMMING.DISTANCE ) return (perm.length)
  if ( dist.id == .CAYLEY.DISTANCE || dist.id == .KENDALL.DISTANCE ) return (perm.length-1)
  stop ("No GMM for the Ulam distance")
}
#' Compute the marginal probability, GMM under the Hamming distance
#' 
#' Compute the marginal probability, GMM under the Hamming distance, 
#' of a distance decomposition vector for which some positions are known and some are not
#'
#' @param h n dimensional distance decomposition vector where h_j = 0 means that $j$ is a fixed point, 
#' h_j = 1 means that $j$ is an unfixed point and otherwise $j$ is not known
#' @param theta n dimensional distance decomposition vector with the dispersion parameters
#' @return The marginal probability
#' @export
#' @examples
#' marginal(c(1,0,1,NA,NA), c(0.1, 0.3, 0.7, 0.1, 1))
#' marginal(c(NA,0,1,NA,NA,0), c(0.1, 0.3, 0.7, 0.1, 0.7, 1))
marginal <- function(h, theta){
  if (length(h) != length(theta))
    stop("Check the length of the h and theta vectors")
  res <- 0
  h[is.na(h)] <- -1
  res<-.C("marginals", as.integer(length(h)), as.integer(.check.distance.name("hamming")), as.integer(h),
          as.numeric(theta), as.double(res))[5]
  return (res[[1]]) 
}

#' Compute the expected distance, MM under the Hamming distance
#' 
#' Compute the expected distance in the MM under the Hamming distance
#'
#' @param theta real dispersion parameter
#' @param perm.length length of the permutation in the considered model
#' @param dist.name optional name of the distance used in the MM. One of: kendall (default), cayley, hamming, ulam
#' @return The expected distance under the MM
#' @export
#' @examples
#' expectation.mm( 1, 7, "kendall" )
#' expectation.mm( 2, 5, "cayley" )
#' expectation.mm( 2, 4, "hamming" )
#' expectation.mm( 1, 6, "ulam" )
expectation.mm <- function(theta, perm.length,  dist.name="kendall"){
  dist_id = .check.distance.name(dist.name, FALSE)
  expec <- rep(0,perm.length)
  expec<-.C("expectation",  as.integer(dist_id), as.integer(.MALLOWS.MODEL), as.integer(perm.length), 
            as.double(theta),as.double(expec))[5]
  return (expec[[1]][1]) 
}

#' Compute the expected distance, GMM under the Hamming distance
#' 
#' Compute the expected distance in the GMM under the Hamming distance
#'
#' @param theta n dimensional real vector with the dispersion parameters
#' @return The expected distance decomposition vector under the GMM
#' @param dist.name optional name of the distance used in the GMM. One of: kendall (default), cayley, hamming
#' @export
#' @examples
#' expectation.gmm(c(0.38, 0.44, 0.1, 0.2, 1, 0.1))
#' expectation.gmm(c(2, 2, 2, 2),"cayley")
#' expectation.gmm(c(0.3, 0.1, 0.5, 0.1),"hamming")
expectation.gmm <- function (theta, dist.name="kendall"){
  dist_id = .check.distance.name(dist.name, TRUE)
  if ( dist_id == .HAMMING.DISTANCE )  perm.length <- length(theta) 
  else  perm.length <- length(theta) + 1
  expec <- rep(0,perm.length)
  expec<-.C("expectation",  as.integer(dist_id), as.integer(.GENERALIZED.MALLOWS.MODEL), as.integer(perm.length), 
           as.double(theta),as.double(expec))[5]
 return (expec[[1]] [1:length(theta)]) 
}

#' Calculate the probability of a permutation in a MM 
#'
#' Calculate the probability of a permutation sigma in a MM of center sigma0, 
#' dispersion parameter theta and under the specified distance
#'
#' @param perm permutation whose probability is asked for
#' @param sigma0 optional central permuation of the MM, by default the identity
#' @param theta dispersion parameter of the MM
#' @param dist.name optional name of the distance used in the MM. One of: kendall (default), cayley, hamming, ulam
#' @return The probability of sigma in the given MM
#' @export
#' @examples
#' data <- matrix(c(1,2,3, 4,1,4,3,2,1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
#' sig<-c(1,2,3,4)
#' log.prob <- apply(data,MARGIN=1,FUN=function(x){log(dmm(x,sig, 1,"cayley"))})
#' sum(log.prob)
#' dmm(c(1,3,2,4), theta=0.1)
#' dmm(c(1,3,2,4), theta=0.1, dist.name="cayley")
#' dmm(c(1,3,2,4), theta=0.1, dist.name="hamming")
#' dmm(c(1,3,2,4), theta=0.1, dist.name="ulam")
dmm <- function(perm, sigma0=identity.permutation(length(perm)), theta, dist.name="kendall"){
if (!is.permutation(perm) || ! is.permutation(sigma0))
    stop("Check parameters")
  dist_id       = .check.distance.name(dist.name)
  if( length(theta) == 1 ) theta = rep( theta, length(perm)  )
  else stop ("This function is for the MM (just one dispersion parameters)")
  return (.proba(perm, sigma0, theta, dist_id))
}

#' Calculate the probability of a permutation in a GMM 
#'
#' Calculate the probability of a permutation sigma in a GMM of center sigma0, 
#' dispersion parameter theta and under the specified distance
#'
#' @param perm permutation whose probability wants to be known
#' @param sigma0 central permuation of the GMM, by default the identity
#' @param theta vector dispersion parameter of the GMM
#' @param dist.name optional name of the distance used in the GMM. One of: kendall (default), cayley, hamming
#' @return The probability of sigma in the given GMM
#' @export
#' @examples
#' data <- matrix(c(1,2,3,4, 1,4,3,2, 1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
#' sig <- c(1,2,3,4)
#' th <- c(0.1, 0.2, 0.3,1)
#' log.prob <- apply(data,MARGIN=1,FUN=function(x){log(dgmm(x,sig, th, "hamming"))})
#' sum(log.prob)
#' dgmm (c(1,2,3,4), theta=c(1,1,1))
#' dgmm (c(1,2,3,4), theta=c(1,1,1), dist.name="cayley")
dgmm <- function(perm, sigma0=identity.permutation(length(perm)), theta, dist.name="kendall"){
  if (!is.permutation(perm) || ! is.permutation(sigma0))
    stop("Check parameters")
  dist_id       = .check.distance.name(dist.name, TRUE)
  mallows.model = .check.theta.length.gmm(dist_id, length(theta), length(perm)  )
  return (.proba(perm, sigma0, theta, dist_id))
}

.proba <- function (perm, sigma0, theta, dist_id){
    res = -1;
  if(.same.length.perms(perm,sigma0) ){
      res<-.C("probability", as.integer(dist_id), as.integer(length(perm)),as.integer(perm),
                as.integer(sigma0), as.numeric(theta), as.double(res))[6]
    }else stop ("The permutations must have the same number of items")
    return( unlist(res) );
}

#' Sample a Mallows Model
#'
#' Generate a sample of n permutations from a Mallows Model (MM). 
#'
#' @param n the number of permutations to be generated
#' @param sigma0 central permuation of the MM
#' @param theta dispersion parameter of the MM
#' @param dist.name optional name of the distance used in the MM. One of: kendall (default), cayley, hamming, ulam
#' @param sampling.method optional name of the sampling algorithm. One of: distances, multistage, gibbs (default)
#' @param disk optional can only be true if using the Distances sampling algorithm for generating under the Ulam distance.
#' Insted of generating the whole set of SYT and count of permutations per distance, it loads the info from a file in the disk
#' @param alert check consistency of the parameters. TRUE by default
#' @return A matrix contaning a sample of permutations from the specified ditribution
#' @export
#' @examples
#' rmm(2,c(1,2,3,4,5),1,"kendall", "distances")
#' rmm(2,c(1,2,3,4,5),1,"cayley", "distances")
#' rmm(2,c(1,2,3,4,5),1,"hamming", "distances")
#' rmm(2,c(1,2,3,4,5),1,"ulam", "distances")
#' rmm(2,c(1,2,3,4,5),1,"kendall", "multistage")
#' rmm(2,c(1,2,3,4,5),1,"cayley", "multistage")
rmm <- function(n, sigma0, theta, dist.name="kendall", sampling.method=NULL, disk=FALSE, alert = TRUE){
  if (! is.permutation(sigma0)) stop("Check parameters")
  num.perms <- n
  perm.length <- length(sigma0)
  if(length(theta) != 1 ) stop("Theta must be one real number") 
  dist_id       = .check.distance.name(dist.name)
  if (is.null(sampling.method)){
    if (dist_id == .ULAM.DISTANCE) sampling.method = 'distances'
    else sampling.method = 'multistage'
  }
  
  if      ( sampling.method == "d"  || sampling.method == "distances")  algorithm_id = 0
  else if ( sampling.method == "m"  || sampling.method == "multistage") algorithm_id = 1
  else if ( sampling.method == "g"  || sampling.method == "gibbs")      algorithm_id = 2
  else {
    if (dist_id == .ULAM.DISTANCE )     
      stop("The argument 'sampling.method' must be one of the following: 'distances', 'gibbs'.")
    else     
      stop("The argument 'sampling.method' must be one of the following: 'multistage', distances', 'gibbs'.")
  }
  

  if (alert){
    if (  algorithm_id == 0 &&  n > 150 ) {
      cat("The distances sampler is reliable for small permutations: 
        For permutations with more than 150 items proceed at your own risk. Do you want to continue? (Y/N)")
      answer = readline()
      if ( answer != 'Y' ) stop("Process cancelled")
    }
    if (algorithm_id == 0 && dist_id == .ULAM.DISTANCE) {
      if ( n > 80 && !disk){
        cat("The memory requirement can be big (it is possible to set the 'disk' to TRUE for a smaller memory req.). Do you want to continue? (Y/N)")
      #if ( n > 100)cat("The process can take . Do you want to continue? (Y/N)")
      answer = readline()
      if ( answer != 'Y' ) stop("Process cancelled")
      }
    }
    
  }

  if ( disk ){
    if ( algorithm_id == 0 && dist_id == .ULAM.DISTANCE){
      dist_id = 4
      if ( ! file.exists(paste('permus_per_dist_',perm.length, sep="")) )
        stop("Generate first the files. Try: ",paste ('generate.aux.files(',perm.length, ')'))
    }
    else stop ("Geneneration form disk can only be used with Ulam and Distances sampling method")
  }
  
  if( dist_id == .ULAM.DISTANCE && algorithm_id == 1 ) 
    stop ("No Multistage algorithm for the Ulam distance, choose the Distances sampling algorithm by setting the parameter <sampling.method='distances' > ")
  
  if (algorithm_id == 0)
    sam <- .Call("distances_sampling", dist_id, perm.length, num.perms, theta)
  else if (algorithm_id == 1 || algorithm_id == 2 ){
    theta <- rep(theta,perm.length)
    sam <- .Call("sampling_multi_gibbs_cayley",dist_id,  perm.length, num.perms, theta, .GENERALIZED.MALLOWS.MODEL, algorithm_id)
  }else stop("Choose one of these algorithms: distances, gibbs, multistage")
  if(!all(sigma0 == identity.permutation(perm.length)) ) 
    sam <- compose(sam, sigma0)
  return(sam)
}

#' Sample a Generalized Mallows Model
#'
#' Generate a sample of n permutations from a Generalized Mallows Model (GMM). 
#'
#' @param n the number of permutations to be generated
#' @param sigma0 central permuation of the GMM
#' @param theta dispersion parameter vector of the GMM
#' @param dist.name optional used name of the distance used in the GMM. One of: kendall (default), cayley, hamming
#' @param sampling.method optional name of the sampling algorithm. One of: multistage, gibbs (default)
#' @return A matrix contaning a sample of permutations from the specified ditribution
#' @export
#' @examples
#' rgmm(2,c(1,2,3,4,5),c(1,1,1,1),"kendall", "multistage")
#' rgmm(2,c(1,2,3,4,5),c(1,1,1,1),"cayley", "multistage")
#' rgmm(2,c(1,2,3,4,5),c(1,1,1,1,1),"hamming", "multistage")
#' rgmm(2,c(1,2,3,4,5),c(1,1,1,1),"cayley", "gibbs")
#' rgmm(2,c(1,2,3,4,5),c(1,1,1,1,1),"hamming", "gibbs")
rgmm <- function(n, sigma0, theta, dist.name="kendall", sampling.method="multistage"){
  if (! is.permutation(sigma0))
    stop("Check parameters")
  num.perms <- n
  perm.length <- length(sigma0)
  theta.length <- length(theta)
  
  dist_id       = .check.distance.name(dist.name, TRUE)
  mallows.model = .check.theta.length.gmm(dist_id, theta.length, perm.length)
  
  if ( sampling.method == "multistage" || sampling.method == "multistage")  algorithm_id = 1
  else if ( sampling.method == "gibbs"      || sampling.method == "gibbs")  algorithm_id = 2
  else stop("Choose one of these sampling.method: multistage, gibbs")
  
  sam <- .Call("sampling_multi_gibbs_cayley",dist_id, perm.length, num.perms, theta, .GENERALIZED.MALLOWS.MODEL, algorithm_id);
  #browser()
  if(!all(sigma0 == identity.permutation(perm.length)) ) 
    sam <- compose(sam, sigma0)
  return(sam)
}

#' Learn a Mallows Model
#'
#' Learn the parameter of the distribution of a sample of n permutations 
#' comming from a Mallows Model (MM). 
#'
#' @param data the matrix with the permutations to estimate
#' @param sigma_0_ini optional the initial guess for the consensus permutation
#' @param dist.name optional the name of the distance used by the model. One of: kendall (default), cayley, hamming, ulam
#' @param estimation optional select the approximated or the exact. One of: approx, exact
#' @param disk optional can only be true if estimating a MM under the Ulam distance.
#' Insted of generating the whole set of SYT and count of permutations per distance, it loads the info from a file in the disk
#' @return A list with the parameters of the estimated distribution: the
#' mode and the dispersion parameter
#' @export
#' @examples
#' data <- matrix(c(1,2,3,4, 1,4,3,2, 1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
#' lmm(data, dist.name="kendall", estimation="approx")
#' lmm(data, dist.name="cayley", estimation="approx")
#' lmm(data, dist.name="cayley", estimation="exact")
#' lmm(data, dist.name="hamming", estimation="exact")
#' lmm(data, dist.name="ulam", estimation="approx")
lmm <- function(data, sigma_0_ini =identity.permutation(dim(data)[2]), dist.name="kendall", estimation="approx", disk=FALSE){
  if ( !is.permutation(data)) stop("Check parameters")
  if ( !.same.length.perms(data[1,],sigma_0_ini)) 
    stop ("Permutation sigma_0_ini must have the same number of items as the permutations in the data")
  num.perms <- dim(data)[1]
  perm.length <- dim(data)[2]
  dist_id = .check.distance.name(dist.name)
  
  if(estimation == "Exact" || estimation == "exact") estim_var = 0
  else if(estimation == "Approx" || estimation == "approx") estim_var = 1
  else stop("Choose one of these estimations: approx, exact")
  
  if ( disk ){
    if ( dist_id == 3){
      dist_id = 4
      if ( ! file.exists(paste('permus_per_dist_',perm.length, sep="")) )
        stop("Generate first the files. Try: ",paste ('generate.aux.files(',perm.length, ')'))
    }
    else stop ("Estimation form disk can only be used with Ulam distance")
  }
  
  if (estim_var == 1 && dist_id == .HAMMING.DISTANCE ) stop ("No approx learning for the Hamming MM, try exact")
  if (estim_var == 0 && (dist_id != .CAYLEY.DISTANCE && dist_id != .HAMMING.DISTANCE) ) 
    stop ("Exact learning only for cayley and Hamming")
  sigma0 <- .Call("consensus", dist_id, data, 0, estim_var, sigma_0_ini)  
  
  the <- .Call("estimate_theta", dist_id, perm.length, num.perms, sigma0, data, 0)
  return(list(mode = sigma0, theta = the[1]));
}

#' Learn a Generalized Mallows Model
#'
#' Learn the parameter of the distribution of a sample of n permutations 
#' comming from a Generalized Mallows Model (GMM). 
#'
#' @param data the matrix with the permutations to estimate
#' @param sigma_0_ini optional the initial guess for the consensus permutation
#' @param dist.name optional name of the distance used by the GMM. One of: kendall (default), cayley, hamming
#' @param estimation optional select the approximated or the exact. One of: approx, exact
#' @return A list with the parameters of the estimated distribution: the
#' mode and the dispersion parameter vector
#' @export
#' @examples
#' data <- matrix(c(1,2,3,4, 1,4,3,2, 1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
#' lgmm(data, dist.name="kendall", estimation="approx")
#' lgmm(data, dist.name="cayley", estimation="approx")
#' lgmm(data, dist.name="cayley", estimation="exact")
#' lgmm(data, dist.name="hamming", estimation="approx")
lgmm <- function(data, sigma_0_ini =identity.permutation(dim(data)[2]), dist.name="kendall", estimation="approx"){
  if (! is.permutation(data)) stop("Check parameters")
  if(! .same.length.perms(data[1,], sigma_0_ini)) 
    stop ("Permutation sigma_0_ini must have the same number of items as the permutations in the data")
  num.perms <- dim(data)[1]
  perm.length <- dim(data)[2]
  dist_id = .check.distance.name(dist.name , TRUE)
  
  if(estimation == "Exact" || estimation == "exact") estim_var = 0
  else if(estimation == "Approx" || estimation == "approx") estim_var = 1
  else stop("Choose one of these estimations: approx, exact")
  
  if (estim_var == 0 && dist_id != .CAYLEY.DISTANCE ) 
    stop ("Exact leraning only for Cayley")
  print("desde r1 ")
  print(sigma_0_ini)
  sigma0 <- .Call("consensus", dist_id, data, 1, estim_var, sigma_0_ini)
  print("desde r 2")
  print(sigma0)
  the <- .Call("estimate_theta", dist_id, perm.length, num.perms, sigma0, data, 1)
  if ( dist_id == 2 ) theta.len = perm.length
  else  theta.len = perm.length - 1
  return(list(mode = sigma0, theta = the[1:theta.len]));
}

#' MLE for theta - Mallows Model
#'
#' Compute the MLE for the dispersion parameter (theta) given a sample of n permutations
#' and a central permutation
#'
#' @param data the matrix with the permutations to estimate
#' @param sigma_0 optional the consensus permutation. If not given it is assumed to be the identity permutation
#' @param dist.name optional the name of the distance used by the model. One of: kendall (default), cayley, hamming, ulam
#' @param disk optional can only be true if estimating a MM under the Ulam distance.
#' Insted of generating the whole set of SYT and count of permutations per distance, it loads the info from a file in the disk
#' @return The MLE for the dispersion parameter
#' @export
#' @examples
#' data <- matrix(c(1,2,3,4, 1,4,3,2, 1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
#' lmm.theta(data, dist.name="kendall")
#' lmm.theta(data, dist.name="cayley")
#' lmm.theta(data, dist.name="cayley", sigma_0=c(1,4,3,2))
#' lmm.theta(data, dist.name="hamming")
#' lmm.theta(data, dist.name="ulam")
lmm.theta <- function(data, sigma_0 =identity.permutation(dim(data)[2]), dist.name="kendall",
                         disk=FALSE){
  if (! is.permutation(data))  stop("Check parameters")
  if(! .same.length.perms(data[1,], sigma_0)) 
    stop ("Permutation sigma_0 must have the same number of items as the permutations in the data")
  num.perms <- dim(data)[1]
  perm.length <- dim(data)[2]
  dist_id = .check.distance.name(dist.name)
  
  if ( disk ){
    if ( dist_id == 3){
      dist_id = 4
      if ( ! file.exists(paste('permus_per_dist_',perm.length, sep="")) )
        stop("Generate first the files. Try: ",paste ('generate.aux.files(',perm.length, ')'))
    }
    else stop ("Estimation form disk can only be used with Ulam distance")
  }
  
  the <- .Call("estimate_theta", dist_id, perm.length, num.perms, sigma_0, data, 0)
  return(theta = the[1]);
}

#' MLE for theta - Generalized Mallows Model
#'
#' Compute the MLE for the dispersion parameter (theta) given a sample of n permutations
#' and a central permutation
#'
#' @param data the matrix with the permutations to estimate
#' @param sigma_0 optional the initial guess for the consensus permutation. If not given it is assumed to be the identity permutation
#' @param dist.name optional name of the distance used by the GMM. One of: kendall (default), cayley, hamming
#' @return The MLE for the dispersion parameter
#' @export
#' @examples
#' data <- matrix(c(1,2,3,4, 1,4,3,2, 1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
#' lgmm.theta(data, dist.name="kendall")
#' lgmm.theta(data, dist.name="cayley")
#' lgmm.theta(data, dist.name="cayley", sigma_0=c(1,4,3,2))
#' lgmm.theta(data, dist.name="hamming")
lgmm.theta <- function(data, sigma_0 =identity.permutation(dim(data)[2]), dist.name="kendall"){
  if (! is.permutation(data))  stop("Check parameters")
  if(! .same.length.perms(data[1,], sigma_0)) 
    stop ("Permutation sigma_0 must have the same number of items as the permutations in the data")
  num.perms <- dim(data)[1]
  perm.length <- dim(data)[2]
  dist_id = .check.distance.name(dist.name , TRUE)
  
  the <- .Call("estimate_theta", dist_id, perm.length, num.perms, sigma_0, data, 1)
  if ( dist_id == 2 ) theta.len = perm.length
  else  theta.len = perm.length - 1
  return(theta = the[1:theta.len]);
}

#' Generates the files for Ulam
#'
#' Generates files for Ulam which are aimed to accelelrate the processes of counting 
#' the number of permutations at ech distance, sampling and learning IFF these operations
#' are going to be computted more than once
#'
#' @param perm.length number of items in the permutations
#' @return Nothing. Only writes in the current folder the axiliary files
#' @export
#' @examples
#' generate.aux.files(4)
generate.aux.files <- function(perm.length){
  .C("save_counts_to_files", as.integer(perm.length))
}

###################   DISTANCE RELATED   ###################

#' Swap two items of a permutation
#'
#' Given a permutation and two position, swap both positions
#'
#' @param perm a permutation
#' @param i position of the permutation
#' @param j position of the permutation
#' @return The permutation in the input in which the two speicfied items have been swapped
#' @export
#' @examples
#' swap(c(1,2,3,4,5),2,5)
swap <- function(perm, i, j){
  if (!.is.permutation.one(perm)) stop("The input parameter perm must be a valid permutation")
  n <- length(perm)
  if ( i < 1 || i > n || j < 1 || j > n ) stop ("Swap only possible at positions 1..",n)
  #if (i == j ) stop ("Set different values for i and j")
  aux <- perm[i]
  perm[i] <- perm[j]
  perm[j] <- aux
  return (perm)
}

#' Inversion operator
#'
#' Given a permutation and a position, swap positions i and i+1
#'
#' @param perm a permutation
#' @param i position of the permutation
#' @return The permutation in the input with an inversion at the specified position
#' @export
#' @examples
#' inversion(c(1,2,3,4,5),2)
inversion <- function(perm, i){
  if (!.is.permutation.one(perm)) stop("The input parameter perm must be a valid permutation")
  n <- length(perm)
  if ( i < 1 || i > n-1 ) stop ("Inversion only possible at positions 1..",n-1)
  aux <- perm[i]
  perm[ i ] <- perm[i+1]
  perm[ i + 1 ] <- aux
  return (perm)
}

#' Insert operator
#'
#' Given a permutation and two positions i, j, move item in position i to position j
#'
#' @param perm a permutation
#' @param i position of the permutation
#' @param j position of the permutation
#' @return The permutation in the input in which th eoperation has been applied
#' @export
#' @examples
#' insert(c(1,2,3,4,5),5,2)
#' insert(c(1,2,3,4,5),2,5)
insert <- function(perm, i, j){
  #remove item at position i and place it after the j-th position 
  if (!.is.permutation.one(perm)) stop("The input parameter perm must be a valid permutation")
  n <- length(perm)
  if(i<1 || i>n || j < 1 || j >n)  stop ("Insertion only possible at positions 1..",n)
  if (i < j)
    sigma <- c(perm[1:j][-i], perm[i], perm[(j+1):n])
  else 
    sigma <- c(perm[1:j], perm[i], perm[(j+1):n][-(i-j)])
  return(sigma)
  item <- perm[ i ]
  mod  <- perm[-i ]
  if(i<j) j <- j-1
  if(j<1)
    sigma <- c(item, mod[(j+1):length(mod)])
  else if((j+1)>length(mod)) 
    sigma <- c(mod[1:j], item)
  else 
    sigma <- c(mod[1:j], item, mod[(j+1):length(mod)])
  return (sigma)
}

#' Decompose a permutation in a set of cycles
#'
#' Factor a given a permutation in the set of independent cycles
#'
#' @param perm a permutation
#' @return The permutation in the input in which the operation has been applied
#' @export
#' @examples
#' perm2cycles(c(1,5,2,3,4))
perm2cycles <- function (perm){
  if (!is.vector(perm) || !is.permutation(perm))
    stop("Check parameters")
  s <- perm
  perm.length <- length(s)
  cont <- 0
  cont_cycle <- 0 
  while (cont < perm.length){
    cont_cycle <- cont_cycle + 1
    current <- 1
    while ( s[ current ] == 0 ) current <- current + 1
    if(cont_cycle == 1) cycles <- list (c(s[ current ]))
    else cycles <- append(cycles, c(s[ current ]))
    repeat{
      cont <- cont + 1
      next_i <- s[ current ]
      s[ current ] <- 0
      current <- next_i
      if( s[ current ] != 0)
        cycles[[cont_cycle]] = c(cycles[[cont_cycle]],  s[ current ])
      else {
        break
      }
    }
  }
  return (cycles)
}

#' Friendly display the cycles
#'
#' Given a list with the cycles of a permutation, displays them in the standard cycle notation
#'
#' @param cy a list with the set of cycles
#' @export
#' @examples
#' cycle2str(perm2cycles(c(1,5,2,3,4)))
cycle2str <- function(cy){
  invisible(lapply(cy, FUN=function(x){
    cat("(") 
    cat(x)
    cat(")")
    return()
  }))
}
#' Get the permutation given the cycles
#'
#' Get the permutation as a vector given the set of cycles in which it factorizes
#'
#' @param cycles a list with the set of disjoint cycles
#' @return The permutation in vector notation
#' @export
#' @examples
#' cycles2perm(perm2cycles(c(1,5,2,3,4)))
cycles2perm <- function(cycles){
  num_cycles <- length(cycles)
  perm <- c()
  for(i in 1:num_cycles){
    cycle_len <- length(cycles[[i]])
    for(j in 1:cycle_len){
      if ( j== cycle_len) perm[cycles[[i]][j]]=cycles[[i]][1]
      else perm[cycles[[i]][j]]=cycles[[i]][j+1]
    }
  }
  if (!.is.permutation.one(perm)) stop("The cycles do not correspond to a valid permutation")
  return (perm)
}

#' Compute the distance between permutations
#'
#' Compute the distance between two given permutations. If only one permutation is 
#' given the other one is assumed to be the identity (1,2,3,....,n)
#' The distance can be kendall, cayley, hamming and ulam
#'
#' @param perm1 a permutation
#' @param perm2 optional a permutation
#' @param dist.name optional. One of: kendall (default), cayley, hamming, ulam
#' @return The distance between the permutations 
#' @export
#' @examples
#' distance(c(1,2,3,5,4))
#' distance(c(1,2,3,5,4), c(1,2,3,5,4))
#' distance(c(1,2,3,5,4), c(1,4,2,3,5), "cayley")
distance<-function(perm1,perm2=identity.permutation(length(perm1)), dist.name="kendall"){
  dist <- 0;
  perm.length<-length(perm1);
  if(! .same.length.perms(perm1,perm2))
    stop ("The permutations must have the same number of items")
  dist_id = .check.distance.name(dist.name)
  if( dist.name == "Hamming" || dist.name == "hamming" ) {
    sigma <- .compose.perms(perm1, inverse.perm(perm2))
    for(i in 1:perm.length)
      if (sigma[i] != i ) dist <- dist + 1
  }else 
    dist<-.C("compute_distance", as.integer(dist_id), as.integer(perm.length), 
             as.integer(perm1), as.integer(perm2), as.integer(dist))[5]
  return(unlist(dist));
}

#' Get the maximum value of the distance ebtween permutations
#'
#' Compute the maximum posible value for the distance between two given permutations.
#' The distance can be kendall, cayley, hamming and ulam
#'
#' @param perm.length number of items in the permutations
#' @param dist.name optional. One of: kendall (default), cayley, hamming, ulam
#' @return The maximum value for the distance between the permutations 
#' @export
#' @examples
#' maxi.dist(4,"cayley")
#' maxi.dist(10,"ulam")
#' maxi.dist(4)
maxi.dist <- function(perm.length, dist.name="kendall"){
  dist_id = .check.distance.name(dist.name)
  if (dist_id == .CAYLEY.DISTANCE || dist_id == .ULAM.DISTANCE ) return (perm.length - 1 )
  if (dist_id == .HAMMING.DISTANCE ) return ( perm.length )
  if (dist_id == .KENDALL.DISTANCE ) return ( perm.length * ( perm.length - 1 ) / 2 ) 
  return (-1)
}

#' Count permutations at a distance
#'
#' Given a distance (kendall, cayley, hamming or ulam), 
#' the number of items in the permutations perm.length and distance value d, 
#' how many permutations are there at distance d from any permutation?
#' It can be used to count the number of derangements and the permutations 
#' with k cycles (Stirling numbers of the first kind)
#'
#' @param perm.length number of items in the permutations
#' @param dist.value the distance
#' @param dist.name optional. One of: kendall (default), cayley, hamming, ulam
#' @param disk optional can only be true if counting the permutations at each Ulam distance.
#' Insted of generating the whole set of SYT and count of permutations per distance, it loads the info from a file in the disk
#' @return The number of permutations at the given distance
#' @export
#' @examples
#' count.perms(4,2,"kendall")
#' count.perms(4,2,"ulam")
#' count.perms(4,2,"hamming")
#' count.perms(4,2,"cayley")
#' # The number of derangements of length 6 is computed as follows
#' len <- 6
#' count.perms(perm.length = len, dist.value = len, dist.name = "h") 
#' # The number of permutations with one cycle is computed as follows
#' num.cycles <- 1 
#' count.perms(perm.length = len, dist.value = len - num.cycles, dist.name = "c") 
count.perms <- function(perm.length, dist.value, dist.name="kendall", disk=FALSE){
  if (dist.value < 0 ) stop("The distance must greater than or equal to 0")
  count <- 0
  dist_id = .check.distance.name(dist.name)
  if ( disk ){
    if (dist_id == 3){
      dist_id = 4
      if ( ! file.exists(paste('permus_per_dist_',perm.length, sep="")) )
        stop("Generate first the files. Try: ",paste ('generate.aux.files(',perm.length, ')'))
    }
    else stop ("Counting form disk can only be used with Ulam distance")
  }
  count <- .C("count_permus_at_dist", as.integer(dist_id), as.integer(perm.length),as.integer( dist.value ),as.double(count))[4]
  return (count[[1]]) ## = unlist(count) 
}

#' Generate a collection of permutations at a given distance
#'
#' Given a number of permutations, the number of items in the permutations,
#' a distance value and a distance name, generate a sample of permutations with 
#' the specified length at the given distance.
#' Can be used to generate derangements and permutations of a given number of cycles
#'
#' @param n number of permutations in the sample
#' @param perm.length number of items in the permutations
#' @param dist.value distance value
#' @param dist.name distance name. One of: kendall (default), cayley, hamming, ulam
#' @return A sample of permutations at the given distance
#' @export
#' @examples
#' rdist(1, 4, 2 ) 
#' rdist(1, 4, 2, "ulam")
#' len <-  3
#' rdist(n = 1, perm.length = len, dist.value = len, "h") #derangement
#' cycles <- 2
#' rdist(n = 1, perm.length = len, dist.value = len - cycles, "c") #permutation with 2 cycles
rdist <- function(n, perm.length, dist.value, dist.name="kendall"){
  num.perms <- n
  dist_id = .check.distance.name(dist.name)
  if ( dist.value < 0 ) stop("The distance must be greater than 0")
  if ( dist.value > maxi.dist(perm.length,dist.name) )
    stop("The distance must be smaller than the largest
         possible ",dist.name," distance for permutations of ",perm.length,
          " items, which is ",maxi.dist(perm.length, dist.name))
  if (dist_id == .HAMMING.DISTANCE && dist.value == 1 ) return (0)
  res<-.Call("get_random_sample_at_dist_d", as.integer(dist_id), 
             as.integer(perm.length),as.integer(num.perms),as.integer( dist.value ))
  return(res)
}

#' Get the decomposition vector
#'
#' Given a permutation and a distance name generate the decomposition vector 
#'
#' @param perm the permutation
#' @param dist.name optional the name of the distance. One of: kendall (default), cayley, hamming
#' @return The distance decomposition vector of the given permutation and distance. For the Kendall distance is the inversion vector
#' @export
#' @examples
#' perm2decomp(c(1,2,4,3,5), "kendall")
#' perm2decomp(c(1,2,4,3,5), "cayley")
#' perm2decomp(c(1,2,4,3,5), "hamming")
perm2decomp <- function(perm, dist.name="kendall"){
  if (!.is.permutation.one(perm)) stop("The input parameter perm must be a valid permutation")
  perm.length<-length(perm)
  dist_id = .check.distance.name(dist.name, TRUE)
  if (dist_id == .HAMMING.DISTANCE) {
    vec <- rep(0,perm.length)
    for (i in 1:perm.length) if (perm[i] != i) vec[i] = 1
    return(vec)
  }
  vec <- c(1:perm.length) #only the first (perm.lenth-1) pos will be returned
  #the last equals 0
  res <- .C("get_altern_repre_for_permu", as.integer(dist_id), as.integer(perm.length), 
            as.integer(perm), as.integer(vec))[4]
  return(unlist(res)[ 1 : perm.length - 1 ] )
}

#' Get a permutation consistent with a decomposition vector
#'
#' Given a distance decomposition vector and a distance name, generate uniformly at random 
#' a permutation consistent with the decomposition vector. 
#'
#' @param vec the permutation
#' @param dist.name optional the name of the distance. One of: kendall (default), cayley, hamming
#' @return The distance decomposition vector of the given permutation and distance
#' @export
#' @examples
#' decomp2perm(c(1,0,1,0,0), "kendall")
#' decomp2perm(c(1,0,1,0,0), "cayley")
#' decomp2perm(c(1,0,1,0,0), "hamming")
decomp2perm <- function(vec, dist.name="kendall"){
  dist_id = .check.distance.name(dist.name, TRUE)
  if (dist_id == .CAYLEY.DISTANCE || dist_id == .HAMMING.DISTANCE )
    if (any(sapply(X = vec, function(x){(x!=0 && x !=1) })))
      stop("The input vector does not correspond to a valid permutation, must be binary")
  if(dist_id == .KENDALL.DISTANCE 
     && any(sapply(X = 1:length(vec), function(x){(vec[x] < 0 ||vec[x] > length(vec)+1-x) }))) #length(sigma)=length(vec)+1
    stop("The input vector does not correspond to a valid permutation")
  perm.length <- length(vec) 
  if (dist_id == 0 || dist_id == 1) {
    perm.length = perm.length + 1
    vec[ perm.length ] = 0 
  }
  sigma <- rep(0,perm.length)
  sigma <- .C("get_permu_given_altern_repre", as.integer(dist_id), 
              as.integer(perm.length), as.integer(vec), as.integer(sigma))[4]
  return(unlist(sigma))
}

###################   DATASET   ###################
#' @name perm.sample.med
#' @docType data
#' @title Sample of permutations
#' @description
#' A rda file containing a sample of permutations
#' @format
#' Each row is a permtuation
NULL
#' @name perm.sample.small
#' @docType data
#' @title Sample of permutations
#' @description
#' A rda file containing a sample of permutations
#' @format
#' Each row is a permtuation
NULL
#' @name data.apa
#' @docType data
#' @title Sample of permutations APA
#' @description
#' A rda file containing a sample of permutations of the American Psychology Association
#' @format
#' Each row is a permtuation
NULL
#' @name data.order
#' @docType data
#' @title Sample of permutations
#' @description
#' A rda file containing a sample of permutations
#' @format
#' Each row is a permtuation
NULL

