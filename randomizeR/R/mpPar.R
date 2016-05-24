###############################################
# --------------------------------------------#
# Class mpPar                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the mpPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatempPar <- function(object) {
  errors <- character()
  mti <- object@mti
  lengthMTI <- length(mti)
  N <- object@N
  ratio <- object@ratio

  if(lengthMTI != 1) {
    msg <- paste("mti has length  ", lengthMTI, ". Should be length one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if (round(mti[1]) != mti) {
    msg <- paste("First element of mti is ", mti, ". Should be an integer.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(mti[1] < 0){
    msg <- "mti must be a positive integer"
    errors <- c(errors, msg)
  }
  
  if(length(N) == 1 && !(N %% sum(ratio) == 0)) {
    msg <- paste("N = ", N, " is not a multiple of sum(ratio) = "
                 , sum(ratio), ".", sep = "")
    errors <- c(errors, msg)
  }

  if(length(ratio) == 1 && length(N) == 1 && sum(ratio) <= N) {
    msg <- paste("Sum of ratio must be smaller than N = ", N, ".",
                 sep = "")
    errors <- c(errors, msg)
  }

  if(length(N) == 1 && length(mti) == 1 &&  !(max(ratio)/sum(ratio) * N >= mti)) {
    msg <- paste("mti has to be smaller than ", max(ratio)/sum(ratio) * N, ".",
                 sep = "")
    errors <- c(errors, msg)
  }

  

  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for mpPar
# --------------------------------------------

# Randomization paramters generic
setClass("mpPar",
         slots = c(mti = "numeric"),
         contains = "randPar",
         validity = validatempPar)


# --------------------------------------------
# Constructor function for mpPar
# --------------------------------------------

#' Representing Maximal Procedure
#' 
#' Represents the Maximal Procedure.
#'
#'
#' @details
#' Fix the total sample size \code{N} and the \code{mti}. Afterwards, the patients
#' are assigend to each treatment arm according to the \code{ratio}.
#' All randomization sequences are equiprobable.
#'
#' @family randomization procedures
#' 
#' @inheritParams overview
#' 
#' @return
#' \code{S4} object of the class \code{mpPar}.
#'
#' @export
#'
#' @references
#' V.W. Berger, A. Ivanova and M.D. Knoll (2003) Minimizing predictability while
#' retaining balance through the use of less restrective randomization
#' procedures. \emph{Statistics in Medicine}, \strong{19}, 3017-28. 
mpPar <- function(N, mti, ratio = c(1, 1), groups = LETTERS[1:2]) {
  new("mpPar", N = N, mti = mti, K = 2, ratio = ratio, groups = groups)
}


# --------------------------------------------
# Sampling algorithm for MP
# --------------------------------------------

# Maximal Procedure
#
# Computes a randomization sequence using a flexible version of the maximal
# procedure. The function permits the number of people to be randomise
# to group A to be any ratio of the total number \code{sum(bc)} of people to
# be allocated. Additionally, a maximum final imbalance \code{finImb=N_A-N_B}
# can be specified.
#
# @inheritParams overview
# 
# @return A vector with the allocation sequence for a clinical trial. 
# It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
# is allocated to treatment A (resp. B).
# 
# @references Salama et al: Efficient generation of constrained block
# allocation sequences. Wiley Online Library (2007). DOI: 10.1002/sim.3014
# 
# @export
getMPRand <- function(S, ratio){ 
  # compute path in graph
  p <- numeric(0)
  for(i in 1:(nrow(S)-1)) {
    p <- c(p, rbinom(1, 1, S[i+1, i-sum(p)]/S[i, i-sum(p)]))
  }
  # check if experimental group is biggest group
  if((which(ratio %in% max(ratio)))[1] == 1) {
    return(p)
  } else {
    return(1-p)
  }
}

# Create matrix for the Maximal Procedure
#
# Creates a matrix with the number of paths from the origin to the endpoint.
# The rows refer to the N steps (i.e. where N is the number of patients to be
# allocated) of the algorithm, where the colums represent the patients to
# allocated to group A (= group 0)
#
# @inheritParams overview
# 
# @return Matrix with number of paths to Well as coordinates.
# 
# @export
createMPMatrix <- function(N, MTI, ratio, FTI = 0) {
  stopifnot(is.numeric(N), length(N)==1,
            round(MTI) == MTI, FTI <= MTI, MTI > 0, 
            all(ratio > 0), all(ratio == round(ratio)),
            MTI <= ratio[1]/sum(ratio)*N, 
            MTI <= ratio[2]/sum(ratio)*N)
  
  bigGrp <- (which(ratio %in% max(ratio)))[1]
  ratio <- ratio[bigGrp]/sum(ratio) # calculate internally with N_E/N
  
  # de facto: number of people in experimental group
  N0 <- round(ratio*N)
  N1 <- N-N0 
  # >= 1 = actual ratio between the groups A and B
  r <- N1/N0 
  
  # corresponds to colums = people in A
  y <- 0:floor(N0+1+MTI/2)
  # corresponds to rows = total number of people at time i
  H <- matrix(0:N, N+1, N+1) 
  
  S <- t(apply(H[ ,1:(N0+2+floor(MTI/2))], 1, feasibleSeqs, y, MTI, r, N0, N1,
             FTI)) 
  
  index <- which(S)
  #indexes of the last line
  lastRow <- (1:length(S))[(1:length(S)) %% (N+1) == 0]
  # take out indexes of last line
  index<-index[!(index %in% lastRow)]
  #rev <- index[(length(index)-1):1]
  rev <- index[length(index):1]
  for(i in rev) S[i] <- countPathsToWell(S, i, N)
  S
}

# Compute feasible sequences for the Maximal Procedure
#
# This function yields the possible paths in the graph of the Maximal
# Procedure. When a node of the graph is reachable this function will
# mark the respective coordinate of the matrix with TRUE.
#
# @inheritParams overview
# @param x coordinate of the row.
# @param y coordinate of the column.
# @param r ratio of N0/N1 = NA/NB > 1.
# @param N0 number of people in group A.
# @param N1 number of people in group B.
# 
# @return TRUE if the node x, y is feasible.
feasibleSeqs <- function(x, y, MTI, r, N0, N1, FTI) {
  (abs((x-y)-r*y) <= MTI) & (x-y >= 0) & (y <= N0+floor(FTI/2)) &
  (abs(x-y) <= N1+floor(FTI/2))
}

# Count feasible sequences for the Maximal Procedure
#
# This function counts the possible paths in the graph of the Maximal
# Procedure. This is done by adding up the respective values of the matrix
# that includes the possible sequences.
#
# @inheritParams overview
# @param x index of possible nodes in S.
# @param N number of people to be allocated.
# 
# @return total number of successors to node x.
# 
countPathsToWell <- function(S, x, N) {
  as.numeric(S[x+1] + S[x+N+2])
}


# --------------------------------------------
# Methods for mpPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "mpPar"),
          function(obj) {
            if(obj@K != 2 || !identical(obj@ratio, c(1,1))) {
              stop("Only possible for K equals 2 and ratio corresponds to c(1,1).")
            }  
            # want to modify ratio for fixed value? or add ratio to mpPar
            allSeqs <- compltSet(obj)
            inside <- apply(allSeqs,1, function(x,MTI) {
              all(abs(cumsum(2*x-1)) <= MTI) & 2*sum(x) == length(x)
              }, MTI = mti(obj))
            new("mpSeq",
                M = allSeqs[inside, ],
                mti = mti(obj),
                N = N(obj),
                K = K(obj),
                ratio = ratio(obj),
                groups = obj@groups)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "mpPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            if(K(obj)>2) stop("MP: K>2 not available.")
            S <- createMPMatrix(N(obj), mti(obj), ratio = ratio(obj))
            new("rMpSeq", 
                M = t(sapply(1:r,function(x) getMPRand(S, ratio(obj)))), 
                mti = mti(obj), 
                N = N(obj),
                K = K(obj),
                ratio = ratio(obj),
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "mpPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            if(K(obj) > 2) stop("MP: K>2 not available.")
            # caution: take care of what group is the bigger one
            S <- createMPMatrix(N(obj), mti(obj), ratio = ratio(obj))
            new("rMpSeq", 
                M = t(getMPRand(S, ratio(obj))), 
                mti = mti(obj), 
                N = N(obj),
                K = K(obj),
                ratio = ratio(obj),
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "mpPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            if(K(obj) > 2) stop("MP: K>2 not available.")
			      S <- createMPMatrix(N(obj), mti(obj), ratio = ratio(obj))
            new("rMpSeq", 
                M = t(sapply(1:r,function(x) getMPRand(S, ratio(obj)))), 
                mti = mti(obj), 
                N = N(obj),
                K = K(obj),
                ratio = ratio(obj),
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "mpPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            if(K(obj) > 2) stop("MP: K>2 not available.")
            # caution: take care of what group is the bigger one
            S <- createMPMatrix(N(obj), mti(obj), ratio = ratio(obj))
            new("rMpSeq", 
                M = t(getMPRand(S, ratio(obj))), 
                mti = mti(obj), 
                N = N(obj),
                K = K(obj),
                ratio = ratio(obj),
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "mpPar"),
          function(obj) {
            paste("MP(", obj@mti, ")", sep = "")
          }
)
