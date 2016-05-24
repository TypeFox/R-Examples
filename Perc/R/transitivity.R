#' calculate transitivity measurements for a matrix
#' 
#' \code{transitivity} calculate transitivity measurements for a matrix 
#' 
#' @param conf an N-by-N conflict matrix whose \code{(i,j)}th element is the number of times \code{i} defeated \code{j}
#' 
#' @param strict a logical vector of length 1 (TRUE or FALSE). It is used in transitivity definition for alpha estimation. 
#' It should be set to TRUE when a transitive triangle is defined as all pathways in the triangle go to the same direction;
#' it should be set to FALSE when a transitive triangle is defined as PRIMARY pathways in the triangle go to the same direction.
#' Strict = FALSE by default.
#'
#' @return A list of four elements.
#' 
#'  \item{transitive}{The number of transitive triangles.}
#'  
#'  \item{intransitive}{The number of intransitive triangles.}
#'  
#'  \item{transitivity}{transitivity, the proportion of transitive triangles.}
#'  
#'  \item{alpha}{The value of alpha corresponding to this value of transitivity.}
#' 
#' 
#' @details \code{transitivity} is calculated as the proportion transitive triangles in the total of transitive and intransitive triangles.
#' transitivity is used to estimate alpha, which is used in turn in imputing information from indirect pathways as to what degree we can trust information from indirect pathways.
#' Greater transitivity is associated with assigning higher weight to information from indirect pathways. 
#' 
#' @seealso \code{\link{countPaths}}, \code{\link{findIDpaths}}, \code{\link{conductance}} 
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # transitivity calculation
#' conftrans <- transitivity(confmatrix, strict = FALSE)
#' conftrans$transitive
#' conftrans$intransitive
#' conftrans$transitivity
#' conftrans$alpha
#' @export

transitivity = function(conf, strict = FALSE){
  
  # making sure conf is of conf.mat
  if (!("conf.mat" %in% class(conf))){
    conf = as.conflictmat(conf)
  }
  
  
  N = nrow(conf)
  
  
  ### These lines set up the transitivity calculation.
  ### We are making a matrix of all the possible sets of three subjects.
  ### We won't really need to refer to this matrix.
  ### The important part of the calculation is at the end.
  
  numrows = 0
  ctr = 0
  for(i in 1:(N-2)){
    ctr = ctr + i
    numrows = numrows + ctr
  }
  
  firstrow = numeric(0)
  
  for(i in 1:(N-2)){
    temp = rep(i, (N-1-i)*(N-i)/2)
    firstrow = c(firstrow, temp)     # first ID
  }
  
  secondrow = numeric(0)
  for(i in 1:(N-2)){ # first row number
    for(j in (i+1):(N-1)){
      temp = rep(j, N-j)
      secondrow = c(secondrow, temp)  # Second ID
    }
  }
  
  thirdrow = numeric(0)
  for(i in 3:N){
    for(j in i:N){
      thirdrow = c(thirdrow, seq(j, N, 1))  # third ID
    }
  }
  
  
  triples = matrix(0, numrows, 3)
  triples[,1] = firstrow
  triples[,2] = secondrow
  triples[,3] = thirdrow
  
  
  ## Here's where the actual transitivity calculation begins.
  
  transitive = 0
  intransitive = 0
  strictTransitive = 0
  tList = matrix(0, 0, 4)
  iList = matrix(0, 0, 4)
  stList = matrix(0, 0, 4)
  for(i in 1:nrow(triples)){
    tA = triples[i,1]                  # first ID
    tB = triples[i,2]                  # second ID
    tC = triples[i,3]                  # third ID
    AB = conf[tA, tB] - conf[tB, tA]  
    AC = conf[tA, tC] - conf[tC, tA]
    BC = conf[tB, tC] - conf[tC, tB]
    ### See if the triangle is transitive...
    if((AC > 0 & BC > 0 & AB != 0) |    # BAC, ABC
       (AB < 0 & AC < 0 & BC != 0) |    # CBA, BCA
       (AB > 0 & BC < 0 & AC != 0)){    # CAB, ACB
      transitive = transitive + 1
      tList = rbind(tList, c(triples[i,], i))
    }
    ### See if the triangle is intransitive...
    if((AB > 0 & BC > 0 & AC < 0) | 
       (AB < 0 & AC > 0 & BC < 0)){
      intransitive = intransitive + 1
      iList = rbind(iList, c(triples[i,], i))
    }
    ABstrict = conf[tA, tB]
    ACstrict = conf[tA, tC]
    BCstrict = conf[tB, tC]
    BAstrict = conf[tB, tA]
    CAstrict = conf[tC, tA]
    CBstrict = conf[tC, tB]
    ### See if the triangle is strict transitive
    if(
      ((BAstrict > 0 | ABstrict > 0) & ACstrict > 0 & BCstrict > 0 & 
       CAstrict == 0 & CBstrict == 0) |
      (BAstrict > 0 & CAstrict > 0 & (BCstrict > 0 | CBstrict > 0) & 
       ABstrict == 0 & ACstrict == 0)|
      (ABstrict > 0 & CBstrict > 0 & (ACstrict > 0 | CAstrict > 0) & 
       BAstrict == 0 &  BCstrict == 0)) {
      strictTransitive = strictTransitive + 1
      stList = rbind(stList, c(triples[i,], i))
    }
  }
  
  if (strict) {
    
    T1 = strictTransitive / (transitive + intransitive)
    strictIntransitive = (transitive + intransitive) - strictTransitive
    alpha = (2 * sqrt(T1) - 1) / (1 - sqrt(T1))
    
    return(list(transitive = strictTransitive, 
                intransitive = strictIntransitive, 
                transitivity = T1, 
                alpha = alpha))
    
  } else {
    ### T1 is the order-1 transitivity.
    T1 = transitive / (transitive + intransitive)
    
    ### From the paper, we estimate alpha as follows.
    alpha = (2 * sqrt(T1) - 1) / (1 - sqrt(T1))
    
    return(list(transitive = transitive, intransitive = intransitive, transitivity = T1, alpha = alpha))
    
  }
}