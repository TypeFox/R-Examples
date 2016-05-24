# Functions analysis

# ' @docType package
# ' ...
# ' @import ggplot2 ellipse grid scales


#' Specific Multiple Correspondence Analysis
#'
#' \code{soc.mca} performs a specific multiple correspondence analysis on a data.frame of factors, where cases are rows and columns are variables.
#' @param active       Defines the active modalities in a data.frame with rows of individuals and columns of factors, without NA's' 
#' @param sup          Defines the supplementary modalities in a data.frame with rows of individuals and columns of factors, without NA's 
#' @param identifier   A single vector containing a single value for each row/individual in x and sup. Typically a name or an id.number.
#' @param passive      A single character vector with the full or partial names of the passive modalities. All names that have a full or partial match will be set as passive.
#' 
#' @return \item{nd}{Number of active dimensions}
#'  \item{n.ind}{The number of active individuals}
#'  \item{n.mod}{The number of active modalities}
#'  \item{eigen}{Eigenvectors}
#'  \item{total.inertia}{The sum of inertia}
#'  \item{adj.inertia}{A matrix with all active dimensions, adjusted and unadjusted inertias. See \link{variance}}
#'  \item{freq.mod}{Frequencies for the active modalities. See \link{add.to.label}}
#'  \item{freq.sup}{Frequencies for the supplementary modalities. See \link{add.to.label}}
#'  \item{ctr.mod}{A matrix with the contribution values of the active modalities per dimension. See \link{contribution}}
#'  \item{ctr.ind}{A matrix with the contribution values of the individuals per dimension.}
#'  \item{cor.mod}{The correlation or quality of each modality per dimension.}
#'  \item{cor.ind}{The correlation or quality of each individual per dimension.}
#'  \item{mass.mod}{The mass of each modality}
#'  \item{coord.mod}{A matrix with the principal coordinates of each active modality per dimension.}
#'  \item{coord.ind}{A matrix with the principal coordinates of each individual per dimension.}
#'  \item{coord.sup}{A matrix with the principal coordinates of each supplementary modality per dimension.}
#'  \item{names.mod}{The names of the active modalities}
#'  \item{names.ind}{The names of the individuals}
#'  \item{names.sup}{The names of the supplementary modalities}
#'  \item{names.passive}{The names of the passive modalities}
#'  \item{modal}{A matrix with the number of modalities per variable and their location}
#'  \item{variable}{A vector with the name of the variable of the active modalities}
#' @name soc.mca
#' @export
#' @references Le Roux, B., og H. Rouanet. 2010. Multiple correspondence analysis. Thousand Oaks: Sage.
#' @author Anton Grau Larsen, University of Copenhagen
#' @author Stefan Bastholm Andrade, University of Copenhagen
#' @author Christoph Ellersgaard, University of Copenhagen
#' @seealso \link{soc.csa}, \link{contribution}
#' @examples
#' # Loads the "taste" dataset included in this package
#' data(taste)
#' # Create a data frame of factors containing all the active variables 
#' taste          <- taste[which(taste$Isup == 'Active'), ]
#'
#' attach(taste)
#' active         <- data.frame(TV, Film, Art, Eat)
#' sup            <- data.frame(Gender, Age, Income)
#' detach(taste)
#' 
#' # Runs the analysis
#' result         <- soc.mca(active, sup)
#' 
#' # Prints the results
#' result
#' 
#' # A specific multiple correspondence analysis
#' # options defines what words or phrases that are looked for in the labels of the active modalities.
#' options(passive = c("Film: CostumeDrama", "TV: Tv-Sport"))
#' soc.mca(active, sup)
#' options(passive = NULL)

soc.mca <- function(active, sup = NULL, identifier = NULL, passive = getOption("passive", default = "Missing")){
  
  active  <- data.frame(lapply(active, factor))               # Turn active variables into factor
  sup     <- data.frame(lapply(sup, factor))                  # Turn sup variables into factor   
  Q       <- ncol(active)                                     # Number of active variables 
  a.r     <- nrow(active)                                     # Number of active rows or the number of individuals
  sup.n   <- sum(unlist(lapply(as.data.frame(sup), nlevels))) # Number of supplementary modalities
  
  if ((nrow(sup) == 0) == TRUE){                                  
    sup             <- matrix(0, nrow = nrow(active), ncol = 2)
    sup[,1:2]       <- cbind(rep(0, nrow(active)), rep(0, nrow(active)))
    colnames(sup)   <- c("No supplementary points defined 1", "No supplementary points defined 2")
    ind.sup         <- sup
  }else{
    ind.sup         <- indicator(sup)  
  }
    
  # Creating the indicatormatrix for active and supplementary variables
  ind.act     <- indicator(active)
  
  # Finding the subset
  sub         <- grepl(paste(passive, collapse = "|"), colnames(ind.act))
  set         <- 1:ncol(ind.act)
  subset      <- set[!sub]
  
  # Finds the amount of variables without passive modalities
  Qm          <- Q
  for (i in seq(Q)){
    lev       <- levels(active[,i])
    pasQ      <- grepl(paste(passive, collapse = "|"), lev)
    if (any(pasQ == TRUE) == TRUE){
      Qm      <- Qm - 1
    }
  }
  
  result       <- subset.ca.indicator(ind.act, ind.sup, subset, Q = Q , Qm = Qm)
  
  if (identical(identifier, NULL) == TRUE){
    identifier  <- 1:nrow(active)
  }
  
  # Names
  result$names.mod        <- colnames(ind.act)[subset]
  result$names.ind        <- as.character(identifier)
  result$names.sup        <- colnames(ind.sup)
  result$names.passive    <- colnames(ind.act)[sub]
  
  # The active indicator matrix
  result$indicator.matrix <- ind.act[,subset]
  
  # List of descriptive values
  # The position and length of the active variables
  varnames      <- colnames(active)
  ml            <- vector()
  for (i in 1:ncol(active)){
    ml          <- c(ml, rep(varnames[i], nlevels(active[, i])))
  }
  ml            <- ml[!sub]
  mm            <- as.matrix(cbind(ml, 1:length(ml)))
  md            <- matrix(, nrow = ncol(active), ncol = 3)
  rownames(md)  <- varnames
  colnames(md)  <- c("Start", "End", "Modalities")
  md            <- as.data.frame(md)
  
  for (i in 1:ncol(active)){
    mr          <- as.numeric(mm[, 2][mm[, 1] == varnames[i]])
    md[i, 1]    <- min(mr)
    md[i, 2]    <- max(mr)
    md[i, 3]    <- length(mr)
  }
  md[, 1]       <- as.numeric(md[, 1])
  md[, 2]       <- as.numeric(md[, 2])
  result$modal  <- md
  
  variable    <- vector()
  for (i in 1:nrow(md)){
    variable  <- c(variable, rep(rownames(md)[i], md[i, 3]))
  }
  result$variable <- variable
  
  # 
  if (identical(sup, NULL) == FALSE){
  
  varnames       <- colnames(sup)
  ml             <- vector()
  for (i in 1:ncol(sup)){
    ml           <- c(ml, rep(varnames[i], nlevels(sup[, i])))
  }
  result$variable.sup <- ml
  }
  result$subset.var   <- Qm
  
  median.standard <- function(result){
    coord.ind     <- result$coord.ind
    coord.median  <- apply(coord.ind, 2, median)  
    dim.ind       <- seq(ncol(coord.ind))[coord.median > 0]
    result        <- invert(result, dim.ind)
    return(result)
  }
  
  result          <- median.standard(result)
  
  class(result)   <- "soc.mca"
  return(result)    
}


# ' Correspondence analysis on a indicator matrix
# ' 
# ' This function is part of the soc.mca function but allows for manipulation of the indicator matrix before analysis.
# ' Most users will not need this function.
# ' 
# ' @param indicator.act   An indicator matrix of all the active modalities (including those that are to be set as passive)
# ' @param indicator.sup   An indicator matrix of the supplementary modalities
# ' @param subset  A vector containing column indices of passive modalities
# ' @param Q       The number of variables
# ' @param Qm      The number of variables without passive modalities
# ' #@export
# ' @return a list of various results. See \link{soc.mca} documentation

subset.ca.indicator <- function(indicator.act, indicator.sup, subset, Q, Qm){
  Z.act     <- indicator.act
  Z.sup     <- indicator.sup
  
  I         <- dim(Z.act)[1]  # Number of individuals
  J         <- dim(Z.act)[2]  # Number of modalities >> Subset
  Q         <- Q              # Number of variables
  
  # Inertias
  P         <- Z.act / sum(Z.act)      
  cm        <- colSums(P)              # Column (modality) mass
  rm        <- rowSums(P)              # Row (individual) mass
  
  diag.cm   <- diag(1/ sqrt(cm))       # This commmand scales badly because it creates a individual X individual matrix - If this number could be obtained differently - for instance - through the Burt matrix - there is a substantial speed gain.
    
  eP        <- rm %*% t(cm)            # Expected distances
  S         <- (P - eP) / sqrt(eP)     # Euclidian distances
  
  # Subsetting
  K         <- length(subset)
  S         <- S[, subset]
  cm        <- cm[subset]
  diag.cm   <- diag.cm[subset, subset]
  
  # Decomposition and eigenvectors
  dec       <- svd(S)                 # Singular decomposition
  eigen     <- dec$d^2                # Eigenvector
    
  # Principal coordinates
  pc.mod    <- diag.cm %*% dec$v %*% diag(dec$d)   # Principal coordinates for modalities
  
  # Fast principal coordinates for individuals
  if (identical(var(rm), 0)){
  sqrm      <- 1/ sqrt(rm)
  pc.ind    <- (sqrm[1] * dec$u) %*% diag(dec$d)    
  }else{
  # Original principal coordinates for individuals
  diag.rm   <- diag(1/ sqrt(rm))  
  pc.ind    <- diag.rm %*% dec$u %*% diag(dec$d)   # Principal coordinates for individuals # This is a slow process, but it scales ok # Anders Holm adjustment  
  }
  
    
  # Fast inertias
  if (identical(var(rm), 0)){
  inr.ind   <- rm[1] * pc.ind^2
  inr.mod   <- diag(cm) %*% pc.mod^2
  }else{
  # Original inertias
  inr.ind   <- diag(rm) %*% pc.ind^2     # Inertia for row (Individuals) (mass x principal coordinates) # This is a slow process and it scales badly - diag(rm) is a individual X individual matrix. It is also sparse - so it might be possible to do it quicker.
  inr.mod   <- diag(cm) %*% pc.mod^2     # Inertia for columns (Modalities)
  }
  
  # Contributions
  ctr.ind   <- t(t(inr.ind) / dec$d^2)   # Contribution for the individuals (inertia / eigenvalue)
  ctr.mod   <- t(t(inr.mod) / dec$d^2)   # Contribution for the modalities
  
  # Squared cosines or correlations
  cor.ind   <- inr.ind/colSums(inr.ind)  # Squared cosines for individuals
  cor.mod   <- inr.mod/rowSums(inr.mod)  # Squared cosines for modalities
    
  # Chi-distances
  
  # Supplementary principal coordinates
  Z.star    <- Z.sup
  I.star    <- dim(Z.sup)[1]
  cs.star   <- apply(Z.sup, 2, sum)
  base      <- Z.star / matrix(rep(cs.star, I.star), nrow = I.star, byrow = TRUE)
  f.s1      <- dec$u * sqrt(eigen) / sqrt(rm)   
  a.s1      <- f.s1 / sqrt(eigen)               
  pc.sup    <- t(base) %*% a.s1
  
  # Supplementary squared cosines or correlations
  #cor.sup <- inr.mod/apply(inr.mod, 1, sum)   # Squared cosines for modalities
  
  # First reduction of dimensionality
  R1.dim        <- K-Qm
  R1.eigen.val  <- eigen[1:R1.dim]
  
  # Second reduction of dimensionality
  Mean.r.eigen  <- mean(R1.eigen.val, na.rm = TRUE)
  R2.eigen.val  <- eigen[eigen >= Mean.r.eigen]
  R2.dim        <- which(R2.eigen.val >= Mean.r.eigen)
  
  # Explained variance
  unadj.var     <- 100*(eigen[1:R1.dim])/sum(eigen[1:R1.dim]) # Unadjusted rates of variance
  sum.adj.var.mod <- (eigen[R2.dim]-Mean.r.eigen)^2
  sum.adj.var   <- sum(sum.adj.var.mod)
  adj.var       <- round(((sum.adj.var.mod/sum.adj.var)*100), digits = 1) # The adjusted rates of variance
  cumpercent    <- cumsum(adj.var)
  adj.inertia   <- cbind(R2.dim, round(eigen[R2.dim], 3), round(unadj.var[R2.dim], 1), adj.var ,cumpercent)
  colnames(adj.inertia) <- c("Dim", "Eigen", "Var" ,"Adj.Var", "Cum%")
  
  freq.mod      <- colSums(Z.act)[subset]
  freq.sup      <- colSums(Z.sup)
  
  # Output
  ca.output <- list(nd        = max(R2.dim),
                    n.ind     = nrow(Z.act),
                    n.mod     = length(subset),
                    eigen     = eigen[R2.dim],
                    total.inertia = sum(eigen[R2.dim]),
                    adj.inertia = adj.inertia,
                    freq.mod  = freq.mod,
                    freq.sup  = freq.sup,
                    ctr.mod   = ctr.mod[, R2.dim],
                    ctr.ind   = ctr.ind[, R2.dim],
                    cor.mod   = cor.mod[, R2.dim],
                    cor.ind   = cor.ind[, R2.dim],
                    mass.mod  = cm, 
                    coord.mod = pc.mod[, R2.dim],
                    coord.ind = pc.ind[, R2.dim],
                    coord.sup = pc.sup[, R2.dim]
                    )
  # Cleanup
  names(ca.output$mass.mod)         <- NULL
  dimnames(ca.output$coord.sup)     <- NULL
  names(ca.output$freq.mod)         <- NULL
  names(ca.output$freq.sup)         <- NULL
  
  return(ca.output)
}


################################### Indicator matrix

#' Indicator matrix
#' 
#' Creates an indicator matrix from a data.frame with questions as columns and individuals as rows.
#' 
#' @param x   a data.frame of factors
#' @param id  a vector defining the labels for the individuals. If id = NULL row number is used.
#' @param ps  the seperator used in the creation of the names of the columns (modalities).
#'
#' @return Returns a indicator matrix
#' @seealso \link{soc.mca}
#' @examples 
#' a  <- rep(c("A","B"), 5)
#' b  <- rep(c("C", "D"), 5)
#' indicator(data.frame(a,b))
#' @export


indicator  <- function(x, id = NULL, ps = ": "){
obj         <- x
I           <- nrow(obj)                                      # Number of individuals
levels.n    <- unlist(lapply(obj, nlevels))
n           <- cumsum(levels.n)                               # Number of modalities for each question
m           <- max(n)                                         # Total number of modalities
Q           <- ncol(obj)                                      # Number of questions
Z           <- matrix(0, nrow = I, ncol = m)                  # Predefinition of the indicatormatrix
newdat      <- lapply(obj, as.numeric)
offset      <- (c(0, n[-length(n)]))
for (i in seq(Q)) Z[seq(I) + (I * (offset[i] + newdat[[i]] - 1))] <- 1 # Indicator matrix
fn          <- rep(names(obj), unlist(lapply(obj, nlevels)))  
ln          <- unlist(lapply(obj, levels))
col.names   <- paste(fn, ln, sep = ps)
colnames(Z) <- col.names

if (identical(id, NULL) == TRUE){
rownames(Z) <- as.character(seq(I))
}else{
rownames(Z) <- id
}    
return(Z)

}

#' Class Specific Multiple Correspondence Analysis
#'
#' \code{soc.csa} performs a class specific multiple correspondence analysis on a data.frame of factors, where cases are rows and columns are variables. Most descriptive and analytical functions that work for \link{soc.mca}, also work for \code{soc.csa}
#' @param object  is a soc.ca class object created with \link{soc.mca}
#' @param class.indicator the row indices of the class specific individuals
#' @param sup          Defines the supplementary modalities in a data.frame with rows of individuals and columns of factors, without NA's 
#' @return \item{nd}{Number of active dimensions}
#'  \item{n.ind}{The number of active individuals}
#'  \item{n.mod}{The number of active modalities}
#'  \item{eigen}{Eigenvectors}
#'  \item{total.inertia}{The sum of inertia}
#'  \item{adj.inertia}{A matrix with all active dimensions, adjusted and unadjusted inertias. See \link{variance}}
#'  \item{freq.mod}{Frequencies for the active modalities. See \link{add.to.label}}
#'  \item{freq.sup}{Frequencies for the supplementary modalities. See \link{add.to.label}}
#'  \item{ctr.mod}{A matrix with the contribution values of the active modalities per dimension. See \link{contribution}}
#'  \item{ctr.ind}{A matrix with the contribution values of the individuals per dimension.}
#'  \item{cor.mod}{The correlation or quality of each modality per dimension.}
#'  \item{cor.ind}{The correlation or quality of each individual per dimension.}
#'  \item{mass.mod}{The mass of each modality}
#'  \item{coord.mod}{A matrix with the principal coordinates of each active modality per dimension.}
#'  \item{coord.ind}{A matrix with the principal coordinates of each individual per dimension.}
#'  \item{coord.sup}{A matrix with the principal coordinates of each supplementary modality per dimension. Notice that the position of the supplementary modalities in class specific analysis is the mean point of the individuals, which is not directly comparable with the cloud of the active modalities.}
#'  \item{indicator.matrix}{A indicator matrix. See \link{indicator}}
#'  \item{names.mod}{The names of the active modalities}
#'  \item{names.ind}{The names of the individuals}
#'  \item{names.sup}{The names of the supplementary modalities}
#'  \item{names.passive}{The names of the passive modalities}
#'  \item{modal}{A matrix with the number of modalities per variable and their location}
#'  \item{variable}{A vector with the name of the variable for each of the active modalities}
#'  \item{variable.sup}{A vector with the name of the variable for each of the supplementary modalities}
#'  \item{original.class.indicator}{The class indicator}
#'  \item{original.result}{The original soc.ca object used for the CSA}
#' @name soc.csa
#' @export
#' @author Anton Grau Larsen, University of Copenhagen
#' @author Stefan Bastholm Andrade, University of Copenhagen
#' @author Christoph Ellersgaard, University of Copenhagen
#' @seealso \link{add.to.label}, \link{contribution}
#' @references Le Roux, B., og H. Rouanet. 2010. Multiple correspondence analysis. Thousand Oaks: Sage.
#' @examples 
#' example(soc.ca)
#' class.age    <- which(taste$Age == '55-64')
#' res.csa      <- soc.csa(result, class.age)
#' res.csa

soc.csa <- function(object, class.indicator, sup = NULL){
  
  
  Z.act   <- object$indicator.matrix                # Original indicator matrix
  Q       <- nlevels(as.factor(object$variable))    # Number of questions
  I       <- nrow(Z.act)                            # Original number of individuals
  
  Z.hat   <- Z.act[class.indicator, ]               # Indicator matrix for the CSA
  i       <- length(class.indicator)                # Number of individuals in CSA
  
  cm      <- apply(Z.hat, 2, sum)
  CM      <- apply(Z.act, 2, sum)
  P       <- Z.hat / sum(Z.hat)      
  cmpc    <- colSums(P)              # Column (modality) mass
  rmpc    <- rowSums(P)  
  
  
  H.hat   <- matrix(, nrow = i, ncol = length(cm))
  for (k in seq(cm)){
    H.hat[,k] <- (1/sqrt(Q)) * (sqrt(I/i)) * (Z.hat[, k]-(cm[k]/i)) * (1/sqrt(CM[k]))
  }
  
  colnames(H.hat)   <- colnames(Z.hat)
  modal.names       <- colnames(Z.hat)
  
  H.svd             <- svd(H.hat)
  dec               <- H.svd
  
  ### Modalitetskoordinater
  csa.m.coord            <- matrix(nrow = nrow(H.svd$v), ncol = ncol(H.svd$v))
  for (ff in 1:length(H.svd$d)){
  csa.m.coord[, ff]      <- H.svd$d[ff] * H.svd$v[, ff] * (1/sqrt(CM )) * (sqrt(Q*I))
  } 
  rownames(csa.m.coord)  <- modal.names
  
  ### Individkoordinater
  csa.i.coord            <- matrix(nrow = nrow(H.svd$u), ncol = ncol(H.svd$u))
  for (ff in 1:length(H.svd$d)){
  csa.i.coord[, ff]      <- sqrt(i) * H.svd$u[, ff] * H.svd$d[ff]
  }
  
  ### Modalitetsbidrag
  csa.m.ctr              <- matrix(nrow = nrow(H.svd$v), ncol = ncol(H.svd$v))
  for (ff in 1:length(H.svd$d)){
  csa.m.ctr[, ff]        <- (((CM/I)/Q) * (csa.m.coord[, ff])^2)/(H.svd$d[ff]^2)
  }
  
  # Eigenvectors
  eigen            <- H.svd$d^2 
  
  ###############################################
  #### Dimension reduction and explained variance
  K                <- ncol(Z.act)
  Qm               <- object$subset.var
  
  # First reduction of dimensionality
  R1.dim           <- K-Qm
  R1.eigen.val     <- eigen[1:R1.dim]
  
  # Second reduction of dimensionality
  Mean.r.eigen     <- mean(R1.eigen.val, na.rm = TRUE)
  R2.eigen.val     <- eigen[eigen >= Mean.r.eigen]
  R2.dim           <- which(R2.eigen.val >= Mean.r.eigen)
  
  # Explained variance
  unadj.var        <- 100*(eigen[1:R1.dim])/sum(eigen[1:R1.dim])               # Unadjusted rates of variance
  sum.adj.var.mod  <- (eigen[R2.dim]-Mean.r.eigen)^2
  sum.adj.var      <- sum(sum.adj.var.mod)
  adj.var          <- round(((sum.adj.var.mod/sum.adj.var)*100), digits = 1)   # The adjusted rates of variance
  cumpercent       <- cumsum(adj.var)
  adj.inertia      <- cbind(R2.dim, round(eigen[R2.dim], 3), round(unadj.var[R2.dim], 1), adj.var ,cumpercent)
  colnames(adj.inertia) <- c("Dim", "Eigen", "Var" ,"Adj.Var", "Cum%")
  
  ###################
  
  freq.mod         <- cm
  
  ###############
  # Principal coordinates
  
  # Principal coordinates for modalities
  diag.cm         <- diag(cmpc)
  pc.mod          <- diag.cm %*% dec$v %*% diag(dec$d)   # Principal coordinates for modalities
  
  # Fast principal coordinates for individuals
  if (identical(var(rmpc), 0)){
    sqrm          <- 1/ sqrt(rmpc)
    pc.ind        <- (sqrm[1] * dec$u) %*% diag(dec$d)    
  }else{
    # Original principal coordinates for individuals
    diag.rm       <- diag(1/ sqrt(rmpc))  
    pc.ind        <- diag.rm %*% dec$u %*% diag(dec$d)   # Principal coordinates for individuals # This is a slow process, but it scales ok # Anders Holm adjustment  
  }
  
  ############
  # Preparing for object object
  names.mod      <- modal.names
  names.ind      <- object$names.ind[class.indicator]
  names.passive  <- object$names.passive
  
  
  ############
  # Supplementary variables
  # They are projected into the cloud of individuals - NOT the cloud of modalities
  
  sup.coord                 <- NULL
  freq.sup                  <- 0
  if(identical(sup,NULL) == FALSE){
    
    sup.ind                 <- indicator(sup)[class.indicator, ]
    sup.coord               <- matrix(nrow = ncol(sup.ind), ncol = ncol(pc.ind))
    rownames(sup.coord)     <- colnames(sup.ind)
    freq.sup                <- colSums(sup.ind)
    
    for (j in 1:ncol(pc.ind)){
      d1                    <- pc.ind[,j]
      Sup.mat               <- matrix(nrow = nrow(sup.ind), ncol = ncol(sup.ind))
      colnames(Sup.mat)     <- colnames(sup.ind)
      Sup.mat[sup.ind == 1] <- 0
      Sup.mat               <- Sup.mat + d1
      sup.coord[, j]        <- apply(Sup.mat, 2, mean, na.rm = TRUE)
    }
  }
  names.sup                 <-  rownames(sup.coord)
  ############# 
  # Result object
  
  csca.result <- list(
                      nd          = object$nd,
                      n.ind       = i, 
                      n.mod       = ncol(Z.hat),
                      eigen       = eigen,
                      total.inertia = sum(eigen[R2.dim]),
                      adj.inertia = adj.inertia,
                      freq.mod    = freq.mod,
                      freq.sup    = freq.sup,
                      ctr.mod     = csa.m.ctr[, R2.dim],
                      #ctr.ind    = ctr.ind[, R2.dim], # Not implemented
                      #cor.mod    = cor.mod[, R2.dim], cor.ind = cor.ind[, R2.dim], # Not implemented
                      mass.mod    = cmpc, # Massen er etchy - se nedenfor
                      mass.ind    = rmpc, # Massen er etchy - se nedenfor
                      coord.mod   = csa.m.coord,
                      #coord.mod  = pc.mod[, R2.dim],
                      coord.ind   = pc.ind[, R2.dim],
                      coord.sup   = sup.coord,
                      coord.mod.standard = csa.m.coord,
                      coord.ind.standard = csa.i.coord,
                      names.mod   = names.mod,
                      names.ind   = names.ind,
                      names.sup   = names.sup,
                      names.passive = names.passive,
                      indicator.matrix = Z.hat,
                      modal       = object$modal,
                      variable    = object$variable,
                      variable.sup = "Not Implemented",
                      subset.var  = object$subset.var,
                      original.result = object,
                      original.class.indicator = class.indicator
  )
  
  median.standard <- function(object){
    
    coord.ind     <- object$coord.ind
    coord.median  <- apply(coord.ind, 2, median)  
    dim.ind       <- seq(ncol(coord.ind))[coord.median > 0]
    object        <- invert(object, dim.ind)
    return(object)
  }
  
  csca.result     <- median.standard(csca.result)
  
  #####################################
  # Class and return
  
  class(csca.result)    <- c("soc.mca", "soc.csa")
  return(csca.result)
  
} 

####################################################
#### Create Quadrant

#' Create categories according to the quadrant position of each individual
#' 
#' Creates a vector from two dimensions from a soc.ca object. Labels are the 
#' cardinal directions with the first designated dimension running East - West.
#' The center category is a circle defined by \code{cut.radius}.
#' 
#' 
#' @param object a soc.ca class object
#' @param dim  the dimensions
#' @param cut.min  Minimum cut value
#' @param cut.max  Maximum cut value
#' @param cut.radius  Radius of the center category
#' @return Returns a character vector with category memberships
#' @seealso \link{soc.mca}
#' @examples 
#' example(soc.ca)
#' create.quadrant(result, dim = c(2, 1))
#' table(create.quadrant(result, dim = c(1, 3), cut.radius = 0.5))
#' @export

create.quadrant <- function(object, dim = c(1,2), cut.min = -0.125, cut.max = 0.125, cut.radius = 0.25){
  
  coord                      <- object$coord.ind
  coord.cut                  <- coord
  coord.cut[coord < cut.min] <- "Min"
  coord.cut[coord > cut.max] <- "Max"
  coord.cut[coord <= cut.max & coord >=cut.min]  <- "Medium"
  
  dim1     <- coord.cut[,dim[1]]
  dim2     <- coord.cut[,dim[2]]
  distance <- sqrt(((coord[, dim[1]]^2) + (coord[, dim[2]]^2)))
  position <- dim1
  
  position[dim1 == "Max" & dim2 == "Max"]        <- "North-East"
  position[dim1 == "Max" & dim2 == "Medium"]     <- "East"
  position[dim1 == "Max" & dim2 == "Min"]        <- "South-East"
  
  position[dim1 == "Medium" & dim2 == "Medium"]  <- "Center"
  position[dim1 == "Medium" & dim2 == "Min"]     <- "South"
  position[dim1 == "Medium" & dim2 == "Max"]     <- "North"
  
  position[dim1 == "Min" & dim2 == "Max"]        <- "North-West"
  position[dim1 == "Min" & dim2 == "Medium"]     <- "West"
  position[dim1 == "Min" & dim2 == "Min"]        <- "South-West"
  
  position[distance < cut.radius]                <- "Center"
  
  return(position)
}


#' Multiple Class Specific Correspondence Analysis on all values in a factor
#' 
#' \code{csa.all} performs a class specific correspondence analysis for each
#' level in a factor variable. Returns a list with soc.csa objects and a list of
#' measures defined by \link{csa.measures}
#' @param object  is a soc.ca class object created with \link{soc.mca}
#' @param variable a factor with the same length and order as the active
#'   variables that created the soc.ca object
#' @param dim is the dimension analyzed
#' @param ... further arguments are directed to \link{csa.measures}
#' @return \item{results}{a list of \link{soc.csa} result objects}
#' @return \item{cor}{a list of correlation matrixes}
#' @return \item{cosines}{a list of matrixes with cosine values}
#' @return \item{angles}{a list of matrixes with cosine angles between
#'   dimensions}
#' @export
#' @examples
#' example(soc.ca)
#' csa.all(result, taste$Age)
#' csa.all(result, taste$Age)$measures
#' @seealso \link{soc.csa}, \link{cor}, \link{csa.measures}

csa.all <- function(object, variable, dim = 1:5, ...){
  lev.variable <- levels(variable)
  result.list     <- list()
  for (i in 1:length(lev.variable)){
    dummy.class         <- which(variable == lev.variable[i])
    result.list[[i]]    <- soc.csa(object, class.indicator = dummy.class)
  }
  
  names(result.list) <- lev.variable
  
  measure.list <- lapply(result.list, csa.measures, format = FALSE, dim = dim, ...)
  
  list(results = result.list, measures = measure.list)
}

#' Add supplementary individuals to a result object
#' 
#' @param object is a soc.ca class object created with \link{soc.mca}
#' @param sup.indicator is a indicator matrix for the supplementary individuals with the same columns as the active variables in object.
#' @param replace if TRUE the coordinates of the active individuals are discarded. If FALSE the coordinates of the supplementary and active individuals are combined. The factor \code{object$supplementary.individuals} marks the supplementary individuals.
#' @return  a soc.ca class object created with \link{soc.mca}
#' @export
#' @examples
#' example(soc.mca)
#' res.pas   <- soc.mca(active, passive = "Costume")
#' res.sup   <- supplementary.individuals(res.pas, sup.indicator = indicator(active))
#' a         <- res.sup$coord.ind[res.sup$supplementary.individuals == "Supplementary",]
#' b         <- res.pas$coord.ind
#' all.equal(as.vector(a), as.vector(b))
#' map.ind(res.sup)

supplementary.individuals <- function(object, sup.indicator, replace = FALSE){
  
  if (length(object$names.passive) > 0) sup.indicator   <- sup.indicator[, -which(colnames(sup.indicator) %in% object$names.passive)]
  
#   sup.ind.dim     <- function(object, sup.indicator, dim){
#     Q             <- length(table(object$variable))
#     yk            <- t(object$coord.mod[, dim] * t(sup.indicator))
#     #Q.ind         <- yk / rowSums(sup.indicator)
#     Q.ind         <- yk / Q
#     out           <- 1/sqrt(object$eigen[dim]) * rowSums(Q.ind)
#     out
#   }
  
  sup.ind.dim     <- function(object, sup.indicator, dim){
    Q             <- length(table(object$variable))
    yk            <- t(object$coord.mod[, dim] * t(sup.indicator))
    pk            <- (colSums(object$indicator.matrix)/object$n.ind) / Q
    sas           <- pk * object$coord.mod[, dim]
    Q.ind         <- yk / Q
    out           <- 1/sqrt(object$eigen[dim]) * (rowSums(Q.ind) - sum(sas))
    out
  }

  
  ndim            <- 1:ncol(object$coord.ind)
  sup.ind.coord   <- sapply(ndim, sup.ind.dim, object = object, sup.indicator = sup.indicator)
  
  if(identical(replace, FALSE)){
  object$coord.ind                 <- rbind(object$coord.ind, sup.ind.coord)
  rownames(object$coord.ind)       <- NULL
  object$supplementary.individuals <- c(rep("Active", object$n.ind), rep("Supplementary", nrow(sup.indicator)))
  object$names.ind                 <- c(object$names.ind, rownames(sup.indicator))
  }
  
  if(identical(replace, TRUE)){
    object$coord.ind                 <- sup.ind.coord
    rownames(object$coord.ind)       <- NULL
    object$names.ind                 <- rownames(sup.indicator)
  }
 object
}
 
# # Test of supplementary individuals
# example(soc.mca)
# 
# # Ingen passive modaliteter
# res.sup   <- supplementary.individuals(result, sup.indicator = indicator(active))
# a         <- res.sup$coord.ind[res.sup$supplementary.individuals == "Supplementary",]
# b         <- result$coord.ind
# all.equal(as.vector(a), as.vector(b))
# map.ind(res.sup)
# 
# # Med passive modaliteter
# res.pas   <- soc.mca(active, passive = "Costume")
# res.sup   <- supplementary.individuals(res.pas, sup.indicator = indicator(active))
# a         <- res.sup$coord.ind[res.sup$supplementary.individuals == "Supplementary",]
# b         <- res.pas$coord.ind
# all.equal(as.vector(a), as.vector(b))
# map.ind(res.sup)
# 
# # Med passive modaliteter 
# res.pas   <- soc.mca(active, passive = "Costume")
# res.sup   <- supplementary.individuals(res.pas, sup.indicator = indicator(active)[1:10,])
# a         <- res.sup$coord.ind[res.sup$supplementary.individuals == "Supplementary",]
# b         <- res.pas$coord.ind[1:10,]
# all.equal(as.vector(a), as.vector(b))
# map.ind(res.sup, point.fill = res.sup$supplementary.individuals, label = T)
# 
# # Med passive modaliteter
# sup.indicator <- indicator(active)[1:10,]
# rownames(sup.indicator) <- paste("Sup", 1:10)
# res.pas   <- soc.mca(active, passive = "Costume")
# res.sup   <- supplementary.individuals(res.pas, sup.indicator = sup.indicator)
# a         <- res.sup$coord.ind[res.sup$supplementary.individuals == "Supplementary",]
# b         <- res.pas$coord.ind[1:10,]
# all.equal(as.vector(a), as.vector(b))
# map.ind(res.sup, point.fill = res.sup$supplementary.individuals, label = T)

