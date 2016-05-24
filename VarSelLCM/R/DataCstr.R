########################################################################################################################
## Declaration et construction des classes S4 relatives aux donnees
########################################################################################################################


########################################################################################################################
## La classe S4 VSLCMdataContinuous est relatives a des donnees continues. Elle possede 5 slots:
## n: nombre d'observations
## d: nombre de variables
## data: matrix ou les colonnes sont numeric et correspondent aux donees
## notNA: matrix of logical valant 1 si la realisation est observee et 0 sinon
## priors: valeur des priors pour chaque variable (en ligne)
########################################################################################################################
setClass(
  Class = "VSLCMdataContinuous", 
  representation = representation(
    n="numeric",
    d="numeric",
    data="matrix",
    notNA="matrix",
    priors="matrix"
  ), 
  prototype = prototype(
    n=numeric(),
    d=numeric(),
    data=matrix(),
    notNA=matrix(),
    priors=matrix()
  )
)



########################################################################################################################
## La fonction VSLCMdata permet de construire un objet de class S4 VSLCMdataContinuous
## de la nature des variables
########################################################################################################################
VSLCMdata <- function(x, redquali=TRUE){
  # Ajout d'un nom de variable si celui-ci est manquant
  if (is.null(colnames(x))) colnames(x) <- paste("X",1:ncol(x), sep="")
  
  n <- nrow(x)
  d <- ncol(x)
  # recherche des indices de variables numeric et factor
  type <- numeric()
  for (j in 1:d) type[j] <- class(x[,j])
  
  if ((all(type %in% c("numeric", "integer"))==FALSE))
    stop("At least one variable is neither numeric nor integer!")
  
  mat <- apply(x, 2, as.numeric)
  # construction des priors
  priors <- matrix(1, d, 4)
  priors[,4] <- 1/100
  priors[,3] <- colMeans(x, na.rm = T)
  colnames(priors) <- c("alpha", "beta", "lambda", "delta")
  ## Pour travailler avec Armadillo on rempli artificellement les NA par 0
  notNA <- (is.na(x)==FALSE)*1
  mat[is.na(mat)] <- 0
  colnames(mat) <-  colnames(x)
  colnames(notNA) <- colnames(x)
  output <-  new("VSLCMdataContinuous", n=n, d=d, data=mat, notNA=notNA, priors=priors)    
  
  return(output)
}