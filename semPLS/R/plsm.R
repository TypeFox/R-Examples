plsm <- function(data, strucmod, measuremod, order=c("generic", "alphabetical"), interactive=FALSE)
{
  if(class(data)!="data.frame")
    stop("The argument 'data' must be of class 'data.frame'.")

  # interactive: using 'edit'
  if(interactive & missing(strucmod)){
    cat("Edit the structural model!\n")
    strucmod <- matrix("Enter!", nrow=1, ncol=2, dimnames=list(NULL, c('source', 'target')))
    strucmod <- edit(strucmod, title="Edit the structural model!")
    cat("Call for structural model:\n",
        paste("strucmod <- matrix(c(",
              paste(as.vector(t(strucmod)), collapse=", "),
              "), ncol=2, byrow=TRUE, dimnames=list(NULL, c('source', 'target')))", sep=""),
        "\n", sep="")
  }
  if(interactive & missing(measuremod)){
    cat("Edit the measurement model!\n")
    measuremod <- matrix("Enter!", nrow=1, ncol=2,
                         dimnames=list(NULL, c('source', 'target')))
    measuremod <- edit(measuremod, title="Edit the measurement model!")
    cat("Call for measurement model:\n",
        paste("measuremod <- matrix(c(",
              paste(as.vector(t(measuremod)), collapse=", "),
              "), ncol=2, byrow=TRUE, dimnames=list(NULL, c('source', 'target')))", sep=""),
        "\n", sep="")
  }

  if(!is.matrix(strucmod)){
    #cat("Choose a .csv file for the structural model!\n")
    strucmod <- as.matrix(read.csv(strucmod))
  }
  if(!is.matrix(measuremod)){
    #cat("Choose a .csv file for the measurement model!\n")
    measuremod <- as.matrix(read.csv(measuremod))
  }
  
  if(ncol(strucmod)!=2 || mode(strucmod)!="character" || class(strucmod)!="matrix")
    stop("The argument 'strucmod' must be a two column character matrix.")
  if(ncol(measuremod)!=2 || mode(measuremod)!="character" || class(measuremod)!="matrix")
    stop("The argument 'measuremod' must be a two column character matrix.")
    
  latent <- unique(as.vector(strucmod))
  if(any(latent %in% colnames(data)))
     stop("The latent variables are not allowed to coincide with names of observed variables.")
  manifest <- sort(setdiff(as.vector(measuremod), latent))

  if(!all(manifest %in% colnames(data)))
     stop("The manifest variables must be contained in the data.")

  order <- match.arg(order)
  # Adjacency matrix D for the structural model
  D <- innerW(strucmod, latent)

  # Ordering of LVs
  if (order=="generic"){
    tmp <- reorder(D)
    latent <- tmp$chain
    strucmod <- tmp$strucmod
  }
  if (order=="alphabetical"){
    latent <- sort(latent)
  }

  # Arranging the rows and columns according to the order of the LVs
  D <- D[latent, latent]

  # build blocks of manifest variables (including 'measurement mode')
  blocks <- block(latent, manifest, measuremod)

  # Ordering of MVs and measuremod
  MVs <- NULL
  mm <- NULL
  for(i in names(blocks)){
    MVs <- append(MVs, blocks[[i]])
    if(attr(blocks[[i]], "mode") == "A"){
      mm <- rbind(mm, (cbind(i, blocks[[i]])))
    }
    if(attr(blocks[[i]], "mode") == "B"){
      mm <- rbind(mm, (cbind(blocks[[i]], i)))
    }
  }
  dimnames(mm) <- dimnames(measuremod)
  measuremod <- mm

  result <- list()
  result$latent <- latent
  result$manifest <- MVs
  result$strucmod <- strucmod
  result$measuremod <- measuremod
  result$D <- D
  result$M <- initM1(model=result)
  result$blocks <- blocks
  result$order <- order
  class(result) <- "plsm"
  if(!connected(model=result)) stop("Structural model must be connected!")
  if(!acyclic(model=result)) stop("Structural model must be acyclic!")
  return(result)
}
