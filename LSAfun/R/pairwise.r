###### Fast pairwise cosine computations

#' @export
#' @importFrom lsa cosine

pairwise <- function(x,y,tvectors=tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    if(length(x) != length(y)){
      stop("x and y must be of same length")
    }
    
    if(length(x[x %in% rownames(tvectors)]) < length(x)){
      cat("Note: not all elements in x were found in rownames(tvectors)\n\n")
    }
    
    if(length(y[y %in% rownames(tvectors)]) < length(y)){
      cat("Note: not all elements in y were found in rownames(tvectors)\n\n")
    }
    
    x[which(!(x %in% rownames(tvectors)))] <- NA
    y[which(!(y %in% rownames(tvectors)))] <- NA
    
    
    if(breakdown==TRUE){x <- breakdown(x)
                        y <- breakdown(y)}
    
    
    vecs1                    <- matrix(nrow=length(x),ncol=ncol(tvectors))
    vecs1[which(!is.na(x)),] <- tvectors[x[!is.na(x)],]
    
    vecs2                    <- matrix(nrow=length(y),ncol=ncol(tvectors))
    vecs2[which(!is.na(y)),] <- tvectors[y[!is.na(y)],]
    
    
    cos <- vector(length=length(x))
    
    for(i in 1:length(x)){
      cos[i] <- as.numeric(cosine(vecs1[i,],vecs2[i,]))
    }
    
    return(cos)
    
  }else{stop("tvectors must be a matrix!")}
  
}