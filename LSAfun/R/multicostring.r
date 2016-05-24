##### Multicostring

#' @export
#' @importFrom lsa cosine
multicostring <- function(x,y,tvectors=tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    if(breakdown==TRUE){
      
      satz1 <- breakdown(x)
      y     <- breakdown(y)
      
    }
    
    if(breakdown==FALSE){
      
      satz1 <- x
      y     <- y
      
    }
    
    ### turn x into a single vector (added vectors of elements)
    
    if(length(satz1) == 1){
      satz1split <- strsplit(satz1,split=" ")[[1]]
    }
    
    if(length(satz1)  > 1){satz1split <- satz1}
    
    used1     <- satz1split[satz1split %in% rownames(tvectors)]
    if(length(used1)==0){(warning("no element of x found in rownames(tvectors)"))
                         return(NA)}
    
    rest1    <- satz1split[!(satz1split %in% rownames(tvectors))]
    
    if(length(used1) >1){satz1vec <- colSums(tvectors[used1,])}
    if(length(used1)==1){satz1vec <- tvectors[used1,]}
    
    ### split up y into single words
    
    
    if(length(y) == 1){
      ysplit <- strsplit(y,split=" ")[[1]]
    }
    
    if(length(y)  > 1){ysplit <- y}
    
    
    used2     <- ysplit[ysplit %in% rownames(tvectors)]
    if(length(used2)==0){(warning("no element of y found in rownames(tvectors)"))
                         return(NA)}
    
    if(length(used2) < length(ysplit)){
      cat("Note: not all elements in y were found in rownames(tvectors)\n\n")}
    
    vecs2 <- matrix(tvectors[used2,],nrow=length(used2))
    
    cosmat <- matrix(nrow=1,ncol=nrow(vecs2))
    rownames(cosmat) <- "expression in x"
    colnames(cosmat) <- used2
    
    
    
    for(j in 1:nrow(vecs2)){
      
      cosmat[1,j] <- cosine(satz1vec,vecs2[j,])
      
    }
    
    
    cosmat
    
  }else{stop("tvectors must be a matrix!")}
  
}
