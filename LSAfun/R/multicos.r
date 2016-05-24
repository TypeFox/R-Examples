##### Vector x Vector Comparison

#' @export
#' @importFrom lsa cosine
multicos <- function(x,y=x,tvectors=tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    ## handle x
    
    if(class(x)=="character"){
      
      if(breakdown==TRUE){x <- breakdown(x)}
      
      if(length(x) == 1){
        xsplit <- strsplit(x,split=" ")[[1]]
      }
      
      if(length(x)  > 1){xsplit <- x}
      
      used1     <- xsplit[xsplit %in% rownames(tvectors)]
      
      if(length(used1)==0){warning("no element of x found in rownames(tvectors)")
                           return(NA)}
      if(length(used1) < length(xsplit)){
        cat("Note: not all elements in x were found in rownames(tvectors)\n\n")}
      
      vecs1 <- matrix(tvectors[used1,],nrow=length(used1))
      
    }
    
    if(class(x)=="numeric"){
      used1 <- "Expression in x"
      vecs1 <- t(as.matrix(x))
    }
    
    ## handle y
    
    if(breakdown==TRUE){y <- breakdown(y)}
    
    if(length(y) == 1){
      ysplit <- strsplit(y,split=" ")[[1]]
    }
    
    if(length(y)  > 1){ysplit <- y}
    
    used2 <- ysplit[ysplit %in% rownames(tvectors)]
    
    if(length(used2)==0){warning("no element of y found in rownames(tvectors)")
                         return(NA)}
    if(length(used2) < length(ysplit)){
      cat("Note: not all elements in y were found in rownames(tvectors)\n\n")}
    
    vecs2 <- matrix(tvectors[used2,],nrow=length(used2))
    
    ## compute cosine matrix
    
    cosmat <- matrix(nrow=nrow(vecs1),ncol=nrow(vecs2))
    rownames(cosmat) <- used1
    colnames(cosmat) <- used2
    
    for(i in 1:nrow(vecs1)){
      
      line <- vecs1[i,]
      
      for(j in 1:nrow(vecs2)){
        
        cosmat[i,j] <- cosine(line,vecs2[j,])
        
      }
    }
    
    return(cosmat)
    
  }else{stop("tvectors must be a matrix!")}
  
}
