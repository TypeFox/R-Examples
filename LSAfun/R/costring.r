##### Sentence Comparisons #####################

#' @export
#' @importFrom lsa cosine 
 
costring <- function(x,y,tvectors=tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    if(breakdown==TRUE){
      
      satz1 <- breakdown(x)
      satz2 <- breakdown(y)
      
    }
    
    if(breakdown==FALSE){
      
      satz1 <- x
      satz2 <- y
      
    }
    
    
    if(length(satz1) == 1){
      satz1split <- strsplit(satz1,split=" ")[[1]]
    }
    if(length(satz2) == 1){
      satz2split <- strsplit(satz2,split=" ")[[1]]
    }
    
    if(length(satz1)  > 1){satz1split <- satz1}
    if(length(satz2)  > 1){satz2split <- satz2}
    
    
    used1     <- satz1split[satz1split %in% rownames(tvectors)]
    if(length(used1)==0){(warning("no element of x found in rownames(tvectors)"))
                         return(NA)}
    
    used2     <- satz2split[satz2split %in% rownames(tvectors)]
    if(length(used2)==0){(warning("no element of y found in rownames(tvectors)"))
                         return(NA)}
    
    rest1    <- satz1split[!(satz1split %in% rownames(tvectors))]
    rest2    <- satz2split[!(satz2split %in% rownames(tvectors))]
    
    if(length(used1) >1){satz1vec <- colSums(tvectors[used1,])}
    if(length(used1)==1){satz1vec <- tvectors[used1,]}
    
    if(length(used2) >1){satz2vec <- colSums(tvectors[used2,])}
    if(length(used2)==1){satz2vec <- tvectors[used2,]}
    
    cos      <- as.numeric(cosine(satz1vec,satz2vec))
    
    #   out <- list(cos=cos,used1=used1,used2=used2,rest1=rest1,rest2=rest2)
    #   print(out)
    
    cos
    
  }else{stop("tvectors must be a matrix!")}
  
}