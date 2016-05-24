##### Generate Targets

#' @export
#' @importFrom lsa cosine 

choose.target <- function(x,lower,upper,n,tvectors=tvectors,
                          breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){ 
    
    allwords <- vector(length=nrow(tvectors))
    
    if(breakdown==TRUE){satz1 <- breakdown(x)}
    if(breakdown==FALSE){satz1 <- x}
    
    satz1split <- strsplit(satz1,split=" ")[[1]]
    
    used1     <- satz1split[satz1split %in% rownames(tvectors)]
    if(length(used1)==0){(warning("no element of x found in rownames(tvectors)"))
                         return(NA)}
    
    if(length(used1) >1){satz1vec <- colSums(tvectors[used1,])}
    if(length(used1)==1){satz1vec <- tvectors[used1,]}
    
    for(i in 1:nrow(tvectors)){
      
      allwords[i] <- cosine(satz1vec,tvectors[i,])
    }     
    
    names(allwords) <- rownames(tvectors)
    
    a <- sample(allwords[allwords >= lower & allwords <= upper])[1:n]
        
    print(a)
    
  }else{warning("tvectors must be a matrix!")}
  
}
