
randomCorrelations <- function(Data1, Data2=NULL, n.items, R=1000, complete.overlap=FALSE, item.overlap=FALSE, trait.overlap=TRUE, Key1=NULL, Key2=NULL){
   
  if(!is.null(Key1) & !is.null(Key2)) item.overlap <- TRUE
  if(is.null(Data2)) Data2 <- Data1
  if(is.null(Key2)) Key2 <- Key1
  if(length(Key1) == length(Key2)) colnames(Data2) <- colnames(Data1)
  if(length(Key1) != length(Key2)) complete.overlap <- FALSE 
  if(FALSE %in% c(unique(Key1) %in% unique(Key2), unique(Key2) %in% unique(Key1))) trait.overlap <- TRUE
  if(complete.overlap==TRUE) trait.overlap <- item.overlap <- TRUE
  
  if(trait.overlap == FALSE & is.null(Key1)) stop("\nKey1 not provided\n\n")
  if((!is.null(Key1) & length(Key1)!= ncol(Data1)) | (!is.null(Key2) & length(Key2)!= ncol(Data2))) stop("\nThe length of Key1 or Key2 does not match the number of items\n\n")
    
  results <- vector()
  
  for(a in 1:R) {
    
    if(trait.overlap==FALSE) {
      K <- as.integer(length(unique(Key1)) / 2)
      draw1 <- sample(unique(Key1), size= length(unique(Key1)) - K)
      draw2 <- sample(unique(Key2)[!unique(Key2) %in% draw1], size=K, replace=T)
      draws <- list(draw1 = draw1, draw2 = draw2)
    }
   
    colnames1 <- colnames(Data1)
    if(trait.overlap==FALSE)  colnames1 <- colnames1[Key1 %in% draws$draw1]
    colnames2 <- colnames(Data2)
    if(trait.overlap==FALSE)  colnames2 <- colnames2[Key2 %in% draws$draw2]
    
    sampler1 <- sample(colnames1, size=n.items)
    sampler2 <- sample(colnames2, size=n.items)
    if(item.overlap==FALSE) sampler2 <- sample(colnames2[!colnames2 %in% sampler1], size=n.items)  
    if(complete.overlap==TRUE) sampler2 <- sampler1
      
    results[a] <- cor(rowMeans(Data1[,sampler1], na.rm=T), rowMeans(Data2[,sampler2], na.rm=T), use="complete.obs")
    
  }
  
  return(results)
  
}


