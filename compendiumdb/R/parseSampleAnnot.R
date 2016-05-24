##Function to split a vector of character strings with sample characteristics in
##a matrix with a column for each tag described in the original sample characteristics
parseSampleAnnot <- function (phenoData, colname)
{
  x <- phenoData
  check <- length(unlist(strsplit(x[1,colname],";;")))
  if(check > 1){
    headers <- character()
    for(j in 1:nrow(x)){
      #cat(".")
      out <- patternCheck(x[j,colname])
      n <- strsplit(unlist(strsplit(out,";;")),":")
      for(y in n){
        y <- unlist(y)
        if(length(y)==2){
          y <- y[seq(1,length(y),by=2)]
          headers=c(headers,y)
        }
      }
    }
    headers <- gsub(" ","_",headers)
    headers <- unique(headers)

    if(length(headers)!=0){
      sampleCharData <- as.data.frame(matrix(nrow=nrow(x),ncol=length(headers)))
      colnames(sampleCharData) <- headers
      
      for(j in 1:nrow(x)){
        #cat(".")
        out <- patternCheck(x[j,colname])
        n <- strsplit(unlist(strsplit(out,";;")),":")
        for(y in n){
          z <- unlist(y)
          label <- gsub(" ","_",z[1])
          ## remove leading space
          sampleCharData[j,label] <- sub("^ ","",z[2])
        }
      }
      rownames(sampleCharData) <- rownames(x)
      ## Two-channel experiment
      if (colname != "samplechar") colnames(sampleCharData) <- paste0(colnames(sampleCharData),"_",colname)			
      phenoData <- cbind(sampleCharData,x[,-which(colnames(x)==colname)])
    }else{phenoData <- x}
  }
  #cat(".")
  as.matrix(phenoData)
}

