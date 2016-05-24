
################################################
##### Generic Summary ##########################

#' @export
#' @importFrom lsa cosine
#' @importFrom lsa textmatrix
#' @importFrom lsa lsa 
genericSummary <- function(text,k,split=c(".","!","?"),min=5,breakdown=FALSE,...){
 
  #### Decompose the document D into individual sentences
  
  sentences <- unlist(strsplit(text,split=split,fixed=T))
  
  if(breakdown==TRUE){sentences <- breakdown(sentences)}
  sentences <- sentences[nchar(sentences) > min]
  
  ## create temporary files
  
  td = tempfile()
  dir.create(td)
  
  ## storing "corpus"
  
  for(i in 1:length(sentences)){
    docname <- paste("sentence",i,".txt",sep="")
    write(sentences[i],file=paste(td,docname,
                                  sep="/"))
  }
  
  ## Construct the terms by sentences matrix A for the
  ## document D
  
  A        <- textmatrix(td,...)
  rownames <- rownames(A)
  colnames <- colnames(A)
  
  A           <- matrix(A,nrow=nrow(A),ncol=ncol(A))
  rownames(A) <- rownames
  colnames(A) <- colnames
  
  ## delete files again
  
  unlink(td,T,T)
  
  ## Perform the SVD on A to obtain the singular value
  ## matrix Sigma, and the right singular vector matrix V^T
  
  Vt <- lsa(A,dims=length(sentences))$dk
  
  ### Vt is V transposed, so the sentences are nor rows instead of columns!
  
  ## Select the k'th right singular vector from matrix V^T
  ## Select the sentence which has the largest index value
  ## with the k'th right singular vector, and include it in
  ## the summary.
  ##
  ## finding the sentence that
  ## has the largest index value with the k'th right singular vector
  ## is equivalent to finding the column vector (for Vt: row vector)  
  ## v_i whose k'th element v_ik is the largest
  ##
  ## 
  
  
  snum <- vector(length=k)
  
  for(i in 1:k){
    snum[i] <- names(Vt[,i][abs(Vt[,i]) == max(abs(Vt[,i]))])
  }
  
  snum <- gsub(snum,pattern="[[:alpha:]]",replacement="")
  snum <- gsub(snum,pattern="[[:punct:]]",replacement="")
  snum <- as.integer(snum)
  
  summary.sentences <- sentences[snum]
  summary.sentences
}
