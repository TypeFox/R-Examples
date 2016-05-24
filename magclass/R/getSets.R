getSets <- function(x,fulldim=TRUE,sep=".") {
  out <- names(dimnames(x))[drop=FALSE]
  if(is.null(out)) return(NULL)
  
  if(fulldim==TRUE){
    tmp<- strsplit(out,split=sep,fixed=TRUE)
    tmp<- lapply(tmp,FUN=function(x){
                            if(length(x)==0) x<-NA
                            return(x)
                          }
    )
    out <- as.vector(unlist(tmp))
  }
  return(out)
}




"getSets<-" <- function(x,fulldim=TRUE,sep=".",value) {
   x <- clean_magpie(x,what="sets")
   if(is.null(value)) return(x)
   if(length(names(dimnames(x)))==0) fulldim <- FALSE  
   if(length(value)==3) fulldim <- FALSE
   if(length(value)==0) fulldim <- FALSE
   if(!fulldim) {
     names(dimnames(x)) <- value
     return(x)
   } else {
     s1 <- getSets(x,fulldim=FALSE)
     s2 <- getSets(x,fulldim=TRUE)
     search_s2 <- paste0("(^|\\.)",s2,"(\\.|$)")
     where <- sapply(search_s2,grep,s1)
     names(where) <- s2
     
     if(length(value)!=length(s2)) stop("Input length does not agree with the number of sets in x!")
     
     for(i in 1:3) {
      s1[i] <- paste(value[where==i],collapse=sep)
     }
     getSets(x,fulldim=FALSE,sep=sep) <- s1
     return(x)
   }  
}
