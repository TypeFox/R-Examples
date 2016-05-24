#<<BEGIN>>
print.mccut <- function(x, lim=c(0.025,0.975), digits=3,...)
#ISALIAS evalmccut
#--------------------------------------------
#
{

  summ <- function(x) c(mean=mean(x,na.rm=TRUE),quantile(x,probs=c(0.5,lim),na.rm=TRUE),Nas=sum(is.na(x)))[c(2,1,3:(length(lim)+3))]
  nbl <- length(x)

  summarize <- function(x){
   if(is.list(x)){
      for(i in 1:length(x)){
        x[[i]] <- summarize(x[[i]])
      }
   }
   else if(inherits(x,"mcnode")) return(x)
   else if(is.numeric(x)){
      dimm <- 1:length(dim(as.array(x)))
      dimm <- dimm[-2]
      x <- apply(x,dimm,summ)

   x <- drop(x)
   if(is.vector(x)) {
    x <- as.matrix(x)
    colnames(x) <- "NoVar"}
   }
   return(x)}

  x <- summarize(x)

  class(x) <- "listof"
  print(x,digits=digits,...)
  return(invisible(x))}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

