#<<BEGIN>>
summary.mccut <- function(object, lim=c(0.025,0.975),...)
#ISALIAS summary.mc
#--------------------------------------------
{
  summ <- function(x) c(mean=mean(x,na.rm=TRUE),quantile(x,probs=c(0.5,lim),na.rm=TRUE),Nas=sum(is.na(x)))[c(2,1,3:(length(lim)+3))]

  quel <- which(sapply(object,inherits,what="summary.mccut"))
  lquel <- length(quel)
  if(lquel == 0) stop("summary.mc was not evaluated in evalmccut : no summary to produce")
  if(lquel > 1) stop("More than one summary.mc was evaluated in evalmccut : impossible to produce a summary")
  
  object <- object[[quel]]
  typen <- sapply(object,"attr",which="type")

  l <- length(object)

  LESSTAT <- function(object,typen){
  
    if(is.list(object)) return(mapply(LESSTAT,object,typen,SIMPLIFY=FALSE))

    if(typen =="0") {
      object <- as.matrix(object[1])
      dimnames(object) = list("NoInc","NoVar")
      return(object)}

    if(typen =="V") {
      object <- t(as.matrix(object[,1,]))
      rownames(object) <- "NoUnc"
      return(object)}

    object <- drop(apply(object,c(1,3),summ))

    if(typen =="U") {
      object <- as.matrix(object)
      colnames(object) <- "NoVar"}
   return(object)
   }

   object <- mapply(LESSTAT,object,typen,SIMPLIFY=FALSE)

  object <- mapply("attr<-",object,"type",typen,SIMPLIFY=FALSE)
	class(object) <- c("summary.mc","listof")
  return(object)
}
