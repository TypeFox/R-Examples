 prep.breaks <- function(sh, breaks){
  nbr <- length(breaks)
  new.br <- breaks
  for(i in 1:nbr){
     br <- breaks[[i]]
   max.br <- max(br, na.rm=TRUE)   
   min.br <- min(br, na.rm=TRUE)   
   len.br <- length(br)
   eps <- min(diff(br),na.rm=TRUE)
   delta <- eps*sh
   if(delta>0) 
    new.br[[i]] <- c(min.br, (br+delta)[-c(len.br)], max.br)
   else	
    new.br[[i]] <- c(min.br, (br+delta)[-c(1)], max.br)
   }
   new.br
 }

shift.cem <- function(obj, data, shifts=NULL, verbose=0, plot=TRUE){
  if(class(obj) != "cem.match")
   stop("obj must be of class `cem.match'")

 if(is.null(shifts))
  return(obj)
 if(is.null(obj$breaks))
  return(obj) 
 nbr <- length(obj$breaks)
 old.br <- obj$breaks
 new.br <- old.br
 tmp <- NULL
 shifts <- c(shifts,-shifts)
 shifts <- sort(unique(shifts))

 
 for(sh in shifts){
  new.br <- prep.breaks(sh, obj$breaks)
  tab <- cem.main(treatment=obj$treatment,data=data, drop=obj$drop,cutpoints=new.br,k2k=obj$k2k)$tab
  if(verbose>1){
   cat(sprintf("\nShift=%f\n", sh))
   print(tab)
  }
  tmp <- rbind(tmp,c(tab[2,], sh))
 }
  idx <- which.max(tmp[,1]+tmp[,2])
  new.br <- prep.breaks(shifts[idx],obj$breaks)
 if(plot){
  plot(tmp[,3], tmp[,1], type="l",col="blue", lty=3, ylim=c(min(tmp[,1:2]), max(tmp[,1]+tmp[,2])),xlab="shift",ylab="matched")
  lines(tmp[,3], tmp[,2], type="l", col="red",lty=3)
  lines(tmp[,3], tmp[,1]+tmp[,2], type="l", col="green")
  abline(h=sum(obj$tab[2,]), lty=3, col="green")
  abline(h=obj$tab[2,1], lty=3, col="blue")
  abline(h=obj$tab[2,2], lty=3, col="red")
  abline(v=shifts[idx], lty=2)
 }
 newobj <- cem.main(treatment=obj$treatment,data=data, drop=obj$drop,cutpoints=new.br,k2k=obj$k2k)
 if(verbose>1){
  cat(sprintf("\nBest shift: %5.3f\nOld Match table:\n", shifts[idx]))
  print(obj$tab) 
  cat(sprintf("\nNew Match table:\n"))
  print(newobj$tab) 
 }
 return( invisible(newobj) )
}
