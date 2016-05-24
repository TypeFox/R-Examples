summary.quint <-
function(object,...,digits=2){
#digits=number of digits at decimal points
if(object$crit=="es"){
  cat("Partitioning criterion: Effect size criterion","\n","\n")}
else{ 
  cat("Partitioning criterion: Difference in treatment means criterion","\n","\n")
  }
         cat("Fit","information:", "\n") 
    if (dim(object$fi)[2]==6){
     cat(c(rep("",22),"Criterion") ,"\n" )
     cat(c(rep("",16),paste(rep("-",15),sep="")),"\n")
        print( round(object$fi[,c(1:4,6)],digits=digits),row.names=FALSE ) }
     else {   
       cat(c(rep("",15),"Criterion") ,"\n" )
     cat(c(rep("",16),paste(rep("-",4),sep="")),"\n")
     print( round(object$fi[,c(1:3)],digits=digits), row.names=FALSE) }
     cat("\n")
     
     object$si<-cbind(object$si[,1:3], splitpoint=round(object$si[,4],digits=digits))
      
      cat("Split information:","\n")
      print(object$si,row.names=TRUE)
      cat("\n")
     cat("Leaf information:","\n")
     print(round(object$li[,c(2:10)],digits=digits))
     }
