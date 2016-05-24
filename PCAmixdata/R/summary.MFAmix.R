summary.MFAmix <-function(object, ...){
  x <- object
  if (!inherits(x, "MFAmix")) 
    stop("use only with \"PCAmix\" objects")
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  n <- x$global.pca$rec$n
  p1 <- x$global.pca$rec$p1
  p <- x$global.pca$rec$p
  p2 <- p-p1
  m<- ncol(x$global.pca$rec$G)
  
  Tab.recap<-data.frame(matrix(0,nrow=length(x$separate.analyses),ncol=3))
  colnames(Tab.recap)<-c("Number of numerical variables","Number of categorical variables","Total number of categories")
  rownames(Tab.recap)<-names(x$separate.analyses)
  for (g in 1:length(x$separate.analyses)){    
    if (!is.null(nrow(x$separate.analyses[[g]]$quanti$coord))){
      Tab.recap[g,1]<-nrow(x$separate.analyses[[g]]$quanti$coord)
    }
    if (!is.null(nrow(x$separate.analyses[[g]]$quali$contrib.pct))){
      Tab.recap[g,2]<-nrow(x$separate.analyses[[g]]$quali$contrib.pct)
      Tab.recap[g,3]<-nrow(x$separate.analyses[[g]]$levels$coord)
    } 
  }
  Tab.recap<-t(Tab.recap)
  
  
  cat("\n")
  cat("\n")
  cat("Data:", "\n")
  cat(paste("   number of observations: ",n),sep=" ") 
  cat("\n")
  if   ((p1!=0)&& (p2==0)) {
    cat(paste("   number of numerical variables: ",p1),sep=" ")  
    cat("\n")
  }
  if   ((p1==0)&& (p2!=0)) {
    cat(paste("   number of categorical variables: ",p2),sep=" ")  
    cat("\n")
    cat(paste("       total number of categories: ",m),sep=" ")  
  }
  if 	((p1!=0)&& (p2!=0)) {
    cat(paste("   number of  variables: ",p),sep=" ")
    cat("\n")
    cat(paste("        number of numerical variables: ",p1),sep=" ")   
    cat("\n")
    cat(paste("        number of categorical variables: ",p2),sep=" ")   
    cat("\n")
    cat(paste("           total number of categories: ",m),sep=" ")  
    
  }
  cat("\n")
  cat("\n")
  cat("\n")
  cat("Summary by groups :","\n")
  print(Tab.recap)
}  