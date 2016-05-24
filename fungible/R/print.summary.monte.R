print.summary.monte<-function(x,digits=3, ...){
z<-x
nclus<-nrow(z$obs.skew)
nvar<-ncol(z$obs.skew)

obs.cor.list<-z$obs.cor
obs.skew.matrix<-obs.kurt.matrix<-matrix(0,nclus,nvar)                                  
                

cat("\nCall monte:","\n")
print(z$call)

cat("\nSeed = ",z$seed,"\n")

cat("\n\nNumber of objects in each group:")
for(i in 1:nclus){
 cat("\nGroup", i,"=",z$clus.size[i])
 }


##---------------------------------------------
##          GROUP CENTROIDS
##---------------------------------------------

cat("\n\n\nGroup Centroids\n")
for(i in 1:nclus){
  cat("Group ",i, "  ")
  cat(round(z$centroids[i,],digits),"\n")
}


##---------------------------------------------
##          WITHIN GROUP VARIANCES
##---------------------------------------------

var.matrix<-matrix(99,nclus,nvar)
cat("\n\nWithin Group Variances\n")
for(i in 1:nclus){
  cat("Group ",i, "  ")
  cat(round(z$var.matrix[i,],digits),"\n")
}

grp.pairs <- expand.grid(1:nclus,1:nclus)
grp.pairs <- grp.pairs[grp.pairs[,1]!=grp.pairs[,2],]
grp.pairs <- grp.pairs[sort.list(grp.pairs[,1]),]

cat("\n\nRatio of Within Group Variances")
for(i in 1:length(grp.pairs[,1])){
   cat("\n\nGroup ",grp.pairs[i,1], "/ Group", grp.pairs[i,2], "  ")
   cat(round( var.matrix[ grp.pairs[i,1], ]/ var.matrix[ grp.pairs[i,2], ], digits))
}   


cat("\n\n\nExpected within group correlations:","\n")

for(i in 1:nclus){
  cat("\nGroup ",i,"\n")
  if(is.null(z$cor.list)) 
    cat("Identity Matrix\n")
  else
  print(round(z$cor.list[[i]],digits))
}


cat("\n\nObserved within group correlations:","\n")
for(i in 1:nclus){
  cat("\nGroup ",i)
  print(round(z$obs.cor[[i]],digits))
}

cat("\nExpected within group indicator skewness:","\n")
  if(is.null(z$skew.list))
    cat("All expected skew = 0\n")
  else {
    ex.sk.mat<-matrix(unlist(z$skew.list),nclus,nvar,byrow=TRUE)    
    dimnames(ex.sk.mat)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
    print(ex.sk.mat)
  }  



cat("\nObserved within group indicator skewness:","\n")
dimnames(z$obs.skew)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
dimnames(z$obs.kurt)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
print(round(z$obs.skew,digits))

cat("\n\nExpected within group indicator kurtosis:","\n")
  if(is.null(z$kurt.list))
     cat("All expected kurtosis = 0\n")
  else{ 
    ex.kt.mat<-matrix(unlist(z$kurt.list),nclus,nvar,byrow=TRUE)
    dimnames(ex.kt.mat)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
    print(ex.kt.mat)
  }  

cat("\nObserved within group indicator kurtosis:","\n")
print(round(z$obs.kurt,digits))



##----------------------------------------------------------
##                   VALIDITIES
##----------------------------------------------------------

     cat("\nObserved indicator validities:","\n")
     if( is.character(z$validities) ) 
         print(z$validities)
     else
         print(round(z$validities,digits)) 
 
 
 ##--------------------------------------------------------
 ##            Total Sample Results
 ##------------------------------------------------------
 if(z$Total.stats){
 cat("\n\nStats for Total Sample\n")
 cat("\nTotal sample correlation matrix\n")
 print(z$Total.cor,digits)
 
 cat("\nTotal sample indicator skewness\n")
 cat( round(z$Total.skew,digits) )
 
 
 cat("\n\nTotal sample indicator kurtosis\n")
 cat(round(z$Total.kurt,digits) )
 cat("\n\n")
 }
 



}
