summary.HiddenF <-
function(object,method="HiddenF",...){
  if(method %in% c("HiddenF","HiddenF","HIDDENF")){
  a<-max(as.numeric(names(table(object$tall$rows))))
  b<-max(as.numeric(names(table(object$tall$cols))))
  ymtx <- matrix(object$tall$y,byrow=T,nrow=a,ncol=b)
  grpassn<-object$config.vector[object$tall$cols==1]
  cat('Number of configurations:',object$cc,'\n')
#  cat('First several pvalues testing for nonadditivity \n  after Bonferroni adjustment:',round(object$pvalue[1:6],4),'... \n\n')
grp1rows <- c(1:a)[grpassn==1]
grp2rows <- c(1:a)[grpassn==0]
grp1means <- apply(ymtx[grp1rows,,drop=FALSE],2,mean)
grp2means <- apply(ymtx[grp2rows,,drop=FALSE],2,mean)
  cat('Minimum adjusted pvalue:',object$adjpvalue,'\n\n')
  cat('Rows in group 1:',grp1rows, '\n')
  cat('Rows in group 2:',grp2rows, '\n\n')
  cat('Column means for grp 1:',grp1means,'\n')
  cat('Column means for grp 2:',grp2means,'\n\n')
  output <- list(group1=grp1rows,group2=grp2rows,grp1means=grp1means,grp2means=grp2means)
#  class(output) <- "summary.HiddenF"
  return(output)
}
# else{stop(cat("The summary function for the HiddenF package is supported \n only for the ACMIF/HiddenF method \n"))}
else{cat("The summary function for the hiddenf package is supported \n only for the ACMIF/HiddenF method \n")}
}
