#' To show the missing values of a dataset
#'@description Plot the dataset with marks where there are missing value. It allows to have a quick idea of the structure of missing values (Missing at Random or not for example).
#'@param X the matrix to analyse (matrix with missing values or correlations matrix)
#'@param what indicates what to plot. If \code{what="correl"} and \code{X} is a correlation matrix then the plot is a correlation plot. Else it shows the missing values positions in the dataset.
#'@param pch for missing, symbol to plot (can set pch="." for large datasets)
#'@export
#'@examples
#'\dontrun{
#'    data<-mtcars
#'    require(CorReg)
#'   datamiss=Terminator(target = data,wrath=0.05)#5% of missing values
#'   showdata(datamiss)#plot positions of the missing values
#'   
#'   #missing values with a structure
#'   datamiss=Terminator(target = data,diag=1)#diag of missing values
#'   showdata(datamiss)#plot positions of the missing values (no full individuals, no full variable)
#'   
#'   opar=par(no.readonly = TRUE)
#'   showdata(X=cor(data),what="correl")
#'   par(opar)

#'
#'}
#'
#'

showdata<-function(X=X,what=c("miss","correl"),pch=7){
   what=what[1]
   if(what=="miss"){
      M=which(is.na(X),arr.ind=T)
      if(nrow(M)>1){
         plot(M[,c(2,1)],pch=pch)
         title("Missing values in the dataset")  
      }else{
         print("No missing values")
      }
   }else if (what=="correl"){
      corrplot(corr=X,addrect=NULL,is.corr=T,method="color",tl.pos="n",diag=F,outline=F)
   }else{
      correl=cor(X[,!is.na(colSums(X)) & apply(X,2,sd)!=0])
      corrplot(corr=correl,addrect=NULL,is.corr=T,method="color",tl.pos="n",diag=F,outline=F)
   }
}
