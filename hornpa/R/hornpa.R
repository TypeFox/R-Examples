#'Perform Horn's Parallel Analysis using Simulated Data
#'
#'User specifies the number of variables in the synthetic dataset, the number of observations, the number of correlation matrices to generate, and the percentile.
#'@description A stand alone function to run a parallel analysis. The program generates a specified number of datasets based on the number of variables entered by the user. Eigenvalues are computed for the each dataset and collected.
#'The mean and the specified percentile (95th is the default) are computed. The output table shows how large eigenvalues can be as a result of 
#'merely using randomly generated datasets. If the user's own dataset has an actual eigenvalue greater than the generated eigenvalue,
#'that lends support to retain that factor. In other words, if the i(th) eigenvalue from the actual data was larger than the percentile of the (i)th eigenvalue generated using
#'randomly generated data, empirical support is provided to retain that factor. Horn's (1965) article suggested using the mean vs. the 95th percentile commonly used today (O'Conner, 2000).
#'
#'If principal axis factoring is specified, squared multiple correlations are substituted in the diagonal
#'of the correlation matrix. 
#'@param k The number of variables in your dataset.
#'@param size The number of observations in your dataset.
#'@param reps How many simulated datasets to create.
#'@param seed Use an optional seed to reproduce the same results each time.
#'@param p Specify the percentile which to evaluate the eigenvalues. Defaults to the 95th percentile.
#'@param method Specify what type of analysis is run (using a principal component analysis or principal axis factoring). Defaults to pca if not specified.
#'@return Mean and the specified percentile eigenvalues based on a user defined set of random datasets. Compare the eigenvalues from your dataset to the eigenvalues generated with the synthetic datasets. Results are comparable to the SPSS and SAS syntax written by O'Connor (2000).
#'@examples hornpa(k=10,size=100,reps=200)
#' hornpa(k=15,size=100,reps=200,method="pa",seed=123)
#' hornpa(k=12,size=200,reps=100,p=.99)
#' @references Horn, J. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 32, 179-185.
#' 
#' O'Connor, B. (2000). SPSS and SAS programs for determining the number of components using parallel analysis and Velicer's MAP test, Behavior Research Methods, Instruments and Computers, 32, 396-402.
#' 
#' 
#'@export
hornpa<-function(k,size=100,reps=100,seed=NA,p=.95,method=c("pca","pa")){

method<-match.arg(method)
dat<-matrix(NA,nrow=size,ncol=k)
ev<-matrix(NA,nrow=reps,ncol=k)
if(is.na(seed)){xx<-1 
                seed<-(abs(rnorm(1))*100)}
xx<-ifelse(exists("xx"),xx,2)
set.seed(seed)
for (r in 1:reps){
for (i in 1:k) {
     dat[,i]<-rnorm(size,0,1)
}

cm<-cor(dat)
if(method=="pa"){
     smc<-1-(1/diag(solve(cm)))
     diag(cm)<-smc
}
v<-eigen(cm)
v$values
ev[r,1:k]<-v$values
}

out<-matrix(rep(1:k),nrow=k,ncol=3)
out[,2]<-colMeans(ev)
for (x in 1:k){
     out[x,3]<-quantile(ev[,x],probs=c(p))
}


cat("\n","Parallel Analysis Results", " \n","\n")
cat("Method:",method,"\n")
cat("Number of variables:",k, "\n")
cat("Sample size:",size,"\n")
cat("Number of correlation matrices:",reps,"\n")
if(xx != 1) {cat("Seed:",seed,"\n")}
cat("Percentile:",p,"\n","\n")


out2<-data.frame(round(out,3))
names(out2)<-c("Factor","Mean",p)
if(method=="pca"){
     names(out2)<-c("Component","Mean",p)
}
cat("Compare your observed eigenvalues from your original dataset to the",p*100,"percentile in the table below generated using random data. If your eigenvalue is greater than the percentile indicated (not the mean), you have support to retain that factor/component.","\n \n")
print(out2,row.names=F)
}


