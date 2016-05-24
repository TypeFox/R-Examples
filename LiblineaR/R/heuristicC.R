#' Fast Heuristics For The Estimation Of the C Constant Of A Support Vector Machine. 
#' 
#' \code{heuristicC} implements a heuristics proposed by Thorsten Joachims in
#' order to make fast estimates of a convenient value for the C constant used by
#' support vector machines. This implementation only works for linear support
#' vector machines.
#' 
#' @return 	A value for the C constant is returned, computed as follows:\cr
#' \eqn{\frac{1}{\frac{1}{n}\sum_{i=1}^{n}\sqrt{G[i,i]}}}{1/(1/n Sum_i=1:n sqrt(G[i,i]))}
#' where
#' \eqn{G=\code{data}\%*\%t(\code{data})}{data \%*\% t(data)}
#' 
#' @param data a nxp data matrix. Each row stands for an example (sample, point)
#'   and each column stands for a dimension (feature, variable)
#' 
#' @references
#' 	\itemize{
#' \item 
#' T. Joachims\cr
#' \emph{SVM light} (2002)\cr
#' \url{http://svmlight.joachims.org}
#' }
#' 
#' @author Thibault Helleputte \email{thibault.helleputte@@dnalytics.com}
#' 
#' @note Classification models usually perform better if each dimension of the
#'   data is first centered and scaled. If data are scaled, it is better to
#'   compute the heuristics on the scaled data as well.
#' 
#' @seealso \code{\link{LiblineaR}}
#' 
#' @examples
#' data(iris)
#' 
#' x=iris[,1:4]
#' y=factor(iris[,5])
#' train=sample(1:dim(iris)[1],100)
#' 
#' xTrain=x[train,]
#' xTest=x[-train,]
#' yTrain=y[train]
#' yTest=y[-train]
#' 
#' # Center and scale data
#' s=scale(xTrain,center=TRUE,scale=TRUE)
#' 
#' # Sparse Logistic Regression
#' t=6
#' 
# Tune the cost parameter of a logistic regression according to the Joachim's heuristics
#' co=heuristicC(s)
#' m=LiblineaR(data=s,labels=yTrain,type=t,cost=co,bias=TRUE,verbose=FALSE)
#' 
#' # Scale the test data
#' s2=scale(xTest,attr(s,"scaled:center"),attr(s,"scaled:scale"))
#' 
#' # Make prediction
#' p=predict(m,s2)
#' 
#' # Display confusion matrix
#' res=table(p$predictions,yTest)
#' print(res)
#' 
#' # Compute Balanced Classification Rate
#' BCR=mean(c(res[1,1]/sum(res[,1]),res[2,2]/sum(res[,2]),res[3,3]/sum(res[,3])))
#' print(BCR)
#' 
#' @keywords classif

heuristicC<-function(data){
	gram=data%*%t(data)
	n=dim(gram)[1]
	kxixi=matrix(ncol=n,nrow=1)
	for(i in 1:n){
		kxixi[1,i]=sqrt(gram[i,i])
	}
	m=mean(kxixi)
	C=1/m
	return(C)
}

