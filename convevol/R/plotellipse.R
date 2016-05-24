#'Plots an ellipse    
#'
#'Plots a minimum ellipse around a set of data  
#'
#'@param ellipse  Gives the parameters of the ellipse - output from the ellipsoidhull functon in cluster.  
#'
#'@details  Routine adapted from a suggestion made on CrossValidated:  http://stats.stackexchange.com/questions/9898/how-to-plot-an-ellipse-from-eigenvalues-and-eigenvectors-in-r
#'
#'@return Nothing - just plots the ellipse.  
#'
#'@import ape 
#'
#'@export
#'
#'@references Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2013).
#'cluster: Cluster Analysis Basics and Extensions. R package version 1.14.4.



plotellipse<-function(ellipse)

{

d<-ellipse$d2
c<-ellipse$loc[1:2]
m<-ellipse$cov[1:2,1:2]

RR<-chol(m)                               	# Cholesky decomposition
angles<-seq(0, 2*pi, length.out=200)          	# angles for ellipse
ell<-(d^1.5)*cbind(cos(angles),sin(angles))%*%RR  	# ellipse scaled with factor 1
ellCtr<-sweep(ell,2,c,"+")               		# center ellipse to the data centroid
lines(ellCtr,type="l",lwd=2,asp=1,col='purple')            	# plot ellipse

}
