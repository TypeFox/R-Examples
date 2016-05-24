isbfReg <-
function(X,Y,epsilon=0.05,K=1,impmin=1/100,favgroups=0,centX=TRUE,centY=TRUE,s=NULL,v=NULL)
{

   D = dim(X)
   n = D[1]
   p = D[2]

   dimension = K*p - K*(K-1)/2

   if (n!=length(Y))
   {
      print("ERROR: The dimension of X should be (n,p) and the dimension of Y (n,1)")
      return(NULL)
   }

   if (is.null(v)) v = var(Y)/2

   if (K>p)
   {
      print("Error: K>p")
      return(NULL)
   }

   if (is.null(s)) s = -sqrt(v)*qnorm(epsilon/(2*dimension))

   if (centX==TRUE) for (cptcent in 1:p) X[,cptcent] = X[,cptcent] - mean(X[,cptcent])
   if (centY==TRUE) Y = Y-mean(Y)

   beta = rep(0,p)
   COV = rep(0,dimension)

   RESULT=.C("isbfRegC",X=as.double(X),Y=as.double(Y),COV=as.double(COV),beta=as.double(beta),s=as.double(s),impmin=as.double(impmin),fgroups=as.double(favgroups),p=as.integer(p),n=as.integer(n),K=as.integer(K),PACKAGE="ISBF")
   return(list(beta=RESULT$beta,s=RESULT$s,impmin=RESULT$impmin,K=RESULT$K))

}

