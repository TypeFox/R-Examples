isbf <-
function(Y,epsilon=0.05,K=1,impmin=1/100,s=NULL,v=NULL)
{

   # DEFINITION DES PARAMETRES DE BASE

   p = length(Y)

   if (is.null(v))
   {
      if (p>=12)
      {
         Ybis = Y
         for (i in 5:(p-4)) Ybis[i] = mean(Y[(i-4):(i+4)])
         v = var(Y-Ybis)
      } else {
         v = var(Y)
      }
   }

   if (K>p)
   {
      print("Error: K>p")
      return(NULL)
   }

   if (is.null(s)) s = -sqrt(v)*qnorm(epsilon/(K*(2*p-K+1)))

   beta = rep(0,p)

   residu = matrix(data=0,nrow=p,ncol=K)
   seuille = residu
   residu[,1] = Y


   RESULT=.C("isbfC",residu=as.double(residu),seuille=as.double(seuille),beta=as.double(beta),s=as.double(s),impmin=as.double(impmin),p=as.integer(p),K=as.integer(K),PACKAGE="ISBF")
   return(list(beta=RESULT$beta,s=RESULT$s,impmin=RESULT$impmin,K=RESULT$K))

}

