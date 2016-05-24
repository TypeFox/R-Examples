scaling <-
function(Y, type="none")
{
   if(type == "pareto")
   {
   	 sdev<-as.matrix(apply(Y,2,sd))
     Y<-sweep(Y,2,sqrt(sdev),"/")             ## Pareto scale Y.
    }
    if(type == "unit")
    {
     Var<-as.matrix(apply(Y,2,var))
     Y<-sweep(Y,2,sqrt(Var),"/")              ## Unit scale Y.
    }
   return(Y)
}

