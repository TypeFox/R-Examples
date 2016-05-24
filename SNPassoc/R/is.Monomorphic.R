`is.Monomorphic` <-
function (x) 
{
   ans<-FALSE   
   if (length(x)==1)
    {
    if (x[1] == "Monomorphic") 
        ans <- TRUE
    }
   else
      ans <- length(table(x)[table(x) > 0]) == 1
   return(ans)
}

