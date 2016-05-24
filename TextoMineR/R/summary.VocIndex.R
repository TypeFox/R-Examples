summary.VocIndex <-
function (object, nword=20, ...) 
{
   res <- object
     if (!inherits(res, "VocIndex")) 
        stop("non convenient object")
     cat("\nVocIndex summary\n")
     cat("\nUse of vacabulary\n")
     Nword1<-min(nword,nrow(res$RegVoc))
	cat("\n",Nword1,"Regular Words\n")
       print(res$RegVoc[c(1:Nword1),])
      Nword2<-min(nword,nrow(res$LocalVoc))
	cat("\n",Nword2,"specialized Words\n")
       print(res$LocalVoc[c(1:Nword2),])
}
