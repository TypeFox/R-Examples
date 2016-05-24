cghISBF <-
function(CGH.Array,chromosome,nucleotide.position,epsilon=0.05,K=1,impmin=1/100,s=NULL,v=NULL)
{
   Result = list(Esti.CopyN=CGH.Array,CGH.Array=CGH.Array,chromosome=chromosome,nucleotide.position=nucleotide.position,FDR=NULL)
   for (i in 1:23)
   {
      Y = CGH.Array[chromosome==i]
      longueur = length(Y)
      if (longueur<K) k = longueur else k = K
      Result$Esti.CopyN[chromosome==i] = isbf(Y,epsilon=0.05,K=k,impmin=impmin,s=s,v=v)$beta
   }
   return(Result)
}

