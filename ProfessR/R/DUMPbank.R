DUMPbank<-function(ofile, QB, sep="\n", append=TRUE)
{
  if(missing(sep)) sep=NULL
  if(missing(append)) append=TRUE
  
#########  dump out an ascii version of an exam question bank

if(append==FALSE)
  {
    gg = deparse(substitute(QB))

    cat(file=ofile,paste(sep="", "############## ", gg) , append=FALSE, sep="\n")
  }

  
  
  for(i in 1:length(QB))
    {
      ###print(i)
      Q1 = NULL
      Q1 = QB[[i]]
      cat(file=ofile,"QUESTION: ", append=TRUE, sep="")
      cat(file=ofile,Q1$Q, append=TRUE, sep="\n")
      for(j in 1:length(  Q1$A))
        {
          if( identical(j, Q1$numANS)  ) { cat(file=ofile,"ANSWER: ", append=TRUE, sep="") }
          cat(file=ofile,Q1$A[j], append=TRUE, sep="\n")
        }
      if(length(Q1$FIG)>0)
        {
         
              cat(file=ofile,"FIG: ", append=TRUE, sep="")
              cat(file=ofile, paste(sep=" ", Q1$FIG$fn ,Q1$FIG$tag)    , append=TRUE, sep="\n")
        
        }
      
      
      if(!is.null(sep)) cat(file=ofile,"", append=TRUE, sep=sep)

    }
#######  example: DUMPbank("dafinal", QBFINAL)
#######  

}

