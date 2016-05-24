`rename.answers` <-
function(Qbank, newnames=letters[1:26] , sep=") " )
  {
    if(missing(sep)) { sep = ") " }
    if(missing(newnames)) { newnames=letters[1:26] }
    

  ####   newnames=as.character(1:10)
    
    Q1 = Qbank
  for(i in 1:length(Qbank))
    {
      z = Qbank[[i]]
      A1 = substring(z$A, 3, 10000)
      
      LENA = length(A1)
      lets = paste(sep="", newnames[1:LENA], sep)
      
     #### paste(sep="", lets,  A1)
      pans = paste(sep="", lets,  A1)

      apans =  paste(sep="", "ANSWER: ",  pans[z$numANS])


          
          Q1[[i]]$A = pans
          Q1[[i]]$a = apans
      
    }
    
    return(Q1)
  }

