`scramble.answers` <-
function(Qbank)
  {
    Q1 = Qbank
  for(i in 1:length(Qbank))
    {
      z = Qbank[[i]]
      A1 = substring(z$A, 3, 10000)
      A2 = paste("\\item {", A1, "}")

      LENA = length(A1)
      lets = paste(sep="", letters[1:LENA], ") ")
     #### paste(sep="", lets,  A1)

      isamp = sample(1:LENA)


      wans = which(isamp==z$numANS)

      pans = paste(sep="", lets,  A1[isamp])

      apans =  paste(sep="", "ANSWER: ",  pans[wans])


      #########  here must be worried about
      #####  answers that include "all of the above"
      #########   or "none of the above"
      ###   or an answer may refer to choices like a) and b) are correct.


      gexclude = c("all of the above",
        "none of the above",
        'None of these are correct',
        'all of the choices are correct',
        'All of the choices are correct',
        'Both choices are correct',
        'None of the choices are correct',
        'Both of the choices are correct',
        'All of these are correct',
        'Neither of these are correct')

      
      gskip = FALSE
      for(j in 1:length(gexclude))
        {
          gtest =grep(gexclude[j],  z$A, ignore.case =TRUE)
          if(length(gtest)>0)
            {
              gskip = TRUE
              break
            }
        }
  
     

      if(gskip)
        {
          
          Q1[[i]]$A = z$A
          Q1[[i]]$a = z$a
          
          
          
        }
      else
        {
          
####print(pans)
#### print(apans)
          
          Q1[[i]]$A = pans
          Q1[[i]]$a = apans
          
        }
    }
    
    return(Q1)
  }

