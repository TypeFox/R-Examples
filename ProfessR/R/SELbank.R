SELbank<-function(QB, N, xclude=NULL)
  {
    ###   make random selection of  question bank
    ch = 1:length(QB)

    
    if(!is.null(xclude) & !identical(xclude, 0) )
      {
       ch =  ch[-xclude]

      }

    if(length(ch)<N) N = length(ch)
    
    ran1 = sample(ch, N )

    KBnew = vector(mode='list')

    for(i in 1:N)
      {
        KBnew[[i]] = QB[[ran1[i]]]
      }
    return(KBnew)
    
  }

