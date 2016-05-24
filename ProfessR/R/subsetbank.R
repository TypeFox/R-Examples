subsetbank<-function(QBANK, sel)
  {
    ### given a question bank and vector of selections
    ###   return selected questions
    V = vector(mode="list")

    for(i in 1:length(sel))
      {
        V[[i]] = QBANK[[sel[i]]]



      }
    return(V)

  }
