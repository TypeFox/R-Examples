gradeSCAN<-function(j, key)
  {
    #######  given a scantron record and the correct key
    #########  cycle through the scores and grade each exam
    ####  return a list of scores for each student
    
    dj = dim(j)
    Nquestions = dj[2] 
    Nstudents = dj[1] 

    
    bscore = rep(0, times=Nstudents)
    for(i in 1:Nstudents)
      {
        L =  length(which(key==j[i,]))
        bscore[i] = L
      }

    return(bscore)
  }

