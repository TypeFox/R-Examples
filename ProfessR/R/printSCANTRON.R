printSCANTRON<-function(B1)
  {
    ####### this fucntion only prints the results of students answers

    studans = B1$studans
    Nams = B1$Nams
    ids =  B1$ids
    Nstudents = dim(studans)[1]
       for(m in 1:Nstudents)
          {
            cat(paste(Nams[m], ids[m]) , sep=" ")
             cat(" ")
            cat(studans[m,], sep=" ")
            cat("\n")
          }
        

  }
