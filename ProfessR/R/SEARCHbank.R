SEARCHbank<-function(gw, y="humidity")
  {
#### search through a test bank and extract out questions
    ## that contain the keyword

    ##  the questions are printed to the screen,
    ####  and the index of the questions are returned
    k = grep(y, gw, ignore.case = TRUE)
    if(length(k)<1)
      {
        print("none")
      }
    else
      {
    for(ik in 1:length(k))
      {
        i = k[ik]
        b = gw[[i]]

        cat("##############\n")
        cat(b$Q, sep="\n")
        for(j in 1:length(b$A))
          {
            cat(b$A[j], sep="\n")
          }
        cat(paste("#####", b$a),sep="\n" )
      }
  }
    return(k)
  }
