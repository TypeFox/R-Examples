editDB<-function(DB, w)
  {
##############   edit the data base by removing or selecting parts
    NDB = DB
    len1 = length(DB$yr)


    lens =   lapply(DB, length)

    w2 = which(lens == len1)

    for(i in 1:length(w2))
      {
        k = DB[[w2[i] ]]
        NDB[[w2[i] ]] = k[w] 
        
      }

    invisible(NDB)


  }
