`readUW.OSTAS` <-
function(OS1)
  {

    OSTAS = vector()
    littlek = 0
    for(j in 1:length(OS1))
      {
        OS2=unlist(strsplit(split=" ", OS1[j]))
        for(m in 2:length(OS2)) { littlek=littlek+1;  OSTAS[littlek]=OS2[m] }
      }
  
    invisible(OSTAS)

  }

