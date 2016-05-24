getKEY<-function(fn)
  {
#####   get the key data from ProfessR generated exam
    key = scan(file=fn, skip=2, what=list(q1=0, j1="", orig=0, j2="", ans=0), flush=TRUE)

    return(key$ans)


  }
