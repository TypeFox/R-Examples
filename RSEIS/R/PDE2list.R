PDE2list<-function(PDF)
  {
    invisible(list(
                   yr= as.vector(unlist(lapply(PDF, getmem, 1))),
                   jd= as.vector(unlist(lapply(PDF, getmem, 2))),
                   mo= as.vector(unlist(lapply(PDF, getmem, 3))),
                   dom= as.vector(unlist(lapply(PDF, getmem, 4))),
                   hr= as.vector(unlist(lapply(PDF, getmem, 5))),
                   mi= as.vector(unlist(lapply(PDF, getmem, 6))),
                   sec= as.vector(unlist(lapply(PDF, getmem, 7))),
                   lat= as.vector(unlist(lapply(PDF, getmem, 8))),
                   lon= as.vector(unlist(lapply(PDF, getmem, 9))),
                   depth= as.vector(unlist(lapply(PDF, getmem, 10))),
                   z= as.vector(unlist(lapply(PDF, getmem, 11))),
                   mag= as.vector(unlist(lapply(PDF, getmem, 12))))
              )


  }
