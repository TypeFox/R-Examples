getPDEscreen<-function(pde = '/home/lees/Site/Santiaguito/pdq.eqs')
  {
    spde  = scan(file=pde, what="", sep="\n")

    kpde = list()
    for(i in 1:length(spde))
      {
        kpde[[i]]  = parse.pde(spde[i])

      }
    invisible(kpde)
  }



