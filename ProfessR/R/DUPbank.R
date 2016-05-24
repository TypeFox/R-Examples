DUPbank<-function(Qbank)
  {

    lens = unlist(lapply(Qbank, "length")) 
    thelens = as.numeric(lens)
    thefiles = names(lens) 

    allqs = vector()
    allind = vector()
    internalind = vector()
    ifile = vector()
    nfile = vector()


    for(i in 1:length(Qbank))
      {
        uq = unlist(Qbank[[i]])
        w = which(names(uq)=="Q")
        nu = seq(from=1, to=length(w))

        allqs =c(allqs, as.vector(uq[w])  )

        internalind = c(internalind, nu)

        allind =  c(allind, w)

        ifile = c( ifile,  rep( thefiles[i], times=length(w)  ))
        nfile = c( nfile,  rep( i, times=length(w)  ))

      }


    wdup = which( duplicated(allqs)  )
    if(length(wdup)<1) { return(NULL) }

#######   cbind( allqs[wdup], ifile[wdup], internalind[wdup] )
    return(list(A=allqs[wdup], F=ifile[wdup], I=internalind[wdup], N=nfile[wdup]  )  )

  }
#########   DQ = DUPbank(Qbank)

#####  source("/home/lees/R_PAX/ProfessR/R/DUPbank.R")

