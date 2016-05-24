COMPbank<-function(Qbank1,  Qbank2)
  {
##########  compare two question banks
    ###  return a vector of indeces
    ###  of questions in Qbank2 that are not in Qbank1 

    unpackbank<-function(QBONE)
      {
        lens1 = unlist(lapply(QBONE, "length")) 
        thelens1 = as.numeric(lens1)
        thefiles1 = names(lens1) 

        allqs1 = vector()
        allind1 = vector()
        internalind1 = vector()
        ifile1 = vector()
        nfile1 = vector()


        for(i in 1:length(QBONE))
          {
            uq = unlist(QBONE[[i]])
            w = which(names(uq)=="Q")
            nu = seq(from=1, to=length(w))

            allqs1 =c(allqs1, as.vector(uq[w])  )

            internalind1 = c(internalind1, nu)

            allind1 =  c(allind1, w)

            ifile1 = c( ifile1,  rep( thefiles1[i], times=length(w)  ))
            nfile1 = c( nfile1,  rep( i, times=length(w)  ))

          }
        return(list(QS=allqs1))
      }



    qb1 = unpackbank(Qbank1)
    qb2 = unpackbank(Qbank2)

#### m12 =   match(qb1$QS, qb2$QS)

    m21  =   match(qb2$QS, qb1$QS)

    unmatched2 =   which(is.na( m21))

    ####  these are the indeces of questions that do not match
    #### i.e. questions in Q2 that are not in Q1
    
####ques = qb2$QS[ which(is.na( m21)) ]


    return( unmatched2 )

  }
#########   
#####  

