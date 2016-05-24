
readSCANTRON<-function(fn="t9543b.raw.csv", nq=50, istart=6)
  {
    ##  B3 = "t9543b.raw.csv"
        PRINT = FALSE
    rb3 = read.csv(fn, header=TRUE, sep=",", quote="\"",  stringsAsFactors =FALSE  )

    

                                        #  first row is the header
    ##  second row is the key information

###   if students mark beyond the end of the exam,
    ##     these rows are read in and recorded.
###      need to get rid of these

    Nams = as.character(rb3$Name)
    Nrows = length(Nams)
    Nams =  Nams[2:Nrows]
    Nstudents =  length(Nams)

    ids = as.character(rb3[2:Nrows,2])

###  need to fix the NA values in the names and IDS


    getcol = istart:(istart+nq-1)

    j = rb3[2:Nrows,6:(6+nq-1)]
    rownames(j)<-NULL
    colnames(j)<-NULL


    
    dj = dim(j)

    K = dj[1]

    studans = matrix(as.numeric(as.matrix(j)), ncol=dj[2], nrow=dj[1])

    if(PRINT)
      {
        for(m in 1:Nstudents)
          {
            cat(paste(Nams[m], ids[m]) , sep=" ")
             cat(" ")
            cat(studans[m,], sep=" ")
            cat("\n")
          }
        
      }  

    key = as.numeric( rb3[1, 6:(6+nq-1)] )

    return(list(Nstudents=Nstudents, Nquestions=nq   ,Nams=Nams, ids=ids, studans=studans, key=key) )
  }

### rawtsts = list.files(path=".", pattern="raw.csv")

###  do1 = readSCANTRON(fn="t9543b.raw.csv") 
###  do1 = readSCANTRON(fn=rawtsts[1], nq=49) 
###  do2 = readSCANTRON(fn=rawtsts[2], nq=49) 
