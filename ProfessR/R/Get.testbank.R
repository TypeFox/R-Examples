`Get.testbank` <-
  function(fn)
  {

##########  this takes an ASCII text file of questions
#########   and returns a list of questions lists.
#########  blank lines are stripped.

    
#########  the key word QUESTION: must appear at the beginning of the line
#########  same for the the key word ANSWER:

    ##  there should be a blank after the keys words but   the program tries
    ##  to compensate for no blanks.

###   there should be no blank lines or lines with NO alphanumeric information:
#######   the program tries to get rid of these also
    
#####     in the file comments are signified by a hash mark # in column 1
    

    ALLQ = scan(file=fn, what="", sep="\n",  quiet=TRUE)

    

    comm = grep("^#", ALLQ)
    
    if(length(comm)>0)
      {
        ALLQ = ALLQ[-comm]

      }

    
    q1 = grep("^QUESTION:", ALLQ)

    
    Qbank = vector(mode="list")

    for(i in 1:length(q1))
      {
        i1 = q1[i]
	if(i<length(q1))
          { i2 = q1[i+1]-1 }	else { i2 = length(ALLQ)	}

        
        ## print(paste(sep=" ", "#####", i, i1, i2))
        ## print(ALLQ[i1:i2])

        
        az1 = substring(ALLQ[i1], 10, 10)

        if(identical(az1, " "))
          {
            quest = substring(ALLQ[i1], 11, nchar(ALLQ[i1]))

          }else{
            
            quest = substring(ALLQ[i1], 10, nchar(ALLQ[i1]))
          }




####################   try to eliminate any line that has no content


        ans1 = ALLQ[(i1+1):i2]

        
        ww=grep("[a-z,A-Z]", ans1)
        ans1 = ans1[ww]
        
        grfig = grep("FIG:", ans1)
        if(length(grfig)>=1)
          {
            fig1 = ans1[grfig]
            ufig = unlist(strsplit(fig1, split=" "))
            
            figname = ufig[2]
            figtag = ufig[3]
            fignewt = list(fn=figname, tag =figtag) 
            ans1 = ans1[-grfig]
            
          }
        else
          {
            fignewt = NULL
          }
        
        a1 = grep("ANSWER:", ans1)
        fa1 = ans1[a1]

        az1 = substring(fa1, 8, 8)

        if(identical(az1, " "))
          {
            a2 = substring(fa1, 9, nchar(fa1))

          }
        else
          {
            
            a2 = substring(fa1, 8, nchar(fa1))
          }




        
###   a2 = substring(fa1, 9, nchar(fa1))
        ans1[a1] = a2

        Qbank[[i]] = list(Q=quest, A=ans1, a=fa1, numANS=a1, FIG=fignewt)

      }

    attr(Qbank, "fn")<-fn

    return(Qbank)



  }

