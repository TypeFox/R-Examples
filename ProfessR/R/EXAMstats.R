EXAMstats<-function(j, key)
  {

###   TRON is the output of reading in the scantron from UNC
###  the input is basically a matrix and a key
    dj = dim(j)
    Nquestions = dj[2] 
    Nstudents = dj[1] 
    key = as.numeric(key)

    bscore = rep(0, times=Nstudents)
    for(i in 1:Nstudents)
      {
        L =  length(which(key==j[i,]))
        bscore[i] = L
      }


    scors = bscore 

    mscor = mean(scors)
    sdscor = sd(scors)
    N = length(scors)   ###  number of students
    K = length(scors)

    QN = Nquestions

    jscors = (scors-mscor)/sdscor

    TOPpct = 0.27

    jnum = round(TOPpct*N)

    hscor =  sort(scors)
    Wtop = which(scors>hscor[N-jnum])
    Wbot = which(scors<hscor[jnum])

    blow = matrix(ncol=5, nrow=Nquestions)
    Correct =  rep(0, Nquestions)
    Desc =  rep(0, Nquestions)
    BiSer =  rep(0, Nquestions)
    PPS =  rep(0,Nquestions )
    QPS =  rep(0, Nquestions)

    for(i in 1:Nquestions)
      {
        g1 = j[,i]
        blow[i,1 ] = length(which(g1==1))
        blow[i,2 ] = length(which(g1==2))
        blow[i,3 ] = length(which(g1==3))
        blow[i,4 ] = length(which(g1==4))
        blow[i,5 ] = length(which(g1==5))

        Correct[i] = blow[i, as.numeric(key[i]) ]


        PX = Correct[i] / N

        PPS[i] = PX
        QPS[i] = 1-PX

        wcorrect = which(g1 == key[i])
        PY = sum(scors[ wcorrect ]) / sum(scors)
        SX =  sqrt(PX * (1-PX))
        BiSer[i] = ((PY - PX) / SX) * (mscor / sdscor)

        Ltop = length(which(g1[Wtop] == key[i]))/length(Wtop)
        Lbot = length(which(g1[Wbot] == key[i]))/length(Wbot)

        Desc[i] = Ltop - Lbot
      }

    JX = (N*Correct-Correct^2)/(N^2)

    kr201 = ((QN)/(QN-1))*(1-(sum(JX)/sdscor^2))

    ######   this is the Kuder–Richardson Formula 20 (KR-20)
    kr20 = (QN)*(1-(sum(PPS*QPS)/sdscor^2))/(QN-1)

    difficulty  = Correct/N

###### Desc = Discrimination index = (Upper Group percent Correct) – (Lower Group percent Correct) 

######  BiSer =  Biserial Correlation Coefficient is: ((PY - PX) / SX) * (Mean / Std)

    
    H = cbind(1:dj[2], Correct, blow, difficulty, Desc, BiSer)



    return(list(H=H, kr20=kr20) )

  }
