eval.stage.gauss<-function(dendat,M,mugrid,siggrid=1,sigeka=TRUE,src="c",
sampstart=FALSE,boost=FALSE,N=60)
{
n<-length(dendat)

  if (src=="R"){
     resu<-stage.gaussR(dendat,M,mugrid,siggrid,sigeka,sampstart)
     return(resu)
  }
  else{
            if (boost) inboost<-1 else inboost<-0
            mu0<-mean(dendat)
            sig0<-sqrt(var(dendat))

            indendat<-c(0,dendat)
            inM<-M
            inmugrid<-c(0,mugrid)
            insiggrid<-c(0,siggrid)
            insigeka<-1
            if (sampstart) insampstart<-1 else insampstart<-0
            dictCard<-length(mugrid)
            dictCardSig<-length(siggrid)
            kg<-.C("stageGauss",
               as.double(mu0),
               as.double(sig0),
               as.double(indendat),
               as.integer(inM),
               as.double(inmugrid),
               as.double(insiggrid),
               as.integer(insigeka),
               as.integer(insampstart),
               as.integer(n),
               as.integer(dictCard),
               as.integer(dictCardSig),
               as.integer(inboost),
               muut = double(inM+1),
               sigit = double(inM+1),
               curmix = double(inM+1))
            sgmuut<-kg$muut[2:(inM+1)]
            sgsigit<-kg$sigit[2:(inM+1)]
            sgcurmix<-kg$curmix[2:(inM+1)]

            minu<-min(sgmuut)-2*max(sgsigit)
            maxi<-max(sgmuut)+2*max(sgsigit)
            support<-c(minu,maxi)
            pcf<-pcf.func("mixt",N,sig=sgsigit,M=sgmuut,p=sgcurmix,
            support=support)

            return(list(value=pcf$value,down=pcf$down,high=pcf$high,
                        support=pcf$support,N=pcf$N,
                        muut=sgmuut,sigit=sgsigit,curmix=sgcurmix))
  }

}

  

