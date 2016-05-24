`mtapspec` <-
function(a, dt, klen=length(a), MTP=NULL)
{
#####  multi-taper spectrum analysis
  ####     Mspec =   mtapspec(a$y,a$dt, klen=4096, MTP=list(kind=2,nwin=5, npi=3,inorm=0)  )   
#####

  if(missing(MTP))
    {
      kind=2;  
      nwin=5;  
      npi=3; 
      inorm=1;
    }
  else
    {
      kind=MTP$kind;  
      nwin=MTP$nwin;  
      npi=MTP$npi; 
      inorm=MTP$inorm ;
    }



  len = length(a)

  if(len<2)
    {
      return(0)

    }


  
  if(missing(klen))
    {
      klen=2*next2(len)
    }
    if(klen<len)
      {
        klen = 2*next2(len)
      }
   



  
  numfreqs = 1+klen/2;
  numfreqtap = numfreqs*nwin;
  nyquist = 0.5/dt;
   df = 2*nyquist/klen;
  freq = df*seq(0,numfreqs-1)
  
  spec1 = rep(0, length=klen )
  dof = rep(0, length=klen )
  Fvalues = rep(0, length=klen )
  ReSpec= rep(0, length= numfreqtap)
  ImSpec= rep(0, length=numfreqtap )

  barf = .C("CALL_Mspec",PACKAGE = "RSEIS",
    as.double(a), 
    as.integer(len),
    as.integer(kind),
    as.integer(nwin) ,
    as.double(npi) ,
    as.integer(inorm) ,
    as.double(dt) ,
    as.double(spec1) ,
    as.double(dof) ,
    as.double(Fvalues) ,
    as.integer(klen) ,
    as.double(ReSpec) ,
    as.double(ImSpec) )

  Ispec=  matrix(unlist(barf[13]), byrow=FALSE, nrow=numfreqs,  ncol=nwin)

  
  Rspec=   matrix(unlist(barf[12]), byrow=FALSE, nrow=numfreqs,  ncol=nwin)

  
  invisible(list(dat=a, dt=dt, spec=unlist(barf[8]), dof=unlist(barf[9]),Fv=unlist(barf[10]),Rspec=Rspec, Ispec=Ispec, freq=freq, df=df, numfreqs=numfreqs, klen=klen,   mtm=list(kind=kind, nwin=nwin, npi=npi, inorm=inorm)))


}

