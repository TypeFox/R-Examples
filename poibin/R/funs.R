ppoibin <-function(kk,pp,method="DFT-CF",wts=NULL)
{

  if(any(pp<0)|any(pp>1))
  {
    stop("invalid values in pp.")
  }

  if(is.null(wts))
  {
    wts=rep(1,length(pp))
  }
  
  switch(method,
  "DFT-CF"={
            mm=length(kk)
            res=double(mm)
            npp=length(pp)
            n=sum(wts)
            avec=double(n+1)
            bvec=double(n+1)
            funcate=1
            ex=0
            tmp=.C("multi_bin_dft_cf",as.double(res),as.integer(kk),
            as.integer(mm),as.integer(n),as.double(pp),as.double(avec),as.double(bvec),
            as.integer(funcate),as.double(ex),as.integer(npp),as.integer(wts),PACKAGE="poibin")
            res=tmp[[1]]

            res[res<0]=0
            res[res>1]=1
            res[kk<0]=0
            res[kk>=sum(wts)]=1
            
           },
  "RF"=    {
              kk1=kk
              kk[kk<0]=0
              pp=rep(pp,wts)
              mm=length(kk)
              res=double(mm)
              n=length(pp)
              mat=rep(0.0,(n+1)*(n+2))
              tmp=.C("multi_bin_bh",as.double(res),as.integer(kk),as.integer(mm),as.integer(n),as.double(pp),as.double(mat),PACKAGE="poibin")
              res=tmp[[1]]
              res[kk1<0]=0
              res[kk1>=sum(wts)]=1
           },
  "RNA"=   {
               pp=rep(pp,wts)
               muk=sum(pp)
               sigmak=sqrt(sum(pp*(1-pp)))
               gammak=sum(pp*(1-pp)*(1-2*pp))
               ind=gammak/(6*sigmak^3)
               kk1=(kk+.5-muk)/sigmak
               vkk.r=pnorm(kk1)+gammak/(6*sigmak^3)*(1-kk1^2)*dnorm(kk1)
               vkk.r[vkk.r<0]=0
               vkk.r[vkk.r>1]=1
               res=vkk.r
           },
  "NA"=    {
              pp=rep(pp,wts)
              muk=sum(pp)
              sigmak=sqrt(sum(pp*(1-pp)))
              gammak=sum(pp*(1-pp)*(1-2*pp))
              kk1=(kk+.5-muk)/sigmak
              res=pnorm(kk1)
           }
  )

  return(res)
}

dpoibin <-function(kk,pp,wts=NULL)
{
 if(any(pp<0)|any(pp>1))
 {
    stop("invalid values in pp.")
 }
  
 if(is.null(wts))
 {
    wts=rep(1,length(pp))
 }
 mm=length(kk)
 res=double(mm)
 npp=length(pp)
 n=sum(wts)
 avec=double(n+1)
 bvec=double(n+1)
 funcate=2
 ex=0
 tmp=.C("multi_bin_dft_cf",as.double(res),as.integer(kk),
 as.integer(mm),as.integer(n),as.double(pp),as.double(avec),as.double(bvec),
 as.integer(funcate),as.double(ex),as.integer(npp),as.integer(wts),PACKAGE="poibin")
 res=tmp[[1]]

 res[res<0]=0
 #res[kk<0|kk>=sum(wts)]=0
 res[kk<0|kk>sum(wts)]=0
 #modified on 13 Feb 2013, reported bug by user.
 return(res)
}

qpoibin <-function(qq,pp,wts=NULL)
{

 if(any(pp<0)|any(pp>1))
 {
    stop("invalid values in pp.")
 }

 if(any(qq<0)|any(qq>1))
 {
    stop("invalid values in qq.")
 }

 if(is.null(wts))
 {
    wts=rep(1,length(pp))
 }
 nn=1
 mm=length(qq)
 res=double(mm)
 npp=length(pp)
 n=sum(wts)
 avec=double(n+1)
 bvec=double(n+1)
 funcate=3
 ex=qq
 tmp=.C("multi_bin_dft_cf",as.double(res),as.integer(nn),
 as.integer(mm),as.integer(n),as.double(pp),as.double(avec),as.double(bvec),
 as.integer(funcate),as.double(ex),as.integer(npp),as.integer(wts),PACKAGE="poibin")
 res=tmp[[1]]
 return(res)
}

rpoibin <-function(m,pp,wts=NULL)
{

 if(any(pp<0)|any(pp>1))
 {
    stop("invalid values in pp.")
 }

 if(is.null(wts))
 {
    wts=rep(1,length(pp))
 }
 qq=runif(m)
 res=qpoibin(qq=qq,pp,wts=wts)
 return(res)
}

