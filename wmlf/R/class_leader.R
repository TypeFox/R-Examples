setClass(
  Class="Leader",
  representation=representation(
    dwtn="list",
    leader="list",
    J="numeric",
    frac="numeric",
    full="logical",
    energy="logical",
    df.mf="data.frame"
    )
)

setMethod(
  f="[",
  signature="Leader",
  definition=function (x, i, j, ..., drop)
  {
    switch(EXPR=i,
           "dwtn"={return(x@dwtn)},
           "leader"={return(x@leaderj)},
           "J"={return(x@J)},
           "frac"={return(x@frac)},
           "full"={return(x@full)},
           "energy"={return(x@energy)},
           "df.mf"={return(x@df.mf)}
           )
  }
)


setMethod("$", "Leader",function(x, name){slot(x,name)})
setGeneric(name="plotfit",def=function(object){standardGeneric("plotfit")})

setMethod(
  f="plotfit",
  signature="Leader",
  definition=function(object) 
  {
    op <- par(no.readonly = TRUE)
    par(mfrow=c(1,1),ask=TRUE)
    df=object@df.mf
    plot(df$j,df$hminj,xlab="j",ylab="Maximum DWT")
    if (object@full==TRUE)
    {
      plot(df$j,df$c1j,xlab="j",ylab="Log Cumulant of order 1")
      plot(df$j,df$c2j,xlab="j",ylab="Log Cumulant of order 2")
    }
    if (object@energy==TRUE)
    {
      plot(df$j,df$t1j,xlab="j",ylab="Log sum of wavelet coefficient")
      plot(df$j,df$t2j,xlab="j",ylab="Log sum of wavelet coefficient squared")
    }
     par(op)
  }
)

setGeneric(
  name="fit",
  def=function(object,scales){
    standardGeneric("fit")
  }
)




setMethod(
  f="fit",
  signature="Leader",
  definition=function(object,scales){
    
    df=object@df.mf[scales,]
    hminj=df$hminj+log2(2^(-scales*object@frac))
    hmin=lm(hminj~scales)
    fit=data.frame(c(hmin$coefficients[2],summary(hmin)$r.squared))
    colnames(fit)="hmin"
    rownames(fit)=c("est","R2")
    if (object@full==T)
    {
      c1=lm(df$c1j~scales)
      c2=lm(df$c2j~scales)
      fit2=data.frame(c(c1$coefficients[2]-object@frac,summary(c1)$r.squared),c(c2$coefficients[2],summary(c2)$r.squared))
      colnames(fit2)=c("c1","c2")
      rownames(fit2)=c("est","R2")
      fit=cbind(fit,fit2)
    }
    if (object@energy==T)
    {
      t1j=df$t1j+log2(2^(-scales*object@frac))
      t2j=df$t2j+log2(2^(-2*scales*object@frac))
      t1=lm(t1j~scales)
      t2=lm(t2j~scales)
      fit3=data.frame(c(t1$coefficients[2],summary(t1)$r.squared),c(t2$coefficients[2],summary(t2)$r.squared))
      colnames(fit3)=c("t1","t2")
      rownames(fit3)=c("est","R2")
      fit=cbind(fit,fit3)
    }
    return(fit)
  }
)


leader<-function(image,wf="d6",frac=0,full=T,energy=F){
  maxj=list()
  J=log2(nrow(image))
  n=switch(EXPR=wf,"d4"={2},"d6"={4},"d8"={6})#pretraitement effet de bord
  Jadmi=sum((2^(2*(J-(1:J)))-(n*2^(1+(J-(1:J))))+n^2)>=25)#pour tenir comptes des effets de bords
  #calcul ondelette
  dwt <- dwt.2d(image,wf,J) 
  dwtn=list()
  for (j in 1:Jadmi){dwtn[[j]]<-2^((-1+frac)*j)*pmax(abs(dwt[paste('HL',j,sep="")][[1]]),abs(dwt[paste('LH',j,sep="")][[1]]),abs(dwt[paste('HH',j,sep="")][[1]]))}
  #calcul leader
  leader=list()
  if (full==TRUE)
  {
    for (j in 1:Jadmi)
    {  
      lead=list()
      for (jj in 1:j)
      {
        pas=2^(j-jj)  
        inter=dwtn[[jj]]
        dim=2^(J-j)
        res= matrix(0,dim,dim)
        for (ix in 1:dim)
        {
          for (iy in 1:dim)
          {
            xmin=max(1,(ix-1)*pas+1-pas)
            xmax=min(dim*pas,(ix-1)*pas+2*pas)
            ymin=max(1,(iy-1)*pas+1-pas)
            ymax=min(dim*pas,(iy-1)*pas+2*pas)
            res[ix,iy]=max(inter[xmin:xmax,ymin:ymax])
          }
        }
        lead[[jj]]=res
      }
      leader[[j]]=do.call(pmax,lead)
    }
  }
  

#On enleve les bords gauche et haut des coeffs et des leaders 
bord=switch(EXPR=wf,"d4"={c(1,2,2)},"d6"={c(2,3,4)},"d8"={c(3,5,6)})
for (j in 1:Jadmi)
{
  tempd=dwtn[[j]]
  if(full==T){templ=leader[[j]]}
  if (j==1){dwtn[[j]]=tempd[-(1:bord[1]),-(1:bord[1])]
            if(full==T){leader[[j]]=templ[-(1:bord[1]),-(1:bord[1])]}}
  if (j==2){dwtn[[j]]=tempd[-(1:bord[2]),-(1:bord[2])]
            if(full==T){leader[[j]]=templ[-(1:bord[2]),-(1:bord[2])]}}
  if (j>=3){dwtn[[j]]=tempd[-(1:bord[3]),-(1:bord[3])]
            if(full==T){leader[[j]]=templ[-(1:bord[3]),-(1:bord[3])]}}
}
#calcul de df.mf
  maxwj<-function(dwtn){unlist(lapply(dwtn,function(x) {log2(max(abs(x)))}))}
  c1<-function(leader){unlist(lapply(leader,function(x) {sum(log2(abs(x[abs(x)>0])))/sum(abs(x)>0)}))}
  c2<-function(leader){unlist(lapply(leader,function(x) {sum((log2(abs(x[abs(x)>0])))^2)/sum(abs(x)>0)}))-c1(leader)^2}
  t1<-function(dwtn){unlist(lapply(dwtn,function(x) {log2(sum(abs(x[abs(x)>0]))/sum(abs(x)>0)  )     }))} 
  t2<-function(dwtn){unlist(lapply(dwtn,function(x) {log2(sum((abs(x[abs(x)>0]))^2)/sum(abs(x)>0)  )     }))} 
  df.mf=data.frame(j=1:Jadmi,hminj=maxwj(dwtn))
  if (full==TRUE)
  {
    df.mf$c1j=c1(leader)
    df.mf$c2j=c2(leader)  
   }
  if (energy==TRUE)
  {
    df.mf$t1j=t1(dwtn)
    df.mf$t2j=t2(dwtn)  
  }
  return(new(Class="Leader",dwtn=dwtn,leader=leader,J=Jadmi,frac=frac,full=full,energy=energy,df.mf=df.mf))
}




