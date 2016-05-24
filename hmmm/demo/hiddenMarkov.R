library(hmmm)

data(drinks)

y<-cbind(drinks$lemon.tea,drinks$orange.juice,drinks$apple.juice)
fm<-c("l-l-l-l")
fmargobs<-marg.list(fm,mflag="m")


# saturated model (fsat<-~lat*tea*ojuice*ajuice is implicit)
# an example to calculate the starting values for the probabilities 

model.obsf0<-hmmm.model(marg=fmargobs,
lev=c(2,3,3,3),names=c("lat","tea","ojuice","ajuice"))

modelsat<-hidden.emfit(y,model.obsf0,y.eps=0.01,maxit=5,
maxiter=2500,norm.diff.conv=0.001,printflag=10)
print(modelsat)

#starting values used in the next models
Ptr<-modelsat$Ptr
Ptobs<-modelsat$Ptobs

## model of constant association among tea, orange and apple juice sales given the latent states

fca<-~lat*tea*ojuice*ajuice-lat:tea:ojuice:ajuice-tea:ojuice:ajuice

model.obsfca<-hmmm.model(marg=fmargobs,
lev=c(2,3,3,3),names=c("lat","tea","ojuice","ajuice"),formula=fca)

modelca<-hidden.emfit(y,model.obsfca,y.eps=0.01,maxit=3,maxiter=2500,printflag=10,
old.tran.p=Ptr,bb=Ptobs)

print(modelca,printflag=TRUE)

## model of independence of tea sales from orange and apple juice sales given the latent states

find<-~lat*tea+lat*ojuice*ajuice 
  
model.obsf<-hmmm.model(marg=fmargobs,
lev=c(2,3,3,3),names=c("lat","tea","ojuice","ajuice"),formula=find)

modelind<-hidden.emfit(y,model.obsf,y.eps=0.01,maxit=5,maxiter=2500,printflag=10,
old.tran.p=Ptr,bb=Ptobs)

print(modelind)

## model of total independence of tea, orange and apple juice sales given the latent states

findtot<-~lat*tea+lat*ojuice+lat*ajuice   

model.obsftot<-hmmm.model(marg=fmargobs,
lev=c(2,3,3,3),names=c("lat","tea","ojuice","ajuice"),formula=findtot)

modelindtot<-hidden.emfit(y,model.obsftot,y.eps=0.01,maxit=5,maxiter=2500,printflag=10,
old.tran.p=Ptr,bb=Ptobs)

print(modelindtot)




