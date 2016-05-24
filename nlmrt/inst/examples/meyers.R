# Meyers Problem from 
# METHODS FOR NON-LINEAR LEAST SQUARES PROBLEMS, 2nd Edition, April 2004
# K. Madsen, H.B. Nielsen, O. Tingleff


formula <- "y ~ b1 * exp(b2/(t+b3))"
isq<-seq(1,16)
tdat<-45+5*isq
udat<-tdat/100
ydat<-c(34780, 28610, 23650, 19630, 16370, 13720, 
      11540, 9744, 8261, 7030, 6005, 5147, 4427, 3820, 3307, 2872)
form2<-"y/1000 ~ z1*exp((10*z2/(u+z3))-13)"
mdata<-data.frame(y=ydat, t=tdat, u=udat)

st1<-c(b1=1, b2=1, b3=1)
stsugg<-c(b1=0.02, b2=4000, b3=250)
require(nlmrt)
# ans<-nlsmnqb(formula, start=st1, trace=TRUE, data=mdata, control=list(watch=TRUE))
amrt1<-nlsmnqb(formula, start=st1, trace=TRUE, data=mdata)
print(amrt1)
anls1<-try(nls(formula, start=st1, trace=TRUE, data=mdata))
print(anls1)

amrt2<-nlsmnqb(formula, start=stsugg, trace=TRUE, data=mdata)
print(amrt2)
anls2<-try(nls(formula, start=stsugg, trace=TRUE, data=mdata))
print(anls2) 

z0<-c(z1=8.85, z2=4, z3=2.5)

amrtz1<-nlsmnqb(form2, start=z0, trace=TRUE, data=mdata)
print(amrtz1)
anlsz1<-try(nls(form2, start=z0, trace=TRUE, data=mdata))
print(anlsz1) 

