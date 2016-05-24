library(hmmm)

data(relpolbirth)
 
y<-getnames(relpolbirth,st=12,sep=";")

names<-c("Rel","Pol","Birth")
  
# variable 1: Religion
# variable 2: Politics
# variable 3: Birthcontrol

#the lower the variable number is the faster the variable sub-script changes in the vectorized table

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Various hierarchical loglinear models: 
# see Table 2.4, pg. 32,
# "Marginal models for dependent, clustered and longitudinal categorical data",
# Bergsma, W., Croon, M. and Hagenaars, J.A.
# Springer, 2009.  
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# (Pol,Rel) _||_ Birth

f<-~Rel*Pol+Birth
model1<-loglin.model(lev=c(3,7,4),formula=f,names=names)

# alternatively 
# model1<-loglin.model(lev=c(3,7,4),int=list(c(1,2),c(3)),names=names)

mod1<-hmmm.mlfit(y,model1)

print(mod1)

# Pol _||_ Birth|Rel

f<-~Rel*Pol+Rel*Birth
model2<-loglin.model(lev=c(3,7,4),formula=f,names=names)

# alternatively 
# model2<-loglin.model(lev=c(3,7,4),int=list(c(1,2),c(1,3)),names=names)

mod2<-hmmm.mlfit(y,model2)

print(mod2)

# Rel _||_ Birth|Pol

f<-~Rel*Pol+Pol*Birth
model3<-loglin.model(lev=c(3,7,4),formula=f,names=names)

# alternatively 
# model3<-loglin.model(lev=c(3,7,4),int=list(c(1,2),c(2,3)),names=names)

mod3<-hmmm.mlfit(y,model3)

print(mod3)
 
# No 3-factor interaction

f<-~Rel*Pol+Rel*Birth+Pol*Birth
model4<-loglin.model(lev=c(3,7,4),formula=f,names=names)

# alternatively 
# model4<-loglin.model(lev=c(3,7,4),int=list(c(1,2),c(1,3),c(2,3)),names=names)

mod4<-hmmm.mlfit(y,model4)

print(mod4)

# Rel _||_ Pol|Birth

f<-~Rel*Birth+Pol*Birth
model5<-loglin.model(lev=c(3,7,4),formula=f,names=names)

# alternatively 
# model5<-loglin.model(lev=c(3,7,4),int=list(c(1,3),c(2,3)),names=names)

mod5<-hmmm.mlfit(y,model5)

print(mod5)


