library("Rfmtool")
env<-fm.Init(3)

#prepare some data

dcol<-3
drow<-10

 mydata<-matrix(runif(dcol*drow,0,1),drow,dcol)

#generate some lambda fuzzy measure
mymeasure<-fm.ConstructLambdaMeasure(c(0.1,0.1,0.9),env)

# calculate choquet integrals of the data

mych<-apply(mydata,1,function(x) fm.Choquet(x, mymeasure$measure, env))

datafit<-cbind(mydata,mych)

#now fit the measure to the data (we should get back the same fuzzy measure we used as a model)

fittedm<-fm.fitting(datafit,env);

mymeasure

fittedm