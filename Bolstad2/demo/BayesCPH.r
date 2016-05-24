AidsSurvival.df<-data(AidsSurvival)
y<-AidsSurvival.df$censor
t<-AidsSurvival.df$time
X<-cbind(AidsSurvival.df$age,AidsSurvival.df$drug)
BayesCPH(y,t,X,plots=TRUE)

