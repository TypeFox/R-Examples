data (fire.df)
fire.df
plot(damage~distance,main="Damage versus distance",data=fire.df)
fire.fit<-lm(damage~distance,data=fire.df)
eovcheck(fire.fit)
summary(fire.fit)
plot(fire.fit,which=1)
normcheck(fire.fit)
ciReg(fire.fit)
fire.predict<-data.frame(c(1,4))
predict20x(fire.fit,fire.predict)

