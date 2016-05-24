data (fruitfly.df)
boxplot(fecundity~strain,main="Boxplots of fecundity by strain",ylab="fecundity",xlab="strain",data=fruitfly.df)
summaryStats(fecundity~strain,data=fruitfly.df)
fruitfly.fit<-lm(fecundity~strain,data=fruitfly.df)
eovcheck(fruitfly.fit)
normcheck(fruitfly.fit)
summary1way(fruitfly.fit)
multipleComp(fruitfly.fit)

