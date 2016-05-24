data (computer.df)
computer.df
computer.df<-within(computer.df,{selfassess1.<-factor(selfassess)})
stripchart(score~selfassess1.,pch=1,vert=T,main="Stripcharts of Selfassess",ylab="score", data=computer.df)
summaryStats(score~selfassess1.,data=computer.df)
computer.fit<-lm(score~selfassess1.,data=computer.df)
eovcheck(computer.fit)
normcheck(computer.fit)
summary1way(computer.fit)
multipleComp(computer.fit)

