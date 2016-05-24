data (rain.df)
rain.df
layout20x(1,2)
stripchart(rain~seed,pch=1,vert=T,ylab="rain",main="Stripcharts of Seeding",data=rain.df)
boxplot(rain~seed,ylab="rain",main="Boxplots of Seeding",data=rain.df)
summaryStats(rain~seed,data=rain.df)
rain.fit<-lm(rain~seed,data=rain.df)
eovcheck(rain.fit)
rain.df<-within(rain.df,{log.rain<-log(rain)})
rain.df[1,]
boxplot(log.rain~seed,ylab="(log) rain",main="Boxplots of (log) Seeding",data=rain.df)
summaryStats(log.rain~seed,data=rain.df)
rain.fit1<-lm(log.rain~seed,data=rain.df)
eovcheck(rain.fit1)
normcheck(rain.fit1)
t.test(log.rain~seed,var.equal=T,data=rain.df)
c.i.<-t.test(log.rain~seed,var.equal=T, data=rain.df)$conf.int
c.i.
exp(c.i.)

