data (nail.df)
boxplot(dry~polish,main="Drying times by polish",xlab="polish",ylab="time in seconds",data=nail.df)
summaryStats(dry~polish,data=nail.df)
nail.fit<-lm(dry~polish,data=nail.df)
eovcheck(nail.fit)
normcheck(nail.fit)
t.test(dry~polish,var.equal=F,data=nail.df)

