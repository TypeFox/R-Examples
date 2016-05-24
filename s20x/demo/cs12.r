data (incomes.df)
boxplot(incomes.df$incomes,main="Boxplot of Incomes")
summaryStats(incomes.df$incomes)
normcheck(incomes.df$incomes)
log.incomes<-log(incomes.df$incomes)
boxplot(log.incomes,main="Boxplot of (log) incomes")
summaryStats(log.incomes)
normcheck(log.incomes)
t.test(log.incomes)
c.i.<-t.test(log.incomes)$conf.int
c.i.
exp(c.i.)

