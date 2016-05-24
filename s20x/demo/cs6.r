data (incomes.df)
incomes.df$incomes
layout20x(1,2)
hist(incomes.df$incomes,main="Histogram of incomes",xlab="incomes")
boxplot(incomes.df$incomes,main="Boxplot of incomes")
summaryStats(incomes.df$incomes)
normcheck(incomes.df$incomes)
t.test(incomes.df$incomes)

