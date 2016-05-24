data (house.df)
house.df$advertised.price
house.df$sell.price
house.df<-within(house.df,{price.difference<-advertised.price-sell.price})
house.df$price.difference
layout20x(1,2)
hist(house.df$price.difference,breaks=10)
boxplot(house.df$price.difference,main="Boxplot of price.difference")
summaryStats(house.df$price.difference)
normcheck(house.df$price.difference)
t.test(house.df$price.difference,alt="greater")

