x <- y <- c(7, 9, NA, NA, 13)
all(is.na(x) == is.na(y)) & all((x == y)[!is.na(x)]) 
x <- factor(c("Huey", "Dewey", "Louie", "Huey"))
y <- c("blue", "red", "green")
x
y[x]
juul.girl <- juul[juul$age >=7 & juul$age < 14 & juul$sex == 2,]
summary(juul.girl)
sapply(1:10, function(i) mean(rexp(20)))
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
