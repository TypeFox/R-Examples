# exemple HM, 24/02/2015

library(planor)
toto <- planor.designkey(factors = c("A","B","C","D","E","F","G","H","J"), nlevels = c(2,2,2,2,2,4,4,4,2), nunits=16, model = ~A+B+C+D+E+F+G+H+J, estimate=~A+B+C+D+E+F+G+H+J)
print(toto)
summary(toto)
