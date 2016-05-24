# Test unit for simulate_dataframe
library(dfexplore)
test <- simulate_dataframe()
test2 <- simulate_dataframe( includeMatrix=T)
str(test)
str(test2)