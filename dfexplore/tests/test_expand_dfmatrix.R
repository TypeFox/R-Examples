# Test unit for dfmatrixToStandard
library(dfexplore)
test2 <- simulate_dataframe( includeMatrix=T)
test3 <- expand_dfmatrix(test2)

# if there is no matrix
aa<-expand_dfmatrix(cars)
all(aa == cars)

# give the position
expand_dfmatrix(test2, matrix_var=9)
