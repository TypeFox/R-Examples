## Create an example data matrix with 50 observations that contains an ID variable for every unit, a dummy variable indicating gender, an age variable (between 18 and 55), a treatment variable and an outcome variable (between 15 and 20). 
set.seed(12345)
id <- seq(1,50,1)
gender <- sample(c(1,2), 50, replace=T)
age <- sample(seq(18,55,1), 50, replace=T)
treat <- sample(c(1,2), 50, replace=T)
out <- treat + sample(seq(15, 20, 1), 50, replace = TRUE)
data <- cbind(id, gender, age, out, treat)

## Check summary statistics for the created data
aggregate(out~treat, data, mean)

## Run invertRIconfInt()
invertRIconfInt(data, outcome.var="out", tr.var="treat", tau.abs.min = -3, tau.abs.max = 3, n.sb.p = 20, id.vars = "id", id.vals = "id", exact.vars = c("gender", "age"), exact.vals = c("gender", "age"))
print("This demo uses 'n.sb.p = 20' to minimize the computational effort. When running invertRIconfInt() outside of demonstrational purposes, 'n.sb.p' should be set to a much higher value.")