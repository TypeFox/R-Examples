library(parfm)
data(kidney)
kidney$sex <- kidney$sex - 1

plot(select.parfm(Surv(time,status) ~ sex + age, 
                   dist=c("exponential", "gompertz", "lognormal"),
                   frailty=c("gamma", "possta"),
                   cluster="id", data=kidney))
