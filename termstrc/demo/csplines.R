oldpar <- par(no.readonly = TRUE)

data(govbonds)

## Cubic splines estimation
cs_res <- estim_cs(govbonds, c("FRANCE"), matrange=c(0,30))
print(cs_res)
summary(cs_res)

plot(cs_res)

## Pricing errors per bond
plot(cs_res,errors="price")

## Remove outliers
cs_res2 <- estim_cs(rm_bond(govbonds, c("FRANCE"),c("FR0000571044", "FR0000571085")), c("FRANCE"), matrange=c(0,30))

summary(cs_res2)

par(oldpar)
