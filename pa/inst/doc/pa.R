### R code from vignette source 'pa.Rnw'

###################################################
### code chunk number 1: pa.Rnw:116-118
###################################################
options(width = 50, digits = 2, scipen = 5)
library(pa)


###################################################
### code chunk number 2: pa.Rnw:121-123
###################################################
data(year)
names(year)


###################################################
### code chunk number 3: pa.Rnw:170-176
###################################################
sample.mat <- rbind(which(year$portfolio == 0 & year$benchmark == 0)[2],
                    which(year$portfolio == 0 & year$benchmark > 0)[2],
                    which(year$portfolio > 0 & year$benchmark == 0)[2],
                    which(year$return > 0 & year$portfolio > 0 & year$benchmark > 0)[500],
                    which(year$barrid == "CANAITC")[11])
year[sample.mat, c(-1, -6, -7, -9: -11, -13) ]


###################################################
### code chunk number 4: pa.Rnw:451-457
###################################################
data(jan)
br.single <- brinson(x = jan, date.var = "date", 
                 cat.var = "sector",
                 bench.weight = "benchmark", 
                 portfolio.weight = "portfolio", 
                 ret.var = "return")


###################################################
### code chunk number 5: pa.Rnw:474-475
###################################################
summary(br.single)


###################################################
### code chunk number 6: pa.Rnw:486-487
###################################################
plot(br.single, var = "sector", type = "return")


###################################################
### code chunk number 7: pa.Rnw:617-623
###################################################
data(quarter)
br.multi <- brinson(quarter, date.var = "date",
                    cat.var = "sector",
                    bench.weight = "benchmark",
                    portfolio.weight = "portfolio",
                    ret.var = "return")


###################################################
### code chunk number 8: pa.Rnw:630-631
###################################################
exposure(br.multi, var = "size")


###################################################
### code chunk number 9: pa.Rnw:641-642
###################################################
returns(br.multi, type = "linking")


###################################################
### code chunk number 10: pa.Rnw:659-660
###################################################
plot(br.multi, type = "return")


###################################################
### code chunk number 11: pa.Rnw:749-759
###################################################
test.df <- data.frame(Return = c(0.3, 0.4, 0.5),
                      Name = c('A', 'B', 'C'),
                      Size = c(1.2, 2, 0.8), 
                      Value = c(3, 2, 1.5),
                      Active_Weight= c(0.5, 0.1, -0.6))
test.df
## model <- lm(Ret ~ Size + Value, data = test.df)
## test.df[,1] %*% test.df[,5]
## model$coefficients[2] * t(test.df[,3]) %*% test.df[,5] ## active exposure of size 
## model$coefficients[3] * t(test.df[,4]) %*% test.df[,5] ## active exposure of value 


###################################################
### code chunk number 12: pa.Rnw:776-783
###################################################
rb.single <- regress(jan, date.var = "date",
                ret.var = "return",
                reg.var = c("sector", "growth",
                  "size"),
                benchmark.weight = "benchmark",
                portfolio.weight = "portfolio")
exposure(rb.single, var = "growth")


###################################################
### code chunk number 13: pa.Rnw:795-796
###################################################
summary(rb.single)


###################################################
### code chunk number 14: pa.Rnw:838-845
###################################################
rb.multi <- regress(quarter, date.var = "date",
                ret.var = "return",
                reg.var = c("sector", "growth", 
                  "size"),
                benchmark.weight = "benchmark",
                portfolio.weight = "portfolio")
rb.multi


###################################################
### code chunk number 15: pa.Rnw:855-856
###################################################
summary(rb.multi)


###################################################
### code chunk number 16: pa.Rnw:870-877
###################################################
rb.multi2 <- regress(year, date.var = "date",
               ret.var = "return",
               reg.var = c("sector", "growth", 
                 "size"),
               benchmark.weight = "benchmark",
               portfolio.weight = "portfolio")
returns(rb.multi2, type = "linking")


###################################################
### code chunk number 17: pa.Rnw:889-890
###################################################
plot(rb.multi2, var = "sector", type = "return")


###################################################
### code chunk number 18: pa.Rnw:980-987
###################################################
data(test)
test.br <- brinson(x = test, date.var = "date", 
                   cat.var = "sector",
                   bench.weight = "benchmark", 
                   portfolio.weight = "portfolio", 
                   ret.var = "return")
returns(test.br)


###################################################
### code chunk number 19: pa.Rnw:998-1005
###################################################
test.reg <- regress(x =test, 
                    date.var = "date", 
                    ret.var = "return", 
                    reg.var = "sector",
                    benchmark.weight = "benchmark", 
                    portfolio.weight = "portfolio") 
returns(test.reg)


###################################################
### code chunk number 20: pa.Rnw:1028-1033
###################################################
lm.test <- lm(return ~ sector - 1, 
              data = test[test$portfolio != 0, ])
lm.test$coefficients 
exposure(br.single, var = "sector")[ ,2] %*% 
  (lm.test$coefficients - test.reg@coefficients)


