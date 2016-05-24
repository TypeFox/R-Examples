"ICC1.CI" <-
function (dv, iv, data, level=.95)
 {
require(multilevel)
 attach(data)
 mod <- aov(dv ~ as.factor(iv), na.action=na.omit)
 detach(data)
 icc <- ICC1(mod)
 tmod <- summary(mod)
 df1 <- tmod[[1]][1,1]
 df2 <- tmod[[1]][2,1]
 Fobs <- tmod[[1]][1,4]
 n <- df2/(df1+1) # k-1
 noma <- 1- level
 Ftabl <- qf(noma/2, df1, df2, lower.tail=F)
 Ftabu <- qf(noma/2, df2, df1, lower.tail=F)
 Fl <- Fobs/Ftabl
 Fu <- Fobs*Ftabu
 lcl <- (Fl-1)/(Fl+n)
 ucl <- (Fu-1)/(Fu+n)
 mat <- data.frame(LCL=lcl, ICC1=icc, UCL=ucl)
 return(mat)
 }

