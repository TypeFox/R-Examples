"distal.med" <-
function (data)
{
yx.lm <- lm (y ~ x, data)# t
yall.lm <- lm (y ~ m2 + m1 + x, data) # c f tp
m2.lm <- lm (m2 ~ m1 + x, data) # b e
m1.lm <- lm (m1 ~x, data)# a

# individual effects
a <- summary(m1.lm)$coef[2,1]
se.a <- summary(m1.lm)$coef[2,2]
t.a <- summary(m1.lm)$coef[2,3]

b <- summary(m2.lm)$coef[2,1]
se.b <- summary(m2.lm)$coef[2,2]
t.b <- summary(m2.lm)$coef[2,3]

c <- summary(yall.lm)$coef[2,1]
se.c <- summary(yall.lm)$coef[2,2]
t.c <- summary(yall.lm)$coef[2,3]

e <- summary(m2.lm)$coef[3,1]
se.e <- summary(m2.lm)$coef[3,2]
t.e <- summary(m2.lm)$coef[3,3]

f <- summary(yall.lm)$coef[3,1]
se.f <- summary(yall.lm)$coef[3,2]
t.f <- summary(yall.lm)$coef[3,3]

t <- summary(yx.lm)$coef[2,1]
se.t <- summary(yx.lm)$coef[2,2]
t.t <- summary(yx.lm)$coef[2,3]

tpr <- summary(yall.lm)$coef[4,1]
se.tpr <- summary(yall.lm)$coef[4,2]
t.tpr <- summary(yall.lm)$coef[4,3]

# covariances of regression coefficients
sbe <- vcov(m2.lm)[3,2]
scf <- vcov(yall.lm)[3,2]

# indirect effects
abc <- a*b*c
se.abc <- se.indirect3 (a,se.a,b,se.b,c,se.c)
t.abc <- abc/se.abc
r.abc <- abc/t

af <- a*f
se.af <- se.indirect2 (a,se.a,f,se.f)
t.af <- af/se.af
r.af <- af/t

ec <- e*c 
se.ec <- se.indirect2 (e,se.e,c,se.c)
t.ec <- ec/se.ec
r.ec <- ec/t

# total indirect effect x -> y
ind.xy <- abc + af + ec

sabcaf <- b*c*f*se.a^2 + a^2*b*scf
sabcec <- a*c^2*sbe + a*b*e*se.c^2
safec  <- a*e*scf
se.ind.xy <- sqrt(se.abc^2 + se.af^2 + se.ec^2 + 2*(sabcaf + sabcec + safec))
t.ind.xy <- ind.xy/se.ind.xy

# indirect effect ratio

med.ratio1 <- (t-tpr)/t
med.ratio2 <- ind.xy/t

# Create a table and return values

rtab <- matrix(nrow=11,ncol=4)
rownames(rtab) <- c("a","b","c","e","f","abc","af","ec","ind.xy","t","t'")
colnames(rtab) <- c("Effect","SE","t-ratio","Med.Ratio")
rtab[1,] <- c((a),(se.a),(t.a),"--")
rtab[2,] <- c((b),(se.b),(t.b),"--")
rtab[3,] <- c((c),(se.c),(t.c),"--")
rtab[4,] <- c((e),(se.e),(t.e),"--")
rtab[5,] <- c((f),(se.f),(t.f),"--")
rtab[6,] <- c((abc),(se.abc),(t.abc),(r.abc))
rtab[7,] <- c((af),(se.af),(t.af),(r.af))
rtab[8,] <- c((ec),(se.ec),(t.ec),(r.ec))
rtab[9,] <- c((ind.xy),(se.ind.xy),(t.ind.xy),(med.ratio2))
rtab[10,] <- c((t),(se.t),(t.t),"--")
rtab[11,] <- c((tpr),(se.tpr),(t.tpr),(med.ratio1))
return(rtab) 

}

