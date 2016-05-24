"proximal.med" <-
function (data)
{

# Models needed to test mediation
yx.lm <- lm (y ~ x, data)# t
yall.lm <- lm (y ~ m + x, data) # b tpr
m1.lm <- lm (m ~ x, data)# a

# individual effects
a <- summary(m1.lm)$coef[2,1]
se.a <- summary(m1.lm)$coef[2,2]
t.a <- summary(m1.lm)$coef[2,3]

b <- summary(yall.lm)$coef[2,1]
se.b <- summary(yall.lm)$coef[2,2]
t.b <- summary(yall.lm)$coef[2,3]

t <- summary(yx.lm)$coef[2,1]
se.t <- summary(yx.lm)$coef[2,2]
t.t <- summary(yx.lm)$coef[2,3]

tpr <- summary(yall.lm)$coef[3,1]
se.tpr <- summary(yall.lm)$coef[3,2]
t.tpr <- summary(yall.lm)$coef[3,3]

# indirect effects
ab <- a*b
se.ab <- se.indirect2 (a,se.a,b,se.b)
t.ab <- ab/se.ab
r.ab <- ab/t
ar.se <- aroian.se.indirect2(a,se.a,b,se.b)
go.se <- goodman.se.indirect2(a,se.a,b,se.b)

# Create a table and return values
rtab <- matrix(nrow=7,ncol=4)
rownames(rtab) <- c("a","b","t","t'","ab","Aroian","Goodman")
colnames(rtab) <- c("Effect","SE","t-ratio","Med.Ratio")
rtab[1,] <- c(a,se.a,t.a,"--")
rtab[2,] <- c(b,se.b,t.b,"--")
rtab[3,] <- c(t,se.t,t.t,"--")
rtab[4,] <- c(tpr,se.tpr,t.tpr,"--")
rtab[5,] <- c(ab,se.ab,t.ab,r.ab)
rtab[6,] <- c("--",ar.se,"--","--")
rtab[7,] <- c("--",go.se,"--","--")
return((rtab))

}

