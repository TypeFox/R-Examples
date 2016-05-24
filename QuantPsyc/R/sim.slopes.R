"sim.slopes" <-
function (mod, z, zsd=1, mcz=FALSE)
{ 
if(!mcz)
{ z <- z - mean(z, na.rm=TRUE)}
else { z <- z }
zhi <- mean(z, na.rm=TRUE) + zsd*sd(z, na.rm=TRUE)
zlo <- mean(z, na.rm=TRUE) - zsd*sd(z, na.rm=TRUE)
zme <- mean(z, na.rm=TRUE)
b0 <- summary(mod)$coef[1,1]
b1 <- summary(mod)$coef[2,1]
b2 <- summary(mod)$coef[3,1]
b3 <- summary(mod)$coef[4,1]

x.zhi <- (b1 + b3*zhi)    # simple slope
x.zlo <- (b1 + b3*zlo)    # simple slope
x.zme <- (b1 + b3*zme)
int.zhi <- (b0 + b2*zhi)  # simple intercept
int.zlo <- (b0 + b2*zlo)  # simple intercept
int.zme <- (b0 + b2*zme)


seb.zhi <- sqrt(vcov(mod)[2,2] + 2*zhi*vcov(mod)[2,4] + zhi^2*vcov(mod)[4,4])
seb.zlo <- sqrt(vcov(mod)[2,2] + 2*zlo*vcov(mod)[2,4] + zlo^2*vcov(mod)[4,4])
seb.zme <- sqrt(vcov(mod)[2,2] + 2*zme*vcov(mod)[2,4] + zme^2*vcov(mod)[4,4])
td <- qt(.975, df = summary(mod)$df[2])

zhi.u <- x.zhi + td*seb.zhi  # upper limit of simple slope
zhi.l <- x.zhi - td*seb.zhi  # lower limit of simple slope
zlo.u <- x.zlo + td*seb.zlo  # upper limit of simple slope
zlo.l <- x.zlo - td*seb.zlo  # lower limit of simple slope
zme.u <- x.zme + td*seb.zme
zme.l <- x.zme - td*seb.zme

# create a matrix with int & slope of y~x at 2 levels of z; se.z; 95% CI 
mat <- matrix(NA,3,5)
colnames(mat) <- c("INT", "Slope", "SE", "LCL", "UCL")
rownames(mat) <- c("at zHigh", "at zMean", "at zLow")
mat[1,] <- c(int.zhi,x.zhi,seb.zhi,zhi.l,zhi.u)
mat[2,] <- c(int.zme,x.zme,seb.zme,zme.l,zme.u)
mat[3,] <- c(int.zlo,x.zlo,seb.zlo,zlo.l,zlo.u)
return(data.frame(mat))

}

