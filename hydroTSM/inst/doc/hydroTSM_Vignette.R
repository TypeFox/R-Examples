### R code from vignette source 'hydroTSM_Vignette.Rnw'

###################################################
### code chunk number 1: hydroTSM_Vignette.Rnw:27-28 (eval = FALSE)
###################################################
## install.packages("hydroTSM")


###################################################
### code chunk number 2: hydroTSM_Vignette.Rnw:34-35 (eval = FALSE)
###################################################
## install.packages("hydroTSM",, "http://rforge.net/", type="source")


###################################################
### code chunk number 3: hydroTSM_Vignette.Rnw:45-46
###################################################
library(hydroTSM)


###################################################
### code chunk number 4: hydroTSM_Vignette.Rnw:51-52
###################################################
data(SanMartinoPPts)


###################################################
### code chunk number 5: hydroTSM_Vignette.Rnw:57-58
###################################################
x <- window(SanMartinoPPts, start=as.Date("1985-01-01"))


###################################################
### code chunk number 6: hydroTSM_Vignette.Rnw:62-63
###################################################
( m <- daily2monthly(x, FUN=sum) )


###################################################
### code chunk number 7: hydroTSM_Vignette.Rnw:67-68
###################################################
dates <- time(x)


###################################################
### code chunk number 8: hydroTSM_Vignette.Rnw:72-73
###################################################
( nyears <- yip(from=start(x), to=end(x), out.type="nmbr" ) )


###################################################
### code chunk number 9: hydroTSM_Vignette.Rnw:85-86
###################################################
smry(x)


###################################################
### code chunk number 10: hydroTSM_Vignette.Rnw:93-95
###################################################
hydroplot(x, var.type="Precipitation", main="at San Martino", 
          pfreq = "dm", from="1987-01-01")


###################################################
### code chunk number 11: hydroTSM_Vignette.Rnw:100-101
###################################################
dwi(x)


###################################################
### code chunk number 12: hydroTSM_Vignette.Rnw:105-106
###################################################
dwi(x, out.unit="mpy")


###################################################
### code chunk number 13: hydroTSM_Vignette.Rnw:112-124
###################################################
# Daily zoo to monthly zoo
m <- daily2monthly(x, FUN=sum, na.rm=TRUE)
     
# Creating a matrix with monthly values per year in each column
M <- matrix(m, ncol=12, byrow=TRUE)
colnames(M) <- month.abb
rownames(M) <- unique(format(time(m), "%Y"))
     
# Plotting the monthly precipitation values
require(lattice)
print(matrixplot(M, ColorRamp="Precipitation", 
           main="Monthly precipitation at San Martino st., [mm/month]"))


###################################################
### code chunk number 14: hydroTSM_Vignette.Rnw:137-138
###################################################
daily2annual(x, FUN=sum, na.rm=TRUE)


###################################################
### code chunk number 15: hydroTSM_Vignette.Rnw:145-146
###################################################
mean( daily2annual(x, FUN=sum, na.rm=TRUE) )


###################################################
### code chunk number 16: hydroTSM_Vignette.Rnw:155-156
###################################################
annualfunction(x, FUN=sum, na.rm=TRUE) / nyears


###################################################
### code chunk number 17: hydroTSM_Vignette.Rnw:171-172
###################################################
monthlyfunction(m, FUN=median, na.rm=TRUE)


###################################################
### code chunk number 18: hydroTSM_Vignette.Rnw:176-177
###################################################
cmonth <- format(time(m), "%b")


###################################################
### code chunk number 19: hydroTSM_Vignette.Rnw:181-182
###################################################
months <- factor(cmonth, levels=unique(cmonth), ordered=TRUE)


###################################################
### code chunk number 20: hydroTSM_Vignette.Rnw:186-188
###################################################
boxplot( coredata(m) ~ months, col="lightblue", main="Monthly Precipitation", 
         ylab="Precipitation, [mm]", xlab="Month")


###################################################
### code chunk number 21: hydroTSM_Vignette.Rnw:203-204
###################################################
seasonalfunction(x, FUN=sum, na.rm=TRUE) / nyears


###################################################
### code chunk number 22: hydroTSM_Vignette.Rnw:208-212
###################################################
( DJF <- dm2seasonal(x, season="DJF", FUN=sum) )
( MAM <- dm2seasonal(m, season="MAM", FUN=sum) )
( JJA <- dm2seasonal(m, season="JJA", FUN=sum) )
( SON <- dm2seasonal(m, season="SON", FUN=sum) )


###################################################
### code chunk number 23: hydroTSM_Vignette.Rnw:217-218
###################################################
hydroplot(x, pfreq="seasonal", FUN=sum, stype="default")


###################################################
### code chunk number 24: hydroTSM_Vignette.Rnw:236-237
###################################################
data(SanMartinoPPts)


###################################################
### code chunk number 25: hydroTSM_Vignette.Rnw:241-242
###################################################
x <- window(SanMartinoPPts, start=as.Date("1988-01-01"))


###################################################
### code chunk number 26: hydroTSM_Vignette.Rnw:246-247
###################################################
hydroplot(x,  ptype="ts", pfreq="o", var.unit="mm")


###################################################
### code chunk number 27: hydroTSM_Vignette.Rnw:258-259
###################################################
( R10mm <- length( x[x>10] ) )


###################################################
### code chunk number 28: hydroTSM_Vignette.Rnw:272-273
###################################################
wet.index <- which(x >= 1)


###################################################
### code chunk number 29: hydroTSM_Vignette.Rnw:278-279
###################################################
( PRwn95 <- quantile(x[wet.index], probs=0.95, na.rm=TRUE) )


###################################################
### code chunk number 30: hydroTSM_Vignette.Rnw:287-288
###################################################
(very.wet.index <- which(x >= PRwn95))


###################################################
### code chunk number 31: hydroTSM_Vignette.Rnw:293-294
###################################################
( R95p <- sum(x[very.wet.index]) )


###################################################
### code chunk number 32: hydroTSM_Vignette.Rnw:310-314
###################################################
x.5max <- rollapply(data=x, width=5, FUN=sum, fill=NA, partial= TRUE, 
                    align="center")

hydroplot(x.5max,  ptype="ts+boxplot", pfreq="o", var.unit="mm")


###################################################
### code chunk number 33: hydroTSM_Vignette.Rnw:318-319
###################################################
(x.5max.annual <- daily2annual(x.5max, FUN=max, na.rm=TRUE))


###################################################
### code chunk number 34: hydroTSM_Vignette.Rnw:335-338
###################################################
sessionInfo()$platform
sessionInfo()$R.version$version.string 
paste("hydroTSM", sessionInfo()$otherPkgs$hydroTSM$Version)


