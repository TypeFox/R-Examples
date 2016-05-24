### R code from vignette source 'paleofire-paper.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: loadPackage
###################################################
#install.packages("paleofire",repo="http://cran.r-project.org")
library(paleofire)


###################################################
### code chunk number 2: selectSites
###################################################
ID <- pfSiteSel(lat>30 & lat<90, long>-100 & long<(-50), 
                date_int<=2500)
length(ID$id_site)


###################################################
### code chunk number 3: selectSites
###################################################
sumID <- summary(ID)
ID <- pfSiteSel(id_site %in% ID$id_site & num_samp>=20)
length(ID$id_site)


###################################################
### code chunk number 4: a
###################################################
plot(ID, zoom="world")


###################################################
### code chunk number 5: fig1
###################################################
plot(ID, zoom="world")


###################################################
### code chunk number 6: TRANSFORM
###################################################
TR1 <- pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"))


###################################################
### code chunk number 7: ADDDATA
###################################################
## Add Ben lake and Small lake data to the 
# analysis (Senici et al., 2013) 
download.file(url="http://blarquez.com/public/data/data_cageo.zip", 
                destfile="data_cageo.zip")
unzip("data_cageo.zip")
mydata=pfAddData(files=c("Ben.csv","Small.csv"), 
                   metadata="metadata.csv", type="CharAnalysis")

## Transform:
TR2 <- pfTransform(ID,add=mydata,BasePeriod=c(200,4000),
                   method=c("MinMax","Box-Cox","MinMax","Z-Score"))
## Delete downloaded files
file.remove(c("Ben.csv","Small.csv","data_cageo.zip","metadata.csv"))


###################################################
### code chunk number 8: composite1
###################################################
COMP1 <- pfComposite(TR2, binning=TRUE, 
                     bins=seq(from=0,to=11000, by=500))


###################################################
### code chunk number 9: composite2
###################################################
COMP2 <- pfCompositeLF(TR2, tarAge=seq(-50,12000,20), 
                       binhw=10, hw=500, nboot=100)


###################################################
### code chunk number 10: plotting
###################################################
par(mfrow=c(2,1))
plot(COMP1,conf=c(0.025,0.975),main="(a)")
plot(COMP2,conf=c(0.05,0.95),main="(b)")


###################################################
### code chunk number 11: fig2
###################################################
par(mfrow=c(2,1))
plot(COMP1,conf=c(0.025,0.975),main="(a)")
plot(COMP2,conf=c(0.05,0.95),main="(b)")


###################################################
### code chunk number 12: circ (eval = FALSE)
###################################################
## circboot <- pfCircular(COMP1, b=NULL, nboot=100,
##                         conf=c(0.005,0.025,0.975,0.995))
## plot(circboot)


