## ------------------------------------------------------------------------
library(hextri)
data(airquality)
airquality$o3group<-with(airquality, cut(Ozone, c(0,18,60,Inf)))

with(na.omit(airquality), 
  hextri(Solar.R, Temp, class=o3group, colours=c("skyblue","grey60","goldenrod"), style="size",
            xlab="Solar radiation (langley)", ylab="Temperature (F)")
)

## ---- fig.height=4, fig.width=12-----------------------------------------
library(lattice)
xyplot(Temp~Solar.R|equal.count(Wind,4), groups=o3group, panel=panel.hextri,
  data=na.omit(airquality), colours=c("royalblue","grey60","goldenrod"),
  strip=strip.custom(var.name="Wind Speed"),
  xlab="Solar Radiation (langley)",ylab="Temperature (F)")

## ------------------------------------------------------------------------
data(NHANES, package="hexbin")
with(na.omit(NHANES[,-8]), hextri(Age,Hemoglobin, class=Sex,colour=c("orange","purple"),
    nbins=20,xlab="Age",ylab="Serum haemoglobin"))

## ------------------------------------------------------------------------
with(na.omit(NHANES[,-8]), hextri(Age,Hemoglobin, class=Sex,colour=c("orange","purple"),
    nbins=20,xlab="Age",ylab="Serum haemoglobin", diffuse=TRUE))

## ------------------------------------------------------------------------
with(na.omit(NHANES[,-8]), hextri(Age,Hemoglobin, class=Sex,colour=c("orange","purple"),
    nbins=20,xlab="Age",ylab="Serum haemoglobin",style="size"))

with(na.omit(NHANES[,-8]), hextri(Age,Hemoglobin, class=Sex,colour=c("orange","purple"),
    nbins=20,xlab="Age",ylab="Serum haemoglobin", diffuse=TRUE,style="size"))

## ------------------------------------------------------------------------
xyplot(Hemoglobin~Age|equal.count(Diet.Iron, 6),groups=Sex, data=na.omit(NHANES[,-8]), 
  colour=c("orange","purple"),panel=panel.hextri, 
  strip=strip.custom(var.name="Dietary iron"),style="size",diffuse=TRUE)

## ------------------------------------------------------------------------
xx<-rnorm(1000)
yy<-rnorm(1000)
cc<-cut(xx*yy,c(-Inf,-.4,0,.4,Inf))

hextri(xx,yy, class=cc, colour=c("#FEEDDE", "#FDD0A2", "#FDAE6B", "#8C2D04"),
   nbins=20, style="size")

hextri(xx,yy, class=cc, colour=c("#FEEDDE", "#8C2D04","#FDD0A2", "#FDAE6B") ,
   nbins=20, style="size")

hextri(xx,yy, class=cc, colour=c( "#8C2D04","#FDD0A2", "#FDAE6B","#FEEDDE") ,
   nbins=20, style="size")

