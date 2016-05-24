### R code from vignette source 'sweetgum.rnw'

###################################################
### code chunk number 1: sweetgum.rnw:4-5
###################################################
options(width=67)


###################################################
### code chunk number 2: Read sweetgum
###################################################
raw.data <- scan("../../data/TX_SGUM2.DAT", 
                 what = "", sep = "\n")
length(raw.data)


###################################################
### code chunk number 3: sweetgum.rnw:33-34
###################################################
raw.data <- raw.data[-c(1:26, 1101)]


###################################################
### code chunk number 4: sweetgum.rnw:43-44
###################################################
metadata <- grep("SWEETGUM", raw.data)


###################################################
### code chunk number 5: sweetgum.rnw:47-49 (eval = FALSE)
###################################################
## metadata <- grep("SWEETGUM", raw.data)
## cbind(metadata, raw.data[metadata])


###################################################
### code chunk number 6: sweetgum.rnw:56-58
###################################################
substr(raw.data[627], 1, 1) <- "4"
substr(raw.data[910], 1, 1) <- "5"


###################################################
### code chunk number 7: sweetgum.rnw:68-73
###################################################
for (i in 1:length(raw.data)) {
  if(substr(raw.data[i], 57, 64) != "SWEETGUM") 
     raw.data[i] <- paste(substr(raw.data[i - 1], 1, 10), 
                          raw.data[i], sep="")
}


###################################################
### code chunk number 8: sweetgum.rnw:80-85
###################################################
tree.data <- raw.data[metadata]
length(tree.data)

sections.data <- raw.data[-metadata]
length(sections.data)


###################################################
### code chunk number 9: sweetgum.rnw:90-103
###################################################
sweetgum <- 
    data.frame(plot = factor(substr(tree.data, 1, 5)),
               tree = substr(tree.data, 6, 10),
               dbh.in = substr(tree.data, 21, 26),
               stump.ht.ft = substr(tree.data, 27, 32),
               height.ft = substr(tree.data, 39, 44))

sections <- 
     data.frame(plot = factor(substr(sections.data, 1, 5)),
                tree = substr(sections.data, 6, 10),
         	meas.ln.ft = substr(sections.data, 11, 16),
         	meas.dob.in = substr(sections.data, 20, 25),
         	meas.dib.in = substr(sections.data, 26, 31))


###################################################
### code chunk number 10: sweetgum.rnw:110-112
###################################################
sapply(sweetgum, class)
sapply(sections, class)


###################################################
### code chunk number 11: sweetgum.rnw:118-122
###################################################
for (i in 3:5) {
   sweetgum[,i] <- as.numeric(as.character(sweetgum[,i]))
   sections[,i] <- as.numeric(as.character(sections[,i]))
}


###################################################
### code chunk number 12: sweetgum.rnw:129-132
###################################################
all.meas <- merge(sweetgum, sections, all = TRUE)
dim(all.meas)
names(all.meas)


###################################################
### code chunk number 13: sweetgum.rnw:137-143
###################################################
all.meas$meas.ht.ft <- with(all.meas, 
                            meas.ln.ft + stump.ht.ft)
all.meas$meas.ht.m <- all.meas$meas.ht.ft / 3.2808399
all.meas$meas.dob.cm <- all.meas$meas.dob.in * 2.54
sweetgum$height.m <- sweetgum$height.ft / 3.2808399
sweetgum$dbh.cm <- sweetgum$dbh.in * 2.54


###################################################
### code chunk number 14: sweetgum.rnw:163-176
###################################################
spline.vol.m3 <- function(hts.m, 
                          ds.cm, 
                          max.ht.m, 
                          min.ht.m = 0) {  
  rs.cm <- c(ds.cm[order(hts.m)] / 2, 0)
  hts.m <- c(hts.m[order(hts.m)], max.ht.m)
  taper <- splinefun(hts.m, rs.cm)
  volume <- integrate(f = function(x) 
                          pi * (taper(pmax(x,0))/100)^2,
                      lower = min.ht.m, 
                      upper = max.ht.m)$value
  return(volume)
}


###################################################
### code chunk number 15: sweetgum.rnw:182-188
###################################################
sweetgum$vol.m3 <- 
  mapply(spline.vol.m3, 
         hts.m = split(all.meas$meas.ht.m, all.meas$tree),
         ds.cm = split(all.meas$meas.dob.cm, all.meas$tree),
         max.ht.m = as.list(sweetgum$height.m),
         min.ht.m = 0.3)


###################################################
### code chunk number 16: fig-sgum-check
###################################################
par(las = 1)
plot(sweetgum$vol.m3, 
  (sweetgum$dbh.cm/200)^2 * pi * sweetgum$height.m / 2,
  ylab = expression(paste("Second-degree paraboloid volume (",
      m^3, ")", sep="")),
  xlab = expression(paste("Integrated spline volume (",
      m^3, ")", sep="")))
abline(0, 1, col="darkgrey")


###################################################
### code chunk number 17: sgum-check
###################################################
par(las = 1)
plot(sweetgum$vol.m3, 
  (sweetgum$dbh.cm/200)^2 * pi * sweetgum$height.m / 2,
  ylab = expression(paste("Second-degree paraboloid volume (",
      m^3, ")", sep="")),
  xlab = expression(paste("Integrated spline volume (",
      m^3, ")", sep="")))
abline(0, 1, col="darkgrey")


