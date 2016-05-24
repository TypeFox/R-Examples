### R code from vignette source 'data_management.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")
if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(lattice, quietly=TRUE)) install.packages("lattice")
if(!require(maptools, quietly=TRUE)) install.packages("maptools")
#if(!require(RODBC, quietly=TRUE)) install.packages("RODBC")
if(!require(colorspace, quietly=TRUE)) install.packages("colorspace")
if(!require(quantreg, quietly=TRUE)) install.packages("quantreg")
lattice.options(default.theme = canonical.theme(color = FALSE))

options(repos="http://www.bioconductor.org/")
if(!require(hexbin, quietly=TRUE)) install.packages("hexbin")

options(width = 65)

rm(list=ls())

source("../../scripts/functions.R")



###################################################
### code chunk number 2: data_management.rnw:122-124
###################################################
fia.plots <- read.table("../../data/fia_plots.csv", sep = ",",
                        header = TRUE, row.names = 1)


###################################################
### code chunk number 3: data_management.rnw:142-143
###################################################
fia.plots <- read.csv("../../data/fia_plots.csv")


###################################################
### code chunk number 4: data_management.rnw:149-150
###################################################
class(fia.plots)


###################################################
### code chunk number 5: data_management.rnw:162-163 (eval = FALSE)
###################################################
## write.table(fia.plots, file = "fia-plots.csv")


###################################################
### code chunk number 6: data_management.rnw:203-210
###################################################
fvs.trees <- 
  read.fwf("../../data/stnd73.fvs", 
           widths = c(4, 3, 6, 1, 3, 4, -3, 3, -7, 1),
           as.is = FALSE, row.names = NULL,
           col.names = c("plot", "tree", "tree.count",
             "history", "species", "dbh", "live.tht",
             "crown.code"))


###################################################
### code chunk number 7: data_management.rnw:219-220
###################################################
library(lattice)


###################################################
### code chunk number 8: fig-fvs-xyplot
###################################################
xyplot(live.tht ~ dbh | species, data = fvs.trees) 


###################################################
### code chunk number 9: data_management.rnw:229-232
###################################################
print(
xyplot(live.tht ~ dbh | species, data = fvs.trees) 
)


###################################################
### code chunk number 10: data_management.rnw:269-271
###################################################
eg <- scan(file = "../../data/scan-example.txt", 
           sep = "\n", what = "")


###################################################
### code chunk number 11: data_management.rnw:282-285
###################################################
n.in <- length(eg)
eg.trees <- eg.plots <- vector(mode = "list", length = n.in)
plot.n <- tree.n <- 1


###################################################
### code chunk number 12: data_management.rnw:298-319
###################################################
for (i in 1 : n.in) {
  chunk <- eg[[i]]
  if (substr(chunk, 1, 4) == "Plot") {
     plot.id <- as.numeric(substr(chunk, 6, 9))
     crew.id <- substr(chunk, 16, 16)
     comments <- ifelse(nchar(chunk) > 17, 
                        substr(chunk, 17, nchar(chunk)),
                        "")
     eg.plots[[plot.n]] <- 
               list(plot.id, crew.id, comments)
     plot.n <- plot.n + 1
   } else {
     tree <- strsplit(chunk, " +")[[1]]
     tree.id <- as.character(tree[1])
     species <- as.character(tree[2])
     dbh.cm <- as.numeric(tree[3])
     eg.trees[[tree.n]] <- 
               list(plot.id, tree.id, species, dbh.cm)
     tree.n <- tree.n + 1
   }
}


###################################################
### code chunk number 13: data_management.rnw:328-331
###################################################
eg.plots <- as.data.frame(do.call(rbind, eg.plots))
names(eg.plots) <- c("plot", "crew", "comments")
eg.plots


###################################################
### code chunk number 14: data_management.rnw:334-337
###################################################
eg.trees <- as.data.frame(do.call(rbind, eg.trees))
names(eg.trees) <- c("plot", "tree", "species", "dbh.cm")
eg.trees


###################################################
### code chunk number 15: data_management.rnw:426-437 (eval = FALSE)
###################################################
## library(RPostgreSQL) 
## drv <- dbDriver("PostgreSQL")
## con <- dbConnect(drv, 
##                  dbname="forestco",
##                  user="hamannj",
##                  host="localhost")
## sql.command <- 
##   sprintf( "select * from plots where plottype = 'fixed';" )
## rs <- dbSendQuery(con, statement = sql.command )
## fixed.plots <- fetch(rs, n = -1)
## dbDisconnect(con)


###################################################
### code chunk number 16: stands
###################################################

library(foreign)
stands <- read.dbf("../../data/stands.dbf")



###################################################
### code chunk number 17: data_management.rnw:582-583
###################################################
names(stands) 


###################################################
### code chunk number 18: fig-stands-hexbin
###################################################
stands.non.zero <- stands[stands$QMD > 0,]
plot(hexbin(stands.non.zero$QMD*2.54 ~ 
            stands.non.zero$TPA*2.47), 
     ylab = "Average Diameter (cm)",
     xlab = "Stem Density (stems/ha)")



###################################################
### code chunk number 19: data_management.rnw:604-605
###################################################
stands.non.zero <- stands[stands$QMD > 0,]
plot(hexbin(stands.non.zero$QMD*2.54 ~ 
            stands.non.zero$TPA*2.47), 
     ylab = "Average Diameter (cm)",
     xlab = "Stem Density (stems/ha)")



###################################################
### code chunk number 20: data_management.rnw:666-671
###################################################

library(maptools)

stands <- readShapePoly("../../data/stands.shp")



###################################################
### code chunk number 21: data_management.rnw:675-679
###################################################

sum( stands$AREA ) / 43560.0




###################################################
### code chunk number 22: data_management.rnw:683-684
###################################################
nrow(stands)


###################################################
### code chunk number 23: data_management.rnw:689-692
###################################################

names(stands)



###################################################
### code chunk number 24: fig-plot-stands
###################################################

plot(stands, axes = TRUE)



###################################################
### code chunk number 25: data_management.rnw:706-707
###################################################

plot(stands, axes = TRUE)



###################################################
### code chunk number 26: data_management.rnw:870-872
###################################################
herbdata <- read.table("../../data/herbdata.txt", 
                       header = TRUE, sep = ",")


###################################################
### code chunk number 27: data_management.rnw:879-880 (eval = FALSE)
###################################################
## str(herbdata)


###################################################
### code chunk number 28: data_management.rnw:882-883
###################################################
str(herbdata, vec.len = 1)


###################################################
### code chunk number 29: data_management.rnw:895-897
###################################################
herbdata$date <- as.POSIXct(strptime(herbdata$date, 
                                     "%m/%d/%Y"))


###################################################
### code chunk number 30: fig-herb-coplot
###################################################
coplot(height ~ dia | treat * rep, type = "p", 
       data = herbdata[herbdata$isalive == 1,], 
       ylab = "Height (cm)", xlab = "Basal Diameter (mm)")


###################################################
### code chunk number 31: data_management.rnw:911-912
###################################################
coplot(height ~ dia | treat * rep, type = "p", 
       data = herbdata[herbdata$isalive == 1,], 
       ylab = "Height (cm)", xlab = "Basal Diameter (mm)")


###################################################
### code chunk number 32: data_management.rnw:974-975
###################################################
head(herbdata[is.na(herbdata$height),])


###################################################
### code chunk number 33: data_management.rnw:983-984
###################################################
table(complete.cases(herbdata), herbdata$isalive)


###################################################
### code chunk number 34: fig-herb-date-coplot
###################################################
coplot(height ~ dia | treat * factor(date),
       data = herbdata[herbdata$isalive == 1,],
       type = "p",
       ylab = "Height (cm)",
       xlab = "Basal Diameter (mm)")


###################################################
### code chunk number 35: data_management.rnw:1020-1021
###################################################
coplot(height ~ dia | treat * factor(date),
       data = herbdata[herbdata$isalive == 1,],
       type = "p",
       ylab = "Height (cm)",
       xlab = "Basal Diameter (mm)")


###################################################
### code chunk number 36: data_management.rnw:1035-1037
###################################################
levels(herbdata$treat) 
levels(herbdata$rep) 


###################################################
### code chunk number 37: data_management.rnw:1041-1042
###################################################
sort(unique(herbdata$date))


###################################################
### code chunk number 38: data_management.rnw:1048-1052
###################################################
bad.index <- herbdata$treat == levels(herbdata$treat)[1] & 
             herbdata$rep == levels(herbdata$rep)[2] & 
             herbdata$date == sort(unique(herbdata$date))[7]
bad.data <- herbdata[bad.index,]


###################################################
### code chunk number 39: data_management.rnw:1055-1056
###################################################
print(head(bad.data))


###################################################
### code chunk number 40: data_management.rnw:1062-1064
###################################################
##bad.data[(nrow(bad.data)-3):nrow(bad.data),]
print(tail(bad.data))


###################################################
### code chunk number 41: data_management.rnw:1075-1076
###################################################
herbdata$dia[bad.index] <- herbdata$dia[bad.index] / 2.54


###################################################
### code chunk number 42: data_management.rnw:1082-1083
###################################################
herbdata$dbh[bad.index] <- herbdata$dbh[bad.index] / 2.54


###################################################
### code chunk number 43: data_management.rnw:1173-1178
###################################################
split.herb <- split(herbdata, herbdata$treat)
class(split.herb)
names(split.herb)
nrow(split.herb$CONTROL)
nrow(split.herb$OUST)


###################################################
### code chunk number 44: data_management.rnw:1183-1184
###################################################
names(split.herb$CONTROL)


###################################################
### code chunk number 45: data_management.rnw:1192-1193
###################################################
lapply(split.herb, nrow)


###################################################
### code chunk number 46: data_management.rnw:1244-1245
###################################################
sort(unique(herbdata$date)) 


###################################################
### code chunk number 47: data_management.rnw:1249-1253
###################################################
herbdata.shorter <- 
     herbdata[herbdata$date == max(herbdata$date), c(1,2,6,8)]
split.herb.shorter <- 
     split(herbdata.shorter, herbdata.shorter$treat)


###################################################
### code chunk number 48: make-st
###################################################
rt <- cbind(herbdata.shorter,  
            dc = cut(herbdata.shorter$dbh, 
            breaks = c(0, 50, 100, 150, 200, 300, 400, 999),
            labels = c("000--050", "050--100", "100--150",
                "150--200", "200--300", "300--400","400+")))

st <- aggregate(x = list(basal.area = pi/(4*10^2) * rt$dbh^2,
                         tht = rt$height, 
                         stems = rep(1, nrow(rt))),
                by = list(treat = rt$treat,
                          diac = rt$dc),
                FUN = sum)
st


###################################################
### code chunk number 49: data_management.rnw:1283-1285
###################################################
st$tht <- st$tht / st$stems / 100
st


###################################################
### code chunk number 50: oust
###################################################
cap <- "OUST herbicide trials."
st <- st[order(st$treat, st$diac),]
st$treat <- as.character(st$treat)
st$diac <- as.character(st$diac)
names(st) <- c("Treatment", "Dia. Class (mm)", 
               "Basal Area ($\\mbox{mm}^2$)", 
               "Mean Total Height (m)", "Stems")


###################################################
### code chunk number 51: data_management.rnw:1354-1355
###################################################
names(split.herb)


###################################################
### code chunk number 52: data_management.rnw:1366-1375
###################################################
areas <- 
   with(herbdata,
        aggregate(x = list(plot.bh.area = pi/400 * dbh^2,
                           plot.bas.area = pi/400 * dia^2),
                  by = list(treat = treat,
                           rep = rep,
                           date = date),
                  FUN = sum))



###################################################
### code chunk number 53: data_management.rnw:1380-1381
###################################################
areas[1:10,]


###################################################
### code chunk number 54: data_management.rnw:1396-1406
###################################################
areas <- 
   with(herbdata,
        aggregate(x = list(plot.bh.area = pi/400 * dbh^2,
                           plot.bas.area = pi/400 * dia^2),
                  by = list(treat = treat,
                           rep = rep,
                           date = date),
                  FUN = sum,
                  na.rm = TRUE))
areas[1:10,]


###################################################
### code chunk number 55: data_management.rnw:1413-1416
###################################################
final.data <- merge(herbdata, areas)
names(final.data)
head(final.data[,c(1,2,3,4,7,10)])


###################################################
### code chunk number 56: data_management.rnw:1434-1449
###################################################
show.cols.with.na <- function(x) {
  ## First, check that object is a data frame
  if (class(x) != "data.frame")      
    stop("x must be a data frame.\n")
  ## Count the missing values by column.
  missing.by.column <- colSums(is.na(x))
  ## Are any missing?
  if (sum(missing.by.column) == 0) {
    cat("No missing values.\n")
  } else {
    ## Only return columns with missing values.
    missing <- which(missing.by.column > 0)
    return(missing.by.column[missing])
  }
}


###################################################
### code chunk number 57: MDD
###################################################

stands <- readShapePoly("../../data/stands.shp", 
                        verbose = FALSE)

perim <- readShapePoly("../../data/boundary.shp", 
                        verbose = FALSE)



###################################################
### code chunk number 58: data_management.rnw:1551-1554
###################################################
stands <- readShapePoly("../../data/stands.shp", 
                        verbose=FALSE)
names(stands)


###################################################
### code chunk number 59: data_management.rnw:1565-1617 (eval = FALSE)
###################################################
## 
## ##stands$id <- 1:nrow(stands)
## 
## ## this can also be done away with...
## 
## blks <- as.numeric(1:length(stands$TAGE))[stands$TAGE == 0]
## stands$TAGE[blks] <- NA
## 
## blks <- as.numeric(1:length(stands$BHAGE))[stands$BHAGE == 0]
## stands$BHAGE[blks] <- NA
## 
## blks <- as.numeric(1:length(stands$DF_SITE))[stands$DF_SITE == 0]
## stands$DF_SITE[blks] <- NA
## 
## blks <- as.numeric(1:length(stands$TPA))[stands$TPA == 0]
## stands$TPA[blks] <- NA
## 
## blks <- as.numeric(1:length(stands$QMD))[stands$QMD == 0]
## stands$QMD[blks] <- NA
## 
## blks <- as.numeric(1:length(stands$BA))[stands$BA == 0]
## stands$BA[blks] <- NA
## 
## blks <- as.numeric(1:length(stands$TOTHT))[stands$TOTHT == 0]
## stands$TOTHT[blks] <- NA
## 
## blks <- as.numeric(1:length(stands$CUBVOL_AC))[stands$CUBVOL_AC == 0]
## stands$CUBVOL_AC[blks] <- NA
## 
## 
## ## blk.rows <- as.numeric(rownames(stands$att.data[stands$att.data$BHAGE == 0,]))
## ## stands$att.data[blk.rows,]$BHAGE <- NA
## 
## ## blk.rows <- as.numeric(rownames(stands$att.data[stands$att.data$DF_SITE == 0,]))
## ## stands$att.data[blk.rows,]$DF_SITE <- NA
## 
## ## blk.rows <- as.numeric(rownames(stands$att.data[stands$att.data$TPA == 0,]))
## ## stands$att.data[blk.rows,]$TPA <- NA
## 
## ## blk.rows <- as.numeric(rownames(stands$att.data[stands$att.data$QMD == 0,]))
## ## stands$att.data[blk.rows,]$QMD <- NA
## 
## ## blk.rows <- as.numeric(rownames(stands$att.data[stands$att.data$BA == 0,]))
## ## stands$att.data[blk.rows,]$BA <- NA
## 
## ## blk.rows <- as.numeric(rownames(stands$att.data[stands$att.data$TOTHT == 0,]))
## ## stands$att.data[blk.rows,]$TOTHT <- NA
## 
## ## blk.rows <- as.numeric(rownames(stands$att.data[stands$att.data$CUBVOL_AC == 0,]))
## ## stands$att.data[blk.rows,]$CUBVOL_AC <- NA
## 
## 


###################################################
### code chunk number 60: data_management.rnw:1623-1626
###################################################

names(stands)



###################################################
### code chunk number 61: data_management.rnw:1655-1658
###################################################

plots <- readShapePoints("../../data/plots.shp")



###################################################
### code chunk number 62: data_management.rnw:1667-1670 (eval = FALSE)
###################################################
## 
## plot(plots, add=TRUE, pch=46)
## 


###################################################
### code chunk number 63: fig-stands-and-plots
###################################################

lev <- as.numeric(stands$ALLOCATION)
fgs <- gray(length(levels(stands$ALLOCATION)):1 / 3) 

plot(stands,
     col=fgs[lev],
     add=FALSE,
     axes=TRUE)
title(paste("McDonald-Dunn Research Forest",
            "Stand Boundaries and Plot Locations",
            sep = "\n"))

legend(1280000, 365000, 
       levels(stands$ALLOCATION)[3:1], 
       fill = fgs[3:1],
       cex = 0.7, 
       title = "Land Allocations")

plot(plots, add=TRUE, pch=46)



###################################################
### code chunk number 64: fig-stands-and-plots
###################################################

lev <- as.numeric(stands$ALLOCATION)
fgs <- gray(length(levels(stands$ALLOCATION)):1 / 3) 

plot(stands,
     col=fgs[lev],
     add=FALSE,
     axes=TRUE)
title(paste("McDonald-Dunn Research Forest",
            "Stand Boundaries and Plot Locations",
            sep = "\n"))

legend(1280000, 365000, 
       levels(stands$ALLOCATION)[3:1], 
       fill = fgs[3:1],
       cex = 0.7, 
       title = "Land Allocations")

plot(plots, add=TRUE, pch=46)



###################################################
### code chunk number 65: read.dbf
###################################################

mdtrees <- read.dbf("../../data/mdtrees.dbf")
head(mdtrees)



###################################################
### code chunk number 66: data_management.rnw:1757-1760
###################################################

mdtrees$EXPF <- NA



###################################################
### code chunk number 67: data_management.rnw:1767-1774
###################################################

mdtrees$EXPF[mdtrees$SUBPLOT == 1] <- 
     20.0 / (0.0054541539 * 
             mdtrees$DBH[mdtrees$SUBPLOT == 1] ^2)
mdtrees$EXPF[mdtrees$SUBPLOT == 2] <- 43560 / (pi * 7.78^2)
mdtrees$EXPF[mdtrees$SUBPLOT == 3] <- 43560 / (pi * 15.56^2)



###################################################
### code chunk number 68: data_management.rnw:1780-1783
###################################################

head(mdtrees[, 3:11])



###################################################
### code chunk number 69: data_management.rnw:1800-1803
###################################################

trees.by.plot <- split(mdtrees, mdtrees$PLOT)



###################################################
### code chunk number 70: data_management.rnw:1829-1859
###################################################
get.plot.sums <- function(trs) {

# /******************************************************/
# /* Bruce, D.  1981.  Consistent height-growth and     */
# /*    growth-rate estimates for remeasured plots.     */
# /*    Forest Science 27:711-725.                      */
# /******************************************************/
  site.index.bruce.1981 <- function(tht, bha) {
  tht * exp(-21.663 * (3.744e-2 - (bha + 8.0)^ -0.809))
  }
  
  not.missing.dbh <- !is.na(trs$DBH)
  bh.idx <- not.missing.dbh & trs$THT > 4.5
  expf.tot <- sum(trs$EXPF)
  expf.bh <- sum(trs$EXPF[not.missing.dbh])
  ba <- sum(0.0054541539 * trs$DBH[not.missing.dbh] ^ 2 * 
        trs$EXPF[not.missing.dbh])
  qmd <- sqrt(ba / expf.bh / 0.0054541539)
  s.trs <- trs[trs$SITETREE == 1 & trs$SPCODE == "DF" & 
               !is.na(trs$THT),]
  nst <- nrow(s.trs)
  site.bar <- 
       ifelse(nst > 0,
              weighted.mean(site.index.bruce.1981(s.trs$THT,
                                                  s.trs$AGE), 
                            s.trs$EXPF),
              NA)
  return(c(nrow(trs), expf.bh, expf.tot, 
           ba, qmd, nst, site.bar))
}


###################################################
### code chunk number 71: data_management.rnw:1866-1870
###################################################

plot.sums  <- 
  data.frame(t(sapply(trees.by.plot, get.plot.sums)))



###################################################
### code chunk number 72: data_management.rnw:1882-1888
###################################################

plot.sums$id <- as.numeric(names(trees.by.plot))
names(plot.sums) <- c("trees","expf.bh","expf.tot",
                      "ba","qmd","nst","site","id")
print(head(plot.sums), digits=3)



###################################################
### code chunk number 73: data_management.rnw:1896-1904
###################################################

plot.id <- as.numeric(as.character(plots$UNIPLOT)) 
plot.centers <- data.frame(cbind(coordinates(plots), plot.id))
names(plot.centers) <- c("x","y","id")

final.plots <- merge(plot.centers, plot.sums, all = TRUE)
print(head(final.plots[,c(1:3,5:10)]), digits = 3)



###################################################
### code chunk number 74: data_management.rnw:1912-1915
###################################################

write.csv( final.plots, "../../data/final-plots.csv")



###################################################
### code chunk number 75: data_management.rnw:1943-1947
###################################################

leusch.ylds <- read.table("../../data/leuschner.txt", 
                          header = TRUE)



###################################################
### code chunk number 76: data_management.rnw:2383-2402
###################################################
Stangle("fia.rnw")
source("fia.R")

Stangle("gutten.rnw")
source("gutten.R")

Stangle("pref.rnw")
source("pref.R")

Stangle("stage.rnw")
source("stage.R")

Stangle("sweetgum.rnw")
source("sweetgum.R")

Stangle("ufc.rnw")
source("ufc.R")
system("rm -fr package-Ch2")
package.skeleton(name = "package-Ch2")


