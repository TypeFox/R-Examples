
###################################################
###
### code chunk: Ch03Init
###
### ARMD Trial: Visual acuity (Section 3.2)
### 
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("nlme")
packageVersion("reshape")
packageVersion("plyr")
packageVersion("lattice")
sessionInfo()

ylims <- c(0, 90) 

###################################################
### code chunk: R3.1
###################################################
data(armd.wide, armd0, package = "nlmeU")       # Data loaded
library(lattice)
armd0.subset <-                                 # Subset
  subset(armd0, as.numeric(subject) %in% seq(1, 240, 10))
xy1 <-                                          # Draft plot 
  xyplot(visual ~ jitter(time) | treat.f, 
         groups = subject, 
         data = armd0.subset, 
         type = "l", lty = 1)  
update(xy1,                                     # Fig. 3.1
       xlab = "Time (in weeks)",
       ylab = "Visual acuity",
       grid = "h")
detach(package:lattice)

###################################################
### code chunk: R3.2
###################################################
table(armd.wide$miss.pat)
with(armd.wide, table(miss.pat))
xtabs(~miss.pat, armd.wide)

###################################################
### code chunk: R3.3a
###################################################
attach(armd0)
flst <- list(time.f, treat.f)                    # "By" factors
(tN <-                                           # Counts
  tapply(visual, flst,
         FUN = function(x) length(x[!is.na(x)])))

###################################################
### code chunk: R3.3b
###################################################
tMn <- tapply(visual, flst, FUN = mean)          # Sample means
tMd <- tapply(visual, flst, FUN = median)        # Sample medians
colnames(res  <- cbind(tN, tMn, tMd))            # Column names
nms1 <- rep(c("P", "A"), 3)
nms2 <- rep(c("n", "Mean", "Mdn"), rep(2, 3))
colnames(res) <- paste(nms1, nms2, sep = ":")    # New column names
res
detach(armd0)

###################################################
### code chunk: R3.4 (Fig. 3.2, no horizontal grid lines)
###################################################
 library(lattice)
 bw1 <-                                            # Draft plot
   bwplot(visual ~ time.f | treat.f, 
          data = armd0)         
 xlims <- c("Base", "4\nwks", "12\nwks", "24\nwks", "52\nwks")
 update(bw1, xlim = xlims, pch = "|")              # Final plot
 detach(package:lattice)


###################################################
### Fig. 3.2 with horizontal grid lines added
###################################################

library(lattice)
myPanel <- function(x, y, subscripts, ...){
  panel.grid(h = -1, v = 0)    
  panel.bwplot(x, y, ...) 
}

bw1 <- bwplot(visual ~ time.f|treat.f,
       data = armd0,
       panel = myPanel, 
       aspect = 1)
xlims <- c("Base", "4\nwks", "12\nwks","24\nwks","52\nwks")
bw1a <- update(bw1, xlim = xlims, ylim = ylims, 
  ylab= "Visual acuity", pch = "|") 
print(bw1a)
detach(package:lattice)

###################################################
### code for figure 3.3 not included
###################################################

###################################################
### code chunk: R3.5a
###################################################
mnt.pat<-                                   # Monotone patterns
  c("----", "---X", "--XX", "-XXX", "XXXX") 
armd.wide.mnt <-                            # Data subset
  subset(armd.wide, miss.pat %in% mnt.pat)
dim(armd.wide.mnt)                          # Number of rows and cols  
levels(armd.wide.mnt$miss.pat)              # Some levels not needed


###################################################
### code chunk: R3.5b
###################################################
armd.wide.mnt1 <- 
  within(armd.wide.mnt, 
         {
           miss.pat <- factor(miss.pat, levels = mnt.pat) 
         })
levels(armd.wide.mnt1$miss.pat)              

###################################################
### code chunk: R3.5c
###################################################
with(armd.wide.mnt1, 
     {
       fl  <- list(treat.f, miss.pat)       # List of "by" factors
       tapply(subject, fl, FUN=function(x) length(x[!is.na(x)]))
     })

###################################################
### code chunk: Fig. 3.4
###################################################
library(lattice)

my.lowerPanel <-  ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
  panel.grid(h = -1, v = -1) 
  panel.xyplot(x, y, ...)
}

my.upperPanel <-  ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
  panel.xyplot(x, y, type = "n", ...)
  corx <- round(cor(x,y,use = "complete.obs"),2)
  abs.corx <- abs(corx)
  cex.value <- 3
  ltext(50,50, corx, cex = abs.corx* cex.value)
}

mySuperPanel <- function(z, subscripts, panel.subscripts,...){
  # z is data frame. Abbreviated variable names on the diagonal of the splom.
  panel.pairs(z, subscripts = subscripts,
              panel.subscripts = panel.subscripts,
              as.matrix = TRUE, 
              upper.panel = "my.upperPanel",
              lower.panel = "my.lowerPanel",
prepanel.limits = function(z) return(c(1, 90))
)}

splom.form <- formula( ~cbind(vis0, vis4, vis12, vis24, vis52))
armd.wide.tmp <- subset(armd.wide, miss.pat == "----",
     select = c(visual0, visual4, visual12, visual24, visual52))
names(armd.wide.tmp) <- c("vis0", "vis4","vis12","vis24","vis52") # Short names

splom.object <- splom(splom.form,
  data = armd.wide.tmp,
  par.settings = list(fontsize = list(points = 4), axis.text = list(cex = 0.9)),
  as.matrix =TRUE,  
  xlab = "",
  superpanel = mySuperPanel
)
print(splom.object)
detach(package:lattice)

###################################################
### code chunk: R3.6
###################################################
visual.x <- subset(armd.wide, select = c(visual0:visual52))
(varx <- var(visual.x, use = "complete.obs")) # Var-cov mtx
print(cor(visual.x, use = "complete.obs"),    # Corr mtx
      digits = 2) 
diag(varx)                  # Var-cov diagonal elements
cov2cor(varx)               # Corr mtx (alternative way) 


###################################################
### code chunk: Cleanup for ARMD study
###################################################
rm(xlims,ylims) 
rm(xy1, flst, tN, tMn, tMd, nms1, nms2, res)
rm(bw1, bw1a)
rm(mnt.pat, armd.wide.mnt, armd.wide.mnt1)
rm(armd.wide, armd.wide.tmp, armd0, armd0.subset)
rm(myPanel, my.lowerPanel, my.upperPanel,  mySuperPanel, splom.object, splom.form)
rm(visual.x, varx)

###################################################
### PRT Study: Muscle fiber specific-force (Section 3.3)
### code chunk: R3.7
###################################################
data(prt.subjects, prt, package = "nlmeU") # Data loaded
with(prt.subjects, tapply(bmi, prt.f, summary))
by(subset(prt.subjects, select = -id), prt.subjects$prt.f, summary)



###################################################
### code chunk: R3.8a
###################################################
fibL <- 
   with(prt, 
        tapply(spec.fo, 
               list(id = id, fiberF = fiber.f, occF = occ.f), 
               length))
dimnms <- dimnames(fibL) 
names(dimnms)         # Shortened names displayed
fibL["5", , ]         # Number of fiber measurements for subject 5
fibL["335", , ]       # Number of fiber measurements for subject 335


###################################################
### code chunk: R3.8b
###################################################
fibM <-
   with(prt, 
        tapply(spec.fo, 
               list(id = id, fiberF = fiber.f, occF = occ.f),
               mean))
fibM["5", , ]

 
###################################################
### code chunk: R3.9a
###################################################

library(reshape)         
idvar <- c("id", "prt.f", "fiber.f", "occ.f")
meas.var <- c("spec.fo", "iso.fo")
prtM <-                                        # Melting data
  melt(prt, id.var = idvar, measure.var = meas.var) 
dim(prtM)
head(prtM, n = 4)                              # First four rows
tail(prtM, n = 4)                              # Last four rows

###################################################
### code chunk: R3.9b
###################################################
prtC <- cast(prtM, fun.aggregate = mean)       # Casting data
names(prtC) 
names(prtC)[5:6] <- c("spec.foMn", "iso.foMn") # Names modified
head(prtC, n = 4)

###################################################
### code chunk: Code for Figs: 3.5, 3.6, 3.7 not shown
###################################################

###################################################
### code chunk: Cleanup for PRT study
###################################################
rm(fibL, fibM, dimnms) 
rm(idvar, meas.var, prtM, prtC)

###################################################
### SII Project: Gain in the math achievement-score (Section 3.4)
### code chunk: R3.10
###################################################

data(SIIdata, package = "nlmeU")
sapply(SIIdata, FUN = function(x) any(is.na(x)))
sum(as.numeric(is.na(SIIdata$mathknow)))
range(SIIdata$mathknow, na.rm = TRUE)

###################################################
### code chunk: R3.11
###################################################
(schlN <- xtabs(~schoolid, SIIdata))  # Number of pupils per school
range(schlN)
xtabs(~schlN) # Distribution of the number of pupils over schools

###################################################
### code chunk: R3.12
###################################################
attach(SIIdata)
##  Long output omitted
##  (mthgM <- by(cbind(mathgain, mathkind), schoolid, colMeans))
detach(SIIdata)

###################################################
### code chunk: R3.13a
###################################################
library(reshape)
idvars <- c("schoolid") 
mvars <- c("classid", "childid")
dtm1 <- melt(SIIdata, id.vars = idvars, measure.vars = mvars)
names(cst1  <- 
  cast(dtm1, 
       fun.aggregate = function(el) length(unique(el))))
names(cst1) <- c("schoolid", "clssn", "schlN")

###################################################
### code chunk: R3.13b
###################################################
mvars <- c("mathgain", "mathkind", "housepov")
dtm2 <- melt(SIIdata, id.vars = idvars, measure.vars = mvars)
names(cst2 <- cast(dtm2, fun.aggregate = mean))
names(cst2) <- c("schoolid", "mthgMn", "mthkMn", "housepov")


###################################################
### code chunk: R3.13c
###################################################
(schlDt  <- merge(cst1, cst2, sort = FALSE))
rm(cst1, cst2)

###################################################
### code chunk: R3.14a
###################################################
summary(schlDt$housepov)

###################################################
### code chunk: R3.14b
###################################################
library(lattice)
xyplot(mthgMn ~ housepov,                         # Fig. 3.8a
       schlDt, type = c("p", "smooth"), grid = TRUE) 
xyplot(mthgMn ~ mthkMn,                           # Fig. 3.8b
       schlDt, type = c("p", "smooth"), grid = TRUE) 

##################################################
## More detailed code for Figs. 3.8a and 3.8b
##################################################

#### Auxiliary step: Setting limits for y-axis
range(SIIdata$mathgain)
ylims <- c(-120, 260)

xyplot(mthgMn ~ housepov,                         # Fig. 3.8a
       schlDt, type = c("p", "smooth"), 
       ylim = ylims, grid = TRUE) 
xyplot(mthgMn ~ mthkMn,                           # Fig. 3.8b
       schlDt, type = c("p", "smooth"), 
       ylim = ylims, grid = TRUE) 

###################################################
### code chunk: R 3.15
###################################################
(clssN <- xtabs(~ classid, SIIdata))
sum(clssN)                           # Total number of pupils
range(clssN)
(clssCnt <- xtabs(~clssN)) # Distribution of no. of pupils/classes
sum(clssCnt)                         # Total number of classes


###################################################
### Preparatory steps for R3.16.
### Steps similar to those shown in R13.3.
###################################################

## See R2.13b how to obtain vector nms2 containing class-specific variables
nms2 <- c("yearstea","mathknow", "mathprep", "classid") 
library(reshape)
(idvars <- c(nms2, "housepov"))
mvars   <- c("childid")
dtm1    <- melt(SIIdata, id.vars = idvars, measure.vars = mvars)
cst1    <- cast(dtm1, 
             classid + housepov + mathknow + mathprep ~ variable, 
             function(el) length(unique(el)))
names(cst1)[5] <- "clssN"

mvars   <- c("mathgain","mathkind")
dtm2    <- melt(SIIdata, id.vars = "classid", measure.vars = mvars)
cst2    <- cast(dtm2, classid ~ variable, mean)
names(cst2) <- c("classid", "mthgMn", "mthkMn")
clssDt  <- merge(cst1, cst2, sort = FALSE)

###################################################
### code chunk: R 3.16
###################################################
head(clssDt)                     # First few records only

###################################################
### code chunk for Figure 3.9
###################################################

library(lattice)

xyplot(mthgMn ~ housepov, data = clssDt,   # Fig. 3.9a
       type = c("p", "smooth"),
       ylim = ylims,
       grid = TRUE)

xyplot(mthgMn ~ mthkMn, data = clssDt, 
       type = c("p", "smooth"),
       ylim = ylims,
       grid = TRUE)

###################################################
### code chunk: R 3.17a
###################################################
auxDt <- merge(SIIdata, clssDt, sort = FALSE)

auxDt2  <- 
  within(auxDt,
         {
          auxL  <- paste(classid, schoolid, sep = "\n:")
          auxL1 <- paste(auxL, clssN, sep = "\n(")    
          auxL2 <- paste(auxL1, ")", sep = "")       
          clssF <-                                 # Factor clssF created
            factor(schoolid:classid, labels = unique(auxL2))
         })
tmpDt <- subset(auxDt2, select = c(classid, schoolid, clssN, clssF))
head(tmpDt, 4)                                    # First four records
tail(tmpDt, 4)                                    # Last four records


###################################################
### code chunk: R 3.17b
###################################################
library(lattice)
dotplot(mathgain ~ clssF,                         # Fig. 3.10a
        subset(auxDt2, schoolid %in% 1:4)) 

xyplot(mathgain ~ housepov, SIIdata,              # Fig. 3.10b
       type = c("p", "smooth"))
detach(package:lattice)

##################################################
### SII Project: Cleanup
##################################################
rm(schlN, mthgM)
rm(idvars, mvars, dtm1, cst1, dtm2, cst2, schlDt)
rm(ylims, clssN, clssCnt, nms2, clssDt)
rm(auxDt, auxDt2, tmpDt)

#################################################
### FCAT Study: Target score (Section 3.6)
#################################################


###################################################
### code chunk: R 3.18a
###################################################
data(fcat, package = "nlmeU")
scM <- with(fcat, tapply(scorec, list(id, target), mean))
scM[c(1, 2, 539), ]


###################################################
### code chunk: R 3.18b
###################################################
library(lattice)
histogram(~scorec|target, 
   data= fcat, breaks = NULL             # Fig. 3.11
)
detach(package:lattice)

##################################################
##### Cleanup and sessionInfo 
##################################################
rm(scM)

sessionInfo()             # At the end of the session
