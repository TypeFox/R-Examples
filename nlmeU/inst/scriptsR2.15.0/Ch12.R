###################################################
### code chunk: Chap12init
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("nlme")
packageVersion("lattice")
packageVersion("reshape")
packageVersion("plyr")

sessionInfo()
require(nlme)    
require(lattice) 


data(armd, package = "nlmeU")

## From Chapter 6
lm1.form <-                   # Fixed effects formula:(6.1)
    formula(visual ~ -1 + visual0 + time.f + treat.f:time.f )

## From Chapter 9
fm9.2 <- gls(lm1.form,                   # R9.2a
       weights = varPower(form = ~time),
       data = armd)                      # varStruct from <delta> group 

##  For plots
xlims <- c("4", "12","24","52wks")
ylims <- c(-4.9, 4.9)

###################################################
### code chunk: R12.1a
###################################################
(Vg1 <- Variogram(fm9.2, form = ~ time | subject))
plot(Vg1, smooth = FALSE, xlab = "Time difference")     # Fig. 12.1a


###################################################
### code chunk: R12.1b
###################################################
(Vg2 <- Variogram(fm9.2, form = ~tp | subject))
plot(Vg2, smooth = FALSE, xlab = "Time Lag")            # Fig. 12.1b



###################################################
### code chunk: R12.2a
###################################################
lm1.form <-                                  # (12.1). See R6.1
   formula(visual ~ -1 + visual0 + time.f + treat.f:time.f )
fm12.1 <-                                    # M12.1
   gls(lm1.form, weights = varPower(form = ~time), 
       correlation = corCompSymm(form = ~1|subject),
       data = armd)

###################################################
### code chunk: R12.2b
###################################################
intervals(fm12.1, which = "var-cov")   # CIs for rho (12.4), delta (12.3), sigma 

###################################################
### code chunk: R12.3a
###################################################
fm12.1vcov <-                          # Ri 
   getVarCov(fm12.1, individual = "2")
nms <- c("4wks", "12wks", "24wks", "52wks")
dnms <- list(nms, nms)                 # Dimnames created
dimnames(fm12.1vcov) <- dnms           # Dimnames assigned
print(fm12.1vcov)
print(cov2cor(fm12.1vcov),             # Ci: (12.4)
      corr = TRUE, stdevs = FALSE) 


###################################################
### code chunk: R12.3b
###################################################
anova(fm9.2, fm12.1)                   # M9.2 nested in M12.1


###################################################
### code chunk: R12.4a
###################################################
fm12.2 <-                            # M12.2 nested in M9.2
   update(fm9.2,                     # (12.5)
          correlation = corAR1(form = ~tp|subject), 
          data = armd)


###################################################
### code chunk: R12.4b
###################################################
intervals(fm12.2, which = "var-cov") # CIs for rho (12.5), delta (12.3), sigma


###################################################
### code chunk: R12.5a
###################################################
fm12.2vcov <- 
   getVarCov(fm12.2, individual = "2")
dimnames(fm12.2vcov) <- dnms
fm12.2vcov                           # Ri matrix 
fm12.2cor <- cov2cor(fm12.2vcov)
print(fm12.2cor, digits = 2,         # Ci:(12.5)
      corr = TRUE, stdevs = FALSE)          


###################################################
### code chunk: R12.5b
###################################################
anova(fm12.1, fm12.2)                # M12.1 versus M12.2


###################################################
### code chunk: R12.6a
###################################################
fm12.3 <-                       # M12.3 from M12.2 
   update(fm12.2, correlation = corSymm(form = ~tp|subject),
          data = armd)


###################################################
### code chunk: R12.6b
###################################################
intervals(fm12.3,               # 95% CIs for rho:(12.6), delta:(12.3), sigma
          which = "var-cov")


###################################################
### code chunk: R12.7
###################################################
fm12.3vcov <-                                    # Ri  
   getVarCov(fm12.3, individual = "2")
dimnames(fm12.3vcov) <- dnms
fm12.3vcov
fm12.3cor <- cov2cor(fm12.3vcov)                 # Ci:(12.6)
print(fm12.3cor, corr = TRUE, stdevs = FALSE)


###################################################
### code chunk: R12.8a
###################################################
anova(fm12.2, fm12.3)              # M12.2 nested in M12.3


###################################################
### code chunk: R12.8b
###################################################
fmA.vc   <-                        # Alternative model
   update(fm12.3, weights = varIdent(form = ~1|time.f))
anova(fm12.3, fmA.vc)              # M12.3 nested in alternative model


###################################################
### code chunk: 12.9a
###################################################
 panel.bwxplot0 <- 
   function(x,y, subscripts, ...)
           {        
            panel.grid(h = -1, v = 0)
            panel.stripplot(x, y, col = "grey", ...)
            panel.bwplot(x, y, pch = "|", ...)
           }
 bwplot(resid(fm12.3) ~ time.f | treat.f,              # Fig. 12.2
        panel = panel.bwxplot0, 
        ylab = "Residuals", data = armd)

###################################################
### code chunk: 12.9b
###################################################
 plot(fm12.3)                                        # Fig. 12.3a
 plot(fm12.3,                                        # Fig. 12.3b
      resid(., type = "p") ~ fitted(.) | time.f)

###################################################
### code chunk:  Fig. 12.3  More detailed syntax
###################################################
 xyplot(resid(fm12.3, type = "p") ~ fitted(fm12.3),           # Fig. 12.3a
   grid = TRUE, ylim = c(-4.5, 4.5), 
   xlab = "", ylab = ""                 # To save some space
 )   
 
 xyplot(resid(fm12.3, type = "p") ~ fitted(fm12.3)| time.f,   # Fig. 12.3b
   grid = TRUE, ylim = c(-4.5, 4.5), 
   xlab = "", ylab = "", 
   data = armd
 )   

###################################################
### code chunk: R12.9b (continued)
###################################################
 stdres.plot <- 
   plot(fm12.3,  resid(., type = "p") ~ jitter(time) | treat.f, 
        id = 0.01, adj = c(-0.3, 0.5 ), grid = FALSE)
 plot(update(stdres.plot,                            # Fig. 12.4
             grid = "h"))


###################################################
### code chunk for Fig. 12.5
###################################################
##### Create auxiliary data

## id for Figure 12.5 (connecting lines for outlying subjects at time = 4)
idq <- 0.02               
id <- armd$subject

residP <- resid(fm12.3, type = "pearson")  # Pearson residuals


attach(armd)

##  Create uid and ix vectors

idx1 <- tp == 1                                    # time = 4 wks 
idx <- (abs(residP) > -qnorm(idq/2)) & idx1        # Logical vector
outliers.idx <- data.frame(subject, time, treat.f, visual, residP, idx)
outliers <- subset(outliers.idx, idx, select = -idx)
nrow(outliers)                                     # Number of outliers
uid <- unique(outliers$subject)
length(uid)                                        # Number of selected subjects
uid
detach(armd)

gin <-  rep(FALSE, length(armd$subject))
gin[id %in% uid] <- TRUE

dt <- data.frame(armd, gin=gin, resid.p = residP)
dtGin <- dt[gin, ]

###################################################
### code chunk for Fig. 12.5
###################################################

myPanel <- function(x, y, subscripts, groups, ...) {
            panel.grid(h = -1, v = 0) 
            gin <- dt$gin
            gins <- gin[subscripts]
            panel.xyplot(x, y)   # All points
            x1 <- x[gins]
            y1 <- y[gins]
            subs1 <- subscripts[gins]
            panel.superpose(x1, y1, subs1, groups, type = "l", lty = "13")
        }

xyplot(resid.p ~ time.f | treat.f, data = dt,
        panel = myPanel,
        subscripts = TRUE,
        groups = subject,
        scales = list(abbreviate = TRUE),
        aspect = 1,
        xlab = "Standardized residuals",
        ylim = ylims)


###################################################
### code chunk for Fig. 12.5 Alternative syntax
###################################################

myPanel <- function(x, y, subscripts, ...) {
           panel.grid(h = -1, v = 0) 
           panel.stripplot(x, y, ...)   
           grps <- dt$subject
           gin <- dt$gin
           gins <- gin[subscripts]
           x1 <- x[gins]
           y1 <- y[gins]
           subs1 <- subscripts[gins]
           panel.superpose(x1, y1, subs1, grps, type = "l", lty = "13")
}

bw.object <- bwplot(resid.p ~ time.f | treat.f, 
        data = dt,
        panel = myPanel,
        ylab = "Standardized residuals", 
        aspect = 1.2 )
update(bw.object, xlim=xlims,
        ylim = ylims )  


###################################################
### code chunk: R12.9c
###################################################
panel.bwxplot <- function(x,y, subscripts, ...){
 panel.grid(h = -1, v = 0)
 bwstats <- tapply(y, x, boxplot.stats)
 outy <-  unlist(lapply(bwstats, FUN = function(el) el$out))
 idx0 <- 1:length(y)
 idx <- y %in% outy
 idx1 <- idx0[idx]
 panel.stripplot(x[-idx1], y[-idx1], jitter.data = TRUE, col = "grey", ...)
 panel.bwplot(x, y, pch = "|", ...)
}

bwplot(                                              # Fig. 12.7
    resid(fm12.3, type = "n") ~ time.f | treat.f,                                 
    panel = panel.bwxplot,                           # User defined panel
    ylab = "Normalized residuals",
    data = armd)                 
qqnorm(fm12.3, ~resid(., type= "n") | time.f)        # Fig. 12.8



###################################################
### code chunk for Fig. 12.6 (splom)
###################################################
r1p <- resid(fm12.3, type = "pearson")
r1n <- resid(fm12.3, type = "normalized")
dtP <- data.frame(armd, r1p, r1n)
library(reshape)
dtPm <- melt(dtP,
     measure.var = c("r1p","r1n"),
     id.var = c("subject","time.f"))
dtPc <- cast(dtPm, subject*variable ~ time.f) #
dtPc <- data.frame(dtPc) 
names(dtPc) <- c("subject","var","P4wks","P12wks","P24wks","P52wks")
range(dtPc$P4wks, na.rm = TRUE)  

myFunL <- function(corx) { 
      ltext(-2.2, 3.2, substitute(paste(rho, corx), list(corx = corx)), cex = 1)
}

myFunU <- function(corx) { 
      ltext(-2.2,-3.9, substitute(paste(rho,corx),list(corx = corx)), cex = 1)
}

my.upperPanel <-   ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
  vr <- dtPc$var == "r1n" 
  subs <- subscripts[vr]         
  x1 <- x[subs]
  y1 <- y[subs]
  panel.abline(h = c(-4, -2, 0, 2, 4), col = "grey", ...)
  panel.abline(v = c(-4, -2, 0, 2, 4), col = "grey", ...)
  panel.xyplot(x1, y1, ...)
  corx <- round(cor(x1, y1, use = "complete.obs"), 2)
  abs.corx <- abs(corx)
  corx <- paste("=", corx, sep = "")
  myFunU(corx)
}

my.lowerPanel <-    ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
  vr <- dtPc$var == "r1p" 
  subs <- subscripts[vr]         
  x1 <- x[subs]
  y1 <- y[subs]
  panel.abline(h = c(-4, -2, 0, 2, 4), col = "grey", ...)
  panel.abline(v = c(-4, -2, 0, 2, 4), col = "grey", ...)
  panel.xyplot(x1, y1, ...)
  corx <- round(cor(x1, y1, use = "complete.obs"), 2)
  abs.corx <- abs(corx)
  corx <- paste("=", corx, sep = "")
  print(corx)
  cex.value <- 2
  myFunL(corx)
}



mySuperPanel <- function(z, subscripts, panel.subscripts, ...){
  panel.pairs(z, subscripts = subscripts,
              panel.subscripts = panel.subscripts,
              as.matrix = TRUE, 
              upper.panel = "my.upperPanel",
              lower.panel = "my.lowerPanel",
        ## simpler syntax: prepanel.limits = function(z) return(c(-4.7,4.7))
              pscales =list(
                P4wks =list(limits=c(-4.7,4.7)),
                P12wks=list(limits=c(-4.7,4.7)),
                P24wks=list(limits=c(-4.7,4.7)),
                P52wks=list(limits=c(-4.7,4.7))  ) )
print(names(z))
}

abbrev.names <- c("vis0", "vis4", "vis12", "vis24", "vis52")
splom.form <- formula(~cbind(P4wks,P12wks,P24wks,P52wks))
splom(splom.form,
  data = dtPc,   #### subset(armd240,miss.pat =="----"),   
  as.matrix = TRUE,  #### varnames = abbrev.names, 
  xlab = "",
  superpanel = mySuperPanel 
)



###################################################
### code chunk for  Fig. 12.8 (alternative syntax)
###################################################
qqPlot <- qqnorm(fm12.3, ~ resid(., type = "n") | time.f)
update(qqPlot, grid = TRUE, xlim = ylims, ylim = ylims, 
       aspect = 1)

###################################################
### code chunk: R12.10
###################################################
anova(update(fm12.3, method = "ML"))                # M12.3



###################################################
### code chunk: R12.11a
###################################################
lm1a.form <- formula (visual ~ visual0 + time.f         # (12.7)
    + treat.f + time.f:treat.f)   
fm12.3a <- update(fm12.3, lm1a.form,          # M12.3a <-  M12.3
    method = "ML", data = armd)
lm2.form <- formula(visual ~ visual0 + time             # (12.8)
    + treat.f + treat.f:time)
fm12.4 <- update(fm12.3, lm2.form,            # M12.4  <- M12.3
    method = "ML", data = armd)
lm3.form <-  update(lm2.form, . ~ . - treat.f:time)     # (12.9)
fm12.5 <- update(fm12.3, lm3.form,            # M12.5 <- M12.3
    method = "ML", data = armd)


###################################################
### code chunk: R12.11b
###################################################
anova(fm12.3a, fm12.4, fm12.5)          # M12.3a,  M12.4, M12.5


###################################################
### code chunk: R12.11c
###################################################
anova(fm12.5) 


###################################################
### code chunk: R12.12
###################################################
fm12.5vcov <- getVarCov(fm12.5,          # R_i
     individual="2")
dimnames(fm12.5vcov) <- dnms             # Dimnames assigned
fm12.5vcov 
fm12.5cor <- cov2cor(fm12.5vcov)         # C_i
print(fm12.5cor, corr=TRUE, stdevs=FALSE)

#### sessionInfo  ####

sessionInfo()                 # with packages attached 
detach(package:nlme)      
detach(package:lattice)
detach(package:reshape)
