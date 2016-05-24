## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE-----------------
library(drLumi)

## ----echo=FALSE, prompt=TRUE, comment=NA, highlight=FALSE, warning=FALSE------
library(xtable)
library(minpack.lm)
library(ggplot2)
library(plyr)
options(width=80)
knitr::opts_chunk$set(echo=FALSE, fig.path='./', cache=TRUE)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
data(mfidata)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(mfidata)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
data(ecdata)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(ecdata)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
datasets <- data_selection(x = mfidata, ecfile = ecdata, 
    byvar.ecfile = c("sample","analyte"),
    backname = "Background0", 
    stanname="Standard",posname = "Controls")

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(datasets$plate_1$background,3)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(datasets$plate_1$standard,3)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(datasets$plate_1$positive,3)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(datasets$plate_1$unknowns,3)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
allanalytes <- scluminex(plateid = "newplate", 
    standard = datasets$plate_1$standard, 
    background = datasets$plate_1$background,
    bkg = "ignore", lfct = c("SSl5","SSl4"), 
    fmfi = "mfi", verbose = FALSE)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, warning=FALSE-------
class(allanalytes)
allanalytes
names(allanalytes$FGF)

## ----echo=FALSE, eval=TRUE, prompt=TRUE, highlight=FALSE----------------------
flog5p <- scluminex(plateid = "plate_1", standard = datasets[[1]]$standard, 
    background = datasets[[1]]$background,
    bkg = "ignore", lfct = "SSl5", 
    fmfi = "mfi", verbose = FALSE)
flog4p <- scluminex(plateid = "plate_1", standard = datasets[[1]]$standard, 
    background = datasets[[1]]$background,
    bkg = "ignore", lfct = "SSl4", fmfi = "mfi", 
    verbose = FALSE)
fexp <- scluminex(plateid = "plate_1", standard = datasets[[1]]$standard, 
    background = datasets[[1]]$background, 
    bkg = "ignore", lfct = "SSexp", fmfi = "mfi", 
    verbose = FALSE)

## ----echo=FALSE, eval=TRUE, prompt=FALSE, highlight=FALSE, fig.show='asis'----
q <- plot(flog5p, subset.list="IL15", size.legend=NA, psize=2.5,  
    color.bkg = NA)
q <- q + ggtitle("5 parameter logistic function") 
q <- q + theme(plot.title = element_text(size = 30))
q

## ----echo=FALSE, eval=TRUE, prompt=FALSE, highlight=FALSE, fig.show='asis'----
q <- plot(flog4p, subset.list="IL15", size.legend=NA, psize=2.5,  
    color.bkg = NA)
q <- q+ggtitle("4 parameter logistic function")  
q <- q + theme(plot.title = element_text(size = 30))
q

## ----echo=FALSE, eval=TRUE, prompt=FALSE, highlight=FALSE, fig.show='asis'----
q <- plot(fexp, subset.list="IL15", size.legend=NA, psize=2.5, 
    color.bkg = NA)
q <- q + ggtitle("Exponential growth") 
q <- q + theme(plot.title = element_text(size = 30))
q

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
fgf <- subset(mfidata, analyte=="FGF" & plate=="plate_1")
dat <- data_selection(fgf, ecdata) 

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE,warning=FALSE--------
ig <- scluminex("plate_1",dat$plate_1$standard, dat$plate_1$background, 
                lfct="SSl4", bkg="ignore", fmfi="mfi", fanalyte="analyte", 
                verbose=FALSE)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
inc <- scluminex("plate_1",dat$plate_1$standard, dat$plate_1$background, 
                lfct="SSl4", bkg="include", fmfi="mfi", fanalyte="analyte", 
                verbose=FALSE)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
sub <- scluminex("plate_1",dat$plate_1$standard, dat$plate_1$background, 
                lfct="SSl4", bkg="subtract", fmfi="mfi", fanalyte="analyte", 
                verbose=FALSE)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
cons <- scluminex("plate_1", dat$plate_1$standard, dat$plate_1$background, 
    lfct="SSl4", bkg="constraint", fmfi="mfi", fanalyte="analyte", 
    verbose=FALSE)

## ----echo=FALSE, eval=FALSE---------------------------------------------------
#  r <- plot(ig, size.legend=NA) + ggtitle("Ignored")
#  s <- plot(inc, size.legend=NA) + ggtitle("Included")
#  u <- plot(sub, size.legend=NA) + ggtitle("Subtracted")
#  v <- plot(cons, size.legend=NA) + ggtitle("Constrained")
#  multiplot(r,s,u,v, cols=2)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
# arguments of the function
args(loq_derivatives)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
der <- loq_derivatives(allanalytes, subset.list="FGF")
der

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
int <- loq_interval(allanalytes, subset.list= "FGF", low.asymp=2, high.asymp=3)
int

## ----loqint, echo=FALSE, warning=FALSE, fig.show='asis', fig.height=4, fig.width=4.5, dev='pdf', fig.align='center' ,fig.cap='Estimation of limits of quantification based on interval method for plate 1, FGF analyte and ignored background'----
par(mar=c(5,4,2,2), cex=0.8)
with(int[[1]]$data, plot(log10_concentration, log10_mfi,
    type="p", main="", axes=TRUE , ylab="MFI", 
    xlab="Concentration", ylim=c(1,3.9), xlim=c(-2,4.8)))
conf <- conf_bands(allanalytes, "FGF", 
    xvalue=seq(min(int[[1]]$data$log10_concentration),
    max(int[[1]]$data$log10_concentration),0.1), interval="prediction")
with(conf,lines(xvalue, log10_mfi))
with(conf,lines(xvalue, log10_mfi.lci, lty=2))
with(conf,lines(xvalue, log10_mfi.uci, lty=2))

model <- allanalytes$FGF$model
ss <- summary(model)
qvalue <- 1-(0.05/2)
tvalue <- qt(qvalue, df = ss$df[2])

coflow <- summary(model)$coef[2,c(1,2,3)]
abline(h=coflow[1], lty=1, col="red")
abline(h=int[[1]]$ul, lty=2, col="red")
abline(v=int[[1]]$lloq, lty=2, col="blue")

cex.lev <- 1
text(2.5, coflow[1] , pos=1, "Lower asymptote coefficient", cex=cex.lev-0.2)
text(3, int[[1]]$ul, pos=3,  "Lower asymptote coefficient upper limit", 
    cex=cex.lev-0.2)
text(int[[1]]$lloq, 2.5, pos=2, "Lower LOQ", cex=cex.lev-0.2)

# text(2, 3,pos=2, "Prediction interval", cex=cex.lev)

coflow <- summary(model)$coef[3,1]
abline(h=coflow[1], lty=1, col="red")
abline(h=int$FGF$ll, lty=2, col="red")
abline(v=int$FGF$uloq, lty=2, col="blue")

text(0, coflow[1], pos=3, "Upper asymptote coefficient", cex=cex.lev-0.2)
text(0, int$FGF$ll, pos=1, "Upper asymptote coefficient lower limit", 
    cex=cex.lev-0.2)
text(int$FGF$uloq,2.5, pos=4, "Higher LOQ", cex=cex.lev-0.2)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
int2 <- loq_interval(allanalytes, subset.list="FGF", high.asymp=3, lowci=1.5)
int2

## ----loqint2, echo=FALSE, warning=FALSE, fig.show='hide', fig.cap='Estimation of limits of quantification based on interval method (fixed lower assymptote value) for plate 1, FGF analyte and ignored background.',fig.show='asis', fig.height=4, fig.width=4.5, dev='pdf', fig.align='center'----
par(mar=c(5,4,2,2), cex=0.8)
with(int[[1]]$data, plot(log10_concentration, log10_mfi,
    type="p", main="", axes=TRUE , ylab="MFI", 
    xlab="Concentration",ylim=c(1,3.9), xlim=c(-1.8,4.4)))
conf <- conf_bands(ig,"FGF", 
    xvalue=seq(min(int[[1]]$data$log10_concentration),
    max(int[[1]]$data$log10_concentration),0.1), 
    interval="prediction")
with(conf,lines(xvalue, log10_mfi))
with(conf,lines(xvalue, log10_mfi.lci, lty=2))
with(conf,lines(xvalue, log10_mfi.uci, lty=2))

model <- allanalytes$FGF$model
ss <- summary(model)
qvalue <- 1-(0.05/2)
tvalue <- qt(qvalue, df = ss$df[2])

coflow <- summary(model)$coef[2,c(1,2,3)]

abline(h=int2[[1]]$ul, lty=2, col="red")
abline(v=int2[[1]]$lloq, lty=2, col="blue")

cex.lev <- 0.8
text(int2[[1]]$lloq, 2, pos=2, "Lower LOQ", cex=cex.lev)
text(2, int2[[1]]$ul, pos=1, "Fixed asymptote", cex=cex.lev)

# text(2, 3,pos=2, "Prediction interval", cex=cex.lev)

coflow <- summary(model)$coef[3,1]
abline(h=coflow[1], lty=1, col="red")
abline(h=int2$FGF$ll, lty=2, col="red")
abline(v=int2$FGF$uloq, lty=2, col="blue")

text(0, coflow[1], pos=3, "Upper asymptote coefficient", cex=cex.lev)
text(0.3, int$FGF$ll, pos=1, "Upper asymptote coefficient lower limit", cex=cex.lev)
text(int$FGF$uloq,2, pos=4, "Higher LOQ", cex=cex.lev)

## ----echo=TRUE,eval=TRUE,prompt=TRUE,comment=NA,highlight=FALSE---------------
cv <- loq_cv(allanalytes, subset.list="FGF", max.cv=0.2)
cv

## ----loqcv, eval=TRUE, echo=FALSE, highlight=FALSE,fig.show='asis', fig.cap='Estimation of limits of quantification based on coefficient of variation method.', fig.height=3, fig.width=3.5, dev='pdf', fig.align='center'----
par(mar=c(2.1, 2.1, 0.1, 2.1), cex=0.7)
cv <- loq_cv(allanalytes, subset.list="FGF", max.cv=0.2,n.cuts=1000)
xax <- c(min(cv[[1]]$data$log10_concentration),
    max(cv[[1]]$data$log10_concentration))
yax <- c(min(cv[[1]]$data$log10_mfi),
    max(cv[[1]]$data$log10_mfi))

conf <- conf_bands(allanalytes,"FGF", 
    xvalue=seq(xax[1], xax[2], 0.1), interval="prediction")

with(conf, plot(xvalue, log10_mfi, type="l", 
    axes=F,xlim=c(xax[1],xax[2]),xlab="",ylab=""))

mtext("Concentration",1, cex=1, line=1)
mtext("Median Fluorescence Intensity",2, cex=1, line=1)

with(cv[[1]]$data, lines(log10_concentration, log10_mfi, type="p"))
abline(v=c(cv[[1]]$lloq, cv[[1]]$uloq),lty=2)

par(new=TRUE)
with(cv[[1]]$cv, plot(log10_concentration.fit, cv, 
    type="l", axes=F, xlab="",ylab="", 
    xlim=c(xax[1],xax[2]),col="red", ylim=c(0,1)))

mtext("Coefficient of Variation", side=4, cex=1, line=1)
box()
abline(v=c(cv[[1]]$lloq, cv[[1]]$uloq), lty=2, col="blue")
abline(h=0.2,lty=2)
legend('topleft', lty=c(1,2,2), col=c("red","black"), 
    legend=c("CV value","CV specified cutoff"), cex=1.3, bg="white")

text(0, 0.65, pos=4, "Lower LOQ", cex=1.15)
text(2.3, 0.65, pos=4, "Higher LOQ", cex=1.15)


## ----eval=TRUE, echo=TRUE, highlight=FALSE, comment=NA, prompt=TRUE-----------
class(der)

## ----eval=TRUE, echo=TRUE, highlight=FALSE, comment=NA, prompt=TRUE-----------
summary(der)

## ----eval=TRUE, echo=TRUE, highlight=FALSE, comment=NA, prompt=TRUE-----------
args(invest)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
invesboot <- invest(ig, "FGF", yvalue = 1.4, ci.method="bootstrap")
invesboot

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
invesdelta <- invest(ig, "FGF", yvalue = 1.4, ci.method="delta")
invesdelta

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
concdf <- subset(datasets$plate_1$positive, analyte=="FGF")

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
est_conc(ig, concdf, fmfi="mfi", dilution=1)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
est_conc(ig, concdf, fmfi="mfi", dilution=2)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
conc_icc_df <- est_conc(allanalytes, datasets$plate_1$positive, 
    fmfi="mfi", dilution=1)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(conc_icc_df)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
icc_positive <- intra_icc(conc_icc_df, id.var=c("sample", "analyte", "plate"), 
    value.var="dil.fitted.conc", type="agreement",model="twoway", 
    unit="single")

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(icc_positive$icc.df)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
names(icc_positive$icc.mod)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
icc_positive$icc.value

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
allanalytes

## ----echo=FALSE, prompt=TRUE, comment=NA, highlight=FALSE, warning=FALSE----------------------------------------------
options(width=120)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, size="small"------------------------------------------------
summary(allanalytes)

## ----echo=FALSE, prompt=TRUE, comment=NA, highlight=FALSE, warning=FALSE------
options(width=80)

## ----echo=TRUE, prompt=TRUE, highlight=FALSE, comment=NA, size="small"--------
as.data.frame(ig)

## ----echo=TRUE, prompt=TRUE, highlight=FALSE, comment=NA----------------------
ss <- summary(ig)
as.data.frame(ss)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center',fig.width=10, fig.height= 10----
plot(allanalytes, type = "scurve", ncol=5, psize=1)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center',  fig.width= 10, fig.height= 10----
plot(allanalytes, type = "residuals", out.limit= 2.5, ncol=5, psize=1)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center', fig.width= 10, fig.height= 10----
plot(allanalytes,  type = "qqplot", ncol=5, psize=1)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center'----
get_outliers(allanalytes, out.limit=2)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center', fig.width= 10, fig.height= 10----
plot(allanalytes, "residuals", out.limit=2, size.text=2.5)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center'----
out <- get_outliers(allanalytes, out.limit=2)
flag.dat <- merge(datasets$plate_1$standard, out, by=c("analyte","well"),all.x=TRUE)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center'----
flag.allanalytes <- scluminex(plateid = "newplate.flag", 
    standard = flag.dat, 
    background = datasets$plate_1$background,
    bkg = "ignore", lfct = c("SSl5","SSl4"), 
    fmfi = "mfi", verbose = FALSE)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE, fig.align='center', fig.width= 10, fig.height=6----
plot(flag.allanalytes, "scurve", 
     subset.list=c("FGF","IL1B", "IL13","EOTAXIN","VEGF","MIG","IL4","IL8"),
     ncol=4, psize=2, size.legend=3)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
imp_path <-  system.file(c("inst","extdata"),"plate1.csv", package="drLumi")
imp <- lum_import(imp_path)
imp

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
expdb <- lum_export(imp)
expdb

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
imp$well_vars
imp$well_vars <- c("Median", "Net MFI")
exp <- lum_export(imp)
head(exp$well)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
imp_path_nozip <- system.file(c("inst","extdata"),"bead_data", 
                              package="drLumi")
bead_nozip <- lum_import(imp_path_nozip)
bead_nozip

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
head(bead_nozip$bead_files)

## ----echo=TRUE, prompt=TRUE, comment=NA, highlight=FALSE----------------------
with(bead_nozip$bead_files, table(well, batch))

