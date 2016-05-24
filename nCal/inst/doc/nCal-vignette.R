## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
# set global chunk options
library(knitr)
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold',dev='pdf',
warning=FALSE,dev.args=list(family="Palatino"), tidy.opts=list(blank=FALSE, width.cutoff=60))
options(formatR.arrow=TRUE,width=90)
#render_listings()

## ----loadnCal, include=FALSE, cache=FALSE, tidy=TRUE, echo=TRUE-------------------------
library(nCal)

## ----Example1data, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE----------------------
set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.2)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), 
    analyte="Test", assay_id="Run 1", sample_id=NA, well_role="Standard", 
    dilution=rep(3**(9:0), each=n.replicate))
# add unknown
dat.unk=rbind(
data.frame(fi=exp(6.75), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=1, well_role="Unknown", dilution=1)
, data.frame(fi=exp(6.70), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=2, well_role="Unknown", dilution=1)
, data.frame(fi=exp(3), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=3, well_role="Unknown", dilution=1)
, data.frame(fi=exp(4.4), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=4, well_role="Unknown", dilution=10)
)
dat=rbind(dat.std, dat.unk)

## ----Example1drm, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE, eval=TRUE, fig.width=8, fig.height=8.5, fig.cap="ncal graphical output, drm fit."----
res.drm = ncal(log(fi)~expected_conc, dat, return.fits = TRUE)

## ----Example1drmres, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE, eval=TRUE, fig.width=8, fig.height=8.5, fig.cap="ncal graphical output, drm fit."----
res.drm 

## ----Example1resultsfit, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE----------------
fit.drm=attr(res.drm, "fits")[[1]]

## ----Example1bcrm, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE, eval=TRUE, fig.width=8, fig.height=8.5, fig.cap="ncal graphical output, bcrm fit."----
res.bcrm = ncal(log(fi)~expected_conc, dat, bcrm.fit=T, return.fits = TRUE, bcrm.model="norm", control.jags=list(n.iter=5e3))
fit.bcrm=attr(res.bcrm, "fits")

## ----Example1bcrmres, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE, eval=TRUE--------
res.bcrm

## ----Example1results, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE-------------------
rbind(cla2gh(coef(fit.drm)), coef(fit.bcrm))
rbind(sqrt(diag(vcov(fit.drm))), sqrt(diag(vcov(fit.bcrm, type="classical"))))

## ----newunknown, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE------------------------
getConc(fit.bcrm, c(5.7,6.3))

## ----Example2bcrm, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE, eval=TRUE-----------

dat=subset(hier.model.ex.2, assay_id %in% paste("Run",1:4))
fit.bcrm=bcrm(log(fi)~expected_conc, dat, error.model="t4", prior="cytokine", n.iter=1e4)

## ----Example2, include=TRUE, cache=FALSE, tidy=TRUE, echo=TRUE, eval=TRUE, fig.width=8, fig.height=8.5, fig.cap="Comparing bcrm fit with drm and Prism fits."----

# parameters from Prism fits
prism.1 = c("c"=1.596,"d"=10.28,"f"=0.7202,"b"=-0.8815,"e"=10^((1.597+1/0.8815*log10(2**(1/0.7202)-1))) )
prism.2 = c("c"=1.350,"d"=11.32,"f"=8.640e+010,"b"=-0.3452,"e"=10^((1.485+1/0.3452*log10(2**(1/8.640e+010)-1))) )
prism.3 = c("c"=1.333,"d"=10.23,"f"=0.7366,"b"=-0.8502,"e"=10^((1.526+1/0.8502*log10(2**(1/0.7366)-1))) )
prism.4 = c("c"=1.580,"d"=10.37,"f"=1.694,"b"=-0.6639,"e"=10^((1.530+1/0.6639*log10(2**(1/1.694)-1))) )
prism=rbind(prism.1,prism.2,prism.3,prism.4)
# start plotting
par(mfrow=c(2,2))
for (i in 1:4) {
    assay.id=paste("Run", i)
    # fit drm model
    fit.drm = drm.fit(log(fi)~expected_conc, data=dat[dat$assay_id==assay.id,], robust="median")
    plot(fit.drm, type="all", col="black", main=assay.id, lty=2)
    plot(get.single.fit(fit.bcrm, assay.id), add=T, log="x", col=1)
    # plot Prism fit
    xx=exp(seq(log(0.51),log(1e4),length=100))
    lines(xx, FivePL.x(xx,prism[i,]), type="l", lty=1, col="darkgray")
    legend(x="bottomright",legend=c("Prism, robust","drm, median","bcrm, t4"),lty=c(1,2,1),col=c("darkgray",1,1),bty="n")
}


