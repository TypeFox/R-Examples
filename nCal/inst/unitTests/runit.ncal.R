library("RUnit")
library("nCal")

test.ncal <- function() {

RNGkind("Mersenne-Twister", "Inversion")
#RNGkind("Marsaglia-Multicarry", "Kinderman-Ramage") 
tolerance.jags=1e-2 # JAGS is not yet reproducible, see http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/6c8c3e6a/
tolerance=1e-1 # drm::drc is not reproducibile across platforms at higher tolerance
# more stringent tolerance for one system to ensure algorithm accuracy
if(file.exists("D:/gDrive/3software/_checkReproducibility")) {
    tolerance.jags=1e-2
    tolerance=1e-6
}


# a dataset without analyte and assay_id column
set.seed(1)
#    print(runif(1))
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.2)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)))

checkException (
    ncal(log(fi)~expected_conc, dat.std, return.fits = TRUE)
, msg="check no analyte or assay_id column")

# add analyte and assay_id column
dat.std$analyte="Test"
dat.std$assay_id="Run 1"

checkEquals (
    nrow(ncal(log(fi)~expected_conc, dat.std, return.fits = TRUE))
, 0, msg="check no analyte or assay_id column")

# add unknown
dat.unk=rbind(
  data.frame(fi=exp(6.75), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=1)
, data.frame(fi=exp(11), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=2)
, data.frame(fi=exp(3),    expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=3)
, data.frame(fi=exp(4.4),  expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=4)
, data.frame(fi=36000,  expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=5)
)
dat.std=cbind(dat.std, sample_id=NA)
dat=rbind(dat.std, dat.unk)

checkException (
    ncal(log(fi)~expected_conc, dat, return.fits = TRUE)
, msg="check no well_role")

# add well_role
dat.std=cbind(dat.std, well_role="Standard")
dat.unk=cbind(dat.unk, well_role="Unknown")
dat=rbind(dat.std, dat.unk)

out=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, additional.plot.func=function() abline(v=10), check.out.of.range=2, verbose=T)

# test two modes of check.out.of.range
checkEqualsNumeric(
    unlist(out[5,c("est.log.conc","se")])
    , c(12.1163415, 15.4486665), tolerance=tolerance)

out=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, additional.plot.func=function() abline(v=10), check.out.of.range=1)

checkEqualsNumeric(
    unlist(out[1,c("est.log.conc","se")])
    , c(3.9388941,    0.1627691), tolerance=tolerance)
checkEqualsNumeric(
    unlist(out[5,c("est.log.conc","se")])
    , c(9.2103404, Inf), tolerance=tolerance)


out.norm=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, bcrm.fit=TRUE, bcrm.model="norm", control.jags=list(n.iter=10, n.adapt=0))

checkEqualsNumeric(
    unlist(out.norm[1,c("est.log.conc","se")])
    , 
    c(3.9400136, 0.1253771)
, tolerance=tolerance.jags)

# weighting
out.w = ncal(fi~expected_conc, dat, return.fits = TRUE, plot.se.profile=TRUE, var.model="power", control.crm.fit=list(max.iter=2), verbose=T)
fit.w=attr(out.w, "fits")[[1]]

checkEqualsNumeric(
    coef(fit.w)
    , 
    c(-1.207390,    76.644395, 31192.072130,   700.299523,     1.148918 )
, tolerance=tolerance)



# test decreasing curve
dat.2 = dat
dat.2$expected_conc=1/dat.2$expected_conc

out.2=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE)

checkTrue(
    coef(attr(out.2, "fits")[[1]])["b"]>0
)

checkEqualsNumeric(
    unlist(out.2[1:3,c("est.log.conc","se")])
    , 
    c(-3.9391352,    -9.9034876,     0.6771702,     0.1624440, Inf, Inf)
, tolerance=tolerance)


# commented out b/c jags fails to fit here
#out.2.norm=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE, bcrm.fit=TRUE, control.jags=list(n.iter=1e1, n.adapt=0), bcrm.model="norm", verbose=FALSE)
#
#checkEqualsNumeric(
#    coef(attr(out.2.norm, "fits"))
#    ,
#    c(4.386389, 10.511270, -4.145641, -1.322411,  2.547013)
#, tolerance=tolerance.jags)
#
#checkEqualsNumeric(
#    unlist(out.2.norm[1:3,c("est.log.conc","se")])
#    , 
#    c(-3.9400136,    -9.9034876,     0.6771702,     0.1253771, Inf, Inf)
#, tolerance=tolerance.jags)


# test 4PL
out.4pl=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, fit.4pl=TRUE, verbose=2)
    
checkTrue(
    is.na(coef(attr(out.4pl, "fits")[[1]])["f"])
)

checkEqualsNumeric(
    unlist(out.4pl[1:3,c("est.log.conc","se")])
    , 
    c(3.9829620, 9.2103404, -1.3703174, 0.1629386, Inf, Inf)
, tolerance=tolerance)


out.4pl.norm=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, fit.4pl=TRUE, bcrm.fit=TRUE, bcrm.model="norm", control.jags=list(n.iter=10, n.adapt=0))

checkEqualsNumeric(
    coef(attr(out.4pl.norm, "fits"))
    , 
    c(4.386389, 10.511270,  4.145641,  1.322411)
, tolerance=tolerance.jags)

# fails to fit
#out.2.4pl.norm=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE, fit.4pl=TRUE, bcrm.fit=TRUE, bcrm.model="norm", control.jags=list(n.iter=10, n.adapt=0))
#
#checkEqualsNumeric(
#    coef(attr(out.2.4pl.norm, "fits"))
#    , 
#    c(4.386389, 10.511270, -4.145641, -1.322411)
#, tolerance=tolerance.jags)


out.2.4pl=ncal(log(fi)~expected_conc, dat.2, return.fits = TRUE, fit.4pl=TRUE)

checkEqualsNumeric(
    coef(attr(out.2.4pl, "fits")[[1]])
    , 
    c(0.86770147,  4.25401765, 10.38578190,  0.01207711)
, tolerance=tolerance)


out.4pl.t4=ncal(log(fi)~expected_conc, dat, return.fits = TRUE, fit.4pl=TRUE, bcrm.fit=TRUE, bcrm.model="t4", control.jags=list(n.iter=10, n.adapt=0))

checkEqualsNumeric(
    coef(attr(out.4pl.t4, "fits"))
    , 
    c(4.386389, 10.511270,  4.145641,  1.322411)
, tolerance=tolerance.jags)


## test read.luminex.xls

dat = read.luminex.xls(paste(system.file(package="nCal")[1],
    '/misc/02-14A22-IgA-Biotin-tiny.xls', sep=""), verbose=FALSE)
out = ncal(log(fi)~expected_conc, dat, return.fits = TRUE, plot.se.profile=FALSE)

checkEqualsNumeric(unlist(out[1:3,c("fi")]), c(15,45,19.33908), tolerance=tolerance)
checkEqualsNumeric(unlist(out[12,c("fi")]), c(183.73622), tolerance=tolerance)


## test weighted LS fit

# GLS-PL fit
fit = crm.fit(formula=fi ~ expected_conc, data=dat.QIL3[dat.QIL3$assay_id=="LMX001",], var.model="power", max.iter=3, verbose=2)
plot(fit, log="xy", type="all")
checkEqualsNumeric(coef(fit), c(-0.6336108,  1121.5930866, 22420.6091038,   497.6121715,     1.5513449), tolerance=tolerance)
checkEqualsNumeric(deviance(fit), 250.9859, tolerance=tolerance)

fit.1=gnls.fit(formula=fi ~ expected_conc, data=dat.QIL3[dat.QIL3$assay_id=="LMX001",], verbose=1) # seems it does not matter optim or nlminb is used
lines5PL(coef(fit.1), xlim=c(.5,1e4), col=2)
checkEqualsNumeric(coef(fit.1), c(1195.540943 ,26798.099196     ,7.059411  ,3589.136932     ,5.384301), tolerance=tolerance)
checkEqualsNumeric(-2*logLik(fit.1), 286.4975, tolerance=tolerance)

# mle fit
fit.2 = crm.fit(formula=fi ~ expected_conc, data=dat.QIL3[dat.QIL3$assay_id=="LMX001",], var.model="power", method="mle", max.iter=10, verbose=2)
lines5PL(coef(fit.2), xlim=c(.5, 1e4), col=2)
checkEqualsNumeric(coef(fit.2), c(1156.076841 ,21883.137557     ,6.739041  ,3705.174635     ,1.868717), tolerance=tolerance)

}
