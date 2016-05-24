library("RUnit")
library("nCal")

test.bcrm <- function() {

RNGkind("Mersenne-Twister", "Inversion")

# JAGS is not yet reproducible, see http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/6c8c3e6a/
tolerance.jags=0.5 
if(file.exists("D:/gDrive/3software/_checkReproducibility")) tolerance.jags=0.2
 

######################################################################################

set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
    assay_id="assay1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
dat.std$replicate=rep(1:2,10)    
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=1, well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=2, well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=3, well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=4, well_role="Unknown", dilution=10, replicate=1)
)
dat=rbind(dat.std, dat.unk)
# second plate
fi=simulate1curve (p.eotaxin[2,], rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
    assay_id="assay2", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
dat.std$replicate=rep(1:2,10)    
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=1, well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=2, well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=3, well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=4, well_role="Unknown", dilution=10, replicate=1)
)
dat=rbind(dat, dat.std, dat.unk)


# check different models
if(file.exists("D:/gDrive/3software/_checkReproducibility")) {

fits = bcrm(log(fi)~expected_conc, dat, parameterization="gh", error.model="mixnorm", prior="cytokine", n.iter=1e4)
par(mfrow=c(1,2)); plot(fits)
checkEqualsNumeric(mean(coef(fits)), 4.514953, tolerance=tolerance.jags)

fits = bcrm(log(fi)~expected_conc, dat, parameterization="gh", error.model="norm", prior="cytokine", n.iter=1e4)
par(mfrow=c(1,2)); plot(fits)
checkEqualsNumeric(mean(coef(fits)),  4.454995, tolerance=tolerance.jags)

fits = bcrm(log(fi)~expected_conc, dat, parameterization="gh", error.model="t4", prior="cytokine", n.iter=1e4, verbose=TRUE)
par(mfrow=c(1,2)); plot(fits)
checkEqualsNumeric(mean(coef(fits)),  4.512197, tolerance=tolerance.jags)

fits = bcrm(log(fi)~expected_conc, dat, parameterization="gh", error.model="lar1", prior="cytokine", n.iter=1e4, verbose=TRUE)
par(mfrow=c(1,2)); plot(fits)
checkEqualsNumeric(mean(coef(fits)),  4.549893, tolerance=tolerance.jags)

fits = bcrm(log(fi)~expected_conc, dat, parameterization="gh", error.model="mix2", prior="cytokine", n.iter=1e4, verbose=TRUE)
par(mfrow=c(1,2)); plot(fits)
checkEqualsNumeric(mean(coef(fits)),  4.454839, tolerance=tolerance.jags)

#fit.drm = drm (log(fi)~expected_conc, data=dat[dat$assay_id=="assay1" & dat$well_role=="Standard",], fct=LL.5())
#summary(fit.drm) # Residual standard error: 0.2832219

}



####################################################################################################

# decreasing curves
set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
p.1=p.eotaxin[1,]
p.1["b"]= -p.1["b"]
fi=simulate1curve (p.1, rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
    assay_id="assay1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=1, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=2, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=3, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=4, well_role="Unknown", dilution=10)
)
dat=rbind(dat.std, dat.unk)
# second plate
p.2=p.eotaxin[2,]
p.2["b"]= -p.2["b"]
fi=simulate1curve (p.2, rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
    assay_id="assay2", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=1, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=2, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=3, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=4, well_role="Unknown", dilution=10)
)
dat=rbind(dat, dat.std, dat.unk)

fits = bcrm(formula=log(fi)~expected_conc, data=dat, parameterization="gh", error.model="norm", prior="cytokine", n.iter=1e1, n.adapt=0, verbose=TRUE)

checkEqualsNumeric(
    mean(fits$median.coef)
    , 4.01789
, tolerance=tolerance.jags)




}
