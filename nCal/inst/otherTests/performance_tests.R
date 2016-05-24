set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.2)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="Test", 
    assay_id="Run 1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate),
    replicate=rep(1:n.replicate, 10))


# why does calling ncal takes much more time than calling drm?

system.time({
    out = ncal(log(fi)~expected_conc, data=dat.std, plot.se.profile=F, plot=T, auto.layout=F, return.fits=T, test.LOD=F, verbose=F, robust="median")
})

system.time({
    fit = drm(log(fi)~expected_conc, data=dat.std, robust="median", fct=LL.5())
    plot(fit)
})
    
