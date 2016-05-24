
test.sampler <- function() {
	set.seed(1)
	data(burkina_faso_females)
	BKFem.Recon.MCMC <- popRecon.sampler(
		n.iter = 10,
		burn.in = 2,
		mean.f = burkina.faso.females$fertility.rates,
		mean.s = burkina.faso.females$survival.proportions,
		mean.g = burkina.faso.females$migration.proportions,
		mean.b = burkina.faso.females$baseline.pop.counts,
		pop.data = burkina.faso.females$census.pop.counts,
		prop.vars = burkina.faso.prop.vars
	)
	stopifnot(dim(BKFem.Recon.MCMC$fert.rate.mcmc)[[1]] == 10)
	stopifnot(dim(BKFem.Recon.MCMC$surv.prop.mcmc)[[1]] == 10)
	stopifnot(dim(BKFem.Recon.MCMC$mig.prop.mcmc)[[1]] == 10)
	stopifnot(dim(BKFem.Recon.MCMC$baseline.count.mcmc)[[1]] == 10)
	stopifnot(dim(BKFem.Recon.MCMC$lx.mcmc)[[1]] == 10)
}