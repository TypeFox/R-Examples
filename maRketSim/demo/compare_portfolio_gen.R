# Comparison of different portfolio generation types
# source("~/Documents/Finances/Analysis/BondFund_vs_BondLadder/maRketSim/demo/compare_portfolio_gen.R")

cat("Meant to be called after demo/public_bond_data.R\n")

nrep=100
target.dur <- 5.2
random.constrained.durs <- rep(NA,nrep)
random.durs <- rep(NA,nrep)

for(i in seq(nrep)) {
	cat("random.constrained - ",i,"\n")
	prt.RC=genPortfolio.bond(
		n=50,
		mkt=TreasuryCurve.hm[[1]],
		dur=fund.dur,
		dur.sd=fund.dur.sd,
		name="Rolling Ladder",
		type="random.constrained"
		)
	random.constrained.durs[i] <- summary(prt.RC)$portfolio.sum$dur
	
}

for(i in seq(nrep)) {
	cat("random - ",i,"\n")
	prt.R=genPortfolio.bond(
		n=50,
		mkt=TreasuryCurve.hm[[1]],
		dur=fund.dur,
		dur.sd=fund.dur.sd,
		name="Rolling Ladder",
		type="random"
		)
	random.durs[i] <- summary(prt.R)$portfolio.sum$dur
}

# Now plot them on the same scale
lmts <- c(min(min(random.durs),min(random.constrained.durs)),max(max(random.durs),max(random.constrained.durs)))
plot(random.constrained.durs~random.durs,xlim=lmts,ylim=lmts)
cat("Constrained is much tighter but at the expense of compromising the true underlying SD\n")
sd.R <- sd(sapply(prt.R$bonds,function(x) summary(x,mkt=TreasuryCurve.hm[[1]])$dur ))
sd.RC <- sd(sapply(prt.RC$bonds,function(x) summary(x,mkt=TreasuryCurve.hm[[1]])$dur ))
cat("SD of most recent random portfolio is",round(sd.R,2),"whereas SD of most recent constrained portfolio is",round(sd.RC,2))
