plot.phtMCMC <- function(x, ...) {
	# Get chain info
	StartEndThin <- attr(x$samples, "mcpar")
	
	# Trace plots
	print(
		ggplot(melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars=NULL)) +
		geom_line(aes_string(x="1:length(value)", y="value")) +
		geom_smooth(aes_string(x="1:length(value)", y="value"), method="glm") +
		#geom_hline(aes(yintercept=value), data=truth, colour="red", linetype="dashed") +
		facet_wrap(~variable, scales="free") +
		#theme_grey(base_family="serif", base_size=11) +
		ggtitle("Parameter Traces") + #, plot.title = theme_text(size=14, face="bold", family="serif")) +
		xlab("Iteration") + ylab("Parameter Value")
	)
	
	# Marginal posterior densities
	print(
		ggplot() +
		geom_density(aes_string(x="value"), melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars=NULL)) +
		#geom_vline(aes(xintercept=value), data=truth, colour="red") +
		facet_wrap(~variable, scales="free") +
		#theme_grey(base_family="serif", base_size=11) +
		ggtitle("Marginal Posterior Densities") + #, plot.title = theme_text(size=14, face="bold", family="serif")) +
		xlab("Parameter Value") + ylab("Density")
	)
}
