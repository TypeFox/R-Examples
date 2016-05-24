summarize_metrics <-
function(v, successes=NULL) {
	if (is.null(successes)) {
		mean_replicates = boot(v, function(y,i) {mean(y[i])}, 1000)
		cis = boot.ci(mean_replicates, conf=0.95, type="perc")$percent[4:5]
		lower = cis[1]
		upper = cis[2]
		return(list(mean=mean(v),median=median(v),lower=lower,upper=upper, num_obs=length(v), total=sum(v)))
	} else {
		cis = prop.test(successes,v)$conf.int
		lower = cis[1]
		upper = cis[2]
		return(list(mean=successes/v,lower=lower,upper=upper, num_obs=v, total=successes))
	}
}
