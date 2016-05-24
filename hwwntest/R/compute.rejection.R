compute.rejection <-
function (ar = NULL, ma = NULL, npow = 100, nom.size = 0.05, 
    ndata = 1024, lapplyfn = lapply, Box.lag = 1, rand.gen=rnorm, 
    hwwn=TRUE, box=TRUE, bartlett=TRUE, d00test=TRUE,
    genwwn=TRUE, hywn=TRUE, hywavwn=TRUE, filter.number=10, family="DaubExPhase", away.from="standard", ...) 
{
    do.all.tests <- function(ndata, ar, ma, Box.lag = 1, rand.gen=rnorm,
	hwwn=TRUE, box=TRUE, bartlett=TRUE, d00test=TRUE,
	genwwn=TRUE, hywn=TRUE, hywavwn=TRUE, filter.number=10, family="DaubExPhase", away.from="standard") {
        x <- arima.sim(n = ndata, model = list(ar = ar, ma = ma),
		rand.gen=rand.gen)

	hywavwn.pval <- hywn.pval <- genwwn.pval <- d00.pval <- gnds.pval <- box.pval <- bartlett.pval <- NULL

	if (hwwn==TRUE)
		gnds.pval <- hwwn.test(x, ...)$p.value
	#gnds_station.pval <- hwwn.test2(x, ..., wavelet.type="station")$p.value
	if (box==TRUE)
		box.pval <- Box.test(x, lag = Box.lag)$p.value

	if (bartlett==TRUE)
		bartlett.pval <- bartlettB.test(x)$p.value

	if (d00test==TRUE)
		d00.pval <- d00.test(x)$p.value

	if (genwwn==TRUE)
		genwwn.pval <- genwwn.test(x, filter.number=filter.number, family=family, away.from=away.from)$p.value

	if (hywn==TRUE)
		hywn.pval <- hywn.test(x)$p.value

	if (hywavwn==TRUE)
		hywavwn.pval <- hywavwn.test(x)$p.value


        answer <- list(gnds.pval = gnds.pval, box.pval = box.pval,
		bartlett.pval=bartlett.pval, 
		d00.pval=d00.pval, genwwn.pval=genwwn.pval, hywn.pval=hywn.pval,
		hywavwn.pval=hywavwn.pval)
        return(answer)
    }

    n.diff.tests <- 0
    if (hwwn==TRUE)
	n.diff.tests <- n.diff.tests + 1
    if (box==TRUE)
	n.diff.tests <- n.diff.tests + 1
    if (bartlett==TRUE)
	n.diff.tests <- n.diff.tests + 1
    if (d00test==TRUE)
	n.diff.tests <- n.diff.tests + 1
    if (genwwn==TRUE)
	n.diff.tests <- n.diff.tests + 1
    if (hywn==TRUE)
	n.diff.tests <- n.diff.tests + 1
    if (hywavwn==TRUE)
	n.diff.tests <- n.diff.tests + 1

    v.in <- as.vector(rep(ndata, npow))
    v.out <- lapplyfn(v.in, do.all.tests, ar = ar, ma = ma, Box.lag = Box.lag,
	rand.gen=rand.gen, hwwn=hwwn, box=box, bartlett=bartlett, d00test=d00test, genwwn=genwwn, hywn=hywn, hywavwn=hywavwn, filter.number=filter.number, family=family,
	away.from=away.from)
    v.out <- matrix(unlist(v.out), ncol = n.diff.tests, byrow = TRUE)

    hywavwn.pow <- hywn.pow <- genwwn.pow <- d00.pow <- gnds.pow <- box.pow <- bartlett.pow <- NULL

    if (hwwn==TRUE)	{
	    gnds.pow <- sum(v.out[, 1] < nom.size)/npow
	    v.out <- matrix(v.out[,-1], ncol=ncol(v.out)-1)
	    }
    if (box==TRUE)	{
	    box.pow <- sum(v.out[, 1] < nom.size)/npow
	    v.out <- matrix(v.out[,-1], ncol=ncol(v.out)-1)
	    }
    if (bartlett==TRUE)	{
	    bartlett.pow <- sum(v.out[, 1] < nom.size)/npow
	    v.out <- matrix(v.out[,-1], ncol=ncol(v.out)-1)
	    }
    if (d00test==TRUE)	{
	    d00.pow <- sum(v.out[, 1] < nom.size)/npow
	    v.out <- matrix(v.out[,-1], ncol=ncol(v.out)-1)
	    }
    if (genwwn==TRUE)	{
	    genwwn.pow <- sum(v.out[, 1] < nom.size)/npow
	    v.out <- matrix(v.out[,-1], ncol=ncol(v.out)-1)
	    }
    if (hywn==TRUE)	{
	    hywn.pow <- sum(v.out[, 1] < nom.size)/npow
	    v.out <- matrix(v.out[,-1], ncol=ncol(v.out)-1)
	    }
    if (hywavwn==TRUE)	{
	    hywavwn.pow <- sum(v.out[, 1] < nom.size)/npow
	    }

    ll <- list(hwwntest.rejprop = gnds.pow, box.rejprop = box.pow,
    	bartlett.rejprop = bartlett.pow,
	d00.pow=d00.pow, genwwn.pow=genwwn.pow, hywn.pow=hywn.pow,
	hywavwn.pow=hywavwn.pow)
    return(ll)
}
