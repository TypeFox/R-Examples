"P.mountscene04" <-
function (filename="mountscene04.pdf") 
{
        if(length(filename)) {
                pdf(file=filename, height=3, width=5)
                par(mar=c(0,0,0,0) + .1, cex.axis=.7, cex.lab=.7)
        } 

	mountscene(seed=11, df=c(seq(4, 40, length=4), 30, 20, 20), levels=7,
		color=c('gray90', 'gray64', 'gray69', 'gray74', 'gray79',
		'gray84', 'gray89'))
	points(runif(500), runif(500), pch=".", col="white", cex=2)

        if(length(filename)) dev.off()
}

