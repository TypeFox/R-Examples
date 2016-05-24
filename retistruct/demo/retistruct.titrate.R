## Load the good raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
r.good <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
r.good <- retistruct.read.markup(r.good)
## Reconstruct
r.good <- retistruct.reconstruct(r.good)

## Titrate to give impression of whether reconstruction is good or bad
t.good <- titrate.reconstructedOutline(r.good)

par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(1.7, 0.5, 0), tcl=-0.2)
with(t.good$dat, plot(phi0d, sqrt.E, type="l", xlab=expression(paste(italic(phi)[0], " (degrees)")), ylab=expression(sqrt(italic(E)[L]))))
with(t.good,
     points(phi0d.orig, dat[which(dat$phi0d==phi0d.orig),"sqrt.E"], col="red"))

## Printing
## dev.copy2pdf(file="retistruct-good-bad.pdf", width=6.83, height=6.83/2)


