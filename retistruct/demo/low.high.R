## Load the good raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
r.low <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
r.low <- retistruct.read.markup(r.low)
## Reconstruct
r.low <- retistruct.reconstruct(r.low)

## Load the bad raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GM182-4/R-CONTRA")
r.high <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
r.high <- retistruct.read.markup(r.high)
## Reconstruct
r.high <- retistruct.reconstruct(r.high)

## Plotting
x11(width=6.83, height=6.83/2)
par(mfrow=c(2, 4))

## Good Retina
par(mar=c(2.2, 2.2, 0.5, 0.5),
    mgp=c(1.1, 0.2, 0),
    tcl=-0.15)
lvsLplot(r.low)
panlabel("A")

par(mar=c(0.5, 0.5, 0.5, 0.5))
flatplot(r.low, strain=TRUE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=FALSE, landmarks=FALSE)
panlabel("B")

flatplot(r.low, strain=FALSE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=TRUE)
panlabel("C")

projection(r.low, datapoints=FALSE, datapoint.means=FALSE, datapoint.contours=FALSE)
panlabel("D")

## Bad Retina
par(mar=c(2.2, 2.2, 0.5, 0.5))
lvsLplot(r.high)
panlabel("E")

par(mar=c(0.5, 0.5, 0.5, 0.5))
flatplot(r.high, strain=TRUE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=FALSE, landmarks=FALSE)
panlabel("F")

flatplot(r.high, strain=FALSE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=TRUE)
panlabel("G")

projection(r.high, datapoints=FALSE, datapoint.means=FALSE, datapoint.contours=FALSE)
panlabel("H")

## Printing
## dev.print(svg, file="fig2-retistruct-low-high.svg", width=6.83, height=6.83/2)

