## Read in data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "orange")
o <- retistruct.read.dataset(dataset)

## Load the human annotation of tears
o <- retistruct.read.markup(o)

## Initial plot
par(mar=c(0.1,0.1,0.1,0.1))
flatplot(o)

## Reconstruct
r <- retistruct.reconstruct(o)

flatplot(r, mesh=FALSE, stitch=FALSE, markup=FALSE)

## Initial plot in 3D space
## plot.retina(p$phi, p$lambda, p$R, m$Tt, m$Rsett)
## dev.print(pdf, file="../figures/orange-outline.pdf", width=5)
## dev.print(png, file="../figures/orange-no-grid.png", width=500)
## dev.print(png, file="../figures/orange-grid.png", width=5)
