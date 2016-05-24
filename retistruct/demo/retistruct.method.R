library(rgl)
## Set up a 3x3 grid for plotting
par(mfrow=c(3, 3))
par(mar=c(0.5, 0.5, 0.5, 0.5))

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GM509/R-CONTRA")
o <- retistruct.read.dataset(dataset)

## Load the human annotation of tears
o <- retistruct.read.markup(o)

## Make this a left eye to help with orientatio of points
o$side="Left"

## Plot of raw data. Axes are reversed to improve comparison with
## polar plot later
flatplot(o, markup=FALSE)
mtext("A", adj=0, font=2, line=-0.9)

## Plot the annotation
flatplot(o, datapoints=FALSE, landmarks=FALSE)
mtext("B", adj=0, font=2, line=-0.9)

## Set up fixed point
o$lambda0 <- 0

## Initial triangulation (with 500 points)
n <- 500
t <- TriangulatedOutline(o, n=n)

## Stitching
s <- StitchedOutline(t)

## Plot triangulation and stitching
flatplot(s, datapoints=FALSE, landmarks=FALSE, markup=FALSE)
mtext("C", adj=0, font=2, line=-0.9)

## Triangulate again, to take into account points added by stitching
m <- TriangulatedOutline(s, n=n,
                         suppress.external.steiner=TRUE)

## Merge the points that have been stitched
m <- mergePointsEdges(m)

## Make a rough projection to a sphere
m <- projectToSphere(m)

## Plot the initial gridlines in 2D
par(mfg=c(3, 3))
flatplot(m, grid=TRUE, 
          datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
          stitch=FALSE, strain=TRUE)
mtext("Dii", adj=0, font=2, line=-0.9)

## Plot the intial projection in 3D
par(mfg=c(2, 3))
plot.new()
mtext("Di", adj=0, font=2, line=-0.9)

sphericalplot(m, strain=TRUE, datapoints=FALSE)
rgl.viewpoint(zoom=0.7)
rgl.postscript("initial-projection.svg", "svg")

## Optimise mapping - this takes a few minutes
alpha <- 8
x0 <- 0.5
r <- solveMappingCart(m, alpha=0, x0=0, nu=1,
                        dtmax=500, maxmove=1E2, tol=2e-7,
                          plot.3d=FALSE)
r <- solveMappingCart(r, alpha=alpha, x0=x0, nu=1,
                        dtmax=500, maxmove=1E2, tol=1e-6,
                        plot.3d=FALSE)
r <- optimiseMapping(r, alpha=alpha, x0=x0, nu=0.5, 
                      plot.3d=FALSE)

## Plot the final projection in 3D and on the grid
par(mfg=c(2, 2))
plot.new()
mtext("Ei", adj=0, font=2, line=-0.9)

sphericalplot(r, strain=TRUE, datapoints=FALSE)
rgl.postscript("final-projection.svg", "svg")

par(mfg=c(3, 2))
flatplot(r, grid=TRUE, 
          datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
          stitch=FALSE, strain=TRUE)
mtext("Eii", adj=0, font=2, line=-0.9)

## Infer locations of datapoints in spherical coordinates
r <- RetinalReconstructedOutline(r)
r <- ReconstructedDataset(r)
r <- RetinalReconstructedDataset(r)

## Plot data in polar coordinates and flattend retina
par(mfg=c(2, 1))
projection(r, datapoints=TRUE, landmarks=TRUE, datapoint.contours=FALSE)
mtext("Fi", adj=0, font=2, line=-0.9)

par(mfg=c(3, 1))
flatplot(r, grid=TRUE, 
          datapoints=TRUE, landmarks=TRUE, mesh=FALSE, markup=FALSE,
          stitch=FALSE)
mtext("Fii", adj=0, font=2, line=-0.9)

## Save to PDF
dev.print(pdf, file="retistruct-method.pdf", width=6.83, height=6.83)
