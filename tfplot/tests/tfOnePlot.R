  require("tfplot")

 Sys.info()

  #x11()
#  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
#             # width=6, height=8, pointsize=10,

  dat <- tframed(matrix(rnorm(300),100,3), list(start=c(1961,1), frequency=12))
  if(dev.cur() == 1) postscript(file=
     tempfile("GraphicsTest2", tmpdir = tempdir(), fileext = "ps"))
 
# plot(dat)
  tfOnePlot(dat)
  seriesNames(dat) <- c("newname 1", "newname 2", "newname 3")
  tfOnePlot(dat)

 
 unlink("Rplots.pdf")
