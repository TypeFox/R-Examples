"divergeZimage" <- function(ob, out=FALSE, 
                            file="divergeZimage.pdf",
                            lin = 1, title = "", center = 0,
                            x2 = vector(), x= vector(),
                            plainmat = FALSE, ylab="wavelength (nm)",
                            xlab = "time (ns)") {
  if(out)
    pdf(file)
  if(!plainmat) {
    x2 <- ob@x2
    x <- ob@x
    psi.df <- ob@psi.df
  }
  else
    psi.df <- ob
  
  xnew <- linloglines(x, mu=0,alpha=lin)
  newlab <- linlogtics(x, mu=0,alpha=lin)
  m <- par("mar")
  par(mar=c(m[1:3],6))
  st <- ceiling( (max(psi.df) - min(psi.df)) / 100 )
  br <- seq(floor(min(psi.df)), ceiling(max(psi.df)),by=st)
  zpt <- which(br == center)
  negbr <- 1:(zpt - 1)
  posbr <- zpt:length(br)
  cl <- c(sequential_hcl(length(negbr)-1, h=260),
          rev(sequential_hcl(length(posbr), h=0)))
  image(xnew, x2, breaks=br, col=cl,
        psi.df, ylab = ylab, xaxt="n",
        main = title, xlab=xlab)
  axis(1, at=newlab[,1], labels=newlab[,2])
  labs <- pretty(br)
  labsp <- vector()
  for(i in 1:length(br)) 
    labsp <- append(labsp, if(br[i] %in% pretty(br)) br[i] else NA) 
  image.plot( xnew, x2, psi.df, ylab = ylab,
             legend.only=TRUE, breaks=br, col=cl,
             lab.breaks=labsp)
  if(out)
    dev.off()
}
