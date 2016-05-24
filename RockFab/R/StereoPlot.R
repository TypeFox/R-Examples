StereoPlot <-
function(my.title = "Stereonet", new = TRUE, pdf.file){
  if(new == TRUE & missing(pdf.file)){
    dev.new(width = 3, height = 3.75, family = 'serif')
  }
  if(missing(pdf.file) == FALSE){
    pdf(file = pdf.file, width = 3, height = 3.75, family = 'serif', useDingbats = FALSE)
  }
  par(mai = c(0, 0, 0, 0), omi = c(0, 0, .5, 0))
  plot(0, 0, pch = '', asp = 1, ann = FALSE, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE)
  lines(c(0, 0), c(1, 1.02), lwd = .5)
  text(0, 1.025, "N", adj = c(.5, 0), cex = .75)
  mtext(my.title, cex = 1.25)
}
