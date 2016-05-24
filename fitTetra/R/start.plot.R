start.plot <-
function(plot.type, plot.filename, ncol=1, nrow=1, 
                       width=6.7, height=10.1) { #A4 with 2 cm margins
  switch (plot.type, 
    pdf = pdf(file    =paste(plot.filename,".pdf",sep=""), width=width, height=height),
    png = png(filename=paste(plot.filename,".png",sep=""), width=width, height=height, units="in", res=96),
    emf = emf(file    =paste(plot.filename,".emf",sep=""), width=width, height=height),
    svg = svg(filename=paste(plot.filename,".svg",sep=""), width=width, height=height)
  )      
  if (ncol>1 || nrow>1)
  layout(matrix(1:(nrow*ncol),byrow=FALSE,ncol=ncol))  #adapt to the number of models tested
}
