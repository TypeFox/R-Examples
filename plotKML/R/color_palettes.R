# Purpose        : Color palettes for visualization of numeric (continuous, binary, ordinary) and categorical variables;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); Pierre Roudier (pierre.roudier@landcare.nz); 
# Status         : pre-alpha
# Note           : a gallery of color palettes is available at [http://plotkml.r-forge.r-project.org/];


# Display default palettes:
display.pal <- function(pal, sel=1:length(pal), names=FALSE) {
  
  if(length(pal)>10) { sel <- 1:10 }
  
  if(names==FALSE){ 
  dev.new(width=2.1, height=length(sel))
	## not needed	
	# fin=c(2.1, length(sel)*.9), 
	par(mfrow=c(length(sel),1), mar=c(1.5,.8,1.5,.5))
	# plot palettes above each other:
	for(j in sel){
	  plot(y=rep(1, length(pal[[j]])), x=1:length(pal[[j]]), axes=FALSE, xlab='', ylab='', pch=15, cex=1.5, col=pal[[j]])
	  mtext(names(pal)[j], cex=.5, side=3)
	}
  } # names == TRUE:
  else {
	sel <- sel[1] # take only the first pallette from the list
	pal.name <- names(pal)[sel]
	pal <- pal[[sel]]
	
	# used to compute plotting region, not figure size
	leg.width <- (max(nchar(names(pal)))*20+150)/100
	leg.height <- length(pal)
	
	par(mar=c(.5,0,1.5,1))
	# plot palette and class names:
	plot(x=rep(1, length(pal)), y=1:length(pal), axes=FALSE, xlab='', ylab='', pch=15, cex=1.5, col=pal, xlim=c(0,.6*leg.width), asp=.6)
	text(x=rep(1, length(pal)), y=1:length(pal), labels=names(pal), cex=.5, pos=4, offset=1)
	mtext(pal.name, cex=.8, side=3)  
  }
}

# end of script;