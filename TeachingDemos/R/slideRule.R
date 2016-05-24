
slideRule <- function( slide=1, rule=1 ) {


  sr.tks <- c( seq(1,2,.1), seq(2.2,3,.2), seq(3.5,10,.5), seq(11,20,1),
	seq(22,30,2), seq(35,100,5) )
  sr.tks2 <- c( 1, 2:10, seq(20,100,10) )
  sr.tl <- c( 1, 2:9, 1, 2:9, 1 )

  op <- par(plt=c(0.03, 0.97, 0.49, 0.51), xpd=TRUE )
  on.exit(par(op))
  plot.new()
  plot.window( xlim=c(0.1, 100), ylim=c(0,1), log='x' )

  axis(3, at=sr.tks, labels=FALSE, tcl=-0.3)
  axis(3, at=sr.tks2, labels=sr.tl, cex.axis=0.4 )

  axis(1, at=sr.tks/slide, labels=FALSE, tcl=-0.3)
  axis(1, at=sr.tks2/slide, labels=sr.tl, cex.axis=0.4, mgp=c(3,.5,0) )
  segments( rule, grconvertY(0.4, from='nfc'), rule, grconvertY(0.6, from='nfc'), col='blue')

 points(1, grconvertY(0.52, from='nfc', to='user'), pch=6)
 points(1/slide, grconvertY(0.48, from='nfc', to='user'), pch=2)

}



slideRule2 <- function( slide=1, rule=1 ) {


  sr.tks <- c( seq(1,2,.1), seq(2.2,3,.2), seq(3.5,10,.5), seq(11,20,1),
	seq(22,30,2), seq(35,100,5) )
  sr.tks2 <- c( 1, 2:10, seq(20,100,10) )
  sr.tl <- c( 1, 2:9, 1, 2:9, 1 )

  op <- par(plt=c(0.03, 0.97, 0.49, 0.51), xpd=TRUE )
  on.exit(par(op))
  plot.new()
  plot.window( xlim=c(0.1, 100), ylim=c(0,1), log='x' )

  axis(3, at=sr.tks, labels=FALSE, tcl=-0.5)
  axis(3, at=sr.tks2, labels=sr.tl, cex.axis=2, line=3 )

  axis(1, at=sr.tks/slide, labels=FALSE, tcl=-0.5)
  axis(1, at=sr.tks2/slide, labels=sr.tl, cex.axis=2, mgp=c(3,.5,0) )
  segments( rule, grconvertY(0.4, from='nfc'), rule, grconvertY(0.6, from='nfc'), col='blue')

 points(1, grconvertY(0.52, from='nfc', to='user'), pch=6)
 points(1/slide, grconvertY(0.48, from='nfc', to='user'), pch=2)

}




TkSlideRule <- function() {

	sl.list <- list( slide=list('slider',init=1, from=0.1, to=9.9, resolution=0.1),
				rule=list('slider',init=1, from=0.1, to=9.9, resolution=0.1))
	tkexamp( slideRule2, sl.list )
}






