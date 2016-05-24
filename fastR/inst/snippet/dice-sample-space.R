dice <- 
	data.frame(red=rep(1:6,each=6),blue=rep(1:6,times=6))
with(dice,dice[(red<blue | blue==5 ) & (red + blue >= 9),])
