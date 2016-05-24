"imagelegend" <-
function(xl,yt,width,nbox,bheight,bgap,col,border=NULL)
{
	x <- c(xl,xl,xl+width,xl+width)
	top <- 0
	bottom <- bheight
	y <- c(yt-bottom,yt-top,yt-top,yt-bottom)
	polygon(x,y,border=border,col=col[1])
	for (i in 2:nbox) {
		top <- top + bheight + bgap
		bottom <- top + bheight
		y <- c(yt-bottom,yt-top,yt-top,yt-bottom)
		polygon(x,y,border=border,col=col[i])
	}
}

