#' @noRd
#' @name precintcon.plot.lorenz.ungrouped
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.plot.lorenz.ungrouped 
#' @title precintcon.plot.lorenz.ungrouped 
#' @description precintcon.plot.lorenz.ungrouped 
#' @usage precintcon.plot.lorenz.ungrouped(x, n, xlab = "SumNi(\%)", ylab = "SumPi(\%)", 
#' axis.text.color = "black", interval, fontsize) 
#' @param x a daily precipitation serie.
#' @param n an arrat if characters containing the name of the datasets.
#' @param xlab the label for the x axis.
#' @param ylab the label for the y axis.
#' @param axis.text.color .
#' @param interval the interval for calculating the Concentration Index.
#' @param fontsize the font size of the graphic.
#' @seealso \code{\link{precintcon.plot.lorenz}} 
#' @keywords precipitation plot lorenz curve
precintcon.plot.lorenz.ungrouped <- 
   function(
      x, 
      n,
		xlab            = "SumNi(%)",
		ylab            = "SumPi(%)",
		axis.text.color = "black",
		interval,
		fontsize
) {
	
	x <- as.precintcon.fd(x, interval)
	
	x <- rbind(data.frame(
		initial.class = 0, 
		final.class   = 0, 
		midpoint      = 0,
		n             = 0,
		sum.n         = 0,
		P             = 0,
		sum.P         = 0,
		p.sum.n       = 0,
		p.sum.P       = 0), x,
		data.frame(
		initial.class = 1, 
		final.class   = 1, 
		midpoint      = 1,
		n             = 1,
		sum.n         = 1,
		P             = 1,
		sum.P         = 1,
		p.sum.n       = 1,
		p.sum.P       = 1))

	x <- cbind(data.frame(dataset = paste(n)), x)

	p <- ggplot(data = x) +  
		 geom_line(data = data.frame(x = seq(0, 1, by=.1)), aes_string(x = "x", y = "x")) +
		 geom_line(aes_string(x = "p.sum.n", y = "p.sum.P")) +
		 geom_ribbon(aes_string(x = "p.sum.n", ymin = "p.sum.P", ymax = "p.sum.n"), alpha = .3) +
		 xlab(xlab) + ylab(ylab) +  facet_grid(. ~ dataset) +
		 theme(text = element_text(size=fontsize), 
				 axis.text = element_text(color = axis.text.color),
				 axis.title.x = element_text(vjust = -.5),
				 axis.title.y = element_text(angle=0, hjust=1)) +
		 scale_x_continuous(labels = percent_format()) +
		 scale_y_continuous(labels = percent_format())
	
	return(p)
}