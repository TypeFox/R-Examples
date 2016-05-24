CEnvelope <- 
function(Profile, LineWidth = 2, ShadeColor = "grey75", BorderColor = "red", ...)
{ 
  if (!is.CommunityProfile(Profile)) 
    stop("Confidence Envelopes require a CommunityProfile.")
  if (!(is.null(Profile$high) | is.null(Profile$low))) {
  	# Shaded polygon (adapted from Didzis Elferts, 
  	# http://stackoverflow.com/questions/14069629/plotting-confidence-intervals)
    graphics::polygon(c(Profile$x, rev(Profile$x)), c(pmax(Profile$low, graphics::par('usr')[3]), pmin(rev(Profile$high), par('usr')[4])), col = ShadeColor, border = FALSE)
  	# Add red lines on borders of polygon
  	graphics::lines(Profile$x, Profile$high, col=BorderColor, lty=2)
  	graphics::lines(Profile$x, Profile$low, col=BorderColor, lty=2)
  }
  graphics::lines(Profile$x, Profile$y, lwd = LineWidth, ...)
}
