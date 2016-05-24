# Compute colors for visual assessment of prognosis
# Author : Sylvain Mareschal <maressyl@gmail.com>
surv.colors <- function(time, event, eventColors=c("#000000", "#CCCCCC"), censColors=c("#FFFFEE", "#FFDD00")) {
	out <- rep(NA, length(event))
	out[ event ] <- rgb(colorRamp(eventColors)(time[ event ] / max(time[ event ], na.rm=TRUE)), maxColorValue=255)
	out[ !event ] <- rgb(colorRamp(censColors)(time[ !event ] / max(time[ !event ], na.rm=TRUE)), maxColorValue=255)
	names(out) <- names(time)
	
	return(out)
}

