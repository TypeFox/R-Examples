f.Rplot <- function(lwd, ylim, L, U, len, pos, est, est.in, use, pch, ...){
##
## BASIC PLOTTING OF EFFECT VALUES AND THEIR CONFIDENCE INTERVALS
##
.use.log <- is.element(seq(along = pos), use)
#
## CHECKING AND FIXING DATA POINTS OUTSIDE ylim. THIS REMOVES BOTHERSOME WARNINGS
## THAT DIDN'T DISAPPEAR WHEN SETTING suppress.graphics.warnings OR OTHER PARAMS.
#
#
.f.segments.rest <- function(x0, y0, x1, y1, use = T, ...){
	## PLOTS ONLY STUFF WITHIN
	.above <- (y0 > ylim[2]) & (y1 > ylim[2])
	.below <- (y0 < ylim[1]) & (y1 < ylim[1])
	.use <- !.above & !.below & use
	
	y0 <- pmin(pmax(ylim[1], y0), ylim[2])
	y1 <- pmax(pmin(ylim[2], y1), ylim[1])
	
	if(sum(.use) == 0) return()
	segments(x0[.use], y0[.use], x1[.use], y1[.use], lwd = lwd, ...)
}
#
## PLOTTING RELATIVE RISKS
# WITH CIs AND LINE-ENDS
## VERTICAL BAR
.f.segments.rest(pos, L, pos, U, use = .use.log, ...)
## HORIZONTAL TOP AND BOTTOM BARS
.f.segments.rest(pos - len, U, pos + len, U, use = .use.log, ...)
.f.segments.rest(pos - len, L, pos + len, L, use = .use.log, ...)
#
## EFFECT
points(pos[est.in], est[est.in], pch = 22, cex = 2, bg = "white", col = "white") # PROVIDES WHITE BACKGROUND FOR PLOTTING CHARACTER
points(pos[est.in], est[est.in], pch = pch, font = 2) # PLOTTING CHARACTER, SUCH AS "s", "d" ETC.
}
