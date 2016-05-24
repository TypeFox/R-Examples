f.windows <- function(n.markers, markers, winlength){
##
## PRODUCE MARKER SEQUENCES (AS A MATRIX) FOR ALL POSSIBLE 
## WINDOWS OF LENGTH (AT MOST) slidelen

if(!missing(markers)) .markers <- markers
else .markers <- 1:n.markers

.n.markers <- length(.markers)

if(winlength > .n.markers) warning('"winlength" is larger than the number of available markers.\n It has been truncated to the right length.', call. = F)
.winlen <- min(.n.markers, winlength) # LENGTH OF SLIDING WINDOW
.nwin <- .n.markers - .winlen + 1 # NUMBER OF SLIDING WINDOWS OF LENGTH (AT MOST) winlen

.slides <- outer(1:.nwin, 0:(.winlen - 1), "+")# CREATE ROWS WITH WINDOW INDICES
.slides <- t(apply(.slides, 1, function(x) .markers[x])) # CREATE ROWS WITH ACTUAL WINDOWS
if(winlength == 1) .slides <- t(.slides)

return(.slides)
}
