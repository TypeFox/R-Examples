approxeq <- function(x, y, tolerance = .Machine$double.eps ^ 0.5, ...) {
# unlike all.equal, this returns a vector of T/F and ignores object properties
#input validation 
if (length(x) != length(y)) warning('x,y lengths differ. Will recycle.')
#don't care about dimensions so long as you're smart about inputs
checkit <- abs(x-y) < tolerance
return(invisible(checkit))
}
