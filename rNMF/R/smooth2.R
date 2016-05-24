
## Smooth a cell with the average of nearby non-outliers.

##  "smooth2" is used in "smoothing2". It takes input i, the index of the row in image.f that needs to be smoothed. It returns the mean of the closest non-outlier pixels. Searching begins with a disk with radius one. The radius increases if the search fails. image.f is a data frame that contains x coordinates, y coordinates, pixel values, logic flags "outs".
##  NOTE: The function needs to be improved by changing the search region from a disk to a square.
 
##  @param i The row number.
##  @param image.f a data.frame. 
##  @param p The number of rows of the image matrix.
##  @param n The number of columns of the image matrix.
 
##  @return The value that is filled into the missing cell.

##  @examples
 ## Load a clean single simulated tumor image.

smooth2 = function(i, image.f, p, n){
    r = 1
    x0 = (i-1) %% p + 1
    y0 = (i-1) %/% p + 1
    neibs = which((image.f$x - x0)^2 + (image.f$y - y0)^2 <= 1)
    neibs.index = which(image.f[neibs,"outs"] == FALSE)
    while(length(neibs.index) == 0){
        r = r + 1
        neibs = which((image.f$x - x0)^2 + (image.f$y - y0)^2 <= r^2)
        neibs.index = which(image.f[neibs,"outs"] == FALSE)
    }
    neibs = neibs[neibs.index]
    return(mean(image.f[neibs,"value"]))
}
