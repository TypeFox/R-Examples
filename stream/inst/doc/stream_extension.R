### R code from vignette source 'stream_extension.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: stream_extension.Rnw:107-108
###################################################
options(width = 75, digits = 3, prompt = 'R> ', scipen = 3)


###################################################
### code chunk number 2: stream_extension.Rnw:197-198
###################################################
library("stream")


###################################################
### code chunk number 3: stream_extension.Rnw:201-216
###################################################
DSD_UniformNoise <- function(d = 2, range = NULL) {
  if(is.null(range)) range <- matrix(c(0, 1), ncol = 2, nrow = d, 
    byrow = TRUE)
  structure(list(description = "Uniform Noise Data Stream", d = d, 
    k = NA_integer_, range = range),
        class = c("DSD_UniformNoise", "DSD_R", "DSD"))
  }
  
get_points.DSD_UniformNoise <- function(x, n = 1, 
  assignment = FALSE, ...) {
    data <- as.data.frame(t(replicate(n, 
      runif(x$d, min = x$range[ , 1], max = x$range[ , 2]))))
    if(assignment) attr(data, "assignment") <- rep(NA_integer_, n)
    data
}


###################################################
### code chunk number 4: dsd_example
###################################################
stream <- DSD_UniformNoise()
stream
plot(stream, main = description(stream))


