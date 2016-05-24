# Spectrogram color vector functions
# based on existing palettes in in grDevices package
# originally by R Development Core Team and contributors worldwide 
gray.1 <- function(n=30) gray(seq(1, 0,length.out=n))
gray.2 <- function(n=30) gray(1-seq(0, 1,length.out=n)^2)
gray.3 <- function(n=30) gray(1-seq(0, 1,length.out=n)^3)
rainbow.1 <- function(n=15) rev(rainbow(n))
topo.1 <- function(n=12) rev(topo.colors(n))
