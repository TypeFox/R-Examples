### R code from vignette source 'hypergeometric.Rnw'

###################################################
### code chunk number 1: set_calculate_from_scratch
###################################################
calculate_from_scratch <- FALSE


###################################################
### code chunk number 2: loadpackages
###################################################
require("hypergeo")
require("elliptic")


###################################################
### code chunk number 3: R_expectation
###################################################
expected <- function(a,b,p){
  Re(
  choose(a+b,b) * p^a * (1-p)^b * (
     p *b/(1+a) * hypergeo(a+b+1,2,a+2,  p) +
  (1-p)*a/(1+b) * hypergeo(a+b+1,2,b+2,1-p) ))
}


###################################################
### code chunk number 4: useit
###################################################
c(expected(8,2,0.8),expected(9,1,0.8))


###################################################
### code chunk number 5: hypergeo_figure_file
###################################################
png("hypergeometric_plot.png",width=800,height=800)


###################################################
### code chunk number 6: wp_figure_plot
###################################################
x <- seq(from=0,to=2,len=200)
y <- seq(from=-1,to=1,len=200)
z <- outer(x,1i*y,"+")
hz <- hypergeo(2,1/2,2/3,z)
par(pty='s')
view(x,y,hz,levels=seq(from=-4,to=4),xlab='Real',ylab='Imag')


###################################################
### code chunk number 7: wp_figure_close
###################################################
null <- dev.off()


