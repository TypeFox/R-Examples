library(rgam)
library(gam)
set.seed(123)
# jiggle data to break ties
x <- ili.visits$week
y <- ili.visits$visits
x <- x + rnorm(length(x), mean=0, sd=.01)

# robust and classical back-fits
#
# the following command takes about 70 mins on an Intel Xeon CPU (3.2GHz)
#
# a <- rgam(x=x, y=y, family='poisson', cv.method='rcv',
#  epsilon=1e-5, alpha=12:40/80, max.it=500)
#
# the optimal is found at alpha = 17/80
a <- rgam(x=x, y=y, family='poisson', cv.method='rcv', epsilon=1e-7, alpha=17/80, max.it=500) #0.2125

b <- gam(y~lo(x, span= 0.2125), family=poisson)

# predicted values
pr.rgam <- predict(a, type='response')
pr.gam <- predict(b, type='response')

# predicted values with the GCV fit
detach('package:gam')
library(mgcv)
pr.gcv <- predict(gam(y~s(x), family=poisson), type='response')


plot(x,y, xlab='Week', ylab='ILI visits', pch=19, col='grey75')
lines(sort(x), pr.rgam[order(x)], lwd=3, col='red')
lines(sort(x), pr.gam[order(x)], lwd=3, col='blue')
lines(sort(x), pr.gcv[order(x)], lwd=3, lty=2, col='black')
legend(2, 27000, legend=c('GAM', 'RGAM', 'MGCV'), lwd=4, col=c('blue', 'red', 'black'))


