### Weighted and robust archetypal analysis: ozone data set
###
### Analysis used in 'Weighted and Robust Archetypal Analysis' by
### Manuel J. A. Eugster and Friedrich Leisch.

library('archetypes')



### Data set:

data('Ozone', package = 'mlbench')

oz <- Ozone[, -c(1, 2, 3, 9)]
oz <- na.omit(oz)
colnames(oz) <- c('OZONE', '500MH', 'WDSP', 'HMDTY', 'STMP',
                  'INVHT', 'PRGRT', 'INVTMP', 'VZBLTY')

oz0 <- scale(oz)



### The three original archetypes:

set.seed(1234)
a.oz <- archetypes(oz0, 3)

parameters(a.oz)
barplot(a.oz, oz0, percentiles = TRUE)

panorama(a.oz, oz0)



### Data set with outliers:

set.seed(1234)
outliers <- t(sapply(runif(5, min = 1.5, max = 2),
                     function(x)
                     x * apply(oz, 2, max) + apply(oz, 2, IQR)))

oz1 <- scale(rbind(oz, outliers))


pairs(oz1)



### Original archetypal algorithm:

set.seed(1234)
a.oz1 <- archetypes(oz1, 3)

parameters(a.oz1)
barplot(a.oz1, oz1, percentiles = TRUE)

panorama(a.oz1, oz1)



### Robust archetypal algorithm:

set.seed(1236)
ra.oz1 <- robustArchetypes(oz1, 3)

parameters(ra.oz1)
barplot(ra.oz1, oz1, percentiles = TRUE)

panorama(a.oz1, oz1)

plot(rss(ra.oz1, type = 'single'), xlab = '', ylab = 'RSS')
plot(weights(ra.oz1, type = 'reweights'), xlab = '', ylab = 'Weight')

