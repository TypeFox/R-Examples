require(ggplot2)

set.seed(0)
x <- sample(0:9, 100, rep=TRUE)

### SPSS like

# Standard QQ-plot
qq <- qqnorm_spss(x, 1)
plot(qq)  
ggplot(qq)

qq <- qqnorm_spss(x, 1, standardize=TRUE)
plot(qq, l.col="red")  
ggplot(qq, line=FALSE)

# Detrended QQ-plot (plottype=2)
plot(qq, plottype=2)  
ggplot(qq, plottype=2)

### R
qqnorm(x, datax=TRUE)
qqline(x, datax=TRUE)


