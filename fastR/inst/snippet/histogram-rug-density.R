x <- rnorm(100)
p <- histogram(~x, type='density',
            panel=function(x,y,...) {
                panel.rug(x,...)
                panel.histogram(x,...)
                panel.mathdensity(
                    dmath=dnorm, args=list(mean=mean(x),sd=sd(x)),
                    lwd=5, col="black", lty=1, alpha=0.5,
                    ...)
            }
     )
print(p)
