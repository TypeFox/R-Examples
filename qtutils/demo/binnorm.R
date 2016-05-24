
## Interactive demo of Normal approximation to binomial

library(qtbase)
library(qtutils)

wtop <- Qt $ QWidget()
ltop <- Qt $ QGridLayout()
wtop $ setLayout(ltop)
wtop $ resize(800, 600)
wtop

pval <- Qt $ QDoubleSpinBox()
pval $ setRange(0, 1)
pval $ setSingleStep(0.01)
pval $ setValue(0.5)

ltop $ addWidget(Qt $ QLabel("Prob"), 0, 0)
ltop $ addWidget(pval, 0, 1)

nval <- Qt $ QSpinBox()
nval $ setRange(1L, 100000L)
nval $ setValue(15L)

ltop $ addWidget(Qt $ QLabel("Size"), 0, 2)
ltop $ addWidget(nval, 0, 3)



plotw <- QT(width = 10, height = 7)
ltop $ addWidget(plotw, 1, 0, 1, 6)


plotApprox <- function(n = nval$value, p = pval$value)
{
    x <- rbinom(1, size = n, prob = p)
    phat <- x/n
    px <- dbinom(0:n, size = n, prob = p)
    plot(0:n/n, n * px, type = "h", col = "darkgrey", lwd = 3,
         ylim = extendrange(c(0, 1.5 * n * max(px))),
         xlab = "x", ylab = "n P(X = x)")
    curve(dnorm(x, mean = p, sd = sqrt(p * (1-p) / n)),
          n = 500, add = TRUE, col = "black")
    curve(dnorm(x, mean = phat, sd = sqrt(phat * (1-phat) / n)),
          n = 500, add = TRUE, col = "red")
    invisible()
}

plotApprox()

timer <- qtimer(150, plotApprox)

button <- Qt$QPushButton("&Start")
ltop $ addWidget(button, 0, 5)

space <- Qt$QSpacerItem(1, 1, Qt$QSizePolicy$MinimumExpanding)
ltop $ addItem(space, 0, 4)

qconnect(button, signal = "clicked", 
         handler = function() {
             cur <- button$text
             switch(cur,
                    "&Start" = {
                        button$text <- "&Stop"
                        timer$start()
                    },
                    "&Stop" = {
                        button$text <- "&Start"
                        timer$stop()
                    })
         })
         


