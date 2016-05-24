## ---- eval=FALSE---------------------------------------------------------
#  install.packages('latex2exp')

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github('stefano-meschiari/latex2exp')

## ------------------------------------------------------------------------
library(latex2exp)

## ---- eval=FALSE---------------------------------------------------------
#  TeX('$\\alpha^\\beta$')

## ---- warning=FALSE------------------------------------------------------
x <- seq(0, 4, length.out=100)
alpha <- 1:5

plot(x, xlim=c(0, 4), ylim=c(0, 10), 
     xlab='x', ylab=TeX('$\\alpha  x^\\alpha$, where $\\alpha \\in 1\\ldots 5$'), 
     type='n', main=TeX('Using $\\LaTeX$ for plotting in base graphics!'))

invisible(sapply(alpha, function(a) lines(x, a*x^a, col=a)))

legend('topleft', legend=TeX(sprintf("$\\alpha = %d$", alpha)), 
       lwd=1, col=alpha)

## ---- warning=FALSE------------------------------------------------------
library(plyr)
library(ggplot2)
x <- seq(0, 4, length.out=100)
alpha <- 1:5
data <- mdply(alpha, function(a, x) data.frame(v=a*x^a, x=x), x)

p <- ggplot(data, aes(x=x, y=v, color=X1)) +
    geom_line() + 
    ylab(TeX('$\\alpha  x^\\alpha$, where $\\alpha \\in 1\\ldots 5$')) +
    ggtitle(TeX('Using $\\LaTeX$ for plotting in ggplot2. I $\\heartsuit$ ggplot!')) +
    coord_cartesian(ylim=c(-1, 10)) +
    guides(color=guide_legend(title=NULL)) +
    scale_color_discrete(labels=lapply(sprintf('$\\alpha = %d$', alpha), TeX)) 
    # Note that ggplot2 legend labels must be lists of expressions, not vectors of expressions

print(p)

## ---- fig.width=10, fig.height=2-----------------------------------------
plot(TeX("A $\\LaTeX$ formula: $\\frac{2hc^2}{\\lambda^5} \\, 
               \\frac{1}{e^{\\frac{hc}{\\lambda k_B T}} - 1}$"), cex=2)

## ---- eval=FALSE---------------------------------------------------------
#  TeX('latexString')

## ---- eval=FALSE---------------------------------------------------------
#  TeX('latexString', output=c('expression', 'character', 'ast'))

## ---- eval=FALSE---------------------------------------------------------
#  latex2exp_examples()

## ---- eval=FALSE---------------------------------------------------------
#  latex2exp_supported(plot=FALSE)

## ---- fig.width=10, fig.height=10----------------------------------------
latex2exp_supported(plot=TRUE)

## ---- fig.width=12, fig.height=7-----------------------------------------
latex2exp_examples()

