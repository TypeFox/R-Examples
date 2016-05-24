library(grImport)
data <- read.table(textConnection('
          name ev.trt n.trt ev.ctrl n.ctrl
1     Auckland     36   532      60    538
2        Block      1    69       5     61
3        Doran      4    81      11     63
4        Gamsu     14   131      20    137
5     Morrison      3    67       7     59
6 Papageorgiou      1    71       7     75
7      Tauesch      8    56      10     71
'
), header=TRUE)

cer <- w.approx.cer(data[["ev.ctrl"]], data[["n.ctrl"]])

sm <- "RR"
if (requireNamespace("meta", quietly = TRUE)) {
    ## Calculate the pooled OR or RR point estimate, we use meta package here
    m <- with(data,
             meta::metabin(ev.trt, n.trt, ev.ctrl, n.ctrl, sm=sm))
    point <- exp(m$TE.random) # meta returns outcomes on the log scale
} else {
    ## Calculated Random Effects RR, using the meta package
    point <- 0.5710092
}

ier <- calc.ier(cer, point, sm)

u <- uplift(ier, cer, F)
plot(u, fig.title="Example", fig.cap="Example from rMeta")


## Test with zero CER
data <- read.table(textConnection('
         name   ai n1i  ci  n2i days
1      Zarate    6  15   0  15     1
2    Murrough   30  47   0  25     1
3      Zarate   13  18   0  18     1
4     Lapidus    8  18   0  18     1
5 Diazgranado   10  18   0  18     1
'), header=TRUE)

# cer <- w.approx.cer(data$ci, data$n2i)

## Draw documentation graphics
## u <- uplift(ier, cer, F)
## pdf("man/figures/green.pdf", 8, 10)
## plot(u, fig.title="Example", fig.cap="Example from rMeta")
## dev.off()
## png("man/figures/green.png", 800, 1000)
## plot(u, fig.title="Example", fig.cap="Example from rMeta")
## dev.off()

## u <- uplift(ier, cer, T)
## pdf("man/figures/red.pdf", 8, 10)
## plot(u, fig.title="Example", fig.cap="Example from rMeta")
## dev.off()
## png("man/figures/red.png", 800, 1000)
## plot(u, fig.title="Example", fig.cap="Example from rMeta")
## dev.off()
