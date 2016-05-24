## ----knitr-setup, include=FALSE,cache=FALSE,echo=FALSE--------------
library(knitr)
opts_chunk$set(comment = NA, size = 'normalsize', prompt = TRUE, highlight = FALSE, 
               cache = TRUE, crop = TRUE, concordance = FALSE, fig.align='center',
               fig.path='paper-figures/Paper-', out.width="0.4\\textwidth", fig.lp = "F:",
               background = "#FFFFFF")
opts_knit$set(out.format = "latex")
knit_hooks$set(crop = hook_pdfcrop)
options(width = 70, prompt = "R> ", continue = "+  ", digits = 3, useFancyQuotes = FALSE)

thm <- knit_theme$get('default')

# Set default font colour (fgcolor) to black
thm$highlight <- "\\definecolor{fgcolor}{rgb}{0, 0, 0}\n
\\newcommand{\\hlnum}[1]{\\textcolor[rgb]{0.686,0.059,0.569}{#1}}%\n
\\newcommand{\\hlstr}[1]{\\textcolor[rgb]{0.192,0.494,0.8}{#1}}%\n
\\newcommand{\\hlcom}[1]{\\textcolor[rgb]{0.678,0.584,0.686}{\\textit{#1}}}%\n
\\newcommand{\\hlopt}[1]{\\textcolor[rgb]{0,0,0}{#1}}%\n
\\newcommand{\\hlstd}[1]{\\textcolor[rgb]{0.345,0.345,0.345}{#1}}%\n
\\newcommand{\\hlkwa}[1]{\\textcolor[rgb]{0.161,0.373,0.58}{\\textbf{#1}}}%\n\\newcommand{\\hlkwb}[1]{\\textcolor[rgb]{0.69,0.353,0.396}{#1}}%\n
\\newcommand{\\hlkwc}[1]{\\textcolor[rgb]{0.333,0.667,0.333}{#1}}%\n
\\newcommand{\\hlkwd}[1]{\\textcolor[rgb]{0.737,0.353,0.396}{\\textbf{#1}}}%"

thm$background <- "#FFFFFF"

thm$fontstyle <- "italic"

knit_theme$set(thm)

## Specific to vignette
options(cl.cores = 1)

## ----spv-example, echo=FALSE, fig.width=5, fig.height=5, fig.cap = "An example of a variance dispersion graph.",out.width="0.35\\textwidth"----
library("Vdgraph")
library("vdg")
data("D310")
set.seed(1)
vdgex <- spv(n = 100000, design = D310, formula = ~.^2, at = TRUE)
plot(vdgex, which = "vdgquantile", tau = c(0, 1))[[1]] + theme(legend.position = "none")

## ----qp-ex, echo = FALSE, fig.width = 6, fig.height=4.5, fig.cap = "An example of a quantile plot, corresponding to the example in Figure~\\ref{F:spv-example}.",out.width="0.35\\textwidth"----
set.seed(1)
qpex <- spv(n = 50000, design = D310, formula = ~.^2, at = TRUE, nr.rad = 6)
my_ecdf <- function(x) {
  xs <- sort(x)
  xun <- unique(xs)
  n <- length(xun)
  prop <- rep(NA, n)
  for (i in seq_along(xun))
    prop[i] <- sum(xs <= xun[i]) / n
  return(cbind(x = xun, y = prop))
}
ds <- formatC(sqrt(rowSums(qpex$sample^2)), format = "f", digits = 3)
lst <- split(qpex$spv[-1], f = factor(ds[-1]))
ecd <- lapply(lst, my_ecdf)
df <- as.data.frame(do.call(rbind, ecd))
df$Radius <- factor(ds[-1])
df$Distance <- sqrt(rowSums(qpex$sample^2))[-1]
ggplot(data = df, mapping = aes(y = x, x = y, group = Radius, colour = Radius)) + geom_line(size = 1) + ylab("SPV Quantile") + xlab("Proportion")

## ----fds-example, echo = FALSE, fig.width = 5.5, fig.height=5.5, out.width="0.35\\textwidth", fig.cap = "An example of and FDS plot, corresponding to Figures~\\ref{F:spv-example} and \\ref{F:qp-ex}.",results='hide'----
set.seed(1)
fdsex <- spv(n = 50000, design = D310, formula = ~.^2)
plot(fdsex, which = "fds", np = 100)

## ----lhs,fig.height=5.5,fig.width=5.5,fig.cap="An example of an LHS of 10 points in a two-dimensional design space.",results='hide'----
library("vdg")
set.seed(8745)
samp <- LHS(n = 10, m = 2, lim = c(-1, 1))
plot(samp, main = "", pty = "s", pch = 16, ylim = c(-1, 1), 
asp = 1, xlab = expression(X[1]), ylab = expression(X[2]))
abline(h = seq(-1, 1, length.out = 10), 
v = seq(-1, 1, length.out = 10), lty = 3, col = "grey")

## ----vign, eval=FALSE-----------------------------------------------
#  vignette(topic = "vdg", package = "vdg")

## ----load-roq-------------------------------------------------------
library("Vdgraph")
data("D416B")
data("D416C")

## ----vdgroq,fig.width=9, fig.height=5.5, results='hide', fig.cap="A VDG for Roquemore's hybrid designs D416B and D416C for a full quadratic model.",out.width="0.7\\textwidth"----
quad4 <- formula( ~ (x1 + x2 + x3 + x4)^2 +  I(x1^2) + I(x2^2) + 
I(x3^2) + I(x4^2))
set.seed(1234)
spv1 <- spv(n = 5000, design = list(D416B = D416B, 
D416C = D416C), formula = quad4)
plot(spv1, which = "vdgboth")

## ----quad4, eval=FALSE----------------------------------------------
#  quad4 <- formula( ~ .^2 +  I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2))

## ----ex1-bothroqfds,fig.width=6, fig.height=5, results='hide',fig.cap="A standard and variance ratio FDS plot for Roquemore's hybrid designs D416B and D416C for a full quadratic model.",fig.show='hold'----
plot(spv1, which = "fds")
plot(spv1, which = "fds", VRFDS = TRUE, np = 100)

## ----vdgroq-theme,fig.width=7, fig.height=4, results='hide', fig.cap="A second version of Figure~\\ref{F:vdgroq}.",out.width="0.5\\textwidth"----
p <- plot(spv1, which = "vdgboth")
p$vdgboth + theme_bw() + theme(panel.grid = element_blank())

## ----make-ccd3------------------------------------------------------
library("rsm") 
ccd3 <- as.data.frame(ccd(basis = 3, n0 = 4, 
alpha = "spherical", oneblock = TRUE))[, 3:5]

## ----algdes-cand,results='hide'-------------------------------------
set.seed(8619) 
cand <- runif_sphere(n = 10000, m = 3)
colnames(cand) <- colnames(ccd3)

## ----algdes-AD------------------------------------------------------
quad3 <- formula( ~ (x1 + x2 + x3)^2 + I(x1^2) + I(x2^2) + I(x3^2))
library("AlgDesign") 
set.seed(3476)
desD <- optFederov(quad3, data = cand, nTrials = 22, criterion = "D")
desA <- optFederov(quad3, data = cand, nTrials = 22, criterion = "A")

## ----ex2-spv,results='hide',fig.show='hide'-------------------------
spv2 <- spv(n = 10000, formula = quad3, 
design = list(CCD = ccd3, D = desD$design, A = desA$design))
plot(spv2, which = 2:3) 

## ----ccdfds,include=FALSE,echo=FALSE,fig.height=5,fig.width=6,fig.show='hide'----
plot(spv2, which = 2)[[1]] + theme(plot.title = element_text(size = 16), legend.position = "none") 

## ----ccdvdg,include=FALSE,echo=FALSE,fig.height=5,fig.width=8,fig.show='hide'----
plot(spv2, which = 3)[[1]] + theme(plot.title = element_text(size = 16)) 

## ----GJdesreg,echo=FALSE,fig.width=5,fig.height=5,fig.cap="The design region and D-optimal design of \\cite{Goo2011}. Some runs are replicated."----
df <- data.frame(Time = c(360, 420, 720, 720, 660, 360, 360), Temperature = c(550, 550, 523, 520, 520, 529, 550))
p <- ggplot(data = df, aes(x = Time, y = Temperature)) + geom_path() + geom_point(data = GJ54)
p

## ----keepfun--------------------------------------------------------
keepfun <- function(x) apply(x >= -1 & x <= 1, 1, all) & 
(x[, 2] <= -1.08 * x[, 1] + 0.28) & (x[, 2] >= -0.36 * x[, 1] - 0.76) 

## ----ex3-for--------------------------------------------------------
cube2 <- formula( ~ (Time + Temperature)^2 + I(Time^2) + 
I(Temperature^2) + I(Time^3) + I(Temperature^3) + 
Time:I(Temperature^2) + I(Time^2):Temperature)
GJmod <- update(cube2, ~ . - I(Time^3) - I(Time^2):Temperature) 

## ----ex3-spv--------------------------------------------------------
spv3 <- spv(n = 10000, design = stdrange(GJ54), type = "lhs", 
formula = list(Cubic = cube2, GoosJones = GJmod), 
keepfun = keepfun) 

## ----vdggj-noplot,eval=FALSE----------------------------------------
#  plot(spv3, which = 1, points.size = 2)

## ----vdggj,fig.width=8,fig.height=4,results='hide',echo=FALSE,fig.cap="VDG for the D-optimal design of \\citet{Goo2011}, for the two models.",out.width="0.5\\textwidth"----
plot(spv3, which = 1, points.size = 2)[[1]] + theme(plot.title = element_text(size = 16)) 

## ----rgl-code,eval=FALSE--------------------------------------------
#  library("rgl")
#  with(spv3$Cubic, plot3d(x = sample[, "Time"],
#  y = sample[, "Temperature"], z = spv))

