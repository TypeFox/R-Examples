#### codeGolliniRougier2015.R,

### This code prepared by 
### Isabella Gollini and Jonty Rougier for the paper
### "Rapidly bounding the probability in the upper tail of the loss distribution"
### arxiv:1507.01853


if (!exists("op"))
  op <- par(no.readonly = TRUE)
graphics.off()
dev.new(width = 6, height = 5)

set.seed(1001)
library(tailloss); 
UShurricane <- ELT(UShurricane) # no compression

######  
#### effect of merging
######

## for formatting summary table of an ELT object

toTable <- function(ELT, file = "", top = 5, bottom = top) {

  fmt <- function(x)
    cbind(ID = as.character(x$ID),
          Rate = formatC(x$Rate, format = "f", digits = 5, width = 10),
          Loss = formatC(x$Loss, format = "f", digits = 0, width = 10))
  
  write.table(rbind(fmt(head(ELT, top)),
                    rep('$\\vdots$',3),
                    fmt(tail(ELT, bottom))),
              file = file,
              sep = ' & ', eol = ' \\\\ \n', quote = FALSE, row.names = FALSE, col.names = FALSE)
}

## sequence of s values in natural units of $

svals <- seq(from = 0e6, to = 40e6, len = 101)

##

par(op)
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), las = 1, cex = 0.9)

dvals <- seq.int(from = -7, to = 0)

EP <- numeric(0); nrows <- numeric(0)
for (d in dvals) {
  ELTcomp <- compressELT(UShurricane, digits = d)
  toTable(ELTcomp, file = sprintf("ELTcompd%i.txt", d))
  defpars <- list(ELT = ELTcomp,
                  s = svals,
                  t = 1,
                  theta = 0)
  defpars$s <- defpars$s * 10^ifelse(is.na(d), 0, d) # local scale
  EP <- cbind(EP, do.call("fMoment", defpars)[, 2])
  nrows <- c(nrows, nrow(ELTcomp))
}

## plot in $m

lty <- seq(along = dvals)
matplot(1e-6 * svals, EP, type = "l",
        lty = lty, col = 1, lwd = 1,
        xlab = "Values for s, $m", ylab = expression(Pr(S >= s)), main = "",
        axes = FALSE)
axis(1, pos = 0); axis(2)
leg <- paste("d = ",
             prettyNum(dvals, width = 2), ": ",
             prettyNum(nrows), sep = "")
legend("topright", legend = leg, pch = NA, lty = lty,
       title = "Number of rows", bty = "n", xpd = NA)

zoombox(1e-6 * svals, EP, x0 = c(20, 40), y0 = c(0, 0.15), y1 = c(0.2, 0.5),
        type = "l", lty = lty, col = 1, lwd = 1)

dev.print(pdf, file = "effectOfCompression.pdf")

###### 
#### Full Experiment
###### 

## Jeffreys confidence interval, see Brown et al (2001)

jeffreysCI <- function(x, n, alpha = 0.05, drop = TRUE) {
  stopifnot(length(n) == 1, x %in% 0:n,
            length(alpha) == 1, 0 < alpha, alpha < 1)
  alpha <- c(alpha / 2, 1 - alpha / 2)
  CI <- outer(x, alpha, function(x, alpha) qbeta(alpha, x + 0.5, n - x + 0.5))
  if (any(x == 0)) CI[x == 0, 1] <- 0
  if (any(x == n)) CI[x == n, 2] <- 1
  colnames(CI) <- paste(formatC(100 * alpha, format = "f", digits = 1), "%", sep = "")
  if (drop && length(x) == 1) CI <- CI[1, ]
  CI
}

## default settings for the arguments of the methods

defpars <- list(ELT = NULL,
                s = seq(from = 5e6, to = 40e6, len = 101),
                t = 1,
                theta = 0,
                cap = Inf)

## other useful arguments that can be changed

dobox <- TRUE       # draw zoombox?
alpha <- 0.05       # 1-confidence level
x0 <- c(30e6, 40e6) # zoombox xlim
stem <- "0"         # for identifying figures & tables for the different experiments

allmethods <- c("Panjer", "MonteCarlo", "Markov", "Cantelli", "Moment", "Chernoff")
names(allmethods) <- allmethods

digits <- -4                                   # these are 
methods <- c("Panjer", "MonteCarlo", "Moment") # just examples

## here's the expression to do one cycle through methods

runCalc <- expression({

  methods <- allmethods[methods] # keep names on, helpful for matching later
  
  ## need to make local modifications
  
  localpars <- defpars
  localpars$ELT <- compressELT(UShurricane, digits = digits)

  ## local scale for $
  
  localpars$s <- localpars$s * 10^ifelse(is.na(digits), 0, digits)
  if (!is.infinite(localpars$cap))
    localpars$cap <- localpars$cap * 10^ifelse(is.na(digits), 0, digits)
  
  ## compute whole curve
  
  compAll <- lapply(methods, function(mm) {
    cat(format(Sys.time(), "%X"), ":", mm, "\n")

    ## very local modifications
    
    if(dobox && mm == "MonteCarlo")
      localpars <- c(localpars, verbose = TRUE)

    st <- system.time(robj <- do.call(sprintf("f%s", mm), localpars)) # calls gc() first
    list(sys = st, res = robj)
  })

  ## produce the picture (B&W)
  
  EP <- sapply(compAll, function(x) x$res[, 2])
  mm <- match(colnames(EP), allmethods) # for line types
  col <- 1

  x <- 1e-6 * defpars$s # measure in millions
  matplot(x, EP, ylim = c(0, 1),
          type = "l", lty = mm, col = col, lwd = 2,
          xlab = "Values for s, $m",
          ylab = substitute(Pr(S[tval] >= s), list(tval = localpars$t)),
          main = "",
          axes = FALSE)
  axis(1, pos = 0); axis(2)

  ## Jeffreys 95% CI for largest value in zoombox (Brown et al, 2001)

  s0 <- x0[2]
  rsam <- attr(compAll$MonteCarlo$res, "rsam")
  n <- length(rsam)
  X <- sum(rsam >= s0 * 10^digits)
  rangeCB <- jeffreysCI(X, n, alpha)
  cat(sprintf("Monte Carlo %.0f%% CI at $%.1fm is (%.3f%%, %.3f%%)\n",
              100 * (1 - alpha), s0 * 1e-6,
              100 * rangeCB[1], 100 * rangeCB[2]))

  leg <- colnames(EP)
  if ("MonteCarlo" %in% leg)
    leg[leg == "MonteCarlo"] <- "Monte Carlo"
  legend("topright", legend = leg, lty = mm, col = col, lwd = 2, pch = NA, bty = "n", xpd = NA)

  if(dobox) local({ # mask change to x0
    
    x0 <- 1e-6 * x0
    y0 <- c(0, .05)
    y1 <- c(.1, .6)
    zoombox(x, EP, x0 = x0, y0 = y0, y1 = y1, type = "l", lty = mm, col = col, lwd = 2)
          
    yCB <- y1[1] + diff(y1) * (rangeCB - y0[1])/diff(y0)
    axis(4, yCB, c("",""), pos = x0[2], col = 1, lwd = 2)
    label <- sprintf("Monte Carlo %.0f%% CI", 100 * (1 - alpha))
    mtext(label, side = 4, cex = 0.7, line = 0, at = sum(yCB)/2, las = 0)
  })

  if (digits == -4) # only draw this case
    dev.print(pdf, file = sprintf("EP%s.pdf", stem))
  
}) # end of expression

## use as, eg,

eval(runCalc)

#######
#######

## expression to wrap runCalc for different values of digits

mymethods <- c("Panjer", "MonteCarlo", "Moment", "Chernoff", "Cantelli", "Markov")

topPanjer <- -3 # In the paper topPanjer <- -2

runDvals <- expression({

  dvals <- c(-4, -3, -2, -1, 0)
  compOut <- vector("list", length(dvals))
  for (i in seq(along = dvals)) {
    digits <- dvals[i]
    cat(sprintf("\n** digits = %.0f **\n", digits))
    methods <- if (digits <= topPanjer)
      mymethods
    else
      setdiff(mymethods, "Panjer") # life's too short!
    eval(runCalc)
    if (!("Panjer" %in% methods))
      compAll$Panjer$sys <- c("Elapsed" = NA)
    compOut[[i]] <- compAll
  }
})

## write a table of timings

timingTable <- function(compOut, dvals, stem, digits = 3) {

  stopifnot(length(compOut) == length(dvals))
  mymethods <- names(compOut[[1]])
  timing <- sapply(compOut, function(x)
                   sapply(mymethods, function(mm)
                          x[[mm]]$sys["elapsed"]))
  colnames(timing) <- paste("$d = ", dvals, "$", sep = "")
  rownames(timing) <- mymethods
  print(timing)
  timing <- formatC(timing, format = "f", digits = digits)
  timing[timing == " NA"] <- ""
  write.table(timing,
              file = sprintf("timings%s.tex", stem),
              quote = FALSE, sep = " & ", eol = " \\\\ \n",
              col.names = NA)
  invisible(timing)
}

#### here we go!

defpars$ELT <- list(NULL)
print(defpars)

set.seed(1001); stem <- "A"

par(op)
par(mar = c(4, 4, 1, 2), mgp = c(2.5, 1, 0), las = 1, cex = 0.9)
eval(runDvals)
timingTable(compOut, dvals, stem = stem)

####### draw a picture for different thetas

tvals <- c(0.1, 0.25, 0.5, 1, 2)
x <- 1e6
xx <- seq(from = 0.01e6, to = 3e6, length = 301)
yy <- sapply(tvals, function(theta) {
  alpha <- 1 / theta^2
  beta <- alpha / x
  dgamma(xx, shape = alpha, rate = beta)
})

matplot(xx, yy, type = "l", col = 1, xlab = "Values for x, $m",
        ylab = expression(paste("Probability density, times ", 10^6)),
        axes = FALSE, bty = "n", lwd = 2)
axis(1, pretty(xx), pretty(xx) * 1e-6)
axis(2, pretty(yy), pretty(yy) * 1e6, cex.axis = 0.9)
legend("topright", legend = prettyNum(tvals), lty = seq(along = tvals), lwd = 2,
       pch = NA,
       title = expression(paste("Values for ", theta)),
       bty = "n", xpd = NA)

dev.print(pdf, file = "thetas.pdf")

####### Same with theta = 0.5, need to increase svals and x0

defpars$theta <- 0.5
defpars$s <- seq(from = 0e6, to = 60e6, len = 101)
x0 <- c(40e6, 60e6)

## local mods hold from now on

topPanjer <- -3
mymethods <- setdiff(mymethods, "Markov")

set.seed(2001); stem <- "B"

par(op)
par(mar = c(4, 4, 1, 2), mgp = c(2.5, 1, 0), las = 1, cex = 0.9)
eval(runDvals)
timingTable(compOut, dvals, stem = stem)

####### now introduce cap at $5m

defpars$cap <- 5e6
defpars$s <- seq(from = 0e6, to = 40e6, length = 101)
x0 <- c(30e6, 40e6)

set.seed(3001); stem <- "C"

par(op)
par(mar = c(4, 4, 1, 2), mgp = c(2.5, 1, 0), las = 1, cex = 0.9)
eval(runDvals)
timingTable(compOut, dvals, stem = stem)

####### and finally run the calculation out to ten years

defpars$t <- 10
defpars$s <- seq(from = 0e6, to = 120e6, length = 101)
x0 <- c(90e6, 120e6)

set.seed(4001); stem <- "D"

par(op)
par(mar = c(4, 4, 1, 2), mgp = c(2.5, 1, 0), las = 1, cex = 0.9)
eval(runDvals)
timingTable(compOut, dvals, stem = stem)
