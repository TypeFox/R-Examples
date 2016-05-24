###################
## Preliminaries ##
###################

## packages
library("coin")
library("strucchange")
library("lattice")

## random seed
rseed <- 20061103

## theme for lattice/trellis graphics
trellis.par.set(theme = canonical.theme(color = FALSE))


######################
## Boston homicides ##
######################

## time series plot
data("BostonHomicide", package = "strucchange")
hom_month <- zoo(as.vector(BostonHomicide$homicides), as.yearmon(as.vector(time(BostonHomicide$homicides))))
hom_year <- aggregate(hom_month, function(x) as.numeric(floor(x)), mean)
plot(hom_month, col = grey(0.7), lwd = 2, ylab = "Number of homicides", xlab = "Time")
lines(hom_year, type = "b")

## monthly data
hom_month <- data.frame(
  log_homicides = log(coredata(hom_month) + 0.5),
  time = as.numeric(time(hom_month)))

## asymptotic unconditional test
sctest(gefp(log_homicides ~ 1, data = hom_month), functional = supLM(0.1))

## asymptotic conditional test
set.seed(rseed)
maxstat_test(log_homicides ~ time, data = hom_month, minprob = 0.1)

## approximate conditional test
set.seed(rseed)
maxstat_test(log_homicides ~ time, data = hom_month, minprob = 0.1, distribution = approximate(9999))

## annual data
hom_year <- data.frame(
  homicides = coredata(hom_year),
  time = time(hom_year))

## asymptotic unconditional test
sctest(gefp(homicides ~ 1, data = hom_year), functional = supLM(0.1))

## asymptotic conditional test
set.seed(rseed)
maxstat_test(homicides ~ time, data = hom_year, minprob = 0.1)

## approximate conditional test
set.seed(rseed)
maxstat_test(homicides ~ time, data = hom_year, minprob = 0.1, distribution = approximate(9999))
## note that the manuscript computes the exact conditional p-value
## by exhaustive search (rather than approximating via 10,000 simulations)


###########################
## Hiring discrimination ##
###########################

## data
hire <- cbind(c(2, 0, 0, 0, 5, 14), c(427, 86, 104, 180, 111, 59))
hire <- data.frame(
  resp = factor(unlist(sapply(1:nrow(hire), function(i) rep(c("female", "male"), hire[i,])))),
  time = as.numeric(rep(1991:1996, rowSums(hire))))

## visualization
set.seed(rseed)
maxstat_test(resp ~ time, data = hire, distribution = approximate(9999), minprob = 0.05)


################
## CO2 reflux ##
################

## data
data("CWD", package = "coin")
lwood <- reshape(CWD, varying = paste("sample", c(2:4,6:8), sep = ""),
  timevar = "tree", direction = "long", v.names =  "sample")
lwood$tree <- factor(lwood$tree)

## visualization
print(xyplot(sample ~ time | tree, data = lwood, type = "b",
  scales = list(y = list(relation = "free")),
  layout = c(3, 2), xlab = "Time", ylab = expression(paste(CO[2], " reflux")),
  as.table = TRUE))

## test
(cwd_mt <- maxstat_test(sample2 + sample3 + sample4 + sample6 + sample7 + sample8 ~ trend, 
  data = CWD, distribution = approximate(1e5)))

## maximally selected statistics with 5% critical value
cwd_st <- statistic(cwd_mt, "standardized")
cwd_st <- data.frame(cwd_st, time = CWD$time[1] + CWD$trend[2:11])
cwd_st <- reshape(cwd_st, varying = paste("sample", c(2:4,6:8), sep = ""),
  timevar = "tree", direction = "long", v.names = "sample")
cwd_st$tree <- factor(cwd_st$tree)
cwd_q <- qperm(cwd_mt, 0.95)
print(xyplot(sample ~ time | tree, data = cwd_st, type = "b",
  panel = function(...) {
    panel.xyplot(...)
    panel.abline(h = c(-cwd_q, cwd_q), col = 2)
    panel.abline(h = 0, col = "gray")
  },
  layout = c(3, 2), xlab = "Time", ylab = "Test statistics", ylim = c(-3.5, 3.5),
  as.table = TRUE))


##################################
## Dow Jones Industrial Average ##
##################################

## data
data("DJIA", package = "strucchange")
djia <- diff(log(DJIA)) * 100
djia_res <- coredata(djia - mean(djia))
djia_trafo <- cbind(Intercept = djia_res, Variance = djia_res^2 - mean(djia_res^2))
djia_time <- time(djia)

## visualization
plot(djia, xlab = "Time", ylab = "Dow Jones stock returns")

## test
set.seed(rseed)
djia_mt <- maxstat_test(djia_trafo ~ as.numeric(djia_time), distribution = approximate(9999), minprob = 0.1)

## maximally selected statistics
apply(abs(statistic(djia_mt, "standardized")), 2, max)

## critical values
qperm(djia_mt, 0.95)

## breakpoint
djia_point <- djia_time[floor(length(djia) * 0.1) + which.max(abs(statistic(djia_mt, "standardized")[,2]))]

## autocorrelation
djia_pre  <- window(djia, end = djia_point)
djia_post <- window(djia, start = djia_point + 1)
ar1 <- function(x, digits = 3) {
  x <- coredata(x)
  x <- x - mean(x)
  c(acf(x, plot = FALSE)$acf[2], acf(x^2, plot = FALSE)$acf[2])
}
ar1(djia_pre)
ar1(djia_post)


#######################
## Economic journals ##
#######################

## data
data("Journals", package = "AER")
Journals$age <- 2000 - Journals$foundingyear
Journals <- Journals[order(Journals$age),]

## model in root node
jour_lm  <- lm(log(subs) ~ log(price/citations), data = Journals)

## test in root node
set.seed(rseed)
(jour_mt <- maxstat_test(estfun(jour_lm) ~ age, data = Journals, distribution = approximate(9999), minprob = 0.1))

## test information
jour_critval <- qperm(jour_mt, 0.95^(1/4))
jour_process <- statistic(jour_mt, "standardized")
jour_process <- zoo(jour_process, as.numeric(sapply(strsplit(rownames(jour_process), "x <= "), tail, 1)))
colnames(jour_process) <- c("Intercept", "Slope")
jour_point <- time(jour_process)[apply(coredata(abs(jour_process)), 2, which.max)]

## visualization of test in root node
mypanel <- function(x, y, subscripts, groups, panel = panel.xyplot,
  col = 1, type = "p", pch = 20, lty = 1, lwd = 1, ...)
{
  col <- rep(as.list(col), length = nlevels(groups))
  type <- rep(as.list(type), length = nlevels(groups))
  pch <- rep(as.list(pch), length = nlevels(groups))
  lty <- rep(as.list(lty), length = nlevels(groups))
  lwd <- rep(as.list(lwd), length = nlevels(groups))

  for(g in 1:nlevels(groups)) {
    idx <- g == groups[subscripts]
    if (any(idx)) panel(x[idx], y[idx], ...,
      col = col[[g]], type = type[[g]], pch = pch[[g]],
      lty = lty[[g]], lwd = lwd[[g]])
    grid::grid.lines(y = grid::unit(0, "native"), gp = grid::gpar(col = "gray"))
    grid::grid.lines(y = grid::unit(jour_critval, "native"), gp = grid::gpar(col = 2))
  }
}
print(xyplot(abs(jour_process), panel = mypanel,
  xlab = "Time", type = "b", ylim = c(0, 6), as.table = TRUE))


### fit node with h(Y) = score and g(X) = maxstat_trafo
ytrf <- function(data) {
    ret <- estfun(lm(data[[1]] ~ data[[2]]))
    attr(ret, "assign") <- 1:2
    ret
}
xtrf <- function(data) trafo(data, numeric_trafo = maxstat_trafo)
mynode <- function(data) {
  vars <- c("society", "pages", "charpp", "age")
  sapply(vars, function(v) {
    f <- as.formula(paste("log(subs) + log(price/citations) ~", v))
    it <- independence_test(f, data = data,
      ytrafo = ytrf, xtrafo = xtrf, distribution = approximate(9999))
    c(statistic(it), 1 - (1 - pvalue(it))^length(vars))
  })
}

## fit all tree elements
jour_node <- factor(Journals$age <= 18, levels = c(TRUE, FALSE), labels = c("Node 2", "Node 3"))
jour_lm2 <- lm(log(subs) ~ log(price/citations), data = Journals[jour_node == "Node 2",])
jour_lm3 <- lm(log(subs) ~ log(price/citations), data = Journals[jour_node == "Node 3",])

## conduct tests in all leaves
set.seed(rseed)
jour_tree <- list(
  mynode(Journals),
  mynode(Journals[jour_node == "Node 2",]),
  mynode(data = Journals[jour_node == "Node 3",])
)

## fitted models
plot(log(subs) ~ log(price/citations), data = Journals,
  xlab = "log(price/citations)", ylab = "log(subscriptions)",
  pch = c(24, 21)[jour_node], bg = hcl(c(0, 240), 50, 70)[jour_node])
abline(coef(jour_lm2), col = hcl(  0, 80, 30), lty = 5, lwd = 1.7)
abline(coef(jour_lm3), col = hcl(240, 80, 30), lty = 1, lwd = 1.7)
legend("bottomleft", c(expression(age > 18), expression(age <= 18)),
  pch = c(19, 17), lty = c(1, 5), col = hcl(c(240, 0), 80, 30), bty = "n")
