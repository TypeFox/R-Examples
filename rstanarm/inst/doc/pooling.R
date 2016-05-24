## ---- SETTINGS-rstan, include=FALSE--------------------------------------
ITER <- 500L
CHAINS <- 2L
CORES <- 1L
SEED <- 12345

## ---- SETTINGS-loo, include=FALSE----------------------------------------
loo.cores <- if (exists("CORES")) CORES else 1L
options(loo.cores = loo.cores)

## ---- SETTINGS-knitr, include=FALSE--------------------------------------
library(knitr)
opts_chunk$set(
  comment = NA, message = FALSE, warning = FALSE, 
  fig.align = 'center', fig.width = 7,fig.height = 5
)

## ---- load-data----------------------------------------------------------
library(rstanarm)
data(bball1970)
bball <- bball1970
print(bball)

## ---- N-K-y--------------------------------------------------------------
# A few quantities we'll use throughout
N <- nrow(bball)
K <- bball$AB
y <- bball$Hits
K_new <- bball$RemainingAB
y_new <- bball$RemainingHits

## ---- create-objects, results="hold"-------------------------------------
batting_avg <- function(x) print(format(round(x, digits = 3), nsmall = 3), quote = FALSE)
player_avgs <- y / K # player avgs through 45 AB
tot_avg <- sum(y) / sum(K) # overall avg through 45 AB

cat("Player averages through 45 at-bats:\n")
batting_avg(player_avgs)
cat("Overall average through 45 at-bats:\n")
batting_avg(tot_avg)

## ---- fig.width=7, fig.height=3, echo=FALSE------------------------------
par(mfrow = c(1,3), las = 1)
p_alpha <- function(alpha) {
  dnorm(alpha, -1, 1)
}
p_theta <- function(theta) {
  dnorm(log(theta) - log1p(-theta), -1, 1) / (theta - theta^2)
}
curve2 <- function(expr, limits, xlab, ...) {
  curve(expr, from = limits[1], to = limits[2], xlab = xlab,
    lwd = 3, bty = "l", ylab = "", cex.lab = 1.5, ...)
}
curve2(p_alpha, c(-3, 1), expression(alpha))
text(x = 0.25, y = 0.35, labels = expression(p(alpha)), cex = 1.5)
curve2(p_theta, c(0, 1), expression(theta), col = "red", ylim = c(0, 2.5))
text(x = 0.575, y = 1.5, labels = expression(p(theta)), cex = 1.5, col = "red")
curve2(p_alpha, c(-3, 1), expression(paste(alpha,", ", theta)), ylim = c(0, 2.5))
curve2(p_theta, c(0,1), col = "red", add = TRUE)
text(x = -1, y = 0.65, labels = expression(p(alpha)), cex = 1.5)
text(x = -0.5, y = 1.5, labels = expression(p(theta)), cex = 1.5, col = "red")

## ---- full-pooling, results="hide"---------------------------------------
SEED <- 101
wi_prior <- normal(-1, 1)  # weakly informative prior on log-odds
fit_pool <- stan_glm(cbind(Hits, AB - Hits) ~ 1, data = bball, family = binomial("logit"),
                     prior_intercept = wi_prior, seed = SEED)

## ---- summary-stats-function---------------------------------------------
invlogit <- plogis  # function(x) 1/(1 + exp(-x))
summary_stats <- function(posterior) {
  x <- invlogit(posterior)  # log-odds -> probabilities
  t(apply(x, 2, quantile, probs = c(0.1, 0.5, 0.9))) 
}

pool <- summary_stats(as.matrix(fit_pool))  # as.matrix extracts the posterior draws
pool <- matrix(pool,  # replicate to give each player the same estimates
               nrow(bball), ncol(pool), byrow = TRUE, 
               dimnames = list(bball$Player, c("10%", "50%", "90%")))
batting_avg(pool)

## ---- no-pooling, results="hide"-----------------------------------------
fit_nopool <- update(fit_pool, formula = . ~ 0 + Player, prior = wi_prior)
nopool <- summary_stats(as.matrix(fit_nopool))
rownames(nopool) <- as.character(bball$Player)
batting_avg(nopool)

## ---- no-pooling-print, echo=FALSE---------------------------------------
batting_avg(nopool)

## ---- partial-pooling, results="hide"------------------------------------
fit_partialpool <- 
  stan_glmer(cbind(Hits, AB - Hits) ~ (1 | Player), data = bball, family = binomial("logit"),
             prior_intercept = wi_prior, seed = SEED)

## ---- partial-pooling-shift-draws----------------------------------------
# shift each player's estimate by intercept (and then drop intercept)
shift_draws <- function(draws) {
  sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
}
alphas <- shift_draws(as.matrix(fit_partialpool))
partialpool <- summary_stats(alphas)
rownames(partialpool) <- as.character(bball$Player)
batting_avg(partialpool)

## ---- plot-observed-vs-estimated-----------------------------------------
library(ggplot2)
models <- c("complete pooling", "no pooling", "partial pooling")
estimates <- rbind(pool, nopool, partialpool)
colnames(estimates) <- c("lb", "median", "ub")
plotdata <- data.frame(estimates, 
                       observed = rep(player_avgs, times = length(models)), 
                       model = rep(models, each = N), 
                       row.names = NULL)

ggplot(plotdata, aes(x = observed, y = median, ymin = lb, ymax = ub)) +
  geom_hline(yintercept = tot_avg, color = "lightpink", size = 0.75) +
  geom_abline(intercept = 0, slope = 1, color = "skyblue") + 
  geom_linerange(color = "gray60", size = 0.75) + 
  geom_point(size = 2.5, shape = 21, fill = "gray30", color = "white", stroke = 0.2) + 
  facet_grid(. ~ model) +
  coord_fixed() +
  scale_x_continuous(breaks = c(0.2, 0.3, 0.4)) +
  labs(x = "Observed Hits / AB", y = "Predicted chance of hit") +
  ggtitle("Posterior Medians and 80% Intervals")

## ---- log_p_new----------------------------------------------------------
newdata <- data.frame(Hits = y_new, AB = K_new, Player = bball$Player)
fits <- list(Pooling = fit_pool, NoPooling = fit_nopool, 
             PartialPooling = fit_partialpool)

# compute log_p_new with each of the models in 'fits'
log_p_new <- sapply(fits, function(x) rowSums(log_lik(x, newdata)))
head(log_p_new)

## ---- log_p_new-mean-----------------------------------------------------
mean_log_p_new <- colMeans(log_p_new)
round(sort(mean_log_p_new, decreasing = TRUE), digits = 1)

## ---- log_sum_exp--------------------------------------------------------
log_sum_exp <- function(u) {
  max_u <- max(u)
  a <- 0
  for (n in 1:length(u)) {
    a <- a + exp(u[n] - max_u)
  }
  max_u + log(a)
}

# Or equivalenty using vectorization
log_sum_exp <- function(u) {
  max_u <- max(u)
  max_u + log(sum(exp(u - max_u)))
}

## ----comment=NA----------------------------------------------------------
M <- nrow(log_p_new) 
new_lps <- -log(M) + apply(log_p_new, 2, log_sum_exp)
round(sort(new_lps, decreasing = TRUE), digits = 1)

## ---- loo----------------------------------------------------------------
compare(loo(fit_partialpool), loo(fit_pool), loo(fit_nopool))

## ---- ppd----------------------------------------------------------------
newdata <- data.frame(Hits = y_new, AB = K_new, Player = bball$Player)
ppd_pool <- posterior_predict(fit_pool, newdata)
ppd_nopool <- posterior_predict(fit_nopool, newdata)
ppd_partialpool <- posterior_predict(fit_partialpool, newdata)
colnames(ppd_pool) <- colnames(ppd_nopool) <- colnames(ppd_partialpool) <- as.character(bball$Player)
colMeans(ppd_partialpool)

## ---- clemente-----------------------------------------------------------
z_1 <- ppd_partialpool[, 1]
clemente_80pct <- (y[1] + quantile(z_1, prob = c(0.1, 0.9))) / (K[1] + K_new[1])
batting_avg(clemente_80pct)

## ---- ppd-stats----------------------------------------------------------
ppd_intervals <- function(x) t(apply(x, 2, quantile, probs = c(0.25, 0.75)))
ppd_summaries <- (1 / K_new) * rbind(ppd_intervals(ppd_pool),
                                     ppd_intervals(ppd_nopool),
                                     ppd_intervals(ppd_partialpool))
df_ppd <- data.frame(player = rep(1:length(y_new), 3),
                     y = rep(y_new / K_new, 3),
                     lb = ppd_summaries[, "25%"],
                     ub = ppd_summaries[, "75%"],
                     model = rep(models, each = length(y_new)))

## ---- plot-ppd-----------------------------------------------------------
ggplot(df_ppd, aes(x=player, y=y, ymin=lb, ymax=ub)) + 
  geom_linerange(color = "gray60", size = 2) + 
  geom_point(size = 2.5, fill = "skyblue2") +
  facet_grid(. ~ model) +
  labs(x = NULL, y = "batting average") + 
  scale_x_continuous(breaks = NULL) +
  ggtitle(expression(
    atop("Posterior Predictions for Batting Average in Remainder of Season",
         atop("50% posterior predictive intervals (gray bars); observed (blue dots)", ""))))

## ---- event-probabilities, results="hold"--------------------------------
draws_partialpool <- shift_draws(as.matrix(fit_partialpool))
thetas_partialpool <- plogis(draws_partialpool)
colnames(thetas_partialpool) <- as.character(bball$Player)
ability_gt_400 <- thetas_partialpool > 0.4
cat("Pr(theta_n >= 0.400 | y)\n")
colMeans(ability_gt_400)[c(1, 5, 10)]

some_gt_350 <- apply(thetas_partialpool, 1, function(x) max(x) > 0.35)
cat("Pr(at least one theta_n >= 0.350 | y)\n")
mean(some_gt_350)

## ---- echo=FALSE---------------------------------------------------------
thetas_pool <- plogis(as.matrix(fit_pool))
thetas_nopool <- plogis(as.matrix(fit_nopool))
some_gt_350_all <- sapply(list(thetas_pool, thetas_nopool, thetas_partialpool), 
                          function(x) apply(x, 1, max) > 0.35)
chance_gt_350 <- round(100 * colMeans(some_gt_350_all))

## ---- ranking------------------------------------------------------------
reverse_rank <- function(x) 1 + length(x) - rank(x) # so lower rank is better
rank <- apply(thetas_partialpool, 1, reverse_rank)
t(apply(rank, 1, quantile, prob = c(0.1, 0.5, 0.9)))

## ---- plot-ranks---------------------------------------------------------
df_rank <- data.frame(name = rep(bball$Player, each = M), 
                      rank = c(t(rank)))

ggplot(df_rank, aes(rank)) +
  stat_count(width = 0.8) +
  facet_wrap(~ name) +
  scale_x_discrete(limits = c(1, 5, 10, 15)) +
  scale_y_discrete(name = "posterior probability", breaks = c(0, 0.1 * M, 0.2 * M),
                   labels = c("0.0", "0.1", "0.2")) + 
  ggtitle("Rankings for Partial Pooling Model")

## ---- plot-best-player---------------------------------------------------
thetas_nopool <- plogis(as.matrix(fit_nopool))
colnames(thetas_nopool) <- as.character(bball$Player)
rank_nopool <- apply(thetas_nopool, 1, reverse_rank)
is_best_nopool <- rowMeans(rank_nopool == 1)
is_best_partialpool <- rowMeans(rank == 1)

df_is_best <- data.frame(unit = rep(bball$Player, 2), 
                         is_best = c(is_best_partialpool, is_best_nopool), 
                         model = rep(c("partial pooling", "no pooling"), each = N))

ggplot(df_is_best, aes(x=unit, y=is_best)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ model) +
  scale_y_continuous(name = "Pr[player is best]") +
  ggtitle("Who is the Best Player?") + 
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))

## ---- plot-ppc-stats-mean, fig.width=3, fig.height=2.5-------------------
pp_check(fit_nopool, check = "test", test = "mean")

## ---- plot-ppc-stats-----------------------------------------------------
tstat_plots <- function(model, stats) {
  lapply(stats, function(stat) {
    graph <- pp_check(model, check = "test", test = stat, 
                      binwidth = .025, seed = SEED) # optional arguments
    graph + xlab(stat) + theme(legend.position = "none")
  })
}
Tstats <- c("mean", "sd", "min", "max")
ppcs_pool <- tstat_plots(fit_pool, Tstats)
ppcs_nopool <- tstat_plots(fit_nopool, Tstats)
ppcs_partialpool <- tstat_plots(fit_partialpool, Tstats)

library(gridExtra)
grid.arrange(
  arrangeGrob(grobs = ppcs_pool, nrow = 1, left = "Pooling"), 
  arrangeGrob(grobs = ppcs_nopool, nrow = 1, left = "No Pooling"),
  arrangeGrob(grobs = ppcs_partialpool, nrow = 1, left = "Partial Pooling"))

## ---- p-value------------------------------------------------------------
yrep <- posterior_predict(fit_nopool, seed = SEED) # seed is optional
Ty <- sd(y)
Tyrep <- apply(yrep, 1, sd)

# tail-area probability
p <- 1 - mean(Tyrep > Ty)
print(p)

## ---- plot-ppc-y-vs-yrep-------------------------------------------------
pp_check(fit_partialpool, check = "distributions", nreps = 15, 
         overlay = FALSE, binwidth = 0.05) +  # optional arguments
  ggtitle("Model: Partial Pooling")

