## ----init, echo=FALSE, message=FALSE-------------------------------------
# Set on_cran to FALSE to get tikzDevice plots and citation
on_cran <- Sys.getenv("USER") != "muelleki"

library(wrswoR)
import::from(tidyr, spread, gather, gather_, .library = .libPaths())
import::from(plyr, ldply, .library = .libPaths())
import::from(dplyr, filter, mutate, select, group_by, summarize, ungroup, tbl_df,
             rename, arrange, mutate_each, funs, .library = .libPaths())
import::from(cluster, daisy, .library = .libPaths())
import::from(magrittr, "%>%", .library = .libPaths())
library(ggplot2)
set.seed(20150710L)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass{jss}\\usepackage{siunitx}\\usepackage{latexsym}")
knitr::opts_chunk$set(cache=TRUE)
if (on_cran) {
  knitr::opts_chunk$set(dev="png")
} else {
  knitr::opts_chunk$set(dev="tikz")
}
knitcitations::cite_options("pandoc")

knit_print.function <- function(x, options) {
  dput(x)
}

## ---- echo=FALSE---------------------------------------------------------
sample_int_R

## ---- echo=FALSE---------------------------------------------------------
sample_int_rank

## ----timing-base, echo=FALSE---------------------------------------------
ggplot_base <- list(
  theme_bw(11)
)

ALGOS <- c("R", "rej", "rank", "crank", "expj")

REMOVE_PROB <- c("linear_double", "linear_half")

combine_prob_mix <- function(prob, mix) {
  prob <- as.character(prob)
  mix <- as.character(mix)
  
  ifelse(prob == "Uniform", prob, paste(prob, mix))
}

timings <- 
  wrswoR.benchmark::timings %>%
  tbl_df %>%
  filter(expr %in% ALGOS) %>%
  filter(!(prob %in% REMOVE_PROB)) %>%
  mutate(expr = factor(expr, levels = ALGOS)) %>%
  mutate(prob = factor(
    prob, levels = c("uniform", "linear", "exp"),
    labels = c("uniform", "linear", "geometric"))) %>% 
  mutate(mix = factor(
    mix, levels = c("asc", "desc", "shuffle"),
    labels = c("$\\nearrow$", "$\\searrow$", "$\\leadsto$"))) %>%
  mutate(prob_mix = kimisc::ofactor(combine_prob_mix(prob, mix)))

BASE <- 1.007
N <- max(timings$n)

## ----break-even-base, echo=FALSE, dependson = "timing-base"--------------
break_even <- 
  wrswoR.benchmark::break_even %>%
  tbl_df %>%
  filter(expr %in% ALGOS) %>%
  filter(!(prob %in% REMOVE_PROB)) %>%
  mutate(expr = factor(expr, levels = ALGOS)) %>%
  mutate(prob = factor(
    prob, levels = c("uniform", "linear", "exp"),
    labels = c("uniform", "linear", "geometric"))) %>% 
  mutate(mix = factor(
    mix, levels = c("asc", "desc", "shuffle"),
    labels = c("$\\nearrow$", "$\\searrow$", "$\\leadsto$"))) %>%
  mutate(prob_mix = kimisc::ofactor(combine_prob_mix(prob, mix)))

## ----ggplot-base, echo=FALSE---------------------------------------------
ggplot_perf_base <-
  ggplot_base %>%
  c(list(
    theme_bw(11),
    scale_color_discrete(name = "Algorithm")
  ))

ggplot_time_base <-
  ggplot_perf_base %>%
  c(list(ylab("Run time (s)")))

ggplot_time_per_item_base <-
  ggplot_base %>%
  c(list(ylab("Run time per element (s)")))

math <- function(x) paste0("$", x, "$")

pow10 <- function(x, in_math = FALSE) {
  x <- ifelse(x %in% 0:1, x, paste0("10^{", log10(x), "}"))
  if (!in_math)
    x <- math(x)
  x
}

pow2 <- function(x) {
  ifelse(x %in% 0:1, x, paste0("$2^{", log2(x), "}$"))
}

percent <- function(x) {
  paste0("\\ensuremath{", round(x * 100, 5), "\\,\\%}")
}

format_expr <- function(x) {
  ifelse(x == "R", "\\proglang{R}", paste0("\\emph{", x, "}"))
}

relabel_expr <- function(x) {
  factor(x, levels = ALGOS, labels = format_expr(ALGOS))
}

arrange_by_func <- . %>%
  mutate(ofunc = factor(func, levels = ALGOS)) %>%
  arrange(ofunc) %>%
  select(-ofunc)

## ----run-time-log, echo=FALSE, fig.height=8, fig.cap="Median run times", dependson=c("timing-base", "ggplot-base")----
timings %>%
  group_by(n, expr, prob_mix, r) %>%
  summarize(median = median(time)) %>%
  ungroup %>%
  mutate(r = paste0("$r = ", r, "$")) %>%
  mutate(expr = relabel_expr(expr)) %>%
  ggplot(aes(x=n, y=median * 1e-9, color=expr)) +
  ggplot_time_base +
  geom_line() +
  scale_x_log10(name = "$n$", breaks = 10 ** c(2,4,6), labels = pow10) +
  scale_y_log10(labels = pow10) +
  facet_grid(prob_mix~r) +
  theme(legend.position="bottom")

## ----run-time-crank-expj, echo=FALSE, fig.cap="Comparison of \\emph{crank} and \\emph{expj} run times", dependson=c("timing-base", "ggplot-base"), fig.height = 3.8----
timings %>%
  filter(expr %in% c("crank", "expj")) %>%
  group_by(n, expr, prob_mix, r) %>%
  summarize(median = median(time)) %>%
  ungroup %>%
  spread(expr, median) %>%
  mutate(crank_vs_expj = crank / expj) %>%
  mutate(r = paste0("$r = ", r, "$")) %>%
  ggplot(aes(x=prob_mix, y=crank_vs_expj, color=n)) +
  scale_x_discrete(name = "Weights distribution") +
  scale_y_log10(name = "Ratio of median run times\n(\\emph{crank} vs. \\emph{expj})", breaks = 2 ** (-1:3)) +
  ggplot_base +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_point() +
  scale_color_continuous(name = "$n$", trans="log", breaks = as.integer(10 ** c(2, 4, 6)), labels = pow10) +
  facet_wrap(~r)

## ----run-time-break-even, echo=FALSE, fig.cap="Comparison of \\proglang{R} and \\emph{expj} run times for linear ascending weights", dependson=c("break-even-base", "ggplot-base")----
break_even %>%
  filter(expr %in% c("R", "expj")) %>%
  filter(prob_mix == "linear $\\nearrow$") %>%
  filter(r %in% c(0.01, 0.1, 1)) %>%
  group_by(n, expr, prob, r) %>%
  summarize(median = median(time)) %>%
  ungroup %>%
  spread(expr, median) %>%
  mutate(R_vs_expj = R / expj) %>%
  ggplot(aes(x=n, y=R_vs_expj, color=factor(r, levels = c("1", "0.1", "0.01")))) +
  ggplot_base +
  geom_line() +
  scale_x_log10(name = "$n$", breaks = c(10, 30, 100, 300, 1000, 3000)) +
  scale_y_log10(name = "Ratio of median run times\n(\\proglang{R} vs. \\emph{expj})", breaks = 2 ** (-1:5)) +
  scale_color_discrete(name = "$r$")

## ----break-even-point, echo = FALSE--------------------------------------
break_even_point <-
  break_even %>%
  filter(expr %in% c("R", "expj")) %>%
  group_by(n, expr, prob, mix, r) %>%
  summarize(median = median(time)) %>%
  ungroup %>%
  spread(expr, median) %>%
  mutate(rate = R / expj)

break_even_point_max <-
  break_even_point %>%
  .[["rate"]] %>%
  min

break_even_point_min_n <-
  break_even_point %>%
  group_by(prob, mix, r) %>%
  summarize(smallest_n_where_expj_is_better = n[tail(which(rate < 1), 1L) + 1L]) %>%
  ungroup

## ----point-break-even, echo=FALSE, fig.cap="Break-even point of \\proglang{R} and \\emph{expj} run times", dependson=c("break-even-base", "ggplot-base", "break-even-point")----
break_even_point_min_n %>% 
  ggplot(aes(x=r, y=smallest_n_where_expj_is_better, color=prob, linetype=mix)) +
  ggplot_base +
  geom_line() +
  scale_x_log10(name = "$r$") +
  scale_y_continuous(name = "$n$", breaks = 1:5 * 100, limits = c(NA, 500)) +
  scale_linetype(name = "Order") +
  scale_color_discrete(name = "Distribution")

## ----def-plot-p-values, echo=FALSE---------------------------------------
n. <- 7
s. <- 4
skew. <- 1.0025
alpha. <- 1.08
N. <- 2 ** 22

p_values_true <- wrswoR.benchmark::p_values_7 %>% filter(N == N., s == s., skew == 1, func == "crank")

p_values_false <- wrswoR.benchmark::p_values_7 %>% filter(N == N., s == s., skew == skew.)

## ----correctness-true, echo=FALSE, dependson="def-plot-p-values", fig.cap="Schweder plot for p-values resulting from comparing the \\emph{crank} and \\proglang{R} implementations", fig.height = 3.7----
p_values_true %>%
  rename(`p-value` = p_value) %>%
  arrange(`p-value`) %>%
  mutate(`Rank` = seq_along(`p-value`)) %>% 
  mutate_each(funs(factor), i, j) %>% 
  ggplot(aes(y = `Rank`, x = `p-value`, color = i, shape = j)) +
  ggplot_base +
  geom_abline(slope = nrow(p_values_true), intercept = 0, linetype = 3) +
  geom_point() +
  coord_fixed(ratio = 1 / nrow(p_values_true)) +
  scale_color_discrete(name = "$i$") +
  scale_shape_discrete(name = "$j$") +
  list()

## ----correctness-false, echo=FALSE, dependson="def-plot-p-values", fig.cap="Schweder plot for p-values resulting from comparing the \\proglang{R} implementation with a skewed version of itself", fig.height = 3.7----
p_values_false %>%
  rename(`p-value` = p_value) %>%
  arrange(`p-value`) %>%
  mutate(`Rank` = seq_along(`p-value`)) %>% 
  mutate_each(funs(factor), i, j) %>% 
  ggplot(aes(y = `Rank`, x = `p-value`, color = i, shape = j)) +
  ggplot_base +
  geom_abline(slope = nrow(p_values_false), intercept = 0, linetype = 3) +
  geom_point() +
  coord_fixed(ratio = 1 / nrow(p_values_false)) +
  scale_color_discrete(name = "$i$") +
  scale_shape_discrete(name = "$j$") +
  list()

## ----comprehensive-base, echo = FALSE, dependson="ggplot-base"-----------
p_value_breakpoints <- c(0, 1e-100, 1e-4, 1e-2, 1e-1, 9e-1, 1)
  
cut_p <- function(x) {
  stopifnot(x >= 0)
  stopifnot(x <= 1)
  x_orig <- x
  x <- kimisc::cut_format(x, p_value_breakpoints, include.lowest = TRUE, format_fun = function(x) ifelse(x >= 0.01, x, pow10(x, in_math = TRUE)), paren = c("$\\left(", "$\\left[", "\\right)$", "\\right]$"))
  x <- factor(as.character(x), levels = rev(levels(x)))
  x
}

comprehensive_base <- function(by) list(
  geom_raster(),
  scale_x_log10(name = "$N$", expand = c(0, 0), breaks = 2 ** seq.int(24, 8, by = -by), labels = pow2),
  scale_fill_brewer(name = "p-value", palette = "YlOrRd", drop = FALSE),
  NULL
)

## ----comprehensive, echo=FALSE, fig.cap="Combining p-values for a comprehensive test", dependson = "comprehensive-base"----
wrswoR.benchmark::p_values_agg_agg %>%
  ungroup %>%
  arrange_by_func %>%
  mutate(facet = kimisc::ofactor(ifelse(skew == 1, paste0("skew ", percent(0)), "\\proglang{R}"))) %>%
  ggplot(aes(x=N, y=kimisc::ofactor(
    ifelse(skew == 1, format_expr(func),
           percent(skew - 1))),
    fill=cut_p(p_value))) +
  ggplot_base +
  comprehensive_base(2) +
  scale_y_discrete(name = "Skew | Implementation", expand = c(0, 0)) +
  facet_wrap(~facet, ncol = 1, scales = "free_y")

## ----comprehensive-true, echo=FALSE, fig.cap="Combined p-values for different values of $n$ and $N$, resulting from comparing each new code to the stock implementation", dependson = "comprehensive-base"----
wrswoR.benchmark::p_values_agg %>%
  ungroup %>% 
  arrange_by_func %>%
  filter(skew == 1 & func != "R") %>%
  mutate(func = kimisc::ofactor(format_expr(func))) %>%
  ggplot(aes(x=N, y=n, fill=cut_p(p_value))) +
  ggplot_base +
  comprehensive_base(4) +
  scale_y_continuous(name = "$n$", expand = c(0, 0)) +
  facet_wrap(~func, nrow = 1)

## ----comprehensive-false, echo=FALSE, fig.cap="Combined p-values for different values of $n$ and $N$, resulting from comparing the stock implementation to a skewed version of itself", dependson = "comprehensive-base"----
wrswoR.benchmark::p_values_agg %>%
  ungroup %>%
  filter(func == "R") %>%
  mutate(skew = kimisc::ofactor(paste0("skew ", percent(skew - 1)))) %>%
  ggplot(aes(x=N, y=n, fill=cut_p(p_value))) +
  ggplot_base +
  comprehensive_base(4) +
  scale_y_continuous(name = "$n$", expand = c(0, 0)) +
  facet_wrap(~skew, ncol = 4)

## ----echo=FALSE, message=FALSE, cache=FALSE------------------------------
# Work around encoding issue on Windows
if (!on_cran) {
  knitcitations::write.bibtex()
} else {
  writeLines(character(), "knitcitations.bib")
}

