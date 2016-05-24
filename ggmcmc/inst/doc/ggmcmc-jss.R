## ----setup, include=FALSE, echo=FALSE, cache=FALSE------------------
library(knitr)
library(Cairo)
options(prompt = "R> ", continue = "+  ", width = 60, useFancyQuotes = FALSE)
options(replace.assign=TRUE, width=70)
opts_chunk$set(fig.path='figure/plot-', fig.align='center', fig.show='hold', pdf=FALSE, dev='cairo_pdf', par=TRUE, out.width='.50\\textwidth', prompt=TRUE, warning=FALSE)
knit_hooks$set(par=function(before, options, envir) {
  if (before && options$fig.show!='none') {
    par(mar=c(4,4,.1,.1), cex.lab=.95, cex.axis=.9,mgp=c(2,.7,0), tcl=-.3)
  }
}, crop=hook_pdfcrop)
knit_theme$set("print")

## ----message=FALSE, warning=FALSE, eval=FALSE-----------------------
#  library("ggmcmc")
#  data("radon")
#  s.radon.short <- radon$s.radon.short

## ----ggs------------------------------------------------------------
S <- ggs(s.radon.short)

## ----ggs_contents---------------------------------------------------
S

## ----echo=TRUE, eval=FALSE------------------------------------------
#  str(S)

## ----echo=FALSE, eval=TRUE------------------------------------------
str(S, width=70, strict.width="cut")

## ----eval=FALSE-----------------------------------------------------
#  ggmcmc(S)

## ----eval=FALSE-----------------------------------------------------
#  ggmcmc(S, file = "model_simple-diag.pdf", param_page = 2)

## ----eval=FALSE-----------------------------------------------------
#  ggmcmc(S, plot = c("density", "running", "caterpillar"))

## ----eval=FALSE-----------------------------------------------------
#  ggs_histogram(S)

## ----histogram_density_part, fig.width=12, fig.height=8, echo=FALSE, out.width='0.8\\textwidth', fig.cap='left) Histogram with the distribution of the posterior values combining all chains by parameter; middle) Density plots.; right) Density plots comparing the whole chain (black) with only the last part (green).', message=FALSE, warning=FALSE----
require(gridExtra)
f1 <- ggs_histogram(S) + ggtitle("ggs_histogram()")
f2 <- ggs_density(S) + ggtitle("ggs_density()")
f3 <- ggs_compare_partial(S) + ggtitle("ggs_compare_partial()")
grid.arrange(f1, f2, f3, ncol = 3)

## ----eval=FALSE-----------------------------------------------------
#  ggs_density(S)

## ----eval=FALSE-----------------------------------------------------
#  ggs_compare_partial(S)

## ----eval=FALSE-----------------------------------------------------
#  ggs_traceplot(S)

## ----eval=FALSE-----------------------------------------------------
#  ggs_running(S)

## ----traceplots_running_means_autocorrelation, fig.width=12, fig.height=8, echo=FALSE, out.width='0.8\\textwidth', fig.cap='left) Traceplots with the time series of the chains; middle) Running means; right) Autocorrelation plots.', message=FALSE, warning=FALSE----
require(gridExtra)
f1 <- ggs_traceplot(S) + ggtitle("ggs_traceplot()")
f2 <- ggs_running(S) + ggtitle("ggs_running()")
f3 <- ggs_autocorrelation(S) + ggtitle("ggs_autocorrelation()")
grid.arrange(f1, f2, f3, ncol = 3)

## ----eval=FALSE-----------------------------------------------------
#  ggs_autocorrelation(S)

## ----eval=FALSE-----------------------------------------------------
#  ggs_crosscorrelation(S)

## ----crosscorrelation, fig.width=10, fig.height=4, out.width='0.6\\textwidth', fig.cap='Tile plots with the crosscorrelations of the parameters in a relative scale (left) or absolute scale (right)', echo=FALSE----
f1 <- ggs_crosscorrelation(S) + ggtitle("Absolute scale")
f2 <- ggs_crosscorrelation(S, absolute_scale = FALSE) + ggtitle("Relative scale")
grid.arrange(f1, f2, ncol = 2)

## ----eval=FALSE-----------------------------------------------------
#  ggs_crosscorrelation(S, absolute_scale = FALSE)

## ----Rhat, eval=FALSE-----------------------------------------------
#  ggs_Rhat(S)

## ----eval=FALSE-----------------------------------------------------
#  ggs_geweke(S)

## ----Rhat_geweke, fig.width=10, fig.height=4, out.width='.8\\textwidth', fig.cap='left) Dotplot with the Potential Scale Reduction Factor. right) Dotplot with the Geweke z-scores. For this model there is no evidence of non-convergence, as $\\hat{R}$ are very close to 1 and z-scores are between -2 and 2.', echo=FALSE----
f1 <- ggs_Rhat(S)
f2 <- ggs_geweke(S)
grid.arrange(f1, f2, ncol = 2)

## -------------------------------------------------------------------
S.full <- ggs(radon$s.radon)

## ----eval=FALSE-----------------------------------------------------
#  ggs_density(S.full, family = "sigma")

## ----eval=FALSE-----------------------------------------------------
#  P <- data.frame(
#    Parameter = c("sigma.alpha", "sigma.beta", "sigma.y"),
#    Label = c("Intercept (sd)", "Covariate (sd)", "Outcome (sd)"))
#  ggs_density(ggs(radon$s.radon, par_labels = P, family = "sigma"))

## ----density_family_par_labels, fig.width=8, fig.height=6, out.width='.80\\textwidth', fig.cap='left) Density plot with a restricted number of parameters controlled by the argument \\code{family}; right) Density plot with the parameters having different labels using the argument \\code{par\\_labels}. Notice that the first two parameters in each column show signals of non-convergence.', echo=FALSE----
f1 <- ggs_density(S.full, family = "sigma") + ggtitle("family")
P <- data.frame(
  Parameter = c("sigma.alpha", "sigma.beta", "sigma.y"),
  Label = c("Intercept (sd)", "Covariate (sd)", "Outcome (sd)"))
f2 <- ggs_density(ggs(radon$s.radon, par_labels = P, family = "sigma")) + ggtitle("par_labels")
grid.arrange(f1, f2, ncol = 2)

## ----caterpillar_preparation----------------------------------------
L.radon.intercepts <- data.frame(
  Parameter = paste("alpha[", radon$counties$id.county, "]", sep = ""),
  Label = radon$counties$County)
head(L.radon.intercepts)
S.full <- ggs(radon$s.radon,
  par_labels = L.radon.intercepts, family = "^alpha")

## ----eval=FALSE-----------------------------------------------------
#  ggs_caterpillar(S.full)

## ----eval=FALSE-----------------------------------------------------
#  Z <- data.frame(
#    Parameter = paste("alpha[", radon$counties$id.county, "]", sep = ""),
#    value = radon$counties$uranium)
#  ggs_caterpillar(ggs(radon$s.radon, family = "^alpha"),
#    X = Z, horizontal = FALSE)

## ----caterpillar, out.width='.80\\textwidth', fig.width=12, fig.height=8, fig.cap='Caterpillar plot for the varying intercepts. Dots represent medians, and thick and thin lines represent  90 and the 95\\% of the Highest Posterior Density regions, respectively. left) Plain caterpillar plot with labels of the parameters, showing Clay and Lake as the counties with highest and lowest $\\alpha$ intercept, respectively; right) Caterpillar plot against continuous values (uranium levels) using the argument \\code{par\\_labels}, suggesting a tendency of higher $\\alpha$ intercepts with increasing uranium levels.', echo=FALSE----
f1 <- ggs_caterpillar(S.full)
Z <- data.frame(
  Parameter = paste("alpha[", radon$counties$id.county, "]", sep = ""),
  value = radon$counties$uranium)
f2 <- ggs_caterpillar(ggs(radon$s.radon, family = "^alpha"), X = Z, horizontal = FALSE)
grid.arrange(f1, f2, ncol = 2)

## ----sample_mu------------------------------------------------------
data("linear")
S.y.rep <- ggs(s.y.rep)
y.observed <- y

## ----eval=FALSE-----------------------------------------------------
#  ggs_ppmean(S.y.rep, outcome = y.observed)

## ----eval=FALSE-----------------------------------------------------
#  ggs_ppsd(S.y.rep, outcome = y.observed)

## ----ppmean_ppsd, out.width='.60\\textwidth', fig.width=8, fig.height=4, fig.cap='Histograms of the posterior predictive distributions of the mean (left) and standard deviation (right) of the replicated datasets, against the vertical lines showing the mean and standard deviations of the data, respectively.', echo=FALSE----
f1 <- ggs_ppmean(S.y.rep, outcome = y.observed)
f2 <- ggs_ppsd(S.y.rep, outcome = y.observed)
grid.arrange(f1, f2, ncol = 2)

## -------------------------------------------------------------------
data("binary")
S.binary <- ggs(s.binary, family="mu")

## ----roc, out.width='0.5\\textwidth', fig.width=6, fig.height=5, fig.cap='ROC (receiver operating characteristic) curve.'----
ggs_rocplot(S.binary, outcome = y.binary)

## ----separation, out.width='0.7\\textwidth', fig.width=6, fig.height=2, fig.cap='Separation plot.'----
ggs_separation(S.binary, outcome = y.binary)

## ----separation_minimalist, out.width='0.3\\textwidth', fig.width=5, fig.height=1, fig.cap='Separation plot. Minimalist version to be used inline.', include=FALSE, echo=TRUE, crop=TRUE----
ggs_separation(S.binary, outcome = y.binary, minimal = TRUE)

## ----eval=FALSE, echo=TRUE------------------------------------------
#  ggs_separation(S.binary, outcome = y.binary, minimal = TRUE)

## ----combination_aesthetics, out.width='0.6\\textwidth', fig.width=10, fig.cap='Combination of the aestheticaly-driven options that complement \\pkg{ggplot2}: \\pkg{ggthemes} and \\pkg{gridExtra}.'----
library("gridExtra")
library("ggthemes")
f1 <- ggs_traceplot(ggs(s, family = "^beta\\[[1234]\\]")) +
  theme_tufte()
f2 <- ggs_density(ggs(s, family = "^beta\\[[1234]\\]")) +
  theme_solarized(light = FALSE)
grid.arrange(f1, f2, ncol = 2, nrow = 1)

## -------------------------------------------------------------------
ci.median <- ci(ggs(radon$s.radon, family = "^alpha|^beta")) %>%
  dplyr::select(Parameter, median)

## -------------------------------------------------------------------
L.radon <- data.frame(
  Parameter = c(
    paste("alpha[", radon$counties$id.county, "]", sep = ""),
    paste("beta[", radon$counties$id.county, "]", sep = "")),
  Label = rep(radon$counties$County, 2),
  Uranium = rep(radon$counties$uranium, 2),
  Location = rep(radon$counties$ns.location, 2),
  Coefficient = gl(2, length(radon$counties$id.county),
    labels = c("Intercept", "Slope")))

## ----message=FALSE--------------------------------------------------
radon.median <- left_join(ci.median, L.radon) %>%
  dplyr::select(Label, Coefficient, median) %>%
  tidyr::spread(Coefficient, median)
head(radon.median)

## ----ggplot_intercepts_slopes, fig.width=4, fig.height=4, out.width='0.25\\textwidth', fig.cap='Scatterplot with the medians of the posteriors of intercepts and slopes for the radon data example.', fig.pos='htbp', crop=TRUE----
ggplot(radon.median, aes(x = Intercept, y = Slope)) + geom_point()

## ----extension_facets_aes, fig.width=12, out.width='0.8\\textwidth', fig.cap='Caterpillar plot of the varying intercepts faceted by North/South location and using county\'s uranium level as color indicator.'----
ggs_caterpillar(ggs(radon$s.radon, par_labels = L.radon, family = "^alpha")) +
  facet_wrap(~ Location, scales = "free") +
  aes(color = Uranium)

