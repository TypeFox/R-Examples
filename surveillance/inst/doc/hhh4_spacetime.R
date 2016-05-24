## ----include = FALSE---------------------------------------------------------------
## load the "cool" package
library("surveillance")

## Compute everything or fetch cached results?
message("Doing computations: ",
        COMPUTE <- !file.exists("hhh4_spacetime-cache.RData"))
if (!COMPUTE) load("hhh4_spacetime-cache.RData", verbose = TRUE)

## ----measlesWeserEms_components, echo=FALSE----------------------------------------
## extract components from measlesWeserEms to reconstruct
data("measlesWeserEms")
counts <- observed(measlesWeserEms)
map <- measlesWeserEms@map
populationFrac <- measlesWeserEms@populationFrac

## ----measlesWeserEms_neighbourhood-------------------------------------------------
weserems_nbOrder <- nbOrder(poly2adjmat(map), maxlag = 10)

## ----measlesWeserEms_construct-----------------------------------------------------
measlesWeserEms <- sts(counts, start = c(2001, 1), frequency = 52,
  population = populationFrac, neighbourhood = weserems_nbOrder, map = map)

## ----measlesWeserEms, fig.cap="Measles infections in the Weser-Ems region, 2001--2002.", fig.subcap=c("Time series of weekly counts.","Disease incidence (per 100\\,000 inhabitants)."), fig.width=5, fig.height=5, out.width="0.47\\linewidth", fig.pos="htb"----
plot(measlesWeserEms, type = observed ~ time)
plot(measlesWeserEms, type = observed ~ unit,
  population = measlesWeserEms@map$POPULATION / 100000,
  labels = list(font = 2), colorkey = list(space = "right"),
  sp.layout = layout.scalebar(measlesWeserEms@map, corner = c(0.05, 0.05),
    scale = 50, labels = c("0", "50 km"), height = 0.03))

## ----measlesWeserEms15, fig.cap=paste("Count time series of the", sum(colSums(observed(measlesWeserEms))>0), "affected districts."), out.width="\\linewidth", fig.width=10, fig.height=6, fig.pos="!h"----
plot(measlesWeserEms, units = which(colSums(observed(measlesWeserEms)) > 0))

## ----measlesWeserEms_animation, eval=FALSE-----------------------------------------
#  animation::saveHTML(
#    animate(measlesWeserEms, tps = 1:52, total.args = list()),
#    title = "Evolution of the measles epidemic in the Weser-Ems region, 2001",
#    ani.width = 500, ani.height = 600)

## ----echo=FALSE, eval=FALSE--------------------------------------------------------
#  ## to perform the following analysis using biweekly aggregated measles counts:
#  measlesWeserEms <- aggregate(measlesWeserEms, by = "time", nfreq = 26)

## ----measlesModel_basic------------------------------------------------------------
measlesModel_basic <- list(
  end = list(f = addSeason2formula(~1 + t, period = measlesWeserEms@freq),
             offset = population(measlesWeserEms)),
  ar = list(f = ~1),
  ne = list(f = ~1, weights = neighbourhood(measlesWeserEms) == 1),
  family = "NegBin1")

## ----measlesFit_basic--------------------------------------------------------------
measlesFit_basic <- hhh4(stsObj = measlesWeserEms, control = measlesModel_basic)

## ----measlesFit_basic_summary------------------------------------------------------
summary(measlesFit_basic, idx2Exp = TRUE, amplitudeShift = TRUE, maxEV = TRUE)

## ----measlesFit_basic_endseason, fig.width=6, fig.height=2.5, out.width=".6\\linewidth", fig.cap="Estimated multiplicative effect of seasonality on the endemic mean.", fig.pos="ht"----
plot(measlesFit_basic, type = "season", components = "end", main = "")

## ----measlesFitted_basic, fig.cap="Fitted components in the initial model \\code{measlesFit\\_basic} for the six districts with more than 20 cases. Dots are only drawn for positive weekly counts.", out.width="\\linewidth", fig.pos="htb"----
districts2plot <- which(colSums(observed(measlesWeserEms)) > 20)
plot(measlesFit_basic, type = "fitted", units = districts2plot, hide0s = TRUE)

## ----------------------------------------------------------------------------------
confint(measlesFit_basic, parm = "overdisp")

## ----measlesFit_basic_Poisson------------------------------------------------------
AIC(measlesFit_basic, update(measlesFit_basic, family = "Poisson"))

## ----Sprop-------------------------------------------------------------------------
Sprop <- matrix(1 - measlesWeserEms@map@data$vacc1.2004,
  nrow = nrow(measlesWeserEms), ncol = ncol(measlesWeserEms), byrow = TRUE)
summary(Sprop[1, ])

## ----SmodelGrid--------------------------------------------------------------------
Soptions <- c("unchanged", "Soffset", "Scovar")
SmodelGrid <- expand.grid(end = Soptions, ar = Soptions)
row.names(SmodelGrid) <- do.call("paste", c(SmodelGrid, list(sep = "|")))

## ----measlesFits_vacc, eval=COMPUTE------------------------------------------------
#  measlesFits_vacc <- apply(X = SmodelGrid, MARGIN = 1, FUN = function (options) {
#    updatecomp <- function (comp, option) switch(option, "unchanged" = list(),
#      "Soffset" = list(offset = comp$offset * Sprop),
#      "Scovar" = list(f = update(comp$f, ~. + log(Sprop))))
#    update(measlesFit_basic,
#      end = updatecomp(measlesFit_basic$control$end, options[1]),
#      ar = updatecomp(measlesFit_basic$control$ar, options[2]),
#      data = list(Sprop = Sprop))
#    })

## ----aics_vacc, eval=COMPUTE-------------------------------------------------------
#  aics_vacc <- do.call(AIC, lapply(names(measlesFits_vacc), as.name),
#    envir = as.environment(measlesFits_vacc))

## ----------------------------------------------------------------------------------
aics_vacc[order(aics_vacc[, "AIC"]), ]

## ----measlesFit_vacc---------------------------------------------------------------
measlesFit_vacc <- update(measlesFit_basic,
  end = list(f = update(formula(measlesFit_basic)$end, ~. + log(Sprop))),
  data = list(Sprop = Sprop))
coef(measlesFit_vacc, se = TRUE)["end.log(Sprop)", ]

## ----------------------------------------------------------------------------------
2^cbind("Estimate" = coef(measlesFit_vacc),
  confint(measlesFit_vacc))["end.log(Sprop)",]

## ----measlesFit_nepop--------------------------------------------------------------
measlesFit_nepop <- update(measlesFit_vacc,
  ne = list(f = ~log(pop)), data = list(pop = population(measlesWeserEms)))

## ----------------------------------------------------------------------------------
cbind("Estimate" = coef(measlesFit_nepop),
  confint(measlesFit_nepop))["ne.log(pop)",]

## ----measlesFit_powerlaw-----------------------------------------------------------
measlesFit_powerlaw <- update(measlesFit_nepop,
  ne = list(weights = W_powerlaw(maxlag = 5)))

## ----------------------------------------------------------------------------------
cbind("Estimate" = coef(measlesFit_powerlaw),
  confint(measlesFit_powerlaw))["neweights.d",]

## ----measlesFit_np-----------------------------------------------------------------
measlesFit_np2 <- update(measlesFit_nepop,
  ne = list(weights = W_np(maxlag = 2)))

## ----measlesFit_neweights, fig.width=5, fig.height=3.5, fig.cap="Estimated weights as a function of adjacency order.", out.width="0.47\\linewidth", fig.subcap=c("Normalized power-law weights.", "Non-normalized weights with 95\\% CIs."), echo=c(1,4)----
library("lattice")
trellis.par.set("reference.line", list(lwd=3, col="gray"))
trellis.par.set("fontsize", list(text=14))
plot(measlesFit_powerlaw, type = "neweights", plotter = stripplot,
  panel = function (...) {panel.stripplot(...); panel.average(...)},
  jitter.data = TRUE, xlab = expression(o[ji]), ylab = expression(w[ji]))
## non-normalized weights (power law and unconstrained second-order weight)
local({
    colPL <- "#0080ff"
    ogrid <- 1:5
    par(mar=c(3.6,4,2.2,2), mgp=c(2.1,0.8,0))
    plot(ogrid, ogrid^-coef(measlesFit_powerlaw)["neweights.d"], col=colPL,
         xlab="Adjacency order", ylab="Non-normalized weight", type="b", lwd=2)
    matlines(t(sapply(ogrid, function (x)
                      x^-confint(measlesFit_powerlaw, parm="neweights.d"))),
             type="l", lty=2, col=colPL)
    w2 <- exp(c(coef(measlesFit_np2)["neweights.d"],
                 confint(measlesFit_np2, parm="neweights.d")))
    lines(ogrid, c(1,w2[1],0,0,0), type="b", pch=19, lwd=2)
    arrows(x0=2, y0=w2[2], y1=w2[3], length=0.1, angle=90, code=3, lty=2)
    legend("topright", col=c(colPL, 1), pch=c(1,19), lwd=2, bty="n", inset=0.1, y.intersp=1.5,
           legend=c("Power-law model", "Second-order model"))
})

## ----------------------------------------------------------------------------------
AIC(measlesFit_nepop, measlesFit_powerlaw, measlesFit_np2)

## ----measlesFit_ri, results="hide"-------------------------------------------------
measlesFit_ri <- update(measlesFit_powerlaw,
  end = list(f = update(formula(measlesFit_powerlaw)$end, ~. + ri() - 1)),
  ar  = list(f = update(formula(measlesFit_powerlaw)$ar,  ~. + ri() - 1)),
  ne  = list(f = update(formula(measlesFit_powerlaw)$ne,  ~. + ri() - 1)))

## ----measlesFit_ri_summary_echo, eval=FALSE----------------------------------------
#  summary(measlesFit_ri, amplitudeShift = TRUE, maxEV = TRUE)

## ----------------------------------------------------------------------------------
head(ranef(measlesFit_ri, tomatrix = TRUE), n = 3)

## ----measlesFit_ri_map, out.width="0.31\\linewidth", fig.width=3.5, fig.height=3.7, fig.pos="htb", fig.cap="Maps of the estimated random intercepts.", fig.subcap=c("Autoregressive $\\alpha_i^{(\\lambda)}$", "Spatio-temporal $\\alpha_i^{(\\phi)}$", "Endemic $\\alpha_i^{(\\nu)}$"), echo=-1----
stopifnot(ranef(measlesFit_ri) > -1.6, ranef(measlesFit_ri) < 1.6)
for (comp in c("ar", "ne", "end")) {
  print(plot(measlesFit_ri, type = "ri", component = comp,
    col.regions = cm.colors(14), labels = list(cex = 0.6),
    at = seq(-1.6, 1.6, length.out = 15)))
}

## ----measlesFitted_ri, out.width="0.93\\linewidth", fig.pos="htb", fig.cap="Fitted components in the random effects model \\code{measlesFit\\_ri} for the six districts with more than 20 cases. Compare to Figure~\\ref{fig:measlesFitted_basic}."----
plot(measlesFit_ri, type = "fitted", units = districts2plot, hide0s = TRUE)

## ----measlesFitted_maps, fig.cap="Maps of the fitted component proportions averaged over all weeks.", fig.pos="hbt", fig.width=10, fig.height=3.7, out.width="0.93\\linewidth"----
plot(measlesFit_ri, type = "maps",
  which = c("epi.own", "epi.neighbours", "endemic"),
  prop = TRUE, labels = list(cex = 0.6))

## ----measlesPreds1-----------------------------------------------------------------
tp <- c(65, 77)
models2compare <- paste0("measlesFit_", c("basic", "powerlaw", "ri"))
measlesPreds1 <- lapply(mget(models2compare), oneStepAhead,
  tp = tp, type = "final")

## ----echo=FALSE--------------------------------------------------------------------
stopifnot(all.equal(measlesPreds1$measlesFit_powerlaw$pred,
                    fitted(measlesFit_powerlaw)[tp[1]:tp[2],],
                    check.attributes = FALSE))

## ----echo=FALSE--------------------------------------------------------------------
stopifnot(all.equal(
    measlesFit_powerlaw$loglikelihood,
    -sum(scores(oneStepAhead(measlesFit_powerlaw, tp = 1, type = "final"),
                which = "logs", individual = TRUE))))

## ----measlesScores1----------------------------------------------------------------
SCORES <- c("logs", "rps", "dss", "ses")
measlesScores1 <- lapply(measlesPreds1, scores, which = SCORES, individual = TRUE)
t(sapply(measlesScores1, colMeans, dims = 2))

## ----measlesPreds2, eval=COMPUTE---------------------------------------------------
#  measlesPreds2 <- lapply(mget(models2compare), oneStepAhead,
#    tp = tp, type = "rolling", which.start = "final")

## ----measlesScores2----------------------------------------------------------------
measlesScores2 <- lapply(measlesPreds2, scores, which = SCORES, individual = TRUE)
t(sapply(measlesScores2, colMeans, dims = 2))

## ----measlesScores_test------------------------------------------------------------
set.seed(321)
sapply(SCORES, function (score) permutationTest(
  measlesScores2$measlesFit_ri[, , score],
  measlesScores2$measlesFit_basic[, , score],
  nPermutation = 999))

## ----measlesPreds2_calibrationTest_echo, eval=FALSE--------------------------------
#  calibrationTest(measlesPreds2[["measlesFit_ri"]], which = "rps")

## ----measlesPreds2_pit, fig.width=8, fig.height=3, out.width="0.93\\linewidth", fig.cap="PIT histograms of competing models to check calibration of the one-week-ahead predictions during the second quarter of 2002.", echo=-1, fig.pos="hbt"----
par(mfrow = sort(n2mfrow(length(measlesPreds2))), mar = c(4.5,4.5,2,1))
for (m in models2compare)
  pit(measlesPreds2[[m]], plot = list(ylim = c(0, 1.25), main = m))

## ----measlesFit_ri_simulate--------------------------------------------------------
(y.start <- observed(measlesWeserEms)[52, ])
measlesSim <- simulate(measlesFit_ri,
  nsim = 100, seed = 1, subset = 53:104, y.start = y.start)

## ----------------------------------------------------------------------------------
summary(colSums(measlesSim, dims = 2))

## ----measlesSim_plot_time, fig.cap="Simulation-based long-term forecast starting from the last week in 2001 (vertical bar on the left), showing the counts aggregated over all districts. The weekly mean of the simulations is represented by dots and the dashed lines correspond to the pointwise 2.5\\% and 97.5\\% quantiles. The actually observed counts are shown in the background.", fig.pos="htb", echo=-1----
par(las = 1, mar = c(4,4,0,0)+.5)
plot(measlesSim, "time", ylim = c(0, 100))

