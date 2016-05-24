## ----include = FALSE---------------------------------------------------------------
## load the "cool" package
library("surveillance")

## Compute everything or fetch cached results?
message("Doing computations: ",
        COMPUTE <- !file.exists("twinSIR-cache.RData"))
if (!COMPUTE) load("twinSIR-cache.RData", verbose = TRUE)

## ----hagelloch.df------------------------------------------------------------------
data("hagelloch")
head(hagelloch.df, n = 5)

## ----hagelloch---------------------------------------------------------------------
hagelloch <- as.epidata(hagelloch.df,
  t0 = 0, tI.col = "tI", tR.col = "tR",
  id.col = "PN", coords.cols = c("x.loc", "y.loc"),
  f = list(household    = function(u) u == 0,
           nothousehold = function(u) u > 0),
  w = list(c1 = function (CL.i, CL.j) CL.i == "1st class" & CL.j == CL.i,
           c2 = function (CL.i, CL.j) CL.i == "2nd class" & CL.j == CL.i),
  keep.cols = c("SEX", "AGE", "CL"))

## ----hagelloch_show, warning=FALSE-------------------------------------------------
head(hagelloch, n = 5)

## ----hagelloch_plot, echo=2, fig.cap="Evolution of the 1861 Hagelloch measles epidemic in terms of the numbers of susceptible, infectious, and recovered children. The bottom \\code{rug} marks the infection times \\code{tI}.", fig.pos="!h"----
par(mar = c(5, 5, 1, 1))
plot(hagelloch, xlab = "Time [days]")

## ----hagelloch_households, fig.cap="Spatial locations of the Hagelloch households. The size of each dot is proportional to the number of children in the household.", fig.pos="ht", echo=-1----
par(mar = c(5, 5, 1, 1))
hagelloch_coords <- summary(hagelloch)$coordinates
plot(hagelloch_coords, xlab = "x [m]", ylab = "y [m]",
  pch = 15, asp = 1, cex = sqrt(multiplicity(hagelloch_coords)))
legend(x = "topleft", pch = 15, legend = c(1, 4, 8), pt.cex = sqrt(c(1, 4, 8)),
  title = "Household size")

## ----hagellochFit, results='hide'--------------------------------------------------
hagellochFit <- twinSIR(~household + c1 + c2 + nothousehold, data = hagelloch)

## ----hagellochFit_summary_echo, eval=FALSE-----------------------------------------
#  set.seed(1)
#  summary(hagellochFit)

## ----hagellochFit_confint----------------------------------------------------------
exp(confint(hagellochFit, parm = "cox(logbaseline)"))

## ----hagellochFit_profile, results='hide', eval=COMPUTE----------------------------
#  prof <- profile(hagellochFit,
#    list(c(match("c1", names(coef(hagellochFit))), NA, NA, 25),
#         c(match("c2", names(coef(hagellochFit))), NA, NA, 25)))

## ----------------------------------------------------------------------------------
prof$ci.hl

## ----hagellochFit_profile_plot, fig.cap="Normalized log-likelihood for $\\alpha_{c1}$ and $\\alpha_{c2}$ when fitting the \\code{twinSIR} model formulated in Equation~\\eqref{eqn:twinSIR:hagelloch} to the Hagelloch data.", fig.pos="ht", fig.height=4.4----
plot(prof)

## ----hagellochFit_plot, echo=2, fig.width=4.5, fig.height=4.5, out.width="0.49\\linewidth", fig.subcap=c("Epidemic proportion.","Transformed residuals."), fig.cap="Diagnostic plots for the \\code{twinSIR} model formulated in Equation~\\ref{eqn:twinSIR:hagelloch}.", fig.pos="htb"----
par(mar = c(5, 5, 1, 1))
plot(hagellochFit, which = "epidemic proportion", xlab = "time [days]")
checkResidualProcess(hagellochFit, plot = 1)

## ----hagellochFit_fstep, results='hide'--------------------------------------------
knots <- c(100, 200)
fstep <- list(
  B1 = function(D) D > 0 & D < knots[1],
  B2 = function(D) D >= knots[1] & D < knots[2],
  B3 = function(D) D >= knots[2])
hagellochFit_fstep <- twinSIR(
  ~household + c1 + c2 + B1 + B2 + B3,
  data = update(hagelloch, f = fstep))

## ----hagellochFit_AIC--------------------------------------------------------------
set.seed(1)
AIC(hagellochFit, hagellochFit_fstep)

