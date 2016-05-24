## ----eval=TRUE,echo=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(fig.width=6, fig.height=6)

library(ggplot2)
library(xtable)
library(dpcR)

size_mod <- -2
cool_theme <- theme(plot.background=element_rect(fill = "transparent",
                                                 colour = "transparent"),
                    panel.grid.major = element_line(colour="lightgrey", linetype = "dashed"),
                    panel.background = element_rect(fill = "white", colour = "black"),
                    legend.background = element_rect(fill="NA"),
                    legend.position = "bottom",
                    axis.text = element_text(size=12 + size_mod),
                    axis.title.x = element_text(size=15 + size_mod, vjust = -0.1), 
                    axis.title.y = element_text(size=15 + size_mod, vjust = 1),
                    strip.text = element_text(size=15 + size_mod, face = "bold"),
                    strip.background = element_rect(fill = "#9ecae1", colour = "black"),
                    legend.text = element_text(size=12 + size_mod), 
                    legend.title = element_text(size=15 + size_mod),
                    plot.title = element_text(size=20 + size_mod))

load("vig.RData")

## ----eval=TRUE-----------------------------------------------------------
# Below we have S4 object
s4 <- sim_adpcr(m = 100, n = 496, times = 100, pos_sums = FALSE, n_panels = 3)
# Is it a dpcr object?
class(s4)
# Yes, it is. Let's see what we have in type slot
slot(s4, "type")
# We can use also shorter notation
s4@type

## ----eval=TRUE-----------------------------------------------------------
# Create single adpcr object. The following code is also true for 
# other objects inhering from dpcr, as ddpcr or qdpcr
single_run <- sim_adpcr(m = 100, n = 765, times = 100, pos_sums = FALSE, n_panels = 1)
two_runs <- bind_dpcr(single_run, single_run)
three_runs <- bind_dpcr(single_run, single_run, single_run)
# It is also possible to bind a list of dpcr objects... 
three_runs_list <- bind_dpcr(list(single_run, single_run, single_run))
# ... which may be useful in do.call statements
dpcr_list <- do.call(bind_dpcr, lapply(5L:10*10, function(n_template)
  sim_adpcr(m = n_template, n = 765, times = 100, pos_sums = FALSE, n_panels = 1)))

## ----eval=TRUE-----------------------------------------------------------
longer_run <- sim_adpcr(m = 10, n = 15, times = 100, pos_sums = FALSE, n_panels = 1)
shorter_run <- sim_adpcr(m = 10, n = 10, times = 100, pos_sums = FALSE, n_panels = 1)
shortest_run <- sim_adpcr(m = 10, n = 5, times = 100, pos_sums = FALSE, n_panels = 1)
# Expect informative message after binding
res <- bind_dpcr(longer_run, shorter_run, shortest_run)
# Print the whole data
slot(res, ".Data")

## ----eval=TRUE-----------------------------------------------------------
five_runs <- sim_adpcr(m = 2, n = 10, times = 100, pos_sums = FALSE, n_panels = 5)
print(five_runs)

# Extract runs by the index
only_first_run <- extract_dpcr(five_runs, 1)
only_first_and_second_run <- extract_dpcr(five_runs, c(1, 2))
# See if proper replicated were extracted
slot(only_first_and_second_run, "replicate")
no_first_run <- extract_dpcr(five_runs, -1)
slot(no_first_run, "replicate")
# Extract runs by the name
run_Experiment1.3 <- extract_dpcr(five_runs, "Experiment1.3")
slot(run_Experiment1.3, "replicate")
run_Experiment1.3and5 <- extract_dpcr(five_runs, c("Experiment1.3", "Experiment1.5"))
slot(run_Experiment1.3and5, "replicate")

## ----eval=TRUE-----------------------------------------------------------
# Generate some data from 15x16 array. Let's presume, that we have results from two plates
sample_runs <- matrix(rpois(480, lambda = 1.5), ncol = 2)
# Check its class - it's a typical R structure
class(sample_runs)
# Save it to adpcr object
adpcr_experiments <- create_dpcr(sample_runs, n = c(240L, 240L), type = "nm", adpcr = TRUE)
class(adpcr_experiments)

## ----eval=TRUE-----------------------------------------------------------
# Create two array dPCR experiments. Mind the difference in the n parameter.
sample_adpcr <- bind_dpcr(sim_adpcr(m = 100, n = 765, times = 100, pos_sums = FALSE, n_panels = 1), 
                          rename_dpcr(sim_adpcr(m = 100, n = 763, times = 100, pos_sums = FALSE, 
                                           n_panels = 1), 
                                 exper = "Experiment2"))

## ----eval=TRUE-----------------------------------------------------------
# It's possible to manipulate data points from dpcr object using all functions that work for matrices
tail(sample_adpcr)

## ----eval=TRUE-----------------------------------------------------------
slot(sample_adpcr, "n")

## ----eval=TRUE-----------------------------------------------------------
# Quickly count positive partitions
colSums(sample_adpcr > 0)
# Baseline fluorescence data
sim_ddpcr(m = 3, n = 2, times = 5, fluo = list(0.1, 0)) - 0.05

## ----eval=TRUE-----------------------------------------------------------
# Inspect all types of data

# Cq
# Load package with qPCR data
library(chipPCR)
qpcr2pp(data = C127EGHP[, 1L:6], type = "ct")

# fluo
sim_ddpcr(m = 3, n = 2, times = 5, fluo = list(0.1, 0)) - 0.05

# nm
sim_adpcr(m = 235, n = 765, times = 100, pos_sums = FALSE, n_panels = 3)

# np
binarize(sim_adpcr(m = 235, n = 765, times = 100, pos_sums = FALSE, n_panels = 3))

# tnp
sim_adpcr(m = 235, n = 765, times = 100, pos_sums = TRUE, n_panels = 3)

## ----eval=TRUE-----------------------------------------------------------
summary(six_panels)

# Save summary data without printing it
summ <- summary(six_panels, print = FALSE)
# Print only the summary table
summ[["summary"]]

# Extract results for Dube's method
summ[["summary"]][summ[["summary"]][["method"]] == "dube", ]

## ----eval=TRUE-----------------------------------------------------------

sample_ddpcr <- sim_ddpcr(m = 3, n = 10, times = 5)
# Standard show method...
show(sample_ddpcr)
# ... which is an equivalent of:
sample_ddpcr
# If you want to see all data points:
slot(sample_ddpcr, ".Data")

## ----eval=TRUE-----------------------------------------------------------
adpcr2panel(six_panels)

## ----eval=TRUE-----------------------------------------------------------
# Remember, you can plot only single panel at once 
plot_panel(extract_dpcr(adpcr_experiments, 1), main = "Experiment 1")

## ----eval=TRUE-----------------------------------------------------------
plot_panel(binarize(extract_dpcr(adpcr_experiments, 1)), main = "Experiment 1")

## ----eval=TRUE-----------------------------------------------------------
# Extract graphical coordinates
panel_data <- plot_panel(extract_dpcr(adpcr_experiments, 1), plot = FALSE)
ggplot_coords <- cbind(panel_data[["ggplot_coords"]], value = as.vector(extract_dpcr(adpcr_experiments, 1)))

# Plot panel using different graphics package
library(ggplot2)
ggplot(ggplot_coords, aes(x = x, y = y, fill = value)) +
  geom_tile()


## ----eval=TRUE-----------------------------------------------------------
# The test_panel function performs a test for each experiment in apdr object.
test_panel(six_panels)

## ----eval=TRUE-----------------------------------------------------------
# Load chiPCR package to access C317.amp data
library(chipPCR)

# Convert data to qdpcr object
qdat <- qpcr2pp(data = C317.amp, type = "np", Cq_range = c(10, 30))

## ----eval=TRUE-----------------------------------------------------------
plot(qdat)

## ----eval=TRUE-----------------------------------------------------------
# Compare experiments using GLM

# 1. Perform test
comp <- test_counts(six_panels)

# 2. See summary of the test
summary(comp)

# 3. Plot results of the test 
plot(comp, aggregate = FALSE)

# 4. Aggregate runs to their groups
plot(comp, aggregate = TRUE)

# 5. Extract coefficients for the further usage
coef(comp)

## ----eval=TRUE-----------------------------------------------------------
#1. Perform multiple test comparison using data from the previous example
comp_ratio <- test_counts(six_panels, model = "ratio")

#2. See summary of the test
summary(comp_ratio)

#3. Plot results of the test 
plot(comp_ratio, aggregate = FALSE)

#4. Aggregate runs to their groups
plot(comp_ratio, aggregate = TRUE)

#5. Extract coefficients for the further usage
coef(comp)

# Compare results of two methods
par(mfrow=c(2,1))
plot(comp, aggregate = FALSE)
title("GLM")
plot(comp_ratio, aggregate = FALSE)
title("Ratio")
par(mfrow=c(1,1))



## ----eval=TRUE,echo=FALSE------------------------------------------------
ggplot(data=madpcr_comp, aes(x = value, fill = method)) +
  geom_density(alpha = 0.3) + 
  scale_fill_discrete("Confidence intervals:") + 
  scale_y_continuous("Density") + 
  scale_x_continuous("Fraction of wrongly assigned experiments") + 
  cool_theme

## ----eval=TRUE,echo=FALSE------------------------------------------------
ggplot(m_coverage2, aes(x = prop, y = value, fill = method)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_y_continuous("Probability coverage") + 
  scale_x_discrete(expression(lambda)) +
  scale_fill_discrete("Confidence intervals:") + 
  geom_hline(yintercept = 0.95, colour = "black", size = 1, linetype = 5) +
  facet_wrap(~ coverage, nrow = 2) + 
  cool_theme

## ----eval=TRUE,echo=FALSE,results="asis"---------------------------------
dat <- as.data.frame(aggregate(value ~ method + coverage, m_coverage2, mean))

colnames(dat) <- c("Method name", "Type of coverage", "Value")
print(xtable(dat), 
      include.rownames = FALSE, type = "html")

