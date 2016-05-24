## ------------------------------------------------------------------------
library(aop)
set.seed(1332142)
bmr_obj <- bootstrap_metaregression(oxybenzone, 15, 1000)

## ------------------------------------------------------------------------
library(ggplot2)
plot_metaregression_spaghetti_plot(bmr_obj, number_to_plot=40)

## ------------------------------------------------------------------------
plot_metaregression_confidence_envelope(bmr_obj, graph_pod = FALSE)

## ------------------------------------------------------------------------
slope_pod <- slope_pod_analysis(bmr_obj, 0.0001, 10, 0.1)
pod_and_threshold <- pod_envelope_analysis(bmr_obj, slope_pod, 10, min(oxybenzone$Concentration), max(oxybenzone$Concentration), 0.1)
plot_metaregression_confidence_envelope(bmr_obj, graph_pod = TRUE, pod = pod_and_threshold$pod, pod_threshold=pod_and_threshold$threshold)

## ------------------------------------------------------------------------
pod_and_threshold

## ------------------------------------------------------------------------
oxybenzone_levels <- c(1.2, 1.5, 2, 2.2, 2.3, 2.3, 2.3, 2.7, 2.8, 3.0, 3.0, 3.5, 3.5, 3.6, 3.7, 3.8, 4.0, 4.1, 4.1, 6.0, 7.1, 8.0, 8.7)

## ----fig.width=4---------------------------------------------------------
oxybenzone_df <- data.frame(levels = oxybenzone_levels, group="Oxybenzone")
r_levels <- data.frame(levels = rpois(2000, var(oxybenzone_levels)), group = "Poisson")
oxybenzone_stuff <- rbind(oxybenzone_df, r_levels)

ggplot(oxybenzone_stuff, aes(x=levels, fill=group)) + 
  geom_density(alpha=.3, binwidth=1)

## ------------------------------------------------------------------------
oxybenzone_g_per_L <- 228.24 * pod_and_threshold[1]
oxybenzone_total_blood <- oxybenzone_g_per_L * 5
administered_population_lower <- NULL
if(qpois(1-0.005, 3.9, lower.tail=FALSE) != 0){
  administered_population_lower <- oxybenzone_total_blood / qpois(1-0.005, 3.9, lower.tail=FALSE)
} else{
  administered_population_lower <- 0
}
administered_population_higher <- oxybenzone_total_blood / qpois(1-0.995, 3.9, lower.tail=FALSE)
administered_population_lower
administered_population_higher

## ------------------------------------------------------------------------
administered_population_1_in_10000 <- oxybenzone_total_blood / qpois(1-0.99995, 3.9, lower.tail=FALSE)

## ------------------------------------------------------------------------
administered_population_higher / 100

## ------------------------------------------------------------------------
administered_population_1_in_10000 / 100

