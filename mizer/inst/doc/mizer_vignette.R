### R code from vignette source 'mizer_vignette.Rnw'

###################################################
### code chunk number 1: mizer_vignette.Rnw:116-117
###################################################
rm(list=ls())


###################################################
### code chunk number 2: mizer_vignette.Rnw:150-151 (eval = FALSE)
###################################################
## install.packages("mizer")


###################################################
### code chunk number 3: mizer_vignette.Rnw:155-156
###################################################
library(mizer)


###################################################
### code chunk number 4: help_demo (eval = FALSE)
###################################################
## help(package="mizer")
## help(mizer)
## help(project)


###################################################
### code chunk number 5: help_plot (eval = FALSE)
###################################################
## help(plot, package="mizer")


###################################################
### code chunk number 6: help_FL (eval = FALSE)
###################################################
## method ? getFeedingLevel


###################################################
### code chunk number 7: help_MP_class (eval = FALSE)
###################################################
## class ? MizerParams


###################################################
### code chunk number 8: help_set_community_model (eval = FALSE)
###################################################
## ?set_community_model


###################################################
### code chunk number 9: demo_comm_model_params
###################################################
params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7)


###################################################
### code chunk number 10: class_comm_params
###################################################
class(params)


###################################################
### code chunk number 11: slots_MP
###################################################
slotNames(params)


###################################################
### code chunk number 12: demo_slot_accessor
###################################################
params@w


###################################################
### code chunk number 13: summary_comm_params
###################################################
summary(params)


###################################################
### code chunk number 14: demo_comm_run
###################################################
sim <- project(params, t_max=50, effort = 0)


###################################################
### code chunk number 15: plot_comm_sim (eval = FALSE)
###################################################
## plot(sim)


###################################################
### code chunk number 16: print_plot_comm_sim
###################################################
plot(sim)


###################################################
### code chunk number 17: get_m2_comm
###################################################
m2 <- getM2(sim)


###################################################
### code chunk number 18: m2_comm_final_time_step
###################################################
m2[51,]


###################################################
### code chunk number 19: print_plot_comm_m2
###################################################
plot(x = sim@params@w, y = m2[51,], log="xy", type="l", xlab = "Size", ylab = "Predation mortality")


###################################################
### code chunk number 20: set_comm_fishing
###################################################
params_knife <- set_community_model(z0 = 0.1, recruitment = 4e7,
    alpha = 0.2, f0 = 0.7, knife_edge_size = 1000)


###################################################
### code chunk number 21: sim_comm_no_fish
###################################################
sim0 <- project(params_knife, effort = 0, t_max = 50)


###################################################
### code chunk number 22: sim_comm_with_fish
###################################################
sim1 <- project(params_knife, effort = 1, t_max = 50)


###################################################
### code chunk number 23: plot_comm_fmort (eval = FALSE)
###################################################
## plot(sim1)


###################################################
### code chunk number 24: print_plot_comm_fmort
###################################################
plot(sim1)


###################################################
### code chunk number 25: 'dim_n'
###################################################
dim(sim0@n)


###################################################
### code chunk number 26: 'relative_comm_abundance'
###################################################
relative_abundance <- sim1@n[51,,] / sim0@n[51,,]


###################################################
### code chunk number 27: plot_relative_comm_abund (eval = FALSE)
###################################################
## plot(x=sim0@params@w, y=relative_abundance, log="x", type="n",
##     xlab = "Size (g)", ylab="Relative abundance")
## lines(x=sim0@params@w, y=relative_abundance)
## lines(x=c(min(sim0@params@w),max(sim0@params@w)), y=c(1,1),lty=2)


###################################################
### code chunk number 28: print_plot_relative_comm_abund
###################################################
plot(x=sim0@params@w, y=relative_abundance, log="x", type="n",
    xlab = "Size (g)", ylab="Relative abundance")
lines(x=sim0@params@w, y=relative_abundance)
lines(x=c(min(sim0@params@w),max(sim0@params@w)), y=c(1,1),lty=2)


###################################################
### code chunk number 29: comm_params_sigma1
###################################################
params_sigma1 <- set_community_model(z0 = 0.1,
     f0 = 0.7, alpha = 0.2, recruitment = 4e7, sigma = 1)


###################################################
### code chunk number 30: comm_project_sigma1
###################################################
sim_sigma1 <- project(params_sigma1, effort = 0, t_max = 50, dt=0.01)


###################################################
### code chunk number 31: plot_comm_biomass_sigma1 (eval = FALSE)
###################################################
## plotBiomass(sim_sigma1)


###################################################
### code chunk number 32: plot_comm_biomass_sigma1
###################################################
plotBiomass(sim_sigma1)


###################################################
### code chunk number 33: help_set_trait_model (eval = FALSE)
###################################################
## ?set_trait_model


###################################################
### code chunk number 34: demo_trait_model_params
###################################################
params <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5)


###################################################
### code chunk number 35: class_trait_params
###################################################
class(params)


###################################################
### code chunk number 36: summary_comm_params
###################################################
summary(params)


###################################################
### code chunk number 37: trait_project_no_fishing
###################################################
sim <- project(params, t_max=75, effort = 0)


###################################################
### code chunk number 38: plot_comm_sim (eval = FALSE)
###################################################
## plot(sim)


###################################################
### code chunk number 39: print_plot_trait_sim
###################################################
plot(sim)


###################################################
### code chunk number 40: set_trait_fishing
###################################################
params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5,
    knife_edge_size = 1000)


###################################################
### code chunk number 41: sim_trait_no_fish
###################################################
sim0 <- project(params_knife, effort = 0, t_max = 75)


###################################################
### code chunk number 42: sim_trait_with_fish
###################################################
sim1 <- project(params_knife, effort = 0.75, t_max = 75)


###################################################
### code chunk number 43: plot_trait_fmort (eval = FALSE)
###################################################
## plot(sim1)


###################################################
### code chunk number 44: print_plot_trait_fmort
###################################################
plot(sim1)


###################################################
### code chunk number 45: trait_dim_n
###################################################
dim(sim0@n)


###################################################
### code chunk number 46: mizer_vignette.Rnw:837-839
###################################################
total_abund0 <- apply(sim0@n[76,,],2,sum)
total_abund1 <- apply(sim1@n[76,,],2,sum)


###################################################
### code chunk number 47: 'relative_comm_abundance'
###################################################
relative_abundance <- total_abund1 / total_abund0


###################################################
### code chunk number 48: plot_relative_comm_abund (eval = FALSE)
###################################################
## plot(x=sim0@params@w, y=relative_abundance, log="xy", type="n", xlab = "Size (g)",
##     ylab="Relative abundance", ylim = c(0.1,10))
## lines(x=sim0@params@w, y=relative_abundance)
## lines(x=c(min(sim0@params@w),max(sim0@params@w)), y=c(1,1),lty=2)


###################################################
### code chunk number 49: show_plot_relative_comm_abund
###################################################
plot(x=sim0@params@w, y=relative_abundance, log="xy", type="n", xlab = "Size (g)",
    ylab="Relative abundance", ylim = c(0.1,10))
lines(x=sim0@params@w, y=relative_abundance)
lines(x=c(min(sim0@params@w),max(sim0@params@w)), y=c(1,1),lty=2)


###################################################
### code chunk number 50: 'get_M2_trait'
###################################################
m2_no_fishing <- getM2(sim0)[76,1,]
m2_with_fishing <- getM2(sim1)[76,1,]


###################################################
### code chunk number 51: plot_relative_trait_m2 (eval = FALSE)
###################################################
## plot(x = sim0@params@w, y = m2_no_fishing, log="x", type="n", xlab = "Size (g)",
##     ylab = "M2")
## lines(x = sim0@params@w, y = m2_no_fishing, lty=2)
## lines(x = sim0@params@w, y = m2_with_fishing)


###################################################
### code chunk number 52: print_plot_relative_trait_m2
###################################################
plot(x = sim0@params@w, y = m2_no_fishing, log="x", type="n", xlab = "Size (g)",
    ylab = "M2")
lines(x = sim0@params@w, y = m2_no_fishing, lty=2)
lines(x = sim0@params@w, y = m2_with_fishing)


###################################################
### code chunk number 53: calc_trait_winf
###################################################
no_sp <- 10
min_w_inf <- 10
max_w_inf <- 1e5
w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)


###################################################
### code chunk number 54: calc_knife_edge_trait
###################################################
knife_edges <- w_inf * 0.05


###################################################
### code chunk number 55: trait_industrial_gear_names
###################################################
other_gears <- w_inf > 500
gear_names <- rep("Industrial", no_sp)
gear_names[other_gears] <- "Other"


###################################################
### code chunk number 56: set_trait_multiple_gears
###################################################
params_multi_gear <- set_trait_model(no_sp = no_sp, min_w_inf = min_w_inf,
    max_w_inf = max_w_inf, knife_edge_size = knife_edges, gear_names = gear_names)


###################################################
### code chunk number 57: looking_at_trait_species_params
###################################################
params_multi_gear@species_params


###################################################
### code chunk number 58: project_multi_gear_trait_single_effort
###################################################
sim_multi_gear <- project(params_multi_gear, t_max = 75, effort = 0.5)


###################################################
### code chunk number 59: plot_multi_gear_trait_single_effort (eval = FALSE)
###################################################
## plot(sim_multi_gear)


###################################################
### code chunk number 60: print_plot_trait_multi_gear_single_effort
###################################################
plot(sim_multi_gear)


###################################################
### code chunk number 61: project_trait_multigear
###################################################
sim_multi_gear <- project(params_multi_gear, t_max = 75,
    effort = c(Industrial = 0.75, Other = 0))


###################################################
### code chunk number 62: plot_multi_gear_trait (eval = FALSE)
###################################################
## plot(sim_multi_gear)


###################################################
### code chunk number 63: print_plot_trait_multi_gear
###################################################
plot(sim_multi_gear)


###################################################
### code chunk number 64: industrial_fishery_simulation
###################################################
sim_industrial0 <- project(params_multi_gear, t_max = 75, effort = 0)
sim_industrial1 <- project(params_multi_gear, t_max = 75,
    effort = c(Industrial = 0.75, Other = 0))
total_abund0 <- apply(sim_industrial0@n[76,,],2,sum)
total_abund1 <- apply(sim_industrial1@n[76,,],2,sum)
relative_abundance <- total_abund1 / total_abund0


###################################################
### code chunk number 65: plot_relative_comm_abund_industrial (eval = FALSE)
###################################################
## plot(x=sim0@params@w, y=relative_abundance, log="xy", type="n", xlab = "Size (g)",
##     ylab="Relative abundance", ylim = c(0.1,10))
## lines(x=sim0@params@w, y=relative_abundance)
## lines(x=c(min(sim0@params@w),max(sim0@params@w)), y=c(1,1),lty=2)


###################################################
### code chunk number 66: show_plot_relative_comm_abund_industrial
###################################################
plot(x=sim0@params@w, y=relative_abundance, log="xy", type="n", xlab = "Size (g)",
    ylab="Relative abundance", ylim = c(0.1,10))
lines(x=sim0@params@w, y=relative_abundance)
lines(x=c(min(sim0@params@w),max(sim0@params@w)), y=c(1,1),lty=2)


###################################################
### code chunk number 67: mizer_vignette.Rnw:1222-1223 (eval = FALSE)
###################################################
## ?knife_edge


###################################################
### code chunk number 68: help_MP_constructor (eval = FALSE)
###################################################
## help(MizerParams)


###################################################
### code chunk number 69: get_location_for_ns_params (eval = FALSE)
###################################################
## system.file("doc/NS_species_params.csv",package="mizer")


###################################################
### code chunk number 70: mizer_vignette.Rnw:1354-1355
###################################################
params_data <- read.csv("NS_species_params.csv")


###################################################
### code chunk number 71: mizer_vignette.Rnw:1361-1362
###################################################
class(params_data)


###################################################
### code chunk number 72: show_simple_params_data
###################################################
params_data


###################################################
### code chunk number 73: first_MP
###################################################
params <- MizerParams(params_data)


###################################################
### code chunk number 74: class_MP
###################################################
class(params)


###################################################
### code chunk number 75: looking_at_params_slot
###################################################
params@species_params


###################################################
### code chunk number 76: summary_params
###################################################
summary(params)


###################################################
### code chunk number 77: MP_200 (eval = FALSE)
###################################################
## params200 <- MizerParams(params_data, no_w=200)
## summary(params200)


###################################################
### code chunk number 78: mizer_vignette.Rnw:1430-1431 (eval = FALSE)
###################################################
## system.file("doc/inter.csv",package="mizer")


###################################################
### code chunk number 79: mizer_vignette.Rnw:1436-1437
###################################################
inter <- read.csv("inter.csv", row.names=1)


###################################################
### code chunk number 80: inter_df_to_matrix
###################################################
inter <- as(inter, "matrix")


###################################################
### code chunk number 81: mizer_vignette.Rnw:1451-1452
###################################################
params <- MizerParams(params_data, interaction = inter)


###################################################
### code chunk number 82: mizer_vignette.Rnw:1472-1477
###################################################
params_data_gears <- params_data
params_data_gears$gear <- c("Industrial","Industrial","Industrial",
		      "Pelagic","Beam","Otter",
		      "Beam","Otter","Beam",
		      "Otter","Otter","Otter")


###################################################
### code chunk number 83: set_params_gears
###################################################
params_gears <- MizerParams(params_data_gears, interaction = inter)


###################################################
### code chunk number 84: remaking_simple_ms_params
###################################################
params <- MizerParams(params_data, interaction = inter)


###################################################
### code chunk number 85: siple_projection_single_effort
###################################################
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)


###################################################
### code chunk number 86: plot_basic_ms_sim (eval = FALSE)
###################################################
## plot(sim)


###################################################
### code chunk number 87: print_plot_basic_ms_sim
###################################################
plot(sim)


###################################################
### code chunk number 88: show_effort
###################################################
head(sim@effort)


###################################################
### code chunk number 89: summary_mizersim (eval = FALSE)
###################################################
## summary(sim)


###################################################
### code chunk number 90: simple_proj_tsave05
###################################################
sim <- project(params_gears, effort = 1, t_max = 10, dt = 0.1, t_save = 0.5)
head(sim@effort)


###################################################
### code chunk number 91: set_effort_vector
###################################################
effort <- c(Industrial = 0, Pelagic = 1, Beam = 0.3, Otter = 0.7)


###################################################
### code chunk number 92: project_effort_vector
###################################################
sim <- project(params_gears, effort = effort, t_max = 10, dt = 1, t_save = 1)
head(sim@effort)


###################################################
### code chunk number 93: plotFMort (eval = FALSE)
###################################################
## plotFMort(sim)


###################################################
### code chunk number 94: print_plotFMort
###################################################
plotFMort(sim)


###################################################
### code chunk number 95: set_empty_effort_array
###################################################
gear_names <- c("Industrial","Pelagic","Beam","Otter")
times <- seq(from = 1, to = 10, by = 1)
effort_array <- array(NA, dim = c(length(times), length(gear_names)),
    dimnames = list(time = times, gear = gear_names))


###################################################
### code chunk number 96: fill_effort_array
###################################################
effort_array[,"Industrial"] <- 0.5
effort_array[,"Pelagic"] <- seq(from = 1, to = 2, length = length(times))
effort_array[,"Beam"] <- seq(from = 1, to = 0, length = length(times))
effort_array[,"Otter"] <- seq(from = 1, to = 0.5, length = length(times))


###################################################
### code chunk number 97: show_effort_array
###################################################
head(effort_array)


###################################################
### code chunk number 98: project_effort_array
###################################################
sim <- project(params_gears,effort=effort_array, dt=0.1, t_save = 1)
head(sim@effort)


###################################################
### code chunk number 99: exploration_simulation
###################################################
sim <- project(params_gears,effort=effort_array, dt=0.1, t_save = 1)


###################################################
### code chunk number 100: show_dim_n
###################################################
dim(sim@n)


###################################################
### code chunk number 101: mizer_vignette.Rnw:1772-1773
###################################################
sim@n[,"Cod",]


###################################################
### code chunk number 102: getSSB_demo
###################################################
ssb <- getSSB(sim)
dim(ssb)
head(ssb)


###################################################
### code chunk number 103: getBiomass_demo
###################################################
biomass <- getBiomass(sim, min_w = 10, max_w = 1000)
head(biomass)


###################################################
### code chunk number 104: getCommunitySlope_demo
###################################################
slope <- getCommunitySlope(sim)
head(slope)


###################################################
### code chunk number 105: getCommunitySlope_with_args
###################################################
dem_species <- c("Dab","Whiting","Sole","Gurnard","Plaice","Haddock",
    "Cod","Saithe")
slope <- getCommunitySlope(sim, min_w = 10, max_w = 5000, 
    species = dem_species)
head(slope)


###################################################
### code chunk number 106: plotBiomass (eval = FALSE)
###################################################
## plotBiomass(sim)


###################################################
### code chunk number 107: print_plotBiomass
###################################################
plotBiomass(sim)


###################################################
### code chunk number 108: plotSpectra_example (eval = FALSE)
###################################################
## plotSpectra(sim, time_range = 5:10, biomass=TRUE)


###################################################
### code chunk number 109: print_plotSpectra_example
###################################################
plotSpectra(sim, time_range = 5:10, biomass=TRUE)


###################################################
### code chunk number 110: demo_summary_plot (eval = FALSE)
###################################################
## plot(sim)


###################################################
### code chunk number 111: plot_example
###################################################
plot(sim)


###################################################
### code chunk number 112: show_load_ns_species_params (eval = FALSE)
###################################################
## params_location <- system.file("doc/NS_species_params.csv",package="mizer")
## params_data <- read.csv(params_location)
## inter_location <- system.file("doc/inter.csv",package="mizer")
## inter <- as(read.csv(inter_location, row.names=1),"matrix")


###################################################
### code chunk number 113: load_ns_species_params
###################################################
params_data <- read.csv("NS_species_params.csv")
inter <- as(read.csv("inter.csv", row.names=1),"matrix")


###################################################
### code chunk number 114: adding_sel_params
###################################################
params_data$sel_func <- "sigmoid_length"
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
    19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
    24.3, 22.9, 43.6)
params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
    3.101, 3.160, 3.173, 3.075)


###################################################
### code chunk number 115: show_load_f_history (eval = FALSE)
###################################################
## f_location <- system.file("doc/NS_f_history.csv",package="mizer")
## f_history <- as(read.csv(f_location, row.names=1), "matrix")


###################################################
### code chunk number 116: load_f_history
###################################################
f_history <- as(read.csv("NS_f_history.csv", row.names=1), "matrix")


###################################################
### code chunk number 117: head_f_history
###################################################
head(f_history)


###################################################
### code chunk number 118: set_catchability
###################################################
params_data$catchability <- as.numeric(f_history["1990",])


###################################################
### code chunk number 119: make_ns_params
###################################################
params <- MizerParams(params_data, inter, kappa = 9.27e10)


###################################################
### code chunk number 120: rescale_effort
###################################################
relative_effort <- sweep(f_history,2,f_history["1990",],"/")
relative_effort[as.character(1988:1992),]


###################################################
### code chunk number 121: add_transient_years
###################################################
initial_effort <- matrix(relative_effort[1,],byrow=TRUE, nrow=100,
    ncol=ncol(relative_effort), dimnames = list(1867:1966))
relative_effort <- rbind(initial_effort,relative_effort)


###################################################
### code chunk number 122: project_ns_model
###################################################
sim <- project(params, effort=relative_effort, dt = 0.5, t_save = 1)


###################################################
### code chunk number 123: plot_ns_biomass (eval = FALSE)
###################################################
## plotBiomass(sim)


###################################################
### code chunk number 124: print_plot_ns_biomass
###################################################
plotBiomass(sim)


###################################################
### code chunk number 125: unexploited_na
###################################################
sim0 <- project(params, effort=0, dt = 0.5, t_save = 1, t_max = 100)


###################################################
### code chunk number 126: ns_comm_ref_point
###################################################
demersal_species <- c("Dab","Whiting","Sole","Gurnard","Plaice","Haddock",
    "Cod","Saithe")
lfi0 <- getProportionOfLargeFish(sim0, species = demersal_species,
    min_w = 10, max_w = 100e3, threshold_l = 40)["100"]
mw0 <- getMeanWeight(sim0, species = demersal_species,
    min_w = 10,max_w = 100e3)["100"]
mmw0 <- getMeanMaxWeight(sim0, species = demersal_species,
    min_w = 10, max_w = 100e3)["100","mmw_biomass"]
slope0 <- getCommunitySlope(sim0, species = demersal_species,
    min_w = 10, max_w = 100e3)["100","slope"]


###################################################
### code chunk number 127: ns_comm_indicators
###################################################
years <- 1967:2010
lfi <- getProportionOfLargeFish(sim, species = demersal_species,
    min_w = 10, max_w = 100e3, threshold_l = 40)[as.character(years)]
mw <- getMeanWeight(sim, species = demersal_species,
    min_w = 10, max_w = 100e3)[as.character(years)]
mmw <- getMeanMaxWeight(sim, species = demersal_species,
    min_w = 10, max_w = 100e3)[as.character(years),"mmw_biomass"]
slope <- getCommunitySlope(sim, species = demersal_species,
    min_w = 10, max_w = 100e3)[as.character(years),"slope"]


###################################################
### code chunk number 128: plot_ns_indicators (eval = FALSE)
###################################################
## library(ggplot2)
## # Simulated data
## community_plot_data <- rbind(
##     data.frame(year = years, measure = "LFI", data = lfi),
##     data.frame(year = years, measure = "Mean Weight", data = mw),
##     data.frame(year = years, measure = "Mean Max Weight", data = mmw),
##     data.frame(year = years, measure = "Slope", data = slope))
## # Unexploited data
## community_unfished_data <- rbind(
##     data.frame(year = years, measure = "LFI", data = lfi0[[1]]),
##     data.frame(year = years, measure = "Mean Weight", data = mw0[[1]]),
##     data.frame(year = years, measure = "Mean Max Weight", data = mmw0[[1]]),
##     data.frame(year = years, measure = "Slope", data = slope0[[1]]))
## # Reference level
## community_reference_level <-
##     data.frame(year=years, measure = "LFI", data = lfi0[[1]] * 0.8)
## # Build up the plot
## p <- ggplot(community_plot_data) + geom_line(aes(x=year, y = data)) +
##     facet_wrap(~measure, scales="free")
## p <- p + geom_line(aes(x=year,y=data), linetype="dashed",
##     data = community_unfished_data)
## p + geom_line(aes(x=year,y=data), linetype="dotted",
##     data = community_reference_level)


###################################################
### code chunk number 129: print_plot_ns_indicators
###################################################
library(ggplot2)
# Simulated data
community_plot_data <- rbind(
    data.frame(year = years, measure = "LFI", data = lfi),
    data.frame(year = years, measure = "Mean Weight", data = mw),
    data.frame(year = years, measure = "Mean Max Weight", data = mmw),
    data.frame(year = years, measure = "Slope", data = slope))
# Unexploited data
community_unfished_data <- rbind(
    data.frame(year = years, measure = "LFI", data = lfi0[[1]]),
    data.frame(year = years, measure = "Mean Weight", data = mw0[[1]]),
    data.frame(year = years, measure = "Mean Max Weight", data = mmw0[[1]]),
    data.frame(year = years, measure = "Slope", data = slope0[[1]]))
# Reference level
community_reference_level <-
    data.frame(year=years, measure = "LFI", data = lfi0[[1]] * 0.8)
# Build up the plot
p <- ggplot(community_plot_data) + geom_line(aes(x=year, y = data)) +
    facet_wrap(~measure, scales="free")
p <- p + geom_line(aes(x=year,y=data), linetype="dashed",
    data = community_unfished_data)
p + geom_line(aes(x=year,y=data), linetype="dotted",
    data = community_reference_level)


###################################################
### code chunk number 130: ns_scenario1_relative_effort
###################################################
scenario1 <- t(array(relative_effort["2010",], dim=c(12,40),
    dimnames=list(NULL,year = 2011:2050)))
scenario1 <- rbind(relative_effort, scenario1)


###################################################
### code chunk number 131: ns_scenario2_relative_effort
###################################################
fmsy <- c(Sprat = 0.2, Sandeel = 0.2, N.pout = 0.2, Herring = 0.25, Dab = 0.2,
    Whiting = 0.2, Sole = 0.22, Gurnard = 0.2, Plaice = 0.25, Haddock = 0.3,
    Cod = 0.19, Saithe = 0.3)
scenario2 <- t(array(fmsy, dim=c(12,40), dimnames=list(NULL,year = 2011:2050)))
scenario2 <- rbind(relative_effort, scenario2)
for (sp in dimnames(scenario2)[[2]]){
    scenario2[as.character(2011:2015),sp] <- scenario2["2010",sp] +
        (((scenario2["2015",sp] - scenario2["2010",sp]) / 5) * 1:5)
}


###################################################
### code chunk number 132: project_ns_future_scenarios
###################################################
sim1 <- project(params, effort = scenario1, dt = 0.5, t_save = 1)
sim2 <- project(params, effort = scenario2, dt = 0.5, t_save = 1)


###################################################
### code chunk number 133: ns_biodiv_ref_point
###################################################
ssb0 <- getSSB(sim0)["100",]


###################################################
### code chunk number 134: projected_ssb_ns (eval = FALSE)
###################################################
## library(reshape2)
## years <- 1967:2050
## ssb1 <- getSSB(sim1)[as.character(years),]
## ssb2 <- getSSB(sim2)[as.character(years),]
## ssb1_df <- melt(ssb1)
## ssb2_df <- melt(ssb2)
## ssb_df <- rbind(
##     cbind(ssb1_df, scenario = "Scenario 1"),
##     cbind(ssb2_df, scenario = "Scenario 2"))
## ssb_unexploited_df <- cbind(expand.grid(
##     sp = names(ssb0),
##     time = 1967:2050),
##     value = as.numeric(ssb0),
##     scenario = "Unexploited")
## ssb_reference_df <- cbind(expand.grid(
##     sp = names(ssb0),
##     time = 1967:2050),
##     value = as.numeric(ssb0*0.1),
##     scenario = "Reference")
## ssb_all_df <- rbind(
##     ssb_df,
##     ssb_unexploited_df,
##     ssb_reference_df)
## p <- ggplot(ssb_all_df) +
##     geom_line(aes(x = time, y = value, colour = scenario)) +
##     facet_wrap(~sp, scales = "free", nrow = 4)
## p + theme(legend.position = "none")


###################################################
### code chunk number 135: print_projected_ssb_ns
###################################################
library(reshape2)
years <- 1967:2050
ssb1 <- getSSB(sim1)[as.character(years),]
ssb2 <- getSSB(sim2)[as.character(years),]
ssb1_df <- melt(ssb1)
ssb2_df <- melt(ssb2)
ssb_df <- rbind(
    cbind(ssb1_df, scenario = "Scenario 1"),
    cbind(ssb2_df, scenario = "Scenario 2"))
ssb_unexploited_df <- cbind(expand.grid(
    sp = names(ssb0),
    time = 1967:2050),
    value = as.numeric(ssb0),
    scenario = "Unexploited")
ssb_reference_df <- cbind(expand.grid(
    sp = names(ssb0),
    time = 1967:2050),
    value = as.numeric(ssb0*0.1),
    scenario = "Reference")
ssb_all_df <- rbind(
    ssb_df,
    ssb_unexploited_df,
    ssb_reference_df)
p <- ggplot(ssb_all_df) +
    geom_line(aes(x = time, y = value, colour = scenario)) +
    facet_wrap(~sp, scales = "free", nrow = 4)
p + theme(legend.position = "none")


