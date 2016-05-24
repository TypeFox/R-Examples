### R code from vignette source 'pauwels2014.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: preliminaries
###################################################
#options(prompt = " ")
#options(continue = " ")
options(width=70)


###################################################
### code chunk number 2: figPlot
###################################################
library(pauwels2014)
data(knobjs)
sapply(	
			 1:length(knobjs),
			 function(k){
				 assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
			 }
)
data(exps)

mean_risks_act_mult <- read_knobjs( sprintf("knobjActMult%s", 1:10) ) 
mean_risks_dream6_mult <- read_knobjs( sprintf("knobjDream6Mult%s", 1:10) )  
mean_risks_rand_mult <- read_knobjs( sprintf("knobjRandMult%s", 1:10) ) 

mean_risks_act_mult <- compute_mean_risks(mean_risks_act_mult, "Bayesian active")
mean_risks_dream6_mult <- compute_mean_risks(mean_risks_dream6_mult, "Dream6")
mean_risks_rand_mult <- compute_mean_risks(mean_risks_rand_mult, "Random")

data_to_plot <- rbind(mean_risks_act_mult,mean_risks_dream6_mult, mean_risks_rand_mult)

ggplot(data = data_to_plot, aes(x=cost, y=risk)) + 
  theme_bw() + 
  facet_grid(.~type) + 
  scale_y_log10() + 
  geom_line(aes(group = chain), alpha = 0.2) + 
  geom_point() + 
  stat_smooth(method = "loess", colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###################################################
### code chunk number 3: fig1
###################################################
library(pauwels2014)
data(knobjs)
sapply(	
			 1:length(knobjs),
			 function(k){
				 assign(names(knobjs)[k], knobjs[[k]], envir = .GlobalEnv)
			 }
)
data(exps)

mean_risks_act_mult <- read_knobjs( sprintf("knobjActMult%s", 1:10) ) 
mean_risks_dream6_mult <- read_knobjs( sprintf("knobjDream6Mult%s", 1:10) )  
mean_risks_rand_mult <- read_knobjs( sprintf("knobjRandMult%s", 1:10) ) 

mean_risks_act_mult <- compute_mean_risks(mean_risks_act_mult, "Bayesian active")
mean_risks_dream6_mult <- compute_mean_risks(mean_risks_dream6_mult, "Dream6")
mean_risks_rand_mult <- compute_mean_risks(mean_risks_rand_mult, "Random")

data_to_plot <- rbind(mean_risks_act_mult,mean_risks_dream6_mult, mean_risks_rand_mult)

ggplot(data = data_to_plot, aes(x=cost, y=risk)) + 
  theme_bw() + 
  facet_grid(.~type) + 
  scale_y_log10() + 
  geom_line(aes(group = chain), alpha = 0.2) + 
  geom_point() + 
  stat_smooth(method = "loess", colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###################################################
### code chunk number 4: preliminaries2
###################################################
rm(list = ls())


###################################################
### code chunk number 5: re-run (eval = FALSE)
###################################################
## ## Load datasets
## data(exps)
## data(experiment_list1)
## data(observables)
## 
## ## Initialize a knowledge list
## knobj <- generate_our_knowledge(transform_params)
## knobj$datas[[1]] <- list(
##   manip = experiment_list1$nothing,
##   data = add_noise(
##    simulate_experiment(
##      knobj$global_parameters$true_params_T, knobj, 
##      experiment_list1$nothing)[
##        knobj$global_parameters$tspan %in% observables[["mrnaLow"]]$reso, 
##        observables[["mrnaLow"]]$obs
##      ]
##   )
## )
## knobj$experiments <- paste("nothing", "mrnaLow")
## 
## ## Update the knowledge list
## knobjActMult1 <- active_design(knobj, sample_function_multi_mod_weight, seed = 1, credits = 5000)


###################################################
### code chunk number 6: perform_simu_one (eval = FALSE)
###################################################
## ## Simulates the active learning experimental design with
## ## A simplified network
## ## Arguments to be passed
## ## active: O or 1
## ## multimodal search: 0 or 1
## ## random seed: integer
## ## Results are saved during the process
## 
## args <- commandArgs(trailingOnly=TRUE)
## if(length(args) != 3){
##   stop("Three argument should be passed: bool active (0 or 1), bool multimod (0 or 1), int seed")
## }
## 
## 
## active <- as.integer(as.numeric(args[1]))
## multi_mod <- as.logical(as.numeric(args[2]))
## seed <- as.integer(args[3])
## 
## 
## ## Define the home directory
## home <- "some_working_directory"
## 
## 
## ## Where to save the results
## file_to_save <- paste(
##   home,
##   "results_subdirectory/",
##   "knobj",
##   c("Rand","Act", "Dream6")[active +1], c("Sing","Mult")[multi_mod +1],
##   seed,
##   sep=""
## )
## 
## 
## ## Load datasets and library
## set.seed(seed)
## library(pauwels2014)
## data(exps)
## data(experiment_list1)
## data(observables)
## 
## 
## ## Initialize a knowledge list
## knobj <- generate_our_knowledge(transform_params)
## knobj$datas[[1]] <- list(
##   manip = experiment_list1$nothing,
##   data = add_noise(
##   simulate_experiment(
##     knobj$global_parameters$true_params_T, 
##     knobj, 
##     experiment_list1$nothing)[
##       knobj$global_parameters$tspan %in% observables[["mrnaLow"]]$reso, 
##       observables[["mrnaLow"]]$obs
##     ]
##   )
## )
## knobj$experiments <- paste("nothing", "mrnaLow")
## 
## 
## ## Chosse a sample function
## if(multi_mod ){
##   sample_function <- sample_function_multi_mod_weight
## }else{
##   sample_function <- sample_function_single_mod
## }
## 
## 
## ## Update the knowledge list with the chosen strategy
## if(active==1){
## 	knobj <- active_design(knobj, sample_function, seed, credits = 5000, file_to_save = file_to_save)
## }
## if(active==0){
## 	knobj <- random_design(knobj, sample_function, exps, seed, credits = 5000, file_to_save = file_to_save)
## }
## if(active==2){
## 	knobj <- dream6_design(knobj, sample_function, exps, seed, credits = 5000, file_to_save = file_to_save)
## }
## 


###################################################
### code chunk number 7: commands
###################################################
home <- "some_home/"
path_to_R <- "some_path/R"

n_seed <- 10
seeds <- 1:n_seed
active <- c(rep(2,n_seed),rep(1,2*n_seed),rep(0,n_seed))
multi_mod <- c(rep(1,2*n_seed),rep(0,n_seed),rep(1,n_seed))

commands <-
  sprintf(
    "%s --vanilla --args '%s' '%s' '%s' < %s",
    path_to_R,
    active,
    multi_mod,
    seeds,
    paste(home, "perform_simu_one.R", 
    sep = ""
  )
)

print(commands)


###################################################
### code chunk number 8: log_prior
###################################################
log_prior


###################################################
### code chunk number 9: log_prior
###################################################
log_likelihood


###################################################
### code chunk number 10: risk_theta_fun
###################################################
risk_theta_fun


###################################################
### code chunk number 11: tr_params
###################################################
knobj <- generate_our_knowledge(transform_params)
knobj$global_parameters$true_params_T
transform_params(knobj$global_parameters$true_params_T)


###################################################
### code chunk number 12: obserevables (eval = FALSE)
###################################################
## superHighRes <-  0:200/2
## highRes <- 0:50 * 2
## lowRes <- 0:25 * 4
## 
## observables <- list(
##   list(name = "p6", obs = c("time", "p6"),reso = superHighRes, cost = 400),
##   list(name = "p7", obs = c("time", "p7"),reso = superHighRes, cost = 400),
##   list(name = "p8", obs = c("time", "p8"),reso = superHighRes, cost = 400),
##   list(name = "mrnaHigh", obs = c("time", "v6_mrna", "v7_mrna", "v8_mrna"), reso = highRes, cost = 1000),
##   list(name = "mrnaLow", obs = c("time", "v6_mrna", "v7_mrna", "v8_mrna"), reso = lowRes, cost = 500)
## )
## 
## names(observables) <- sapply(observables, FUN = function(x){x$name})


###################################################
### code chunk number 13: experiments (eval = FALSE)
###################################################
## delete_gene6 <- function(theta, init){
##   theta[names(theta) %in% c("pro6_strength","rbs6_strength")] <- 0
##   init[names(init) == "p6"] <- 0
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 800
##   res
## }
## 
## delete_gene7 <- function(theta, init){
##   theta[names(theta) %in% c("pro7_strength","rbs8_strength")] <- 0
##   init[names(init) == "p7"] <- 0
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 800
##   res
## }
## 
## delete_gene8 <- function(theta, init){
##   theta[names(theta) %in% c("pro9_strength","rbs7_strength")] <- 0
##   init[names(init) == "p8"] <- 0
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 800
##   res
## }
## 
## 
## 
## knockdown_gene6 <- function(theta, init){
##   theta[names(theta) %in% "mrna6_degradation_rate"] <- 10
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 350
##   res
## }
## 
## knockdown_gene7 <- function(theta, init){
##   theta[names(theta) %in% "mrna7_degradation_rate"] <- 10
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 350
##   res
## }
## 
## knockdown_gene8 <- function(theta, init){
##   theta[names(theta) %in% "mrna8_degradation_rate"] <- 10
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 350
##   res
## }
## 
## 
## 
## decrease_rbs_gene6 <- function(theta, init){
##   theta[names(theta) %in% "rbs6_strength"] <- theta[names(theta) %in% "rbs6_strength"] / 10
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 450
##   res
## }
## 
## decrease_rbs_gene7 <- function(theta, init){
##   theta[names(theta) %in% "rbs8_strength"] <- theta[names(theta) %in% "rbs8_strength"] / 10
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 450
##   res
## }
## 
## decrease_rbs_gene8 <- function(theta, init){
##   theta[names(theta) %in% "rbs7_strength"] <- theta[names(theta) %in% "rbs7_strength"] / 10
##   res <- c()
##   res$theta <- theta
##   res$initial_conditions <- init
##   res$cost <- 450
##   res
## }
## 
## nothing <- function(theta, init){
##   list(theta = theta, initial_conditions = init, cost = 0)
## }
## 
## experiment_list1 <- c()
## experiment_list1$delete_gene6 <- delete_gene6
## experiment_list1$delete_gene7 <- delete_gene7
## experiment_list1$delete_gene8 <- delete_gene8
## experiment_list1$knockdown_gene6 <- knockdown_gene6
## experiment_list1$knockdown_gene7 <- knockdown_gene7
## experiment_list1$knockdown_gene8 <- knockdown_gene8
## experiment_list1$decrease_rbs_gene6 <- decrease_rbs_gene6
## experiment_list1$decrease_rbs_gene7 <- decrease_rbs_gene7
## experiment_list1$decrease_rbs_gene8 <- decrease_rbs_gene8
## experiment_list1$nothing <- nothing
## 
## exps <- lapply(names(experiment_list1), function(name){
##     t(sapply(observables, function(obs){
##         cbind(obs$name, obs$cost + experiment_list1[[name]](0,0)$cost, name)
##       }
##     ))
##   }
## )
## 
## exps <- Reduce(exps, f = function(s,t){rbind(s,t)}) 
## 
## exps <- data.frame(Measurement = exps[,1], Cost = as.numeric(exps[,2]), exp = exps[,3])


###################################################
### code chunk number 14: add_noise
###################################################
add_noise


