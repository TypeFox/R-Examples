library(PhaseType)

set.seed(34752076)

# I) Basic test of function
x <- c(1.45353415045187, 1.85349532001349, 2.01084961814576, 0.505725921290172, 
1.56252630012213, 3.41158665930278, 1.52674487509487, 4.3428662377235, 
8.03208018151311, 2.41746547476986, 0.38828086509283, 2.61513815012196, 
3.39148865480856, 1.82705817807965, 1.42090953713845, 0.851438991331866, 
0.0178808867191894, 0.632198596390046, 0.959910259815998, 1.83344199966323)

# Prior on starting state
dirpi <- c(1, 0, 0)
# Define the structure of the Phase-type generator
T <- matrix(c(0,"R","R",0,"F",0,0,0,"F",0,0,0,0,"F","F",0), 4)
# Gamma prior: shape hyperparameters (one per model parameter)
nu <- list("R"=180, "F"=24)
# Gamma prior: reciprocal scale hyperparameters (one per model parameter)
zeta <- c("R"=16,"F"=16)
# Perform 6 MCMC iterations
print(phtMCMC2(x, T, dirpi, nu, zeta, 20), digits=2)
