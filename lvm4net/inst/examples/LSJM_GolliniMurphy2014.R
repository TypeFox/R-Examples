####
## This code accompany the paper
## Gollini I., and Murphy T.B. (2014)
## " Joint Modelling of Multiple Network Views "
## Journal of Computational and Graphical Statistics
## arXiv:1301.3759
####
## author: Gollini I., <isabella.gollini@bristol.ac.uk>
## date: Sept 26th 2014

library(lvm4net)

D <- 2

########################
########################
#######  Girls Datasets #######
########################
########################

library(RSiena)  # load the datasets

Y1 <- s501
Y2 <- s502
Y3 <- s503

Nb <- nrow(Y1)

Y123 <- list(Y1 = Y1, Y2 = Y2, Y3 = Y3)

########################
########## LSM ##########
########################

modLSM1 <- lsm(Y1, D)
modLSM2 <- lsm(Y2, D)
modLSM3 <- lsm(Y3, D)

########################
########## LSJM ##########
########################

modLSJM123 <- lsjm(Y123, D) # It takes ~ 40 seconds

########################
#######match rotations #######
########################

Z123 <- modLSJM123$EZ

Z1 <- rotXtoY(modLSM1$lsmEZ, Z123)$X
Z2 <- rotXtoY(modLSM2$lsmEZ, Z123)$X
Z3 <- rotXtoY(modLSM3$lsmEZ, Z123)$X

XYlimb <- range(Z123, Z1, Z2, Z3)
namesb <- paste('Network ', 1:3, sep ='')

colPl <- rainbow(Nb,alpha=.8)

########################
######## LSM PLOTS ########
########################

par(mfrow = c(1,3))

plotY(Y1, Ndata = 1, EZ = Z1, VZ = modLSM1$lsmVZ, main = namesb[1], xlim = XYlimb, ylim = XYlimb, colPl = colPl)
plotY(Y2, Ndata = 1, EZ = Z2, VZ = modLSM2$lsmVZ, main = namesb[2], xlim = XYlimb, ylim = XYlimb, colPl = colPl)
plotY(Y3, Ndata = 1, EZ = Z3, VZ = modLSM3$lsmVZ, main = namesb[3], xlim = XYlimb, ylim = XYlimb, colPl = colPl)

Zlsm <- list()
Zlsm[[1]] <- Z1
Zlsm[[2]] <- Z2
Zlsm[[3]] <- Z3

bpbLSM <- boxroc(Y123, 
	EZ = Zlsm,
	xi = c(modLSM1$xiT, modLSM2$xiT, modLSM3$xiT), 
	Lroc = 150, 
	ROC = TRUE, 
	BOXPLOT = TRUE,
	labelsPlot = namesb
	)


########################
######## LSJM PLOTS ########
########################

plot(modLSJM123, Y123, drawCB = TRUE, plotZtilde = TRUE, colPl = colPl)
plot(modLSJM123, Y123, drawCB = TRUE, colPl = colPl, main = 'Multiple networks')
	
bpbLSJM <- boxroc(Y123, 
	EZ = modLSJM123$lsmEZ,	
	xi = modLSJM123$xiT, 
	Lroc = 150, 
	ROC = TRUE, 
	BOXPLOT = TRUE,
	labelsPlot = namesb
	)


########################
########################
#######  PPI Datasets ########
########################
########################

data(PPInet) # provided in the package

Yg <- PPIgen
Yp <- PPIphy

Na <- nrow(Yg)

Ygp <- list()

Ygp[[1]] <- Yg
Ygp[[2]] <- Yp

########################
########## LSM ##########
########################

modLSMg <- lsm(Yg, D) 
modLSMp <- lsm(Yp, D) 

########################
########## LSJM ##########
########################

modLSJMgp <- lsjm(Ygp, D) # It takes ~ 2 minutes

########################
#######match rotations #######
########################
Zgp <- modLSJMgp$EZ

Zg <- rotXtoY(modLSMg$lsmEZ, Zgp)$X
Zp <- rotXtoY(modLSMp$lsmEZ, Zgp)$X

XYlima <- range(Zgp, Zg, Zp)
namesa <- c('Genetic Interactions', 'Physical Interactions')

colPl <- rainbow(Na,alpha=.8)

########################
######## LSM PLOTS ########
########################
par(mfrow = c(1,2))

plotY(Yg, Ndata = 1, EZ = Zg, VZ = modLSMg$lsmVZ, main = namesa[1], xlim = XYlima, ylim = XYlima, colPl = colPl)
plotY(Yp, Ndata = 1, EZ = Zp, VZ = modLSMp$lsmVZ, main = namesa[2], xlim = XYlima, ylim = XYlima, colPl = colPl)

Zlsmgp <- list()
Zlsmgp[[1]] <- Zg
Zlsmgp[[2]] <- Zp

bpaLSM <- boxroc(Ygp, 
	EZ = Zlsmgp,
	xiT = c(modLSMg$xiT, modLSMp$xiT), 
	Lroc = 150, 
	ROC = TRUE, 
	BOXPLOT = TRUE,
	labelsPlot = namesa
	)

########################
######## LSJM  PLOTS #######
########################

plot(modLSJMgp, Ygp, drawCB = TRUE, plotZtilde = TRUE, mainZtilde = namesa, colPl = colPl)
plot(modLSJMgp, Ygp, drawCB = TRUE, colPl = colPl)
	
bpaLSJM <- boxroc(Ygp, 
	EZ = modLSJMgp$lsmEZ,	
	xiT = modLSJMgp$xiT, 
	Lroc = 150, 
	ROC = TRUE, 
	BOXPLOT = TRUE,
	labelsPlot = namesa
	)

########################
########################
#### Supplementary  Material ####
########################
########################


##########################
#  Remove unconnected nodes from Yg  #
##########################

cond <- seq(1:nrow(Yg))[rowSums(Yg) + colSums(Yg) == 0]
Yg2 <- Yg[-cond, -cond]

direct <- !all(Yg2 == t(Yg2))
nYg2 <-network(Yg2, directed = direct)

timings <- numeric(3)
names(timings) <-  c('latentnet', 'VBLPCM', 'lvm4net')

########################
### Euclidean distance -- mcmc ###
########################

library(latentnet)

timings[1] <- system.time(euclmcmc <- ergmm(nYg2 ~ euclidean(d = D)))[[3]]

#######################
# Euclidean distance -- variational  #
#######################

library(VBLPCM)

timings[2] <- system.time(euclvar<- vblpcmfit(vblpcmstart(nYg2, G = 1, d = D), maxiter = 300))[[3]]

##########################
# (Euclidean distance)^2 -- variational  #
##########################

library(lvm4net)

timings[3] <- system.time(eucl2var <- lsm(Yg2, D) )[[3]]

########################
#######match rotations #######
########################

Z <- eucl2var$lsmEZ
Zm <- rotXtoY(euclmcmc$mkl$Z,Z)$X
Zv <- rotXtoY(euclvar$V_z,Z)$X

N <- nrow(Yg2)
Att<- rainbow(N,alpha=.8)

########################
######### LSM plots ########
########################

par(mfrow = c(1,3))

plotY(Yg2, Ndata = 1, EZ = Zm, main = 'latentnet', font.main = 1, cex.main = 1, colPl = Att)
plotY(Yg2, Ndata = 1, EZ = Zv, main = 'VBLPCM', font.main = 1, cex.main = 1, colPl = Att)
plotY(Yg2, Ndata = 1, EZ = Z, main = 'lvm4net', font.main = 1, cex.main = 1, colPl = Att)


########################
######### ROC plots ########
########################

Zmv <- list()
Zmv[[1]] <- Zm
Zmv[[2]] <- Zv
Zmv[[3]] <- Z


a <- boxroc(Yg2, 
	EZ = Zmv,
	xiT = c(euclmcmc$mkl$beta, euclvar$V_xi_e, eucl2var$xiT), 
	Lroc = 150, 
	ROC = TRUE, 
	BOXPLOT = TRUE,
	labelsPlot = c('latentnet', 'VBLPCM', 'lvm4net'),
	powdist = c(1, 1, 2),
	main =''
	)
res <- cbind(a$AUC, timings)
colnames(res) <- c('AUC', 'Timings in Sec')
res	




