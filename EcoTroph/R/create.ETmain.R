create.ETmain <- function (ecopath, smooth_type = NULL, sigmaLN_cst = NULL, pas = NULL, shift = NULL, smooth_param = NULL) 
{
B <- apply(biomass <- Transpose(tab_smooth <- create.smooth(tab_input = ecopath, smooth_type, sigmaLN_cst, pas, shift, smooth_param), ecopath, "biomass"), 1, sum)
B_acc <- apply(biomass_acc <- sweep(biomass, 2, ecopath$accessibility, FUN = "*"), 1, sum)
P <- apply(flowP <- sweep(biomass, 2, ecopath$prod, FUN = "*"), 1, sum)
P_acc <- apply(flowP_acc <- sweep(flowP, 2, ecopath$accessibility, FUN = "*"), 1, sum)
Kin <- P/B
Kin_acc <- P_acc/B_acc

Y <- list()
somme_pecheries <- biomass
somme_pecheries[] <- 0
for (pecheries in colnames(ecopath)[grep("catch", colnames(ecopath))]) 
{
Y[[paste(pecheries)]] <- Transpose(tab_smooth, ecopath, pecheries)
somme_pecheries <- somme_pecheries + Y[[paste(pecheries)]]
}
Y_tot <- apply(somme_pecheries, 1, sum)

F_loss <- Y_tot/P
F_loss_acc <- Y_tot/P_acc
V <- as.numeric(rownames(biomass))[-1] - as.numeric(rownames(biomass))[-nrow(biomass)]
N_loss <- c(log(P[-length(P)]/P[-1])/V - F_loss[-length(F_loss)], NA)
N_loss[1]=log(P[1]/P[2]*V[2])-F_loss[1]

Fish_mort <- Y_tot/B
Fish_mort_acc <- Fish_mort/(B_acc/B)
Selec <- B_acc/B
Time <- cumsum(c(0, V/Kin[-length(Kin)]))
N_loss_acc <- c(log(P_acc[-length(P_acc)]/P_acc[-1])/V - F_loss_acc[-length(F_loss_acc)], NA)
N_loss_acc[1]=log(P_acc[1]/P_acc[2]*V[2])-F_loss_acc[1]
ET_Main <- cbind(B, B_acc, P, P_acc, Kin, Kin_acc, Y_tot, F_loss, F_loss_acc, N_loss, Fish_mort, Fish_mort_acc, Selec, Time, N_loss_acc)
retour<-list(ET_Main = as.data.frame(ET_Main), biomass = biomass, biomass_acc = biomass_acc, prod = flowP, prod_acc = flowP_acc, tab_smooth = tab_smooth,Y=Y)

class(retour)<-"ETmain"
return(retour)
}
