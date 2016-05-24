#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%
#%%%%%%%  Potts parameter estimation 
#%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cd ../../data1/etude/etude25/2891/

require(MRIaggr)
require(snowfall)

multipar(legend = FALSE)

#### parametrization ####
cpus <- 2

mode <- 2

seq_n <- switch(mode,
                "1" = c(10,20,30,50,75,100,175),
                "2" = c(10,20,30)
)
n.n <- length(seq_n)

seq_rho <- switch(mode,
                  "1" = seq(0,10,length.out = 15),
                  "2" = seq(0,10,length.out = 4)
)
n.rho <- length(seq_rho)

n.rep <- switch(mode,
                "1" = 250,
                "2" = 10
) 

iter_max <- switch(mode,
                "1" = 1000,
                "2" = "auto"
)

# number of groups
G <- 3 

#### initialisation ####
distband_SR <- sqrt(2) + 0.001

rho_hat <- array(NA,dim = c(n.rep,n.rho,n.n),dimnames = list(1:n.rep,seq_rho,seq_n))
rho_bias <- array(NA,dim = c(n.rep,n.rho,n.n),dimnames = list(1:n.rep,seq_rho,seq_n))
rho_rbias <- array(NA,dim = c(n.rep,n.rho,n.n),dimnames = list(1:n.rep,seq_rho,seq_n))
Usimul <- array(NA,dim = c(n.rep,n.rho,n.n),dimnames = list(1:n.rep,seq_rho,seq_n))

set.seed(10)


#### fct ####

estimateN_rho <- function(iter_n,trace=TRUE, save = TRUE){

  Mrho_hat <- matrix(NA, nrow = n.rep , ncol = n.rho,dimnames = list(1:n.rep,seq_rho))
  Mrho_bias <- matrix(NA, nrow = n.rep , ncol = n.rho,dimnames = list(1:n.rep,seq_rho))
  Mrho_rbias <- matrix(NA, nrow = n.rep , ncol = n.rho,dimnames = list(1:n.rep,seq_rho))
  MUsimul <- matrix(NA, nrow = n.rep , ncol = n.rho,dimnames = list(1:n.rep,seq_rho))
  
  if (trace == TRUE) {
  cat("****************************************** \n") 
  cat(iter_n,") field n=",seq_n[iter_n],"\n") 
  }
  
  ### field ####
  n_px <- seq_n[iter_n]
  n <- (n_px*G) ^ 2 # nombre de pixels total
  
  M_coords <- matrix(1,nrow = n_px*G,ncol = n_px*G)
  coords <- which(M_coords == 1,arr.ind = T)
  
  #### neighborhood matrix ####
  resW <- calcW(as.data.frame(coords),range = distband_SR,upper = NULL, row.norm = TRUE ,calcBlockW = TRUE)
  site_order <- unlist(resW$blocks$ls_groups) - 1
    
  #### simulation and estimation ####
  
  for (iter_rho in 1:n.rho) {
    
    rho <- seq_rho[iter_rho]
	if(iter_max=="auto"){	
	iter_max.sample <- if(rho<3.5){100} else if(rho<4){1000} else{2500}	
	}else{
	iter_max.sample <- iter_max
	}
    if (trace == TRUE) {cat("theoric rho :",rho,"\n")}
    
    for (iter_rep in 1:n.rep) {
      if (trace == TRUE) {cat("*")}
      
      sample <- simulPottsFast_cpp(W_i = resW$W@i, W_p = resW$W@p, W_x = resW$W@x, 
                         site_order = site_order, 
                         sample = t(stats::rmultinom(n, size = 1, prob = rep(1/G,G))), 
                         rho = rho, n = n, p = G, iter_nb = iter_max.sample)
      
      Mrho_hat[iter_rep,iter_rho] <- rhoMF(Y = sample,W_SR = resW$W)
      
      MUsimul[iter_rep,iter_rho] <- sum(sample*(resW$W %*% sample))      #sum(sapply(1:G,function(g){ sum( as.numeric(sample.vec == g) * resW$W %*% as.numeric(sample.vec == g) ) }))/n
  
    }
    if (trace == TRUE) {cat("\n \n")}
  }
  
  Mrho_bias <- Mrho_hat - matrix(seq_rho, byrow = T,ncol = n.rho, nrow = n.rep)
  Mrho_rbias <- Mrho_hat * matrix(1/seq_rho, byrow = T,ncol = n.rho, nrow = n.rep)
  
  if (any(Mrho_hat > 30)) {
    Mrho_bias[Mrho_hat > 30] <- NA
    Mrho_rbias[Mrho_hat > 30] <- NA
  }
  
  #### export ####
  
  if (save == TRUE) {
  save(Mrho_hat,file = paste("MRIaggrTests-rho_hat(tempo",iter_n,").RData",sep = ""))
  save(Mrho_bias,file = paste("MRIaggrTests-rho_bias(tempo",iter_n,").RData",sep = ""))
  save(Mrho_rbias,file = paste("MRIaggrTests-rho_rbias(tempo",iter_n,").RData",sep = ""))
  save(MUsimul,file = paste("MRIaggrTests-Usimul(tempo",iter_n,").RData",sep = ""))
  }
  
  return(list(Mrho_hat = Mrho_hat,
              Mrho_bias = Mrho_bias,
              Mrho_rbias = Mrho_rbias,
              MUsimul = MUsimul,
              n = n_px))
}

#### loop ###
snowfall::sfInit(parallel = TRUE, cpus = cpus)
snowfall::sfExport("estimateN_rho",
                   "distband_SR","G","iter_max",
                   "seq_n","n.n","n.rep","seq_rho","n.rho")
sfLibrary( MRIaggr )

system.time(
listRes <- snowfall::sfClusterApplyLB(seq(n.n,1,-1),function(iter_n){
  estimateN_rho(iter_n,trace = TRUE, save = FALSE)
})
)
snowfall::sfStop()

## computation time 
# (n=75) (6*5.5 + 1*50 + 9*112)*250/3600/24
# (n=100) (6*10 + 1*79 + 9*208)*250/3600/24

#### export ####
for (iter_n in 1:n.n) {
 index_ls <-  which(unlist(lapply(listRes,"[[","n")) == seq_n[iter_n])
 
 rho_hat[,,iter_n] <- listRes[[index_ls]]$Mrho_hat 
 rho_bias[,,iter_n] <- listRes[[index_ls]]$Mrho_bias
 rho_rbias[,,iter_n] <- listRes[[index_ls]]$Mrho_rbias
 Usimul[,,iter_n] <- listRes[[index_ls]]$MUsimul
}



save(rho_hat,file = paste("MRIaggrTests-rho_hat.RData",sep = ""))
save(rho_bias,file = paste("MRIaggrTests-rho_bias.RData",sep = ""))
save(rho_rbias,file = paste("MRIaggrTests-rho_rbias.RData",sep = ""))
save(Usimul,file = paste("MRIaggrTests-Usimul.RData",sep = ""))
