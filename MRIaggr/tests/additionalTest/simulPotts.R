#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%
#%%%%%%%  Potts simulation
#%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


require(MRIaggr)

#### parametrization ####
# spatial field
n <- 25
G <- 3
coords <- data.frame(which(matrix(0,nrow = n*G,ncol = n*G) == 0,arr.ind = TRUE), 1)
optionsMRIaggr(legend = FALSE, quantiles.legend = FALSE,axes = FALSE,num.main = FALSE)

# neighborhood matrix
resW <- calcW(coords,range = sqrt(2),row.norm = TRUE,calcBlockW = TRUE)
W <- resW$W
site_order <- unlist(resW$blocks$ls_groups) - 1

# parameters
seq_rho <- seq(0.5,8,by = 0.5)
seq_iter_max <- c(10,50,100,500,1000,2500,5000,7500,10000)
n.seq_rho <- length(seq_rho)
n.seq_iter_max <- length(seq_iter_max)
n.rep <- 10

# storage
res_V <- list()
res_U <- array(NA,dim = c(n.seq_rho,n.seq_iter_max,n.rep),
              dimnames = list(seq_rho,seq_iter_max))

#### loop - Local Potts model ####

for (iter_rho in 1:n.seq_rho) {
  set.seed(10)
  rho <- seq_rho[iter_rho]
  cat(rho," : ")
  
  res_V[[iter_rho]] <- matrix(NA,nrow = (n*G) ^ 2,ncol = n.seq_iter_max)
  colnames(res_V[[iter_rho]]) <- seq_iter_max
  
  for (iter_rep in 1:n.rep) {
    cat("*")
    
    res_tempo <- simulPotts(W, G = 3,rho = rho, iter_max = seq_iter_max[1], site_order = site_order,
                            verbose = FALSE, fast = TRUE)$simulation
    
    res_V[[iter_rho]][,1]  <- apply(res_tempo, 1, which.max)
    
    for (iter_simul in 2:n.seq_iter_max) {
      
      res_tempo <-  simulPotts(W, G = 3, rho = rho, initialization = res_tempo, iter_max = seq_iter_max[iter_simul], 
                               site_order = site_order, verbose = FALSE, fast = TRUE)$simulation
      
      res_V[[iter_rho]][,iter_simul]  <- apply(res_tempo, 1, which.max)
    }
    
    res_U[iter_rho,,iter_rep] <- apply(res_V[[iter_rho]],2,function(x){
      sum(sapply(1:G,function(g){ mean( as.numeric(x == g) * W %*% as.numeric(x == g) ) }))
    })
  }
  cat("\n")
}
names(res_V) <- seq_rho

#### export ####
save(res_V,file = "res_V.RData")
save(res_U,file = "res_U.RData")

# path_data <- "E:/Creation_package/Package_MRIaggr/MRIaggr/inst/Data_Tests"
# load(file.path(path_data,"res_U.RData"))
# load(file.path(path_data,"res_V.RData"))

### display ####
# require(lattice)
# require(reshape2)
# require(ggplot2)

# multiplot(coords, res_V[["0.5"]][,"1000"], main = "rho 0.5")
# multiplot(coords, res_V[["3"]][,"1000"], main = "rho 3")
# multiplot(coords, res_V[["4"]][,"1000"], main = "rho 4")
# multiplot(coords, res_V[["5"]][,"1000"], main = "rho 5")
# multiplot(coords, res_V[["6"]][,"1000"], main = "rho 6")

### potential
# res_medianU <- apply(res_U,c(1,2),median)
# res_005U <- apply(res_U,c(1,2),quantile,probs = 0.05)
# res_095U <- apply(res_U,c(1,2),quantile,probs = 0.95)

# dfU <- reshape(data.frame(rho = rownames(res_medianU),res_medianU), 
        # varying = list(2:(n.seq_iter_max + 1)), 
        # timevar = "iter_max", times = seq_iter_max , v.names = "U",
        # idvar = "rho", direction = "long")
        
# dfU$Uinf <- reshape(data.frame(rho = rownames(res_005U),res_005U), 
               # varying = list(2:(n.seq_iter_max + 1)), 
               # timevar = "iter_max", times = seq_iter_max , v.names = "Uinf",
               # idvar = "rho", direction = "long")$Uinf

# dfU$Usup <- reshape(data.frame(rho = rownames(res_095U),res_095U), 
               # varying = list(2:(n.seq_iter_max + 1)), 
               # timevar = "iter_max", times = seq_iter_max , v.names = "Usup",
               # idvar = "rho", direction = "long")$Usup

# rownames(dfU) <- 1:nrow(dfU)

# dfU$iter_max.num <- as.numeric(as.factor(dfU$iter_max))
# dfU$rho.num <- as.numeric(as.factor(dfU$rho))

# Ulist.axis <- list(arrows = FALSE,  
                   # x = list(at = seq(1,n.seq_iter_max,2), lab = seq_iter_max[seq(1,n.seq_iter_max,2)],cex =  1),
                   # y = list(at = seq(1,n.seq_rho,2), lab = seq_rho[seq(1,n.seq_rho,2)],cex = 1) , 
                   # z = list(at = seq(0,1,0.2),  lab = seq(0,1,0.2),cex = 1)
# )

## 3D plot
# lattice::wireframe(U ~ iter_max.num + rho.num,data = dfU,
                   # scales = Ulist.axis,
                   # xlab = list("Nombre d iterations",rot = 45),
                   # ylab = list(expression(rho),rot = 0),
                   # zlab = list("Potentiel",rot = 90),
                   # zlim = c(0,1),
                   # at = seq(0,1,0.05), col.regions = terrain.colors(20),
                   # drape = TRUE, colorkey = TRUE,
                   # screen = list(z = 60, x = -70, y = 2.5))

## boxplot
# boxplot(t(res_U[1,,]),ylim = c(min(res_U),1))
# for (iter_rho in 2:nrow(res_U)) {
  # boxplot(t(res_U[iter_rho,,]),add = TRUE)
# }

## lines
# ggU <- ggplot2::ggplot(data = dfU,aes(x = as.factor(iter_max),
                             # y = U, 
                             # group = as.factor(rho),
                             # col = rho))
# ggU <- ggU + ggplot2::geom_point() + ggplot2::geom_line()
# ggU <- ggU + ggplot2::labs(title = "", x = "number of iterations", y = "local potential")
# ggU <- ggU + ggplot2::geom_errorbar(aes(ymin = Uinf, ymax = Usup, col = rho),
                            # width = 0.1)
# ggU

# ggU + ggplot2::coord_cartesian(ylim = c(0.9,1))
# 2500 iteration if rho> 3.5

