### Spot rates functions

## Nelson/Siegel spot rate function
spr_ns <- function(beta, m){
  (beta[1] + beta[2]*((1-exp(-m/beta[4]))/(m/beta[4])) +
   beta[3]*(((1-exp(-m/beta[4]))/(m/beta[4]))-exp(-m/beta[4])))
}

## Svensson spot rate function 
spr_sv <- function(beta, m){
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
   beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
   beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6])))
}

## Adjusted Svensson spot rate function
spr_asv <- function(beta, m){
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
   beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
   beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-(2*m)/beta[6])))
}

## Diebold/Li spot rate function
spr_dl <- function(beta,m,lambda){
  (beta[1] + beta[2]*((1-exp(-m*lambda))/(m*lambda))+
   beta[3]*(((1-exp(-m*lambda))/(m*lambda))-exp(-m*lambda)))
}

## Spot rate wrapper function
spotrates <- function(method,beta,m,lambda = 0.0609*12){ 
  switch(method,
         "ns" = spr_ns(beta,m),
         "sv" = spr_sv(beta,m),
         "asv"= spr_asv(beta,m),
         "dl" = spr_dl(beta,m,lambda))
}

### Forward rate functions

## Nelson/Siegel forward rate function
fwr_ns <- function(beta,m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
   beta[3]*(m/beta[4]*exp(-m/beta[4])))
}

## Svensson forward rate function
fwr_sv <- function(beta, m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
   beta[3] *m/beta[4]*exp(-m/beta[4]) +
   beta[5] *m/beta[6]*exp(-m/beta[6]))
}

## Adjusted Svensson forward rate function
fwr_asv <- function(beta, m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
   beta[3] *m/beta[4]*exp(-m/beta[4]) +
   beta[5] *(exp(-m/beta[6])+(2*m/beta[6] -1)*exp(-2*m/beta[6])))
}

## Diebold/Li forward rate function
fwr_dl <- function(beta, m,lambda) {
  (beta[1] + beta[2]*exp(-m*lambda)
   + beta[3]*(m*lambda*exp(-m*lambda)))
}

## Forward rates wrapper function
forwardrates <- function(method,beta,m,lambda){
  switch(method,
         "ns" = fwr_ns(beta,m),
         "sv" = fwr_sv(beta,m),
         "asv"= fwr_asv(beta,m),
         "dl"= fwr_dl(beta,m,lambda))
}

## Implied foreward rates calculation
impl_fwr <- function(m,s) {
  impl_fwr <- c(s[1],(s[-1]*m[-1] - s[-length(s)]*m[-length(m)])/(diff(m)))
  impl_fwr[1] <- impl_fwr[2]
  impl_fwr	
}

get_paramnames <- function(method){
  names <- c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")
  switch(method,"ns"= names[1:4],"sv"=names,"asv"=names,"dl"=names[1:3])
}

get_realnames <- function(method){
  switch(method,"dl"="Diebold/Li","ns"="Nelson/Siegel","sv"="Svensson","asv"="Adj. Svensson")
}

### Loss function for parametric methods
get_objfct <- function(method) {
  objfct <- switch(method,
                   "dl" = objfct_dl,
                   "ns" = objfct_ns,
                   "sv" = objfct_sv,
                   "asv" = objfct_asv)
}

### Gradient of loss function for parametric methods
get_grad_objfct <- function(method) {
  grad_objfct <- switch(method,
                   "dl" = grad_dl,
                   "ns" = grad_ns,
                   "sv" = grad_sv,
                   "asv" = grad_asv)
}

### Diebold/Li loss function for yields
objfct_dl <- function(beta, lambda, m, y)
      {
        sum((y - spr_dl(beta,m, lambda))^2)
      }

### Nelson/Siegel loss function for yields
objfct_ns <- function(beta, m, y)
      {
        sum((y - spr_ns(beta,m))^2)
      }

### Nelson/Siegel grid loss function for yields
objfct_ns_grid <- function(beta, tau, m, y)
      {
        betans <- c(beta, tau)
        sum((y - spr_ns(betans,m))^2)
      }

### Svensson loss function for yields
objfct_sv <- function(beta, m, y)
      {
        sum((y - spr_sv(beta,m))^2)
      }


### Svensson grid loss function for yields
objfct_sv_grid <- function(beta, tau,  m, y)
      {
        betasv <- c(beta[1:3], tau[1], beta[4], tau[2])
        sum((y - spr_sv(betasv,m))^2)
      }

### Adjusted Svensson loss function for yields
objfct_asv <- function(beta, m, y)
      {
        sum((y - spr_asv(beta,m))^2)
      }


### Adjusted Svensson grid loss function for yields
objfct_asv_grid <- function(beta, tau, m, y)
      {
        betasv <- c(beta[1:3], tau[1], beta[4], tau[2])
        sum((y - spr_asv(betasv,m))^2)
      }

### Constraints for constrOptim()

get_constraints <- function(method, tauconstr) {

  ## tauconstr = c(upper, lower, gridsize, distance)
  
  ## Diebold/Li
  
  if (method == "dl") {
    ui <- rbind(c(1,0,0),               # beta0 > 0
                c(1,1,0))               # beta0 + beta1 > 0
    ci <- c(0,0)
   }
  
  ## Nelson/Siegel
  
  if (method == "ns") {
    ui <- rbind(c(1,0,0,0),             # beta0 > 0
                c(1,1,0,0),             # beta0 + beta1 > 0
                c(0,0,0,1),             # tau1 > tauconstr[1]
                c(0,0,0,-1))            # tau1 <= tauconstr[2]
    ci <- c(0,0,tauconstr[1],-tauconstr[2])
  }

  ## Svensson
  
  if (method == "sv") {
     ui <- rbind(c(1,0,0,0,0,0),        # beta0 > 0
                 c(1,1,0,0,0,0),        # beta0 + beta1 > 0
                 c(0,0,0,1,0,0),        # tau1 > tauconstr[1]
                 c(0,0,0,-1,0,0),       # tau1 <= tauconstr[2]
                 c(0,0,0,0,0,1),        # tau2 > tauconstr[1]
                 c(0,0,0,0,0,-1),       # tau1 <= tauconstr[2]
                 c(0,0,0,-1,0,1))       # tau2 - tau1 > tauconstr[4]
     ci <- c(0,0,tauconstr[1],-tauconstr[2],tauconstr[1],-tauconstr[2],tauconstr[4]) 
   }

## Adjusted Svensson
  
  if (method == "asv") {
     ui <- rbind(c(1,0,0,0,0,0),        # beta0 > 0
                 c(1,1,0,0,0,0),        # beta0 + beta1 > 0
                 c(0,0,0,1,0,0),        # tau1 > tauconstr[1]
                 c(0,0,0,-1,0,0),       # tau1 <= tauconstr[2]
                 c(0,0,0,0,0,1),        # tau2 > tauconstr[1]
                 c(0,0,0,0,0,-1),       # tau1 <= tauconstr[2]
                 c(0,0,0,-1,0,1))       # tau2 - tau1 > 0
     ci <- c(0,0,tauconstr[1],-tauconstr[2],tauconstr[1],-tauconstr[2],0) 
   }
  
  constraints <- list(ui = ui, ci = ci)
  constraints
}

### Loss function for parametric methods
get_objfct_bonds <- function(method) {
  objfct <- switch(method,
                   "dl" = objfct_dl_bonds,
                   "ns" = objfct_ns_bonds,
                   "sv" = objfct_sv_bonds,
                   "asv" = objfct_asv_bonds)
}

### Gradient of loss function for parametric methods
get_grad_objfct_bonds <- function(method) {
  grad_objfct <- switch(method,
                   "dl" = grad_dl_bonds,
                   "ns" = grad_ns_bonds,
                   "sv" = grad_sv_bonds,
                   "asv" = grad_asv_bonds)
}

### Diebold/Li loss function for bonds 
objfct_dl_bonds <- function(beta, lambda, m, cf, w, p) {
      phat <- bond_prices("dl",beta,m,cf, lambda)$bond_prices
      sum(w*((p - phat)^2))
    }

### Nelson/Siegel loss function for bonds 
objfct_ns_bonds <- function(beta, m, cf, w, p) {
  .Call("objfct_ns_bonds_Cpp", beta, m, cf, w, p)
    }

### Nelson/Siegel grid loss function for bonds

objfct_ns_bonds_grid <- function(beta, tau, m, cf, w, p) {
   .Call("objfct_ns_bonds_gridCpp", beta, tau, m, cf, w, p)
    }

### Svensson loss function for bonds 
objfct_sv_bonds <- function(beta, m, cf, w, p) {
  .Call("objfct_sv_bonds_Cpp", beta, m, cf, w, p)
    }

### Svensson grid loss function for bonds
objfct_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
  .Call("objfct_sv_bonds_gridCpp", beta, tau, m, cf, w, p)
    }

### Adjusted Svensson loss function for bonds 
objfct_asv_bonds <- function(beta, m, cf, w, p) {
   .Call("objfct_asv_bonds_Cpp", beta, m, cf, w, p)
    }

### Adjusted Svensson grid loss function for bonds
objfct_asv_bonds_grid <- function(beta, tau, m, cf, w, p) {
  .Call("objfct_asv_bonds_gridCpp", beta, tau, m, cf, w, p)
}




