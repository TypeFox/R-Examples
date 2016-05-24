SPF.BinBin <- function(x){ 

r_1_1_all <- r_m1_1_all <- r_0_1_all <- r_1_0_all <- r_m1_0_all <- r_0_0_all <- r_1_m1_all <- r_m1_m1_all <- r_0_m1_all <- Monotonicity_all <- NULL 
                         
Object <- x                                 

if (Object$Monotonicity == "No"){

pi_m1_m1 <- Object$Pi.Vectors$Pi_1010 
pi_1_m1 <- Object$Pi.Vectors$Pi_0110
pi_m1_1 <- Object$Pi.Vectors$Pi_1001
pi_1_1 <- Object$Pi.Vectors$Pi_0101
pi_0_0 <- Object$Pi.Vectors$Pi_0000 + Object$Pi.Vectors$Pi_0011 + Object$Pi.Vectors$Pi_1100 + Object$Pi.Vectors$Pi_1111
pi_0_m1 <- Object$Pi.Vectors$Pi_0010 + Object$Pi.Vectors$Pi_1110 
pi_m1_0 <- Object$Pi.Vectors$Pi_1000 + Object$Pi.Vectors$Pi_1011
pi_0_1 <- Object$Pi.Vectors$Pi_0001 + Object$Pi.Vectors$Pi_1101
pi_1_0 <- Object$Pi.Vectors$Pi_0100 + Object$Pi.Vectors$Pi_0111 

pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  

r_1_1 <- pi_1_1 / pi_1_Delta_S 
r_m1_1 <- pi_m1_1 / pi_1_Delta_S
r_0_1 <- pi_0_1 / pi_1_Delta_S 

r_1_0 <- pi_1_0 / pi_0_Delta_S 
r_m1_0 <- pi_m1_0 / pi_0_Delta_S
r_0_0 <- pi_0_0 / pi_0_Delta_S 

r_1_m1 <- pi_1_m1 / pi_m1_Delta_S 
r_m1_m1 <- pi_m1_m1 / pi_m1_Delta_S
r_0_m1 <- pi_0_m1 / pi_m1_Delta_S 

r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))

Monotonicity_all <- 
  cbind(Monotonicity_all, t(rep("No", times=length(pi_0_0))))
}

if (Object$Monotonicity == "Surr.Endp"){
  
  pi_m1_m1 <- 0 
  pi_1_m1 <- 0
  pi_m1_1 <- Object$Pi.Vectors$Pi_1001
  pi_1_1 <- Object$Pi.Vectors$Pi_0101
  pi_0_0 <- Object$Pi.Vectors$Pi_0000 + Object$Pi.Vectors$Pi_0011 + Object$Pi.Vectors$Pi_1100 + Object$Pi.Vectors$Pi_1111
  pi_0_m1 <- 0
  pi_m1_0 <- Object$Pi.Vectors$Pi_1000 + Object$Pi.Vectors$Pi_1011
  pi_0_1 <- Object$Pi.Vectors$Pi_0001 + Object$Pi.Vectors$Pi_1101
  pi_1_0 <- Object$Pi.Vectors$Pi_0100 + Object$Pi.Vectors$Pi_0111 
  
  pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
  pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
  pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  
  
  r_1_1 <- pi_1_1 / pi_1_Delta_S 
  r_m1_1 <- pi_m1_1 / pi_1_Delta_S
  r_0_1 <- pi_0_1 / pi_1_Delta_S 
  
  r_1_0 <- pi_1_0 / pi_0_Delta_S 
  r_m1_0 <- pi_m1_0 / pi_0_Delta_S
  r_0_0 <- pi_0_0 / pi_0_Delta_S 
  
  r_1_m1 <- rep(NA, times=length(r_1_0))
  r_m1_m1 <- rep(NA, times=length(r_1_0))
  r_0_m1 <- rep(NA, times=length(r_1_0))
  
  r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
  r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
  r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
  r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
  r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
  r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
  r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
  r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
  r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))
  
  Monotonicity_all <- 
    cbind(Monotonicity_all, t(rep("Surr.Endp", times=length(pi_0_0))))
}

if (Object$Monotonicity == "True.Endp"){
  
  pi_m1_m1 <- 0 
  pi_1_m1 <- Object$Pi.Vectors$Pi_0110
  pi_m1_1 <- 0
  pi_1_1 <- Object$Pi.Vectors$Pi_0101
  pi_0_0 <- Object$Pi.Vectors$Pi_0000 + Object$Pi.Vectors$Pi_0011 + Object$Pi.Vectors$Pi_1100 + Object$Pi.Vectors$Pi_1111
  pi_0_m1 <- Object$Pi.Vectors$Pi_0010 + Object$Pi.Vectors$Pi_1110 
  pi_m1_0 <- 0
  pi_0_1 <- Object$Pi.Vectors$Pi_0001 + Object$Pi.Vectors$Pi_1101
  pi_1_0 <- Object$Pi.Vectors$Pi_0100 + Object$Pi.Vectors$Pi_0111 
  
  pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
  pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
  pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  
  
  r_1_1 <- pi_1_1 / pi_1_Delta_S 
  r_m1_1 <- rep(NA, times=length(r_1_1))
  r_0_1 <- pi_0_1 / pi_1_Delta_S 
  
  r_1_0 <- pi_1_0 / pi_0_Delta_S 
  r_m1_0 <- rep(NA, times=length(r_1_1)) 
  r_0_0 <- pi_0_0 / pi_0_Delta_S 
  
  r_1_m1 <- pi_1_m1 / pi_m1_Delta_S 
  r_m1_m1 <- rep(NA, times=length(r_1_1))
  r_0_m1 <- pi_0_m1 / pi_m1_Delta_S 
  
  r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
  r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
  r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
  r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
  r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
  r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
  r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
  r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
  r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))
  
  Monotonicity_all <- 
    cbind(Monotonicity_all, t(rep("True.Endp", times=length(pi_0_0))))
}


if (Object$Monotonicity == "Surr.True.Endp"){
  
  pi_m1_m1 <- 0
  pi_1_m1 <- 0
  pi_m1_1 <- 0
  pi_1_1 <- Object$Pi.Vectors$Pi_0101
  pi_0_0 <- Object$Pi.Vectors$Pi_0000 + Object$Pi.Vectors$Pi_0011 + Object$Pi.Vectors$Pi_1100 + Object$Pi.Vectors$Pi_1111
  pi_0_m1 <- 0 
  pi_m1_0 <- 0
  pi_0_1 <- Object$Pi.Vectors$Pi_0001 + Object$Pi.Vectors$Pi_1101
  pi_1_0 <- Object$Pi.Vectors$Pi_0100 + Object$Pi.Vectors$Pi_0111 
  
  pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
  pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
  pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  
  
  r_1_1 <- pi_1_1 / pi_1_Delta_S 
  r_m1_1 <- rep(NA, times=length(r_1_1))
  r_0_1 <- pi_0_1 / pi_1_Delta_S 
  
  r_1_0 <- pi_1_0 / pi_0_Delta_S 
  r_m1_0 <- rep(NA, times=length(r_1_1))
  r_0_0 <- pi_0_0 / pi_0_Delta_S 
  
  r_1_m1 <- rep(NA, times=length(r_1_1)) 
  r_m1_m1 <- rep(NA, times=length(r_1_1))
  r_0_m1 <- rep(NA, times=length(r_1_1)) 
  
  r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
  r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
  r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
  r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
  r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
  r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
  r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
  r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
  r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))
  
  Monotonicity_all <- 
    cbind(Monotonicity_all, t(rep("Surr.True.Endp", times=length(pi_0_0))))
}


if (Object$Monotonicity == "General"){
  
  Object_temp <- Object$Pi.Vectors
  
  Object_No <- subset(Object_temp, Object_temp$Monotonicity %in% "No")
  Object_Surr <- subset(Object_temp, Object_temp$Monotonicity %in% "Surr")
  Object_True <- subset(Object_temp, Object_temp$Monotonicity %in% "True")
  Object_SurrTrue <- subset(Object_temp, Object_temp$Monotonicity %in% "SurrTrue")
  
  if (exists("Object_No") == TRUE) {
    
    pi_m1_m1 <- Object_No$Pi_1010 
    pi_1_m1 <- Object_No$Pi_0110
    pi_m1_1 <- Object_No$Pi_1001
    pi_1_1 <- Object_No$Pi_0101
    pi_0_0 <- Object_No$Pi_0000 + Object_No$Pi_0011 + Object_No$Pi_1100 + Object_No$Pi_1111
    pi_0_m1 <- Object_No$Pi_0010 + Object_No$Pi_1110 
    pi_m1_0 <- Object_No$Pi_1000 + Object_No$Pi_1011
    pi_0_1 <- Object_No$Pi_0001 + Object_No$Pi_1101
    pi_1_0 <- Object_No$Pi_0100 + Object_No$Pi_0111 
    
    pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
    pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
    pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  
    
    r_1_1 <- pi_1_1 / pi_1_Delta_S 
    r_m1_1 <- pi_m1_1 / pi_1_Delta_S
    r_0_1 <- pi_0_1 / pi_1_Delta_S 
    
    r_1_0 <- pi_1_0 / pi_0_Delta_S 
    r_m1_0 <- pi_m1_0 / pi_0_Delta_S
    r_0_0 <- pi_0_0 / pi_0_Delta_S 
    
    r_1_m1 <- pi_1_m1 / pi_m1_Delta_S 
    r_m1_m1 <- pi_m1_m1 / pi_m1_Delta_S
    r_0_m1 <- pi_0_m1 / pi_m1_Delta_S 
    
    r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
    r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
    r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
    r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
    r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
    r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
    r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
    r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
    r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))
    
    Monotonicity_all <- 
      cbind(Monotonicity_all, t(rep("No", times=length(pi_0_0))))
  }
  
  if (exists("Object_Surr") == TRUE) {
    
    pi_m1_m1 <- 0 
    pi_1_m1 <- 0
    pi_m1_1 <- Object_Surr$Pi_1001
    pi_1_1 <- Object_Surr$Pi_0101
    pi_0_0 <- Object_Surr$Pi_0000 + Object_Surr$Pi_0011 + Object_Surr$Pi_1100 + Object_Surr$Pi_1111
    pi_0_m1 <- 0
    pi_m1_0 <- Object_Surr$Pi_1000 + Object_Surr$Pi_1011
    pi_0_1 <- Object_Surr$Pi_0001 + Object_Surr$Pi_1101
    pi_1_0 <- Object_Surr$Pi_0100 + Object_Surr$Pi_0111 
    
    pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
    pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
    pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  
    
    r_1_1 <- pi_1_1 / pi_1_Delta_S 
    r_m1_1 <- pi_m1_1 / pi_1_Delta_S
    r_0_1 <- pi_0_1 / pi_1_Delta_S 
    
    r_1_0 <- pi_1_0 / pi_0_Delta_S 
    r_m1_0 <- pi_m1_0 / pi_0_Delta_S
    r_0_0 <- pi_0_0 / pi_0_Delta_S 
    
    r_1_m1 <- rep(NA, times=length(r_1_0))
    r_m1_m1 <- rep(NA, times=length(r_1_0))
    r_0_m1 <- rep(NA, times=length(r_1_0))
    
    r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
    r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
    r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
    r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
    r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
    r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
    r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
    r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
    r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))
    
    Monotonicity_all <- 
      cbind(Monotonicity_all, t(rep("Surr.Endp", times=length(pi_0_0))))
  }
  
  if (exists("Object_True") == TRUE) {
    
    pi_m1_m1 <- 0 
    pi_1_m1 <- Object_True$Pi_0110
    pi_m1_1 <- 0
    pi_1_1 <- Object_True$Pi_0101
    pi_0_0 <- Object_True$Pi_0000 + Object_True$Pi_0011 + Object_True$Pi_1100 + Object_True$Pi_1111
    pi_0_m1 <- Object_True$Pi_0010 + Object_True$Pi_1110 
    pi_m1_0 <- 0
    pi_0_1 <- Object_True$Pi_0001 + Object_True$Pi_1101
    pi_1_0 <- Object_True$Pi_0100 + Object_True$Pi_0111 
    
    pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
    pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
    pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  
    
    r_1_1 <- pi_1_1 / pi_1_Delta_S 
    r_m1_1 <- rep(NA, times=length(r_1_1))
    r_0_1 <- pi_0_1 / pi_1_Delta_S 
    
    r_1_0 <- pi_1_0 / pi_0_Delta_S 
    r_m1_0 <- rep(NA, times=length(r_1_1))
    r_0_0 <- pi_0_0 / pi_0_Delta_S 
    
    r_1_m1 <- pi_1_m1 / pi_m1_Delta_S 
    r_m1_m1 <- rep(NA, times=length(r_1_1))
    r_0_m1 <- pi_0_m1 / pi_m1_Delta_S 
    
    r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
    r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
    r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
    r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
    r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
    r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
    r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
    r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
    r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))
    
    Monotonicity_all <- 
      cbind(Monotonicity_all, t(rep("True.Endp", times=length(pi_0_0))))
  }
  
  
  if (exists("Object_SurrTrue") == TRUE) {
    
    pi_m1_m1 <- 0
    pi_1_m1 <- 0
    pi_m1_1 <- 0
    pi_1_1 <- Object_SurrTrue$Pi_0101
    pi_0_0 <- Object_SurrTrue$Pi_0000 + Object_SurrTrue$Pi_0011 + Object_SurrTrue$Pi_1100 + Object_SurrTrue$Pi_1111
    pi_0_m1 <- 0 
    pi_m1_0 <- 0
    pi_0_1 <- Object_SurrTrue$Pi_0001 + Object_SurrTrue$Pi_1101
    pi_1_0 <- Object_SurrTrue$Pi_0100 + Object_SurrTrue$Pi_0111 
    
    pi_m1_Delta_S <- pi_m1_m1 + pi_1_m1 + pi_0_m1
    pi_0_Delta_S <- pi_0_0 + pi_m1_0 + pi_1_0
    pi_1_Delta_S <- pi_m1_1 + pi_1_1 + pi_0_1  
    
    r_1_1 <- pi_1_1 / pi_1_Delta_S 
    r_m1_1 <- rep(NA, times=length(r_1_1))
    r_0_1 <- pi_0_1 / pi_1_Delta_S 
    
    r_1_0 <- pi_1_0 / pi_0_Delta_S 
    r_m1_0 <- rep(NA, times=length(r_1_1))
    r_0_0 <- pi_0_0 / pi_0_Delta_S 
    
    r_1_m1 <- rep(NA, times=length(r_1_1)) 
    r_m1_m1 <- rep(NA, times=length(r_1_1))
    r_0_m1 <- rep(NA, times=length(r_1_1)) 
    
    r_1_1_all <- cbind(r_1_1_all, t(r_1_1))
    r_m1_1_all <- cbind(r_m1_1_all, t(r_m1_1))
    r_0_1_all <- cbind(r_0_1_all, t(r_0_1))
    r_1_0_all <- cbind(r_1_0_all, t(r_1_0))
    r_m1_0_all <- cbind(r_m1_0_all, t(r_m1_0))
    r_0_0_all <- cbind(r_0_0_all, t(r_0_0))
    r_1_m1_all <- cbind(r_1_m1_all, t(r_1_m1))
    r_m1_m1_all <- cbind(r_m1_m1_all, t(r_m1_m1))
    r_0_m1_all <- cbind(r_0_m1_all, t(r_0_m1))
    
    Monotonicity_all <- 
      cbind(Monotonicity_all, t(rep("Surr.True.Endp", times=length(pi_0_0))))
  }
  
}


fit <- 
list(r_1_1=as.numeric(r_1_1_all), r_min1_1=as.numeric(r_m1_1_all), r_0_1=as.numeric(r_0_1_all),
     r_1_0=as.numeric(r_1_0_all), r_min1_0=as.numeric(r_m1_0_all), r_0_0=as.numeric(r_0_0_all),
     r_1_min1=as.numeric(r_1_m1_all), r_min1_min1=as.numeric(r_m1_m1_all), r_0_min1=as.numeric(r_0_m1_all),   
     Monotonicity=as.factor(Monotonicity_all), Call=match.call())      
                                 
class(fit) <- "SPF.BinBin"
fit
                                 
} 