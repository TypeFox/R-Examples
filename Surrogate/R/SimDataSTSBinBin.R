Sim.Data.STSBinBin <- function(Monotonicity=c("No"), 
                               N.Total=2000, Seed=sample(1:1000, size=1)){
  
  
  if (Monotonicity=="No"){
    set.seed(seed=Seed)
    Pi_s <- 
      RandVec(a=0, b=1, s=1, n=16, m=1)  
    Pi_0000 <- Pi_s$RandVecOutput[1]  
    Pi_0100 <- Pi_s$RandVecOutput[2]  
    Pi_0010 <- Pi_s$RandVecOutput[3]  
    Pi_0001 <- Pi_s$RandVecOutput[4]  
    Pi_0101 <- Pi_s$RandVecOutput[5]  
    Pi_1000 <- Pi_s$RandVecOutput[6]  
    Pi_1010 <- Pi_s$RandVecOutput[7]  
    Pi_1001 <- Pi_s$RandVecOutput[8]  
    Pi_1110 <- Pi_s$RandVecOutput[9]  
    Pi_1101 <- Pi_s$RandVecOutput[10]  
    Pi_1011 <- Pi_s$RandVecOutput[11]  
    Pi_1111 <- Pi_s$RandVecOutput[12]  
    Pi_0110 <- Pi_s$RandVecOutput[13]  
    Pi_0011 <- Pi_s$RandVecOutput[14]  
    Pi_0111 <- Pi_s$RandVecOutput[15]  
    Pi_1100 <- Pi_s$RandVecOutput[16] 
    Pi_s_all <- cbind(Pi_0000, Pi_0100, Pi_0010, Pi_0001, Pi_0101, Pi_1000, Pi_1010, Pi_1001, Pi_1110, Pi_1101, 
                      Pi_1011, Pi_1111, Pi_0110, Pi_0011, Pi_0111, Pi_1100)
  }
  
  if (Monotonicity=="True.Endp"){
    set.seed(seed=Seed)
    Pi_s <- RandVec(a=0, b=1, s=1, n=12, m=1)  
    Pi_0000 <- Pi_s$RandVecOutput[1]  
    Pi_0100 <- Pi_s$RandVecOutput[2]  
    Pi_0010 <- Pi_s$RandVecOutput[3]  
    Pi_0001 <- Pi_s$RandVecOutput[4]  
    Pi_0101 <- Pi_s$RandVecOutput[5]  
    Pi_1000 <- c(0)  
    Pi_1010 <- c(0)   
    Pi_1001 <- c(0)   
    Pi_1110 <- Pi_s$RandVecOutput[6]  
    Pi_1101 <- Pi_s$RandVecOutput[7]  
    Pi_1011 <- c(0) 
    Pi_1111 <- Pi_s$RandVecOutput[8]  
    Pi_0110 <- Pi_s$RandVecOutput[9]  
    Pi_0011 <- Pi_s$RandVecOutput[10]  
    Pi_0111 <- Pi_s$RandVecOutput[11]  
    Pi_1100 <- Pi_s$RandVecOutput[12]   
    Pi_s_all <- cbind(Pi_0000, Pi_0100, Pi_0010, Pi_0001, Pi_0101, Pi_1000, Pi_1010, Pi_1001, Pi_1110, Pi_1101, 
                      Pi_1011, Pi_1111, Pi_0110, Pi_0011, Pi_0111, Pi_1100)
  }  
  
  
  if (Monotonicity=="Surr.Endp"){
    set.seed(seed=Seed)
    Pi_s <- RandVec(a=0, b=1, s=1, n=12, m=1)  
    Pi_0000 <- Pi_s$RandVecOutput[1]  
    Pi_0100 <- Pi_s$RandVecOutput[2]  
    Pi_0010 <- c(0)  
    Pi_0001 <- Pi_s$RandVecOutput[3]  
    Pi_0101 <- Pi_s$RandVecOutput[4]  
    Pi_1000 <- Pi_s$RandVecOutput[5]  
    Pi_1010 <- c(0)  
    Pi_1001 <- Pi_s$RandVecOutput[6]  
    Pi_1110 <- c(0)  
    Pi_1101 <- Pi_s$RandVecOutput[7]  
    Pi_1011 <- Pi_s$RandVecOutput[8]  
    Pi_1111 <- Pi_s$RandVecOutput[9]  
    Pi_0110 <- c(0)  
    Pi_0011 <- Pi_s$RandVecOutput[10]  
    Pi_0111 <- Pi_s$RandVecOutput[11]  
    Pi_1100 <- Pi_s$RandVecOutput[12] 
    Pi_s_all <- cbind(Pi_0000, Pi_0100, Pi_0010, Pi_0001, Pi_0101, Pi_1000, Pi_1010, Pi_1001, Pi_1110, Pi_1101, 
                      Pi_1011, Pi_1111, Pi_0110, Pi_0011, Pi_0111, Pi_1100)
    
  }  
  
  if (Monotonicity=="Surr.True.Endp"){
    set.seed(seed=Seed)
    Pi_s <- RandVec(a=0, b=1, s=1, n=9, m=1)  
    Pi_0000 <- Pi_s$RandVecOutput[1]  
    Pi_0100 <- Pi_s$RandVecOutput[2]  
    Pi_0010 <- c(0)  
    Pi_0001 <- Pi_s$RandVecOutput[3]  
    Pi_0101 <- Pi_s$RandVecOutput[4]  
    Pi_1000 <- c(0)  
    Pi_1010 <- c(0)  
    Pi_1001 <- c(0)  
    Pi_1110 <- c(0)  
    Pi_1101 <- Pi_s$RandVecOutput[5]  
    Pi_1011 <- c(0)  
    Pi_1111 <- Pi_s$RandVecOutput[6]  
    Pi_0110 <- c(0)  
    Pi_0011 <- Pi_s$RandVecOutput[7]  
    Pi_0111 <- Pi_s$RandVecOutput[8]  
    Pi_1100 <- Pi_s$RandVecOutput[9] 
    Pi_s_all <- cbind(Pi_0000, Pi_0100, Pi_0010, Pi_0001, Pi_0101, Pi_1000, Pi_1010, Pi_1001, Pi_1110, Pi_1101, 
                      Pi_1011, Pi_1111, Pi_0110, Pi_0011, Pi_0111, Pi_1100)
    
  }  
  
  Pi_0000_ma <- matrix(rep(c(0, 0, 0, 0), round(Pi_0000*N.Total)), ncol=4, byrow=T)
  Pi_0100_ma <- matrix(rep(c(0, 1, 0, 0), round(Pi_0100*N.Total)), ncol=4, byrow=T)
  Pi_0010_ma <- matrix(rep(c(0, 0, 1, 0), round(Pi_0010*N.Total)), ncol=4, byrow=T)
  Pi_0001_ma <- matrix(rep(c(0, 0, 0, 1), round(Pi_0001*N.Total)), ncol=4, byrow=T)
  Pi_0101_ma <- matrix(rep(c(0, 1, 0, 1), round(Pi_0101*N.Total)), ncol=4, byrow=T)
  Pi_1000_ma <- matrix(rep(c(1, 0, 0, 0), round(Pi_1000*N.Total)), ncol=4, byrow=T)
  Pi_1010_ma <- matrix(rep(c(1, 0, 1, 0), round(Pi_1010*N.Total)), ncol=4, byrow=T)
  Pi_1001_ma <- matrix(rep(c(1, 0, 0, 1), round(Pi_1001*N.Total)), ncol=4, byrow=T)
  Pi_1110_ma <- matrix(rep(c(1, 1, 1, 0), round(Pi_1110*N.Total)), ncol=4, byrow=T)
  Pi_1101_ma <- matrix(rep(c(1, 1, 0, 1), round(Pi_1101*N.Total)), ncol=4, byrow=T)
  Pi_1011_ma <- matrix(rep(c(1, 0, 1, 1), round(Pi_1011*N.Total)), ncol=4, byrow=T)
  Pi_1111_ma <- matrix(rep(c(1, 1, 1, 1), round(Pi_1111*N.Total)), ncol=4, byrow=T)
  Pi_0110_ma <- matrix(rep(c(0, 1, 1, 0), round(Pi_0110*N.Total)), ncol=4, byrow=T)
  Pi_0011_ma <- matrix(rep(c(0, 0, 1, 1), round(Pi_0011*N.Total)), ncol=4, byrow=T)
  Pi_0111_ma <- matrix(rep(c(0, 1, 1, 1), round(Pi_0111*N.Total)), ncol=4, byrow=T)
  Pi_1100_ma <- matrix(rep(c(1, 1, 0, 0), round(Pi_1100*N.Total)), ncol=4, byrow=T)
  
  mat <- 
    data.frame(rbind(Pi_0000_ma, Pi_0100_ma, Pi_0010_ma, Pi_0001_ma, Pi_0101_ma, Pi_1000_ma, Pi_1010_ma, 
                     Pi_1001_ma, Pi_1110_ma, Pi_1101_ma, Pi_1011_ma, Pi_1111_ma, Pi_0110_ma, Pi_0011_ma, Pi_0111_ma, Pi_1100_ma))
  
  colnames(mat) <- c("T0", "T1", "S0", "S1")
  
  set.seed(Seed) 
  Z <- rbinom(dim(mat)[1], 1, 0.5)  
  Z[Z==0] <- c(-1)
  mat <- cbind(mat, Z)
  
  mat_obs <- data.frame(matrix(NA, nrow=dim(mat)[1], ncol=3))
  colnames(mat_obs) <- c("T", "S", "Z")
  for (i in 1: dim(mat)[1]){    
    if (mat$Z[i]==-1) {
      mat_obs$T[i] <- mat$T0[i]
      mat_obs$S[i] <- mat$S0[i]
      mat_obs$Z[i] <- mat$Z[i]
    }
    if (mat$Z[i]==1) {
      mat_obs$T[i] <- mat$T1[i]
      mat_obs$S[i] <- mat$S1[i]
      mat_obs$Z[i] <- mat$Z[i]
    }
  }
  
  colnames(mat_obs) <- c("T", "S", "Z")
  
  if ((dim(mat_obs)[1]) != N.Total) {
    cat("\nNOTE: The number of patients requested in the function call equals ", N.Total, ", but the actual number of generated ", sep="")
    cat("\nobservations was ", (dim(mat)[1]), " (due to rounding).", sep="")
  } 
  
  Data.STSBinBin_Obs <- Data.STSBinBin_Counterfactuals <- NULL
  Data.STSBinBin_Obs <<- mat_obs 
  Data.STSBinBin_Counterfactuals <<- mat
  
  # Marginal probs
  pi1_1_ <- Pi_s_all[7]+Pi_s_all[9]+Pi_s_all[11]+Pi_s_all[12]
  pi1_0_ <- Pi_s_all[6]+Pi_s_all[8]+Pi_s_all[10]+Pi_s_all[16] 
  pi_1_1 <- Pi_s_all[5]+Pi_s_all[10]+Pi_s_all[12]+Pi_s_all[15] 
  pi_1_0 <- Pi_s_all[2]+Pi_s_all[9]+Pi_s_all[13]+Pi_s_all[16]
  pi0_1_ <- Pi_s_all[3]+Pi_s_all[13]+Pi_s_all[14]+Pi_s_all[15]
  pi_0_1 <- Pi_s_all[4]+Pi_s_all[8]+Pi_s_all[11]+Pi_s_all[14]
  Pi_Marginals <- cbind(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1)
  
  # True pars 
  Pi_s <- Pi_s_all
  mat1 <- Pi_s[7] 
  mat2 <- Pi_s[3] + Pi_s[9] 
  mat3 <- Pi_s[13] 
  mat4 <- Pi_s[6] + Pi_s[11] 
  mat5 <- Pi_s[1] + Pi_s[14] + Pi_s[16] + Pi_s[12] 
  mat6 <- Pi_s[2] + Pi_s[15]
  mat7 <- Pi_s[8]
  mat8 <- Pi_s[4] + Pi_s[10] 
  mat9 <- Pi_s[5]
  
  Delta_c_mat <-
    matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)
  
  #C3
#  pi_a <- mat1+mat5+mat9
#  pi_e <- ((mat1+mat4+mat7)*(mat1+mat2+mat3))+((mat2+mat5+mat8)*(mat4+mat5+mat6))+((mat3+mat6+mat9)*(mat7+mat8+mat9))
#  kappa <- (pi_a - pi_e)/(1-pi_e)
#  pi_max <- 
#    min(mat1+mat4+mat7, mat1+mat2+mat3) + min(mat2+mat5+mat8, mat4+mat5+mat6) + min(mat3+mat6+mat9, mat7+mat8+mat9)
#  k_max <- (pi_max - pi_e)/(1-pi_e)
#  C3 <- kappa/k_max
  
  #R2_H
  sum_S_min1 <- mat1+mat2+mat3
  sum_S_0 <- mat4+mat5+mat6
  sum_S_1 <- mat7+mat8+mat9
  
  sum_T_min1 <- mat1+mat4+mat7
  sum_T_0 <- mat2+mat5+mat8
  sum_T_1 <- mat3+mat6+mat9
  
  
  if (Monotonicity=="No"){
    
    I_Delta_T_Delta_S <- 
      (mat1*log2(mat1/(sum_S_min1*sum_T_min1)))+(mat2*log2(mat2/(sum_S_min1*sum_T_0)))+(mat3*log2(mat3/(sum_S_min1*sum_T_1)))+
      (mat4*log2(mat4/(sum_S_0*sum_T_min1)))+(mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
      (mat7*log2(mat7/(sum_S_1*sum_T_min1)))+(mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
    
    H_Delta_T <-  
      -(((mat1+mat4+mat7)*log2(mat1+mat4+mat7))+ 
          ((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+
          ((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
    H_Delta_S <-   
      -(((mat1+mat2+mat3)*log2(mat1+mat2+mat3))+ 
          ((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+
          ((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
    R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
  }
  
  if (Monotonicity=="True.Endp"){
    
    I_Delta_T_Delta_S <- 
      0+(mat2*log2(mat2/(sum_S_min1*sum_T_0)))+(mat3*log2(mat3/(sum_S_min1*sum_T_1)))+
      0+(mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
      0+(mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
    
    H_Delta_T <-  
      -(0+((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
    H_Delta_S <-   
      -(((mat1+mat2+mat3)*log2(mat1+mat2+mat3))+ 
          ((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+
          ((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
    R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
  }
  
  if (Monotonicity=="Surr.Endp"){
    
    I_Delta_T_Delta_S <- 
      0+(mat4*log2(mat4/(sum_S_0*sum_T_min1)))+(mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
      (mat7*log2(mat7/(sum_S_1*sum_T_min1)))+(mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
    
    H_Delta_T <-  
      -(((mat1+mat4+mat7)*log2(mat1+mat4+mat7))+ 
          ((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+
          ((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
    H_Delta_S <-   
      -(0+((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+
          ((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
    R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
  }
  
  if (Monotonicity=="Surr.True.Endp"){
    
    I_Delta_T_Delta_S <- 
      (mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
      (mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
    
    H_Delta_T <-  
      -(((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
    H_Delta_S <-   
      -(((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
    R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
  }
  
  # Association potential outcomes true (theta_T) and surrogate (theta_S)
  pi_T_00 <- Pi_s[1] + Pi_s[3] + Pi_s[4] + Pi_s[14] 
  pi_T_01 <- Pi_s[2] + Pi_s[5] + Pi_s[13] + Pi_s[15]
  pi_T_10 <- Pi_s[6] + Pi_s[7] + Pi_s[8] + Pi_s[11]
  pi_T_11 <- Pi_s[9] + Pi_s[10] + Pi_s[12] + Pi_s[16]
  
  pi_S_00 <- Pi_s[1] + Pi_s[2] + Pi_s[6] + Pi_s[16]  
  pi_S_01 <- Pi_s[4] + Pi_s[5] + Pi_s[8] + Pi_s[10] 
  pi_S_10 <- Pi_s[3] + Pi_s[7] + Pi_s[9] + Pi_s[13]   
  pi_S_11 <- Pi_s[11] + Pi_s[12] + Pi_s[14] + Pi_s[15]
  
  theta_T <- (pi_T_00 * pi_T_11)/(pi_T_10 * pi_T_01)
  theta_S <- (pi_S_00 * pi_S_11)/(pi_S_10 * pi_S_01)
  
  fit <- list(Data.STSBinBin.Obs=mat_obs, Data.STSBinBin.Counter=mat, Vector_Pi=Pi_s_all, Pi_Marginals=Pi_Marginals, 
              True.R2_H=R2_H, True.Theta_T=theta_T, True.Theta_S=theta_S)
  
  class(fit) <- "Sim.Data.STSBinBin"
  fit
  
}
