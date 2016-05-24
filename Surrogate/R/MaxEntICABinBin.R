MaxEntICABinBin <- function(pi1_1_, pi1_0_, pi_1_1,
                         pi_1_0, pi0_1_, pi_0_1, Method="BFGS", 
                         Fitted.ICA=NULL){
  
vector_b <- matrix(data=c(1, pi1_1_, pi1_0_, pi_1_1,
                          pi_1_0, pi0_1_, pi_0_1), ncol=1)

A_r <- matrix(data=c(1, 0, 0, 0, 0, 0, 0,
                     1, 0, 0, 0, 1, 0, 0,
                     1, 0, 0, 0, 0, 1, 0,
                     1, 0, 0, 0, 0, 0, 1,
                     1, 0, 0, 1, 0, 0, 0,
                     1, 0, 1, 0, 0, 0, 0,
                     1, 1, 0, 0, 0, 0, 0), ncol=7)
A_f <- matrix(data=c(1, 0, 1, 0, 0, 0, 1,
                     1, 1, 0, 0, 1, 0, 0,
                     1, 0, 1, 1, 0, 0, 0,
                     1, 1, 0, 0, 0, 0, 1,
                     1, 1, 0, 1, 0, 0, 0,
                     1, 0, 0, 0, 1, 1, 0,
                     1, 0, 0, 0, 0, 1, 1,
                     1, 0, 0, 1, 0, 1, 0,
                     1, 0, 1, 0, 1, 0, 0), ncol=9)
A_mat <- A <- cbind(A_r, A_f)


min_fun <- function(par){
  1/16 *   
    ((exp(t(A[,1])%*%par-1) - t(vector_b) %*% par)  + 
       (exp(t(A[,2])%*%par-1) - t(vector_b) %*% par)  + 
       (exp(t(A[,3])%*%par-1)  - t(vector_b) %*% par)+ 
       (exp(t(A[,4])%*%par-1) - t(vector_b) %*% par)  + 
       (exp(t(A[,5])%*%par-1) - t(vector_b) %*% par)  + 
       (exp(t(A[,6])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,7])%*%par-1) - t(vector_b) %*% par)  + 
       (exp(t(A[,8])%*%par-1) - t(vector_b) %*% par)  + 
       (exp(t(A[,9])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,10])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,11])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,12])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,13])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,14])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,15])%*%par-1) - t(vector_b) %*% par) + 
       (exp(t(A[,16])%*%par)-1) - t(vector_b) %*% par) 
}



if (Method=="BFGS"){

res <- optim(par=rep(-2500, times=7), min_fun, control=c(maxit=2e9, reltol=1e-350), 
             method="BFGS")  

if (res$convergence != 0){  
  cat("\nWarning: the optimization algorithm (FFGS) did not converge. The results may not be trustworty. \n\n")
} 

p_i_all <- NULL
for (i in 1: 16){
  p_i_part <- 
    1/16 * (exp((t(A_mat[,i]) %*% res$par) - 1))
  p_i_all <- cbind(p_i_all, p_i_part)
}

pi <- p_i_all_BFGS <- p_i_all
pi[pi==0] <- 1e-20
mat1 <- pi[7] 
mat2 <- pi[3] + pi[9] 
mat3 <- pi[13] 
mat4 <- pi[6] + pi[11] 
mat5 <- pi[1] + pi[14] + pi[16] + pi[12] 
mat6 <- pi[2] + pi[15]
mat7 <- pi[8]
mat8 <- pi[4] + pi[10] 
mat9 <- pi[5]

Delta_c_mat <-
  matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)

#R2_H
sum_S_min1 <- mat1+mat2+mat3
sum_S_0 <- mat4+mat5+mat6
sum_S_1 <- mat7+mat8+mat9

sum_T_min1 <- mat1+mat4+mat7
sum_T_0 <- mat2+mat5+mat8
sum_T_1 <- mat3+mat6+mat9

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
R2_H_BFGS <- R2_H  #OUTPUT
Pi.Vector.BFGS <- data.frame(p_i_all)
names(Pi.Vector.BFGS) <- c("p*_0000", "p*_0100", "p*_0010", "p*_0001", "p*_0101", "p*_1000", "p*_1010",
                           "p*_1001", "p*_1110",
                           "p*_1101",  "p*_1011", "p*_1111", "p*_0110", "p*_0011", "p*_0111", "p*_1100")
p_i_all_orig <- data.frame(p_i_all)
names(p_i_all_orig) <- c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010",
                           "Pi_1001", "Pi_1110",
                           "Pi_1101",  "Pi_1011", "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100")
p_i_all <- Pi.Vector.BFGS

pi_m1_m1 <- p_i_all_orig$Pi_1010 
pi_1_m1 <- p_i_all_orig$Pi_0110
pi_m1_1 <- p_i_all_orig$Pi_1001
pi_1_1 <- p_i_all_orig$Pi_0101
pi_0_0 <- p_i_all_orig$Pi_0000 + p_i_all_orig$Pi_0011 + p_i_all_orig$Pi_1100 + p_i_all_orig$Pi_1111
pi_0_m1 <- p_i_all_orig$Pi_0010 + p_i_all_orig$Pi_1110 
pi_m1_0 <- p_i_all_orig$Pi_1000 + p_i_all_orig$Pi_1011
pi_0_1 <- p_i_all_orig$Pi_0001 + p_i_all_orig$Pi_1101
pi_1_0 <- p_i_all_orig$Pi_0100 + p_i_all_orig$Pi_0111 

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

    } # einde BFGS



if (Method=="CG"){
  
res <- optim(par=rep(-25, times=7), min_fun, control=c(maxit=10000), method="CG") 

if (res$convergence != 0){  # must be 0; als 1 iteration limit reached
  cat("\nWarning: the optimization algorithm (CG) did not converge. The results may not be trustworty. \n\n")
} 

p_i_all <- NULL
for (i in 1: 16){
  p_i_part <- 
    1/16 * (exp((t(A_mat[,i]) %*% res$par) - 1))
  p_i_all <- cbind(p_i_all, p_i_part)
}

pi <- p_i_all
pi[pi==0] <- 1e-20
mat1 <- pi[7] 
mat2 <- pi[3] + pi[9] 
mat3 <- pi[13] 
mat4 <- pi[6] + pi[11] 
mat5 <- pi[1] + pi[14] + pi[16] + pi[12] 
mat6 <- pi[2] + pi[15]
mat7 <- pi[8]
mat8 <- pi[4] + pi[10] 
mat9 <- pi[5]

Delta_c_mat <-
  matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)

#R2_H
sum_S_min1 <- mat1+mat2+mat3
sum_S_0 <- mat4+mat5+mat6
sum_S_1 <- mat7+mat8+mat9

sum_T_min1 <- mat1+mat4+mat7
sum_T_0 <- mat2+mat5+mat8
sum_T_1 <- mat3+mat6+mat9

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
R2_H_CG <- R2_H  #OUTPUT
Pi.Vector.CG <- data.frame(p_i_all)
names(Pi.Vector.CG) <- names(p_i_all) <- c("p*_0000", "p*_0100", "p*_0010", "p*_0001", "p*_0101", "p*_1000", "p*_1010",
                         "p*_1001", "p*_1110",
                         "p*_1101",  "p*_1011", "p*_1111", "p*_0110", "p*_0011", "p*_0111", "p*_1100")
p_i_all_orig <- data.frame(p_i_all)
names(p_i_all_orig) <- c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010",
                         "Pi_1001", "Pi_1110",
                         "Pi_1101",  "Pi_1011", "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100")
p_i_all <- Pi.Vector.CG

pi_m1_m1 <- p_i_all_orig$Pi_1010 
pi_1_m1 <- p_i_all_orig$Pi_0110
pi_m1_1 <- p_i_all_orig$Pi_1001
pi_1_1 <- p_i_all_orig$Pi_0101
pi_0_0 <- p_i_all_orig$Pi_0000 + p_i_all_orig$Pi_0011 + p_i_all_orig$Pi_1100 + p_i_all_orig$Pi_1111
pi_0_m1 <- p_i_all_orig$Pi_0010 + p_i_all_orig$Pi_1110 
pi_m1_0 <- p_i_all_orig$Pi_1000 + p_i_all_orig$Pi_1011
pi_0_1 <- p_i_all_orig$Pi_0001 + p_i_all_orig$Pi_1101
pi_1_0 <- p_i_all_orig$Pi_0100 + p_i_all_orig$Pi_0111 

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

   }  #einde CG


if (Method=="MD"){
  
  res <- optim(par=rep(-2500, times=7), min_fun, control=c(maxit=2e9, reltol=1e-350), method="BFGS")  
  p_i_all <- NULL
  for (i in 1: 16){
    p_i_part <- 
      1/16 * (exp((t(A_mat[,i]) %*% res$par) - 1))
    p_i_all <- cbind(p_i_all, p_i_part)
  }
  pi <- p_i_all_BFGS <- p_i_all

sq_dev_all <- NULL
for (k in 1: dim(Fitted.ICA$Pi.Vectors)[1]){
  sq_dev <- sum((p_i_all_BFGS - Fitted.ICA$Pi.Vectors[k,][1:16])**2)
  sq_dev_all <- rbind(sq_dev_all, cbind(k, sq_dev))
}
sq_dev_all <- data.frame(sq_dev_all)
sorted <- sq_dev_all[order(sq_dev_all$sq_dev),]
min_sq <- sorted[1, 1]
Pi.Vector.Min.SQ <- Fitted.ICA$Pi.Vectors[min_sq,][1:16]
pi <- Pi.Vector.Min.SQ
mat1 <- pi[7] 
mat2 <- pi[3] + pi[9] 
mat3 <- pi[13] 
mat4 <- pi[6] + pi[11] 
mat5 <- pi[1] + pi[14] + pi[16] + pi[12] 
mat6 <- pi[2] + pi[15]
mat7 <- pi[8]
mat8 <- pi[4] + pi[10] 
mat9 <- pi[5]

Delta_c_mat <-
  matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)

#R2_H
sum_S_min1 <- mat1+mat2+mat3
sum_S_0 <- mat4+mat5+mat6
sum_S_1 <- mat7+mat8+mat9

sum_T_min1 <- mat1+mat4+mat7
sum_T_0 <- mat2+mat5+mat8
sum_T_1 <- mat3+mat6+mat9

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

R2_H_Min.Diff <- as.numeric(R2_H)  #OUTPUT
Pi.Vector.Min.Diff <- p_i_all <- data.frame(Pi.Vector.Min.SQ, row.names = " ")
names(Pi.Vector.Min.Diff) <- names(p_i_all) <- c("p*_0000", "p*_0100", "p*_0010", "p*_0001", "p*_0101", "p*_1000", "p*_1010",
                               "p*_1001", "p*_1110",
                               "p*_1101",  "p*_1011", "p*_1111", "p*_0110", "p*_0011", "p*_0111", "p*_1100")

p_i_all_orig <- data.frame(Pi.Vector.Min.Diff)
names(p_i_all_orig) <- c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010",
                         "Pi_1001", "Pi_1110",
                         "Pi_1101",  "Pi_1011", "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100")
p_i_all <- p_i_all_orig

pi_m1_m1 <- p_i_all_orig$Pi_1010 
pi_1_m1 <- p_i_all_orig$Pi_0110
pi_m1_1 <- p_i_all_orig$Pi_1001
pi_1_1 <- p_i_all_orig$Pi_0101
pi_0_0 <- p_i_all_orig$Pi_0000 + p_i_all_orig$Pi_0011 + p_i_all_orig$Pi_1100 + p_i_all_orig$Pi_1111
pi_0_m1 <- p_i_all_orig$Pi_0010 + p_i_all_orig$Pi_1110 
pi_m1_0 <- p_i_all_orig$Pi_1000 + p_i_all_orig$Pi_1011
pi_0_1 <- p_i_all_orig$Pi_0001 + p_i_all_orig$Pi_1101
pi_1_0 <- p_i_all_orig$Pi_0100 + p_i_all_orig$Pi_0111 

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

     }  #einde MD

#H_max
h_max_part <- h_max_all <- NULL
for (z in 1: 16){
  h_max_part <- as.numeric(p_i_all[z]) * log(as.numeric(p_i_all[z]))
  h_max_all <- cbind(h_max_all, h_max_part)
}
H_max <- -sum(h_max_all)


fit <- 
  list(R2_H=as.numeric(R2_H), Vector_p = p_i_all, H_max=H_max, Call=match.call()) #,
#r_1_1=as.numeric(r_1_1), r_min1_1=as.numeric(r_m1_1), r_0_1=as.numeric(r_0_1),
#r_1_0=as.numeric(r_1_0), r_min1_0=as.numeric(r_m1_0), r_0_0=as.numeric(r_0_0),
#r_1_min1=as.numeric(r_1_m1), r_min1_min1=as.numeric(r_m1_m1), r_0_min1=as.numeric(r_0_m1))

class(fit) <- "MaxEntICA.BinBin"
fit

} 


