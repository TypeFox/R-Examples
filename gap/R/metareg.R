# 23-7-2009 MRC-Epid JHZ

metareg <- function(data,N,verbose="Y",prefixb="b",prefixse="se")
{
   M <- dim(data)[1]
   eps <- .Machine$double.eps
   B <- SE <- W <- OK <- WC <- matrix(0,M,N)
   for (j in 1:M)
   {
       for (i in 1:N)
       {
           B[j,i] <- data[paste(prefixb,i,sep="")][j,]
           SE[j,i] <- data[paste(prefixse,i,sep="")][j,]
       }
   }
   K <- BW <- SW <- SWW <- QW <- DL <- P_HETER <- I2 <- rep(0,M)
   BETA_F <- SE_F <- Z_F <- P_F <- BETA_R <- SE_R <- Z_R <- P_R <- rep(0,M)
   for (j in 1:M)
   {
       for (i in 1:N)
       {
           OK[j,i] <- 0
           S1 <- as.numeric(B[j,i])
           S2 <- as.numeric(SE[j,i])
           if (!is.na(S1+S2))
           {
              OK[j,i] <- 1
              T <- max(S2,eps)
              K[j] <- K[j] + 1
              W[j,i] <- 1/T/T
              BW[j] <- BW[j] + S1*W[j,i]
              SW[j] <- SW[j] + W[j,i]
              SWW[j] <- SWW[j]+W[j,i]^2
           }
       }
       if (K[j]>1)
       {
          BETA_F[j] <- BW[j]/SW[j]
          SE_F[j] <- sqrt(1/SW[j])
          Z_F[j] <- BETA_F[j]/SE_F[j]
          P_F[j] <- 2*pnorm(-abs(Z_F[j]))
      }
      for (i in 1:N) if(OK[j,i]==1) QW[j] <- QW[j]+W[j,i]*(B[j,i]-BETA_F[j])^2
      if (K[j]>1)
      {
         P_HETER[j] <- pchisq(QW[j],K[j]-1,lower.tail=FALSE)
         DL[j] <- max(0,(QW[j]-(K[j]-1))/(SW[j]-SWW[j]/SW[j]))
         bwr <- swc <- 0
         for (i in 1:N)
         {
             if (OK[j,i]==1)
             {
                WC[j,i] <- 1/(1/W[j,i]+DL[j])
                bwr <- bwr + B[j,i]*WC[j,i]
                swc <- swc + WC[j,i]
             }
         }
         BETA_R[j] <- bwr/swc
         SE_R[j] <- sqrt(1/swc)
         Z_R[j] <- BETA_R[j]/SE_R[j]
         P_R[j] <- 2*pnorm(-abs(Z_R[j]))
         I2[j] <- (QW[j]-K[j]+1)/QW[j]
         if (I2[j]<0) I2[j] <- 0
      }
      if (toupper(verbose)=="Y")
         cat ("\nMeta-analysis of", N, "studies:\n\n",
           "p_f=", P_F[j], "\n",
           "p_r=", P_R[j], "\n",
           "beta_f=", BETA_F[j], "\n",
           "beta_r=", BETA_R[j], "\n",
           "se_f=",SE_F[j], "\n",
           "se_r=", SE_R[j], "\n",
           "z_f=", Z_F[j], "\n",
           "z_r=", Z_R[j], "\n",
           "p_heter=", P_HETER[j], "\n",
           "i2=", I2[j], "\n",
           "k=", K[j], "\n",
           "eps=",eps, "\n")
   }
   if (toupper(verbose)=="Y")
      cat ("\nwhere\n\n",
           "p_f=P value (fixed effects model)  \n",
           "p_r=P value (random effects model) \n",
           "beta_f=regression coefficient      \n",
           "beta_r=regression coefficient      \n",
           "se_f=standard error                \n",
           "se_r=standard error                \n",
           "z_f=z value                        \n",
           "z_r=z value                        \n",
           "p_heter=heterogeneity test p value \n",
           "i2=I^2                             \n",
           "k=No of tests used                 \n",
           "eps=smallest double-precision number\n")
   invisible(data.frame(beta_f=BETA_F,se_f=SE_F,z_f=Z_F,beta_r=BETA_R,se_r=SE_R,p_heter=P_HETER,i2=I2,k=K,eps=eps))
}
