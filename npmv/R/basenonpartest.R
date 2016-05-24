basenonpartest <-
function(frame,group,vars,tests=c(1,1,1,1)){

   if(!is.factor(frame[,group]))
   {
      frame[,group] <- factor(frame[,group])
   }

   # Orders data by group
   o <- order(frame[,group])
   frame <- frame[o,]

   # Compute sample sizes per group
   p <- length(vars)
   a <- length(levels(frame[,group]))
   N <- length(frame[,1])
   ssize <- array(NA,a)
   lims <- matrix(NA,2,a)
   for(i in 1:a){
      ssize[i] <- length(frame[frame[,group]==levels(frame[,group])[i],1])
      lims[1,i] <- min(which(frame[,group]==levels(frame[,group])[i]))
      lims[2,i] <- max(which(frame[,group]==levels(frame[,group])[i]))
   }

   # Sets up R matrix
   Rmat <- matrix(NA,N,p)

   for(j in 1:p){
      Rmat[,j] <- rank(frame[,vars[j]],ties.method="average")
   }

   # Manipulating R

   Rbars <- matrix(NA,a,p)
   for(i in 1:a){
      for(j in 1:p){
         Rbars[i,j] <- mean(Rmat[(lims[1,i]:lims[2,i]),j])
      }
   }

   Rtilda <- (1/a)*colSums(Rbars)
   Rbarovr <- (1/N)*colSums(Rmat)

   H1 <- matrix(0,p,p)
   H2 <- matrix(0,p,p)
   G1 <- matrix(0,p,p)
   G2 <- matrix(0,p,p)
   G3 <- matrix(0,p,p)
   for(i in 1:a){
      H1 <- H1 + ssize[i]*(Rbars[i,] - Rbarovr)%*%t(Rbars[i,] -Rbarovr)
      H2 <- H2 + (Rbars[i,] - Rtilda)%*%t(Rbars[i,] - Rtilda)
      for(j in 1:ssize[i]){
         G1 <- G1 + (((Rmat[(lims[1,i]+j-1),(1:p)])-Rbars[i,])%*%t((Rmat[(lims[1,i]+j-1),(1:p)])-Rbars[i,]))
         G2 <- G2 + (1-(ssize[i]/N))*(1/(ssize[i]-1))*(((Rmat[(lims[1,i]+j-1),(1:p)])-Rbars[i,])%*%t((Rmat[(lims[1,i]+j-1),(1:p)])-Rbars[i,]))
         G3 <- G3 + (1/(ssize[i]*(ssize[i]-1)))*(((Rmat[(lims[1,i]+j-1),(1:p)])-Rbars[i,])%*%t((Rmat[(lims[1,i]+j-1),(1:p)])-Rbars[i,]))
      }
   }
   H1 <- (1/(a-1))*H1
   H2 <- (1/(a-1))*H2
   G1 <- (1/(N-a))*G1
   G2 <- (1/(a-1))*G2
   G3 <- (1/a)*G3

   # ANOVA TYPE TEST
   fhat <- ((a-1)*sum(diag(G3))^2)/sum(diag(G3%*%G3))
   fhat0 <- ((a^2)/((a-1)*sum(1/(ssize-1))))*fhat
   Fanova <- sum(diag(H2))/sum(diag(G3))
   pvalanova <- 1 - pf(Fanova,fhat,fhat0)
   df1anova<-fhat
   df2anova<-fhat0

   #McKeon approximation for the Lawley Hotelling Test
   if(tests[2]==1){
      U <- sum(diag(((a-1)*H1) %*% solve((N-a)*G1)))
      K <- p*(a-1)
      B <- ((N-p-2)*(N-a-1))/((N-a-p)*(N-a-p-3))
      D <- 4 + (p*(a-1) +2)/(B-1)
      g <- (p*(a-1)*(D-2))/((N-a-p-1)*D)
      FLH <- U/g
      pvalLH <- 1- pf(FLH,K,D)
      df1LH<-K
      df2LH<-D
   }else{
      FLH <- NA
      pvalLH <- NA
      df1LH<-NA
      df2LH<-NA
   }

   #Muller approximation for the Bartlett-Nanda-Pillai Test
   if(tests[3]==1){
      V     <- sum(diag(((a-1)*H1)%*%solve(((a-1)*H1)+(N-a)*G1)))
      gamma <- min(c((a-1),p))
      nu1   <- ((p*(a-1))/(gamma*(N-1)))*  (((gamma*(N-a+gamma-p)*(N+2)*(N-1))/((N-a)*(N-p)))-2)
      nu2   <- ((N-a+gamma-p)/N) * (( (gamma*(N-a+gamma-p)*(N+2)*(N-1))/((N-a)*(N-p)))-2)
      FBNP  <- ((V/gamma)/nu1)/((1-(V/gamma))/nu2)
      pvalBNP <- 1-pf(FBNP,nu1,nu2)
      df1BNP<-nu1
      df2BNP<-nu2
   }else{
      FBNP  <- NA
      pvalBNP <- NA
      df1BNP<-NA
      df2BNP<-NA
   }

   #Wilk's Lambda Test
   if(tests[4]==1){
      lambda=det((N-a)*G1 )/det( (N-a)*G1+(a-1)*H1 )
      I=diag(p)
      r_WL=(N-a)-(p-(a-1)+1)/2
      u_WL=(p*(a-1)-2)/4
      t_d_WL=p*p+(a-1)*(a-1)-5
      if (t_d_WL > 0) {
         t_WL=sqrt( (p^2*(a-1)^2-4) / t_d_WL )
      }else {
         t_WL = 1
      }
      df1_WLF=p*(a-1)
      df2_WLF=r_WL*t_WL-2*u_WL
      WLF=(1-lambda^(1/t_WL))/lambda^(1/t_WL)*df2_WLF/df1_WLF
      pvalWL=1-pf(WLF,df1_WLF,df2_WLF)
      df1WL<-df1_WLF
      df2WL<-df2_WLF
   }else{
      WLF=NA
      pvalWL=NA
      df1WL<-NA
      df2WL<-NA
   }

   out <- list('Fanova'=Fanova,'FWL'=WLF,'pvalWL'=pvalWL,'fhat'=fhat,'fhat0'=fhat0,'pvalanova'=pvalanova,'FLH'=FLH,'pvalLH'=pvalLH,'FBNP'=FBNP,'pvalBNP'=pvalBNP,
               'df1anova'=df1anova,'df2anova'=df2anova,'df1LH'=df1LH,'df2LH'=df2LH,'df1BNP'=df1BNP,'df2BNP'=df2BNP,'df1WL'=df1WL,'df2WL'=df2WL)
   return(out)

}
