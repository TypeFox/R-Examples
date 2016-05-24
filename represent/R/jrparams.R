jrparams <-
function(BLOCK.1, BLOCK.2, ncomp=min(c(dim(BLOCK.1), dim(BLOCK.2))), Cscrit=0.6, Rscrit=0.6) {  #Default 'ncomp' setting based on assumption that # of PCs retained by prcomp() function by default depends on 'tol' rather than on actual rank
  #Arguments 'Cscrit' and 'Rscrit' represent the the values of 'C*' and 'R*' below which 'C*' and 'R*' are deemed to be 'unacceptable' by the user

  #Number of objects in both data blocks to be compared
  n1 <- nrow(BLOCK.1)
  n2 <- nrow(BLOCK.2)
  
  
  ## I.) Comparison of the directions of the two data sets ('§2.3')
  #Carry out two separate PCA, one on the training set and one on the test
  #set, both centered independently
  
  model.block1 <- prcomp(BLOCK.1, center=T, scale.=F)  #'ncomp1' and 'ncomp2' might be incorporated somewhere around here as well, but it seems impossible to do this in the call to prcomp() rightaway
  model.block2 <- prcomp(BLOCK.2, center=T, scale.=F)
  
  
  #Retrieve loading vectors
  loads.block1 <- model.block1$rotation
  loads.block2 <- model.block2$rotation
  
  #Retrieve eigenvalues
  eigs.block1 <- as.matrix((model.block1$sdev)^2)
  eigs.block2 <- as.matrix((model.block2$sdev)^2)

  
  #Only keep PCs deemed significant
  loads.block1 <- as.matrix(loads.block1[, 1:ncomp])
  loads.block2 <- as.matrix(loads.block2[, 1:ncomp])
  #
  eigs.block1 <- as.matrix(eigs.block1[1:ncomp, ])
  eigs.block2 <- as.matrix(eigs.block2[1:ncomp, ])
  
  
  #Give all eigenvectors same sign (see Jouan-Rimbaud Chemom Intell Lab Syst
  #1998 p. 131 right hand column)
  refvect <- matrix(nrow=1, ncol=nrow(loads.block1), data=1) / norm(matrix(nrow=1, ncol=nrow(loads.block1), data=1)) #The 'ref'erence vector [1 1 1 ... 1]/norm([1 1 1 ... 1]) as on page 131 of Jouan-Rimbaud Chemom Intell Lab Syst 1998, right hand column
  #We assume that both datasets have the same number of variables, and
  #therefore also the same number of loadings. Therefore we can use 'refvect'
  #for both datasets
  #
  #First for 'BLOCK.1'
  refloads.block1 <- matrix(nrow=nrow(loads.block1), ncol=ncol(loads.block1), data=NA) #Initialization; array with loadings vectors (eigenvectors) all having the same sign
  for (pc in 1:ncol(loads.block1)) {
      refproduct <- refvect%*%loads.block1[,pc] #Scalar product
      if (refproduct<0) {
          refloads.block1[,pc] <- -loads.block1[,pc] #Use ' - eigenvector'
      } else {
          refloads.block1[,pc] <- loads.block1[,pc]
      }
  }
  #
  #Same for 'BLOCK.2'
  refloads.block2 <- matrix(nrow=nrow(loads.block2), ncol=ncol(loads.block2), data=NA) #Initialization; array with loadings vectors (eigenvectors) all having the same sign
  for (pc in 1:ncol(loads.block2)) {
      refproduct <- refvect%*%loads.block2[,pc] #Scalar product
      if (refproduct<0) {
          refloads.block2[,pc] <- -loads.block2[,pc] #Use ' - eigenvector'
      } else {
          refloads.block2[,pc] <- loads.block2[,pc]
      }
  }
  
  
  JR.RESULTS <- matrix(nrow=6, ncol=ncol(loads.block1), data=NA) #Matrix to be filled with the results of analyses as proposed by Jouan-Rimbaud in 1998 article: rows [1:2]='P' and 'P*'; rows [3:4]='C' and 'C*'; rows [5:6]='R' and 'R*'. Instead of size(loads.block1,2) we could have taken size(loads.block2,2) as well here for the number of columns (provided that the same number of PCs is kept for both data sets)
  rownames(JR.RESULTS) <- c("P", "P*", "C", "C*", "R", "R*")
  colnamevect <- "1PC"
  if (ncomp > 1) {
    colnamevect <- c(colnamevect, paste(2:ncomp, "PCs", sep=""))
  }
  colnames(JR.RESULTS) <- colnamevect
  #
  for (p in 1:ncol(loads.block1)) {
  
      d1 <- matrix(nrow=nrow(refloads.block1),ncol=1, data=0) #Average direction vector for training set (see Jouan-Rimbaud Chemom Intell Lab Syst 1998 p. 131 eq. 2a)
      for (pc in 1:p) {
          d1 <- d1 + eigs.block1[pc,1]*refloads.block1[,pc]
      }
      #
      d2 <- matrix(nrow=nrow(refloads.block2),ncol=1, data=0) #Average direction vector for test set (see Jouan-Rimbaud Chemom Intell Lab Syst 1998 p. 131 eq. 2b)
      for (pc in 1:p) {
          d2 <- d2 + eigs.block2[pc,1]*refloads.block2[,pc]
      }
  
      P <- abs((t(d1)%*%d2) / sqrt(t(d1)%*%d1%*%t(d2)%*%d2)) #Jouan-Rimbaud Chemom Intell Lab Syst 1998 p. 131 eq. 3
      JR.RESULTS[1,p] <- P
      #
      Ps <- acos(P)*180/pi #Jouan-Rimbaud Chemom Intell Lab Syst 1998 p. 131 eq. 4
      Ps.scaled <- 1-Ps/45 #== 'P*' from page 131 of Jouan-Rimbaud Chemom Intell Lab Syst 1998 onwards
      JR.RESULTS[2,p] <- Ps.scaled
  }
  
  
  ## II.) Comparison of the variance-covariance matrices ('§2.4')
  #Alternative PCA needed, project test set in model set up on basis of
  #training set
  
  #1.) Perform PCA on training set 'BLOCK.1' only, apply meancentering as
  #preprocessing step. Save PCA model as 'model.block1' to ML workspace.
  #2.) Perform PCA on test set 'BLOCK.2' only, using 'model.block1' as PCA
  #model. Save Prediction (!) as 'pred.block2' to <R> workspace.
  
  #model.block1 <- prcomp(BLOCK.1, center=T, scale.=F) #Could be superfluous; this PCA has already been computed for 'P' parameter value computation above
  #
  pred.block2 <- predict(object=model.block1, newdata=BLOCK.2) #Apply PCA based on 'BLOCK.1' to 'BLOCK.2' data. Note that predict() by default does meancentering on 'newdata'

  #loads.block1 <- model.block1$rotation #Could be superfluous; PCA loadings for 'BLOCK.1' have already been computed for 'P' and 'P*' parameter value computation above
  
  
  #Combine scores after PCA using same model on training and test data
  #To compute scores based on original data matrix and PCA loadings matrix: see Smilde et al Metabolomics 2009
  scores.block1 <- as.matrix(scale(BLOCK.1, center=T, scale=F))%*%loads.block1  #Could be superfluous as well
  scores.block2 <- pred.block2[, 1:ncomp] #Not superfluous; new scores, for 'significant' PCs only
  #
  scores <- rbind(scores.block1, scores.block2)
  
  #Weight scores according to % explained variance
  percvar <- as.matrix(eigs.block1/sum(eigs.block1))
  
  wscores <- matrix(nrow=nrow(scores), ncol=ncol(scores), data=NA) #Array with weighted scores
  for (pc in 1:ncol(scores)) { #For each PC considered in current Box's test
      wscores[,pc] <- scores[,pc]/percvar[pc,1]
  }
  
 
  groups <- as.matrix(c(rep(1, n1), rep(2, n2))) #Define groupings, necessary for MBoxtest function
  #
  for (p in 1:ncol(loads.block1)) {
  
      GSCORES <- cbind(groups, wscores[,1:p])  #Weighted scores
      #
      MB.OUT <- MBoxtest(GSCORES, nmanvars=ncol(BLOCK.1))
      MB <- MB.OUT$MB #Box's M-statistic
      Sp <- MB.OUT$Sp #Pooled covariance matrix
  
  
      #Compute parameter 'C' (REF Jouan-Rimbaud et al Chemom Intell Lab Syst
      #1998, p. 132)
      
      cparam <- exp(-MB/(n1+n2-2)) #'cparam' should be in the range [0,1]
      #(REF Jouan-Rimbaud et al Chemom Intell Lab Syst 1998, p. 132). Formally it
      #is safer to parametrize the number of groups in 'deno' as done here using
      #the parameter 'g', than to set the number of groups to 2 as done in the
      #article of Jouan-Rimbaud et al.
      if (cparam<0) {
          cparam <- 0
      }
      JR.RESULTS[3,p] <- cparam
  
      #Compute 'C*' parameter value (REF Jouan-Rimbaud et al Chemom Intell Lab
      #Syst 1998, p. 132):
      #1.) Compute critical value Mcrit = X2(p*(p+1)*(q-1)/2) [REF Michie et al
      #(1994) p. 113 ("This statistic has an asymptotic ... distribution")] for
      #the M statistic
      q <- 2 #Number of groups; let's assume that this always equals 2 for our applications
      
      Mcrit <- qchisq(p=0.95, df=(p*(p+1)*(q-1)/2), lower.tail=T)      
      #
      #2.) Compute 'C*' parameter value (REF Jouan-Rimbaud et al Chemom Intell
      #Lab Syst 1998, p. 132):
      #Cs <- 1-(.4/Mcrit)*MB
      Cs <- 1-((1-Cscrit)/Mcrit)*MB
      if (Cs<0) {
          Cs <- 0
      }
      JR.RESULTS[4,p] <- Cs
  
  
  ## III.) Comparison of the data set centroids ('§2.5')
      #Mahalanobis distance computation
      #[D2,CI,prob,power] <- Mahal(wscores(:,1:p),groups) #Using weighted scores; see De Maesschalck et al, Chemom Intell Lab Syst 2000
      D2 <- JRsMahaldist(GSCORES)$Ds  #Using weighted scores; see De Maesschalck et al, Chemom Intell Lab Syst 2000

      #F parameter computation (see Jouan-Rimbaud Chemom Intell Lab Syst 1998 p.
      #132)
      #
      numer <- (n1*n2*(n1+n2-p-1))
      denom <- (p*(n1+n2)*(n1+n2-2))
      # F=numer/denom*D2; #Their eq. 11; value of 'F' is not used furthermore
  
      #Calculate critical value for the Mahalanobis distance for significance
      #level of 0.05
      #x <- finv(0.95,p,(n1+n2-p-1)) #Find a value that should exceed 95% of the
      #samples from an F distribution with 'p' degrees of freedom in the
      #numerator and (n1+n2-p-1) degrees of freedom in the denominator (see help
      #on function 'finv'). Of course, the numerator and the denominator are
      #specific for given numbers of objects and variables (n1, n2, p), causing
      #D2crit to be dependent on these numbers of objects and variables as well.
      #
      x <- qf(p=0.95, df1=p, df2=(n1+n2-p-1))
      #
      D2crit <- x/(numer/denom) #Their eq. 11 as well
  
  
      #R parameter computation (see Jouan-Rimbaud Chemom Intell Lab Syst 1998 p.
      #132)
      if (D2 <= D2crit) {
          R <- 1-(D2/D2crit)
      } else {
          R <- 0
      }
      #
      JR.RESULTS[5,p] <- R
      #
      #Rs <- (-.4/D2crit)*D2+1 #'R*' parameter in Jouan-Rimbaud 1998 p. 133
      Rs <- (-(1-Rscrit)/D2crit)*D2+1 #'R*' parameter in Jouan-Rimbaud 1998 p. 133
      if (Rs<0) {
          Rs <- 0
      }
      #
      JR.RESULTS[6,p] <- Rs
  }


  #Prepare output
  #OUT <- list()
  #OUT$JR.RESULTS <- JR.RESULTS
  OUT <- JR.RESULTS
  #
  return(OUT)

}
