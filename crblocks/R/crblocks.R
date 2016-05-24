print.crblocks_output <- function(x,...) {
#######################################################################
# Function to format our output nicely                                #
#######################################################################
 Nsigfigs=4
 cat("\n")
 if (!is.null(x$Sstatistic)){ ### format output from catrandstat()
  cat(paste(" Statistic   dof   data value   chi^2 p-value\n"))
  cat(paste("  S          ",(x$Nproducts-1)*(x$Ncategories-1),"   ",signif(x$Sstatistic,Nsigfigs),"      ",signif(x$Schi2pvalue,Nsigfigs),"\n"))
  cat(paste("  M          ",x$Ncategories-1,"   ",signif(x$Mstatistic,Nsigfigs),"      ",signif(x$Mchi2pvalue,Nsigfigs),"\n"))
  cat(paste("  L^2         1    ",signif(x$L2statistic,Nsigfigs),"      ",signif(x$L2chi2pvalue,Nsigfigs),"\n"))
 } else if (!is.null(x$Smontecarlo)){ ### format output from catrandpvalue()
  cat(paste(" Statistic   dof   data value   chi^2 p-value   Simulated p-value\n"))
  cat(paste("  S          ",(x$Nproducts-1)*(x$Ncategories-1),"   ",signif(x$Sdata,Nsigfigs),"      ",signif(x$Schi2pvalue,Nsigfigs),"       ",signif(x$Spvalue,Nsigfigs),"\n"))
  cat(paste("  M          ",x$Ncategories-1,"   ",signif(x$Mdata,Nsigfigs),"      ",signif(x$Mchi2pvalue,Nsigfigs),"      ",signif(x$Mpvalue,Nsigfigs),"\n"))
  cat(paste("  L^2         1    ",signif(x$L2data,Nsigfigs),"      ",signif(x$L2chi2pvalue,Nsigfigs),"      ",signif(x$L2pvalue,Nsigfigs),"\n"))
 } else if (!is.null(x$Spermute)){ ### format output from catrandpvaluepermute()
  cat(paste(" Statistic   dof   data value   chi^2 p-value   Simulated p-value\n"))
  cat(paste("  S          ",(x$Nproducts-1)*(x$Ncategories-1),"   ",signif(x$Sdata,Nsigfigs),"      ",signif(x$Schi2pvalue,Nsigfigs),"       ",signif(x$Spvalue,Nsigfigs),"\n"))
  cat(paste("  M          ",x$Ncategories-1,"   ",signif(x$Mdata,Nsigfigs),"      ",signif(x$Mchi2pvalue,Nsigfigs),"      ",signif(x$Mpvalue,Nsigfigs),"\n"))
  cat(paste("  L^2         1    ",signif(x$L2data,Nsigfigs),"      ",signif(x$L2chi2pvalue,Nsigfigs),"      ",signif(x$L2pvalue,Nsigfigs),"\n"))
 }
 cat("\n")
}




catrandstat <- function(rawdata){
#######################################################################
# Function to compute the statistic                                   #
#######################################################################

### Extract the meta variables from the data:
 ### Njudges is the number of rows:
  Njudges <- nrow(rawdata)
 ### Nproducts is the number of columns:
  Nproducts <- ncol(rawdata)

### Extract the category labels that are present:
 if (is.matrix(rawdata)){
  ### for the Monte Carlo data (which is a matrix):
   categories <- sort(unique(sort(rawdata)))
 } else {
  ### for the data read from file (which is a data frame):
   categories <- sort(unique(stack(rawdata)[,1]))
 }
 ### Count how many we have:
 Ncategories <- length(categories)

 ### Construct the product-X-categories counts:
  catCounts <- matrix(NaN,Nproducts,Ncategories) # initialise
  for (i in 1:Nproducts){
   for (j in 1:Ncategories){
    catCounts[i,j] <- length(which(rawdata[,i]==categories[j]))
   }
  }
 catColSums <- apply(catCounts,2,sum)

 ### Construct the d vectors (as columns of a matrix:):
  d <- matrix(NaN,Ncategories-1,Nproducts) # initialise
  for (i in 1:Nproducts){
   d[,i]=(Nproducts*catCounts[i,])[1:(Ncategories-1)]-catColSums[1:(Ncategories-1)]
 }

 ### Construct the judges-X-categories counts (judgeCatCounts):
  U <- matrix(NaN,Njudges,Ncategories) # initialise
  for (i in 1:Njudges){
   for (j in 1:Ncategories){
    U[i,j] <- length(which(rawdata[i,]==categories[j]))
   }
  }
  judgeCatCountsFull=U
  U=U[,1:Ncategories-1] # for the statistic we only want the first Ncategories-1 columns

 ### Construct the covariance matrix, V:
  V <- matrix(NaN,Ncategories-1,Ncategories-1) # initialise
  for (i in 1:Ncategories-1){
   for (j in 1:Ncategories-1){
    if (i==j){
     V[i,j] <- Nproducts*sum(U[,i])-crossprod(U[,i])
    } else {
     V[i,j] <- -crossprod(U[,i],U[,j])
    }
   }
  }
 ### find the inverse of V: (first test if it is singular (ie. det(V)==0) (or almost))
  if (det(V)>1e-2){
   ### not singular? okay, compute the inverse:
    Vinverse <- qr.solve(V)

   ### Compute the S statistic:
    S = 0  # initialise
    for (i in 1:Nproducts){
     S = S + (Nproducts-1)/Nproducts*crossprod(t(crossprod(d[,i],Vinverse)),d[,i])
    }
    Schi2pvalue = pchisq(S,(Nproducts-1)*(Ncategories-1),lower.tail=FALSE);
  } else {
   ### flag the singular matrix by returning NaN
    S = NaN
    Vinverse <- NaN
  }

  ### Compute the M and L^2 statistics:
   w=sum(apply(rawdata,1,var))
   Abar=apply(rawdata,2,mean)
   ### obtain the linear contrasts:
   if (Nproducts%%2){ ### odd number of products:
    # eg. [-1 0 -1]
    lambda = -floor(Nproducts/2):floor(Nproducts/2)
   } else {
    # eg. [-3 -1 1 3]
    lambda = 2*(1:Nproducts)-Nproducts-1
   }
   M = (Njudges^2/w)*sum((Abar-mean(Abar))^2)
   L2 = ((Njudges/sqrt(w))*sum(lambda*Abar)/sqrt(sum(lambda^2)))^2
   Mchi2pvalue = pchisq(M,(Nproducts-1),lower.tail=FALSE);
   L2chi2pvalue = pchisq(L2,1,lower.tail=FALSE);


#######################################################################
# Return some variables of interest                                   #
#######################################################################
  returnoutput = list(
  Njudges=Njudges,
  Nproducts=Nproducts,
  rawdata=rawdata,
  categories=categories,
  Ncategories=Ncategories,
  catCounts=catCounts,
  judgeCatCounts=judgeCatCountsFull,
  Sstatistic=S,
  Mstatistic=M,
  L2statistic=L2,
  Schi2pvalue = Schi2pvalue,
  Mchi2pvalue = Mchi2pvalue,
  L2chi2pvalue = L2chi2pvalue
 )

 ### specify our print method to display the output nicely:
  class(returnoutput) <- c("crblocks_output",class(returnoutput))
 ### print (using the method determined by R):
  returnoutput

} ### end of function










catrandpvalue <- function(datafilename,Nrepeats){
#######################################################################
# Function to compute the p-value for the data                        #
#######################################################################

 ### Check if the data file exists:
  if(!file.exists(datafilename))
   stop('File \'',datafilename,'\' not found')

 ### Open the file read-only:
  inputfile  <-  file(datafilename,"r")
 ### Read the data:
 ### (note that comments (starting with #) are allowed in the data file)
  rawdata  <-  read.table(inputfile, header=FALSE)
 ### Close the file:
  close(inputfile)

 #######################################################################
 # Compute the statistic for the data                                  #
 #######################################################################
  dataoutput <- catrandstat(rawdata)
  Sdata <- dataoutput$Sstatistic
  Mdata <- dataoutput$Mstatistic
  L2data <- dataoutput$L2statistic
  Schi2pvalue <- dataoutput$Schi2pvalue
  Mchi2pvalue <- dataoutput$Mchi2pvalue
  L2chi2pvalue <- dataoutput$L2chi2pvalue

 #######################################################################
 # Run the Monte Carlo simulations to compute a p-value                #
 #######################################################################
 ### Firstly, we generate the Monte Carlo data:
  Ngenerated=0
  montecarlodata=array(NaN,c(dataoutput$Njudges,dataoutput$Nproducts,Nrepeats))
  for (i in 1:dataoutput$Njudges){
   indx = 1:Nrepeats # initialise
   while (length(indx)){
    Ngenerated=Ngenerated+length(indx) # count how many we actully generate
    montecarlodata[i,,indx]=dataoutput$categories[apply(rmultinom(length(indx)*dataoutput$Nproducts,1,dataoutput$judgeCatCounts[i,])==1,2,which)]
    ### we have the data, now loop over them and see which are tied, so that we can replace them (those in 'indx'):
    indx=indx[which(apply(apply(apply(montecarlodata[i,,indx],2,sort),2,diff),2,sum)==0)]
    if (length(indx)==1){ # if there is only one element in indx, add another since R doesn't handle the multidimensional matrix indexing well with only one element
     if (indx[1]==1){
      indx[2]=2
     } else {
      indx[2]=1
     }
    }
   }
  }

 ### Now we can compute the statistic for each Monte Carlo data set:
  Smontecarlo=matrix(NaN,Nrepeats,1) # initialise
  Mmontecarlo=matrix(NaN,Nrepeats,1) # initialise
  L2montecarlo=matrix(NaN,Nrepeats,1) # initialise
  for (n in 1:Nrepeats){
   montecarlooutput = catrandstat(montecarlodata[,,n])
   Smontecarlo[n] = montecarlooutput$Sstatistic
   Mmontecarlo[n] = montecarlooutput$Mstatistic
   L2montecarlo[n] = montecarlooutput$L2statistic
  }

Spvalue=length(Smontecarlo[Smontecarlo>c(Sdata)])/length(Smontecarlo)
Mpvalue=length(Mmontecarlo[Mmontecarlo>c(Mdata)])/length(Mmontecarlo)
L2pvalue=length(L2montecarlo[L2montecarlo>c(L2data)])/length(L2montecarlo)

#######################################################################
# Return some variables of interest                                   #
#######################################################################
 returnoutput = list(
  rawdata=rawdata,
  Nproducts=dataoutput$Nproducts,
  Ncategories=dataoutput$Ncategories,
  Njudges=dataoutput$Njudges,
  Ngenerated=Ngenerated,
  Sdata=Sdata,
  Mdata=Mdata,
  L2data=L2data,
  Smontecarlo=Smontecarlo,
  Mmontecarlo=Mmontecarlo,
  L2montecarlo=L2montecarlo,
  Spvalue=Spvalue,
  Mpvalue=Mpvalue,
  L2pvalue=L2pvalue,
  Schi2pvalue=Schi2pvalue,
  Mchi2pvalue=Mchi2pvalue,
  L2chi2pvalue=L2chi2pvalue
 )

 ### specify our print method to display the output nicely:
  class(returnoutput) <- c("crblocks_output",class(returnoutput))
 ### print (using the method determined by R):
  returnoutput

} ### end of function







catrandpvaluepermute <- function(datafilename,Nrepeats){
#######################################################################
# Function to compute the p-value for the data using permutation test #
#######################################################################

 ### Check if the data file exists:
  if(!file.exists(datafilename))
   stop('File \'',datafilename,'\' not found')

 ### Open the file read-only:
  inputfile  <-  file(datafilename,"r")
 ### Read the data:
 ### (note that comments (starting with #) are allowed in the data file)
  rawdata  <-  read.table(inputfile, header=FALSE)
 ### Close the file:
  close(inputfile)

 #######################################################################
 # Compute the statistic for the data                                  #
 #######################################################################
  dataoutput <- catrandstat(rawdata)
  Sdata = dataoutput$Sstatistic
  Mdata = dataoutput$Mstatistic
  L2data = dataoutput$L2statistic
  Schi2pvalue <- dataoutput$Schi2pvalue
  Mchi2pvalue <- dataoutput$Mchi2pvalue
  L2chi2pvalue <- dataoutput$L2chi2pvalue

 ### Compute the statistic for each permuted data set:
  Spermute=matrix(NaN,Nrepeats,1) # initialise
  Mpermute=matrix(NaN,Nrepeats,1) # initialise
  L2permute=matrix(NaN,Nrepeats,1) # initialise
  for (n in 1:Nrepeats){
   permdata=t(apply(rawdata,1,sample)) ### this permutes each row separately (but transposes too)
   permuteoutput = catrandstat(permdata)
   Spermute[n] = permuteoutput$Sstatistic
   Mpermute[n] = permuteoutput$Mstatistic
   L2permute[n] = permuteoutput$L2statistic
  }

 Spvalue=length(Spermute[Spermute>c(Sdata)])/length(Spermute)
 Mpvalue=length(Mpermute[Mpermute>c(Mdata)])/length(Mpermute)
 L2pvalue=length(L2permute[L2permute>c(L2data)])/length(L2permute)

#######################################################################
# Return some variables of interest                                   #
#######################################################################
 returnoutput = list(
  rawdata=rawdata,
  Nproducts=dataoutput$Nproducts,
  Ncategories=dataoutput$Ncategories,
  Njudges=dataoutput$Njudges,
  Sdata=Sdata,
  Mdata=Mdata,
  L2data=L2data,
  Spermute=Spermute,
  Mpermute=Mpermute,
  L2permute=L2permute,
  Spvalue=Spvalue,
  Mpvalue=Mpvalue,
  L2pvalue=L2pvalue,
  Schi2pvalue=Schi2pvalue,
  Mchi2pvalue=Mchi2pvalue,
  L2chi2pvalue=L2chi2pvalue
 )

 ### specify our print method to display the output nicely:
  class(returnoutput) <- c("crblocks_output",class(returnoutput))
 ### print (using the method determined by R):
  returnoutput

} ### end of function
