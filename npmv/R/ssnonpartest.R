#Function to perform the subset algorithm.  Algorithm is comprised of 3 steps.  Step 1: Global Hypothesis Test.  Step 2: Test all a-1 (or p-1) subsets.  
# Step 3: Test all remaing subsets per closed multiple testing principle


ssnonpartest <- function(formula,data,alpha=.05,test=c(0,0,0,1),factors.and.variables=FALSE){

##########################################################################################
#Sets up data from formula, and other checks before subset algorthim can be performed
##########################################################################################

#Checks to see if formula
  if(!is(formula,"formula")){
    return('Error: Please give a formula')
  }  

#Checks to ensure only one test is requested
  if(sum(test)!=1){
    return('Error:Please specify a single test')
  }
  
#Creates the data frame
  formula=Formula(formula)
  frame=model.frame(formula,data=data)
  
   
#Checks for missing data
   if(sum(is.na(frame))>0)
   {
      return('Error: Missing Data')
   }
  
#Assigns group variable and response variables
  groupvar.location=length(frame[1,])
  groupvar=names(frame)[groupvar.location]
  vars=names(frame)[1:(groupvar.location-1)]

#Changes factor levels to a factor if not already
   if(!is.factor(frame[,groupvar]))
   {
      frame[,groupvar] <- factor(frame[,groupvar])
   }
   
#Gives levels of factor
   levels=levels(frame[,groupvar])

# Orders data by group
   o <- order(frame[,groupvar])
   frame <- frame[o,]

# Compute a-number of factors and --number of variables
   p <- length(vars)
   a <- length(levels(frame[,groupvar]))

#Defines logical variable that tells when to exit
  exit=FALSE

#Checks to see if R Matrix is singular, if so returns warning and chances to ANOVA type statistic
if(test[1]!=1){
  # Compute sample sizes per group
  N <- length(frame[,1])
  ssize <- array(NA,a)
  lims <- matrix(NA,2,a)
  for(i in 1:a){
    ssize[i] <- length(frame[frame[,groupvar]==levels(frame[,groupvar])[i],1])
    lims[1,i] <- min(which(frame[,groupvar]==levels(frame[,groupvar])[i]))
    lims[2,i] <- max(which(frame[,groupvar]==levels(frame[,groupvar])[i]))
  }

  if(sum(ssize<2)>0){return('Error: Each group must have sample size of at least 2')}
  
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
  
  if(det(H1)==0 | det(H2)==0 | det(G1)==0 | det(G2)==0 | det(G3)==0 && test[1]!=1){
    cat('Rank Matrix is Singular, only ANOVA test can be calculated \n')
    test=c(1,0,0,0)
  }else{
  
  H1 <- (1/(a-1))*H1
  H2 <- (1/(a-1))*H2
  G1 <- (1/(N-a))*G1
  G2 <- (1/(a-1))*G2
  G3 <- (1/a)*G3

  if(det(H1)==0 | det(H2)==0 | det(G1)==0 | det(G2)==0 | det(G3)==0 && test[1]!=1){
    test=c(1,0,0,0)
    cat('Rank Matrix is Singular, only ANOVA test can be calculated \n')
    }
  }
}  

##################################
#Output of which statistic is used
##################################
if(test[1]==1){cat('\nThe ANOVA type statistic will be used in the following test \n')}
if(test[2]==1){cat('\nThe Lawley Hotelling type (McKeon\'s F approximation) statistic will be used in the following test \n')}
if(test[3]==1){cat('\nThe  Bartlett-Nanda-Pillai type (Muller\'s F approximation) statistic will be used in the following test \n')}
if(test[4]==1){cat('\nThe Wilks\' Lambda type statistic will be used in the following test \n')}

#######################
#Global Hypothesis Test
####################### 
 
   base <- basenonpartest(frame,groupvar,vars,tests=test)
   if (test[1]==1){testpval=base$pvalanova}
   if (test[2]==1){testpval=base$pvalLH}
   if (test[3]==1){testpval=base$pvalBNP}
   if (test[4]==1){testpval=base$pvalWL}
   
   if( testpval < alpha) {cat('The Global Hypothesis is significant, subset algorithm will continue \n')
   } else {return ('The Global Hypothesis is not significant, subset algorithm will not continue')}


#######################
#Algorithm for when p>a
#######################
  
if(p>a || factors.and.variables==TRUE){  #Only runs if user wants to check subsets of factor levels
  
#Since the Global hypothesis is significant outputs first subset, which is the subset of all factor levels
cat('\n~Performing the Subset Algorithm based on Factor levels~\nThe Hypothesis of equality between factor levels ', levels, 'is rejected \n')

#Exit if only 2 factor levels are being tested
if (length(levels)<= 2 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
if (length(levels)<= 2 && factors.and.variables==TRUE){
  cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
  exit=TRUE}

# Step 2: Subset Algorithm for Testing Factor Levels with 1 factor removed, returns matrix where the rows are subsets(size a-1) that are significant

#Creates new matrix of all possible intersections and finds which are signficant, creates new matrix of significant subsets
if(exit==FALSE){

  #Test each of the a-1 subgroups, creates list of a-1 groups which are significant, the p-value is multiplied by a/(a-1) if a>3
  step2subsets=vector("list",a)
  for(i in 1:a)
    {
    subsetframe <-subset(frame, frame[,groupvar] != levels[i])
    subsetframe<- droplevels(subsetframe)
    groupvarsub=names(subsetframe)[groupvar.location]
    base <- basenonpartest(subsetframe,groupvarsub,vars,tests=test)
      if (test[1]==1){testpval=base$pvalanova}
      if (test[2]==1){testpval=base$pvalLH}
      if (test[3]==1){testpval=base$pvalBNP}
      if (test[4]==1){testpval=base$pvalWL}
   
      #p-value is multiplied by a/k if a>3 , where k=a-1 is the dimension of the subset being tested
    k=a-1
    if(a>3){
            if( testpval*a/k < alpha) {step2subsets[[i]]=levels(subsetframe[,groupvarsub])}
            else{step2subsets[[i]]=NA}
            if( testpval*a/k < alpha) {cat('The Hypothesis of equality between factor levels ', siglevels=levels[-i], 'is rejected  \n')}
          }else{
            if( testpval < alpha) {step2subsets[[i]]=levels(subsetframe[,groupvarsub])}
            else{step2subsets[[i]]=NA}
            if( testpval < alpha) {cat('The Hypothesis of equality between factor levels ', siglevels=levels[-i], 'is rejected  \n')} 
          }
    
    }

      step2subsets=step2subsets[!is.na(step2subsets)] #Step 2 subsets will be of length a-1

  #Checks to see if there are step2subs and if there is only one step 2 subset in either case exit is set to TRUE
    if (length(step2subsets)<= 1 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    if (length(step2subsets)<= 1 && factors.and.variables==TRUE){
    cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
    exit=TRUE}
 
}

#Checks to see if the significant subsets are of length 2, if length 2 or less, if length 2 or less function is done checking factor levels
if(exit==FALSE){
    if (length(step2subsets[[1]])<= 2 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    if (length(step2subsets[[1]])<= 2 && factors.and.variables==TRUE){
        cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
        exit=TRUE}
}

#Creates a new list of all the possible intersections of the previous significant subsets, step2subsets
if(exit==FALSE){
  step2subsetcount=length(step2subsets)
  num.intersections=((step2subsetcount-1)*step2subsetcount)/2
  newsubsets=vector("list",num.intersections)
  k=1   
   for(i in 1:(step2subsetcount-1)){
     h=i+1
      for(j in h:step2subsetcount){
        newsubsets[[k]]=intersect(step2subsets[[i]],step2subsets[[j]])
        k=k+1
      }
   }  
  
  
  #Creates a list of significant subsets and non-significant subsets
  newsubsetcount=length(newsubsets)
  sigfactorsubsets=vector("list",newsubsetcount)
  nonsigfactorsubsets=vector("list",newsubsetcount)
  for(i in 1:newsubsetcount)
   {
     subsetstotest=as.factor(newsubsets[[i]])
     subsetframe <-subset(frame, frame[,groupvar] %in% subsetstotest)
     subsetframe<- droplevels(subsetframe)
     groupvarsub=names(subsetframe)[groupvar.location]
     base <- basenonpartest(subsetframe,groupvarsub,vars,tests=test)
     if (test[1]==1){testpval=base$pvalanova}
     if (test[2]==1){testpval=base$pvalLH}
     if (test[3]==1){testpval=base$pvalBNP}
     if (test[4]==1){testpval=base$pvalWL}
     
     #p-value is multiplied by a/k if a>3 , where k is the dimension of the subset being tested
     k=length(subsetstotest)
     if(a>3){
            if( testpval*(a/k) >= alpha) {nonsigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else{nonsigfactorsubsets[[i]]=NA}
            if( testpval*(a/k) < alpha) {sigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else{sigfactorsubsets[[i]]=NA}
            if( testpval*(a/k) < alpha) {cat('The Hypothesis of equality between factor levels ', newsubsets[[i]], 'is rejected \n')}
     }else{
            if( testpval >= alpha) {nonsigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else{nonsigfactorsubsets[[i]]=NA}
            if( testpval < alpha) {sigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else{sigfactorsubsets[[i]]=NA}
            if( testpval< alpha) {cat('The Hypothesis of equality between factor levels ', newsubsets[[i]], 'is rejected \n')}
     }
  
  }
   nonsigfactorsubsets=nonsigfactorsubsets[!is.na(nonsigfactorsubsets)]
   sigfactorsubsets=sigfactorsubsets[!is.na(sigfactorsubsets)]
  
    #Checks to see if there are step2subs and if there is only one significant subset in either case exit is set to TRUE 
    if (length(sigfactorsubsets)<= 1 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    if (length(sigfactorsubsets)<= 1 && factors.and.variables==TRUE){
      cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
      exit=TRUE}
}

#Checks to see if the significant subsets are of length 2, if length 2 or less if length 2 or less function is done checking factor levels
if(exit==FALSE){
   if (length(sigfactorsubsets[[1]])<= 2 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
   if (length(sigfactorsubsets[[1]])<= 2 && factors.and.variables==TRUE){
     cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
     exit=TRUE}
}


#Step 3: Reads in the significant and non significant subsets, intersects the significant sets, checks to see if intersection is 
# non-significant matrix, for all intersections that are not in non-significant matrix creates a new matrix of subsets to check. 
#Creates new matrix of significant and non significant sets
   
#Creates a new list of all the possible intersections of the previous significant subsets

number.elements=a-2 #This is the number of elements currently in the subsets, the 2 is from removing elements in previous steps

for(l in 1:(a-4)) {  #We only run from 1:(a-4) because we are checking subsets of size a-2 or less, an stop at subsets of size 2

  if(exit==FALSE){  
    newsubsetcount=length(sigfactorsubsets)
    rows=((newsubsetcount-1)*newsubsetcount)/2
    newsubsets=vector("list",rows)
    k=1   
      for(i in 1:(newsubsetcount-1)){
      h=i+1
        for(j in h:newsubsetcount){
        newsubsets[[k]]=intersect(sigfactorsubsets[[i]],sigfactorsubsets[[j]])
        k=k+1
        }
      } 
    newsubsets=unique(newsubsets)

    #Checks to see if intsections are in non-significant subsets, assigns NA to subsets that are in non-significant sets
    if(length(nonsigfactorsubsets)>0){
      for(i in 1:length(nonsigfactorsubsets))
      {
        for(j in 1:length(newsubsets))
        {
          if(sum(!(newsubsets[[j]] %in% nonsigfactorsubsets[[i]]))==0){newsubsets[[j]]=NA}
        }
      }
    }  
    
    #Checks to see if there are new subsets to check
    if (length(newsubsets)==1 && is.na(newsubsets) && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    if (length(newsubsets)==1 && is.na(newsubsets) && factors.and.variables==TRUE){cat('All appropriate subsets using factor levels have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
        exit=TRUE}
    
    if(exit==FALSE){
    #Removes duplicates
    newsubsets=unique(newsubsets)

    
    #Removes subsets of size less than current size (1:a-2)
    for(i in 1:length(newsubsets))
    {
      if(length(newsubsets[[i]])<(number.elements-1)){newsubsets[[i]]=NA}
    }
    
    # Removes NA
    newsubsets=newsubsets[!is.na(newsubsets)]

    #Checks to see if there are new subsets to check
    if (length(newsubsets)==0  && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    if (length(newsubsets)==0 && factors.and.variables==TRUE){cat('All appropriate subsets using factor levels have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
        exit=TRUE}
    
    if(exit==FALSE){
    #Creates a list of significant subsets and non-significant subsets
    newsubsetcount=length(newsubsets)
    sigfactorsubsets=vector('list',newsubsetcount)
    nonsigfactorsubsets.new=vector('list',newsubsetcount)

    for(i in 1:newsubsetcount)
    {
     subsetstotest=as.factor(newsubsets[[i]])
     subsetframe <-subset(frame, frame[,groupvar] %in% subsetstotest)
     subsetframe<- droplevels(subsetframe)
     groupvarsub=names(subsetframe)[groupvar.location]
     base <- basenonpartest(subsetframe,groupvarsub,vars,tests=test)
     if (test[1]==1){testpval=base$pvalanova}
     if (test[2]==1){testpval=base$pvalLH}
     if (test[3]==1){testpval=base$pvalBNP}
     if (test[4]==1){testpval=base$pvalWL}
     
     #p-value is multiplied by a/k, if a>3 where k is the dimension of the subset being tested
     k=length(subsetstotest)
     if(a>3){
            if( testpval*(a/k)  >= alpha) {nonsigfactorsubsets.new[[i]]=levels(subsetframe[,groupvarsub])}else {nonsigfactorsubsets.new[[i]]=NA}
            if( testpval*(a/k)  < alpha) {sigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else {sigfactorsubsets[[i]]=NA}
            if( testpval*(a/k)  < alpha) {cat('The Hypothesis of equality between factor levels ', newsubsets[[i]], 'is rejected \n')}
     }else{
       if( testpval  >= alpha) {nonsigfactorsubsets.new[[i]]=levels(subsetframe[,groupvarsub])}else {nonsigfactorsubsets.new[[i]]=NA}
       if( testpval  < alpha) {sigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else {sigfactorsubsets[[i]]=NA}
       if( testpval  < alpha) {cat('The Hypothesis of equality between factor levels ', newsubsets[[i]], 'is rejected \n')}
       
     }
    }
   
    nonsigfactorsubsets.new=nonsigfactorsubsets.new[!is.na(nonsigfactorsubsets.new)]
    sigfactorsubsets=sigfactorsubsets[!is.na(sigfactorsubsets)]
    
    #Combine previous non-significant subsets with new
    nonsigfactorsubsets=c(nonsigfactorsubsets,nonsigfactorsubsets.new)

    #Checks to see if there are step2subs and if there is only one significant subset in either case exit is set to TRUE
    if (length(sigfactorsubsets)<= 1 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    if (length(sigfactorsubsets)<= 1 && factors.and.variables==TRUE){
      cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
      exit=TRUE}
    }
  }
}

#Checks to see if the significant subsets are of length 2, if length 2 or less if length 2 or less function is done checking factor levels
if(exit==FALSE){
if (length(sigfactorsubsets[[1]])<= 2 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
if (length(sigfactorsubsets[[1]])<= 2 && factors.and.variables==TRUE){
  cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
  exit==TRUE}
}

number.elements=number.elements-1
}
}

  
########################
#Algorithm for when p<=a
########################

if(p<=a ||factors.and.variables==TRUE){
cat('\n~Performing the Subset Algorithm based on Response Variables~ \n The Hypothesis of equality using response variables ', vars, 'is rejected \n')
if (length(vars)<= 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}

  # Step 2: Subset Algorithm for Testing Different response variabless with 1 variable removed, returns matrix where the rows are subsets(size p-1) that are significant
  #Creates new matrix of all possible intersections and finds which are signficant, creates new matrix of significant subsets
  
  #Create multiplier for multiple testing procedure
 
  #Test each of the p-1 subgroups, creates list of p-1 groups which are significant, for the multiple testing procedure the p-value must be 
  #multiplied by p choose # in subset, p-1, so in this case the multiplier is p
  step2subsets=vector("list",p)
  for(i in 1:p)
  {
    base <- basenonpartest(frame,groupvar,vars[-i],tests=test)
    if (test[1]==1){testpval=base$pvalanova}
    if (test[2]==1){testpval=base$pvalLH}
    if (test[3]==1){testpval=base$pvalBNP}
    if (test[4]==1){testpval=base$pvalWL}
    
    
    if( testpval*p < alpha) {step2subsets[[i]]=vars[-i]}
    else{step2subsets[[i]]=NA}
    if( testpval*p < alpha) {cat('The Hypothesis of equality using response variables ', sigvariables=vars[-i], 'is rejected \n')}
  }
  
  step2subsets=step2subsets[!is.na(step2subsets)]
  
 
  if (length(step2subsets)<= 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  if (length(step2subsets[[1]])<= 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  
#Creates a new list of all the possible intersections of the previous significant subsets
  
  step2subsetcount=length(step2subsets)
  num.intersections=((step2subsetcount-1)*step2subsetcount)/2
  newsubsets=vector("list",num.intersections)
  k=1   
  for(i in 1:(step2subsetcount-1)){
    h=i+1
    for(j in h:step2subsetcount){
      newsubsets[[k]]=intersect(step2subsets[[i]],step2subsets[[j]])
      k=k+1
    }
  }  
  
  
#Creates a list of significant subsets and non-significant subsets, subsets are now of length p-2,  p-value is multiplied by the number of test performded
  newsubsetcount=length(newsubsets)
  sig.variable.subsets=vector("list",newsubsetcount)
  nonsig.variable.subsets=vector("list",newsubsetcount)
  multiplier=newsubsetcount  
  
  for(i in 1:newsubsetcount)
  {

    base <- basenonpartest(frame,groupvar,vars=newsubsets[[i]],tests=test)
    if (test[1]==1){testpval=base$pvalanova}
    if (test[2]==1){testpval=base$pvalLH}
    if (test[3]==1){testpval=base$pvalBNP}
    if (test[4]==1){testpval=base$pvalWL}
    
    if( testpval*multiplier >= alpha) {nonsig.variable.subsets[[i]]=newsubsets[[i]]}else{nonsig.variable.subsets[[i]]=NA}
    if( testpval*multiplier < alpha) {sig.variable.subsets[[i]]=newsubsets[[i]]}else{sig.variable.subsets[[i]]=NA}
    if( testpval*multiplier < alpha) {cat('The Hypothesis of equality using response variables ', newsubsets[[i]], 'is rejected  \n')}
  }

  nonsig.variable.subsets=nonsig.variable.subsets[!is.na(nonsig.variable.subsets)]
  sig.variable.subsets=sig.variable.subsets[!is.na(sig.variable.subsets)]
  if (length(sig.variable.subsets)<= 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  if (length(sig.variable.subsets[[1]])<= 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  
  #Step 3: Reads in the significant and non significant subsets, intersects the significant sets, checks to see if intersection is 
  # non-significant list, for all intersections that are not in non-significant matrix creates a new matrix of subsets to check. 
  #Creates new list of significant and non significant sets
  
  #Creates a new list of all the possible intersections of the previous significant subsets
  
  number.elements=p-2
  if (number.elements== 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  
  for(l in 1:(p-3)) {  
    
    newsubsetcount=length(sig.variable.subsets)
    
    rows=((newsubsetcount-1)*newsubsetcount)/2
    newsubsets=vector("list",rows)
    k=1   
    for(i in 1:(newsubsetcount-1)){
      h=i+1
      for(j in h:newsubsetcount){
        newsubsets[[k]]=intersect(sig.variable.subsets[[i]],sig.variable.subsets[[j]])
        k=k+1
      }
    } 
    newsubsets=unique(newsubsets)
    
    
    #Checks to see if intsections are in non-significant subsets, assigns NA to subsets that are in non-significant sets
    if(length(nonsig.variable.subsets)>0){
      for(i in 1:length(nonsig.variable.subsets))
      {
        for(j in 1:length(newsubsets))
        {
          if(sum(!(newsubsets[[j]] %in% nonsig.variable.subsets[[i]]))==0){newsubsets[[j]]=NA}
        }
      }
    }  
    #Removes duplicates
    newsubsets=unique(newsubsets)
    
    #Checks to see if there are new subsets to check
    if (length(newsubsets)==1 && is.na(newsubsets)){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    
    #Removes subsets of size less than current size (1:p-2)
    for(i in 1:length(newsubsets))
    {
      if(length(newsubsets[[i]])<(number.elements-1)){newsubsets[[i]]=NA}
    }
    # Removes NA
    newsubsets=newsubsets[!is.na(newsubsets)]
       
    #Creates a list of significant subsets and non-significant subsets
    newsubsetcount=length(newsubsets)
    sig.variable.subsets=vector('list',newsubsetcount)
    nonsig.variable.subsets.new=vector('list',newsubsetcount)
    multiplier=newsubsetcount   #The multiplier is the number of test peformed
    for(i in 1:newsubsetcount)
    {
     
      base <- basenonpartest(frame,groupvar,vars=newsubsets[[i]],tests=test)
      if (test[1]==1){testpval=base$pvalanova}
      if (test[2]==1){testpval=base$pvalLH}
      if (test[3]==1){testpval=base$pvalBNP}
      if (test[4]==1){testpval=base$pvalWL}
      
      if( testpval*multiplier >= alpha) {nonsig.variable.subsets.new[[i]]=newsubsets[[i]]}else {nonsig.variable.subsets.new[[i]]=NA}
      if( testpval*multiplier < alpha) {sig.variable.subsets[[i]]=newsubsets[[i]]}else {sig.variable.subsets[[i]]=NA}
      if( testpval*multiplier < alpha) {cat('The Hypothesis of equality using response variables ', newsubsets[[i]], 'is rejected \n')}
    }
    
    nonsig.variable.subsets.new=nonsig.variable.subsets.new[!is.na(nonsig.variable.subsets.new)]
    
    #Combine previous non-significant subsets with new
    nonsig.variable.subsets=c(nonsig.variable.subsets,nonsig.variable.subsets.new)
    
    sig.variable.subsets=sig.variable.subsets[!is.na(sig.variable.subsets)]
    
    #Removes Duplicates
    sig.variable.subsets=unique(sig.variable.subsets)
    
    # Removes NA
    sig.variable.subsets=sig.variable.subsets[!is.na(sig.variable.subsets)]
    
    number.elements=number.elements-1
    if (length(sig.variable.subsets)<= 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
    if (length(sig.variable.subsets[[1]])<= 1){return(cat('All appropriate subsets using response variables have been checked using a multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  } 
  
  
  
}
}   

