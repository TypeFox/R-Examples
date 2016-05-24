#######
##Confusion matrix
#######
confusionMatrix <- function (actual, predicted){
    actual<-as.vector(actual)
    predicted<-as.vector(predicted)

    t<-table(predicted,actual)
    if(which.max(t)!=1)
    {
      t<-t[apply(t,2,function(x) which.max(x)),]
    }
    return(t)
  
}


###########################################
#The function belows returns the purity 
#values for  the classes returned by the
#EM algorithm
##########################################

matchCluster<-function(actual,predicted)
{
  actual<-as.vector(actual)
  predicted<-as.vector(predicted)
  #This function returns the corresponding label for each component returned
  #by the algorithm.
  
  
  #numSpecies <- table(actual)
  
  right_predict<-list() #keeps track of the number of data correctly classified
  purity<-0
  species <-unique(actual)
  # species<-species[1:length(levels(as.factor(predicted)))]
  components<-rep(NA,length(species)) #Keeps track of the species corresponding to what label
  componentTaken<-list() #This list ensures that there is uniqueness between components and labels
  for(sp in species)
  {
    rownames<-which(actual==sp) #get the data that has this species
    component<-predicted[rownames]
    #gets the component that appears most in the classification for the specie
    tab<-table(component)
    major_component<-max(tab) #gets the highest component occurence
    component_name<-which.max(tab) #gets the corresponding component
    label<-names(component_name)
    
    
    #ensure uniqueness of component values and labels
    
    while(length(component)!=0&&label %in% componentTaken)
    {
      
      componentlabel<-component[component==label]
      if((length(component)-length(componentlabel))>0)
      {
        component<-component[component!=label]
        tab<-table(component)
        if(dim(tab)!=0)
        {  major_component<-max(tab)  
           component_name<-which.max(tab) #gets the corresponding component
        }
        
        
      }
      else
      {
        break
      }
      
      label<-names(component_name)
      
    }
    
    label<-as.integer(label)
    purity<-purity + major_component
    componentTaken[[label]]<-label
    
    if(is.na(components[[label]])==TRUE)
    {
      components[[label]]<-sp  #places the corresponding actual label at the index of the predicted label
    }
    
    
    
    
    
    right_predict[[label]]<-major_component #places the number of correctly classified data at the index of the predicted label
    
  }
  unallocated <- which(is.na(components)) #The components with no label
  unidentified <- which(!species %in% components) #species not in components
  alloc_size<-length(unallocated)
  
  if(alloc_size>0)
  {
    for(k in 1:alloc_size)
    {
      components[[unallocated[[k]] ]]<- species[[unidentified[[k]] ]]
    }
  }
  
  
  components<-unlist(components) #changes the list of lists to a single list
  componentData<-data.frame(species=species,labels=components)
  
  
  
  return(components)
  
}



##############################################################################
#The function belows returns the Normalized Mutual Information of the algorithm
##############################################################################
normalizedMI <- function(trueLabel,predictedLabel)
{
  trueLabel<-as.vector(trueLabel)
  predictedLabel<-as.vector(predictedLabel)
  
  if(length(trueLabel)!=length(predictedLabel))
  {
    stop("The actual and predicted values must have the same size")
  }
  size <-length(trueLabel)
  truedataLabels<-unique(trueLabel)
  predicteddataLabels<-unique(predictedLabel)
  
  
  #combine the vectors so as to get the priorXY
  combine<-cbind(trueLabel,predictedLabel)
  overallsum<-0
  entropyX<-0
  entropyY<-0
  for(i in 1:length(truedataLabels))
  {
    X <-truedataLabels[i]
    priorX <- (length(trueLabel[trueLabel==X]))/size
    entropyX <- entropyX + (priorX * log(priorX))
    innersum<-0
    for(j in 1:length(predicteddataLabels))
    {
      Y<-predicteddataLabels[j]
      priorY <- (length(predictedLabel[predictedLabel==Y]))/size
      
      
      priorXY <- (nrow(subset(combine,trueLabel==X&predictedLabel==Y)))/size
      
      if(priorXY==0)
      {
        innersum <-innersum+0
      }
      else
      {
        innersum <-innersum+ (priorXY*(log(priorXY)-(log(priorX)+log(priorY))))
      }
      
      
      
      
      if(i==1)
      {
        entropyY <- entropyY + (priorY*log(priorY))
      }
      
    }
    
    overallsum<-overallsum+innersum
    
  }
  entropyX = entropyX * (-1)
  entropyY = entropyY * (-1)
  
  
  normalized<- overallsum/{{entropyX+entropyY}/2}
  return(normalized)
}

####################################################################################
# A function to calculate the False Positive Rate, and False Negative Rate of the
# classifier.
###################################################################################

errorRate<-function(actual,predicted,beta=1)
{
  #This function calculates the false positive and false negative rate
  #The false positive rate is the number of data that are in the same class but
  #different clusters, while the false positive rate is the number of data
  #that are in different classes but the same cluster.
  
  actual <- as.vector(actual)
  predicted<-as.vector(predicted)
  predicted<-as.factor(predicted)
  #level<-as.vector(levels(predicted))
  level<-unique(predicted)
  classify<-data.frame(actual=actual,predicted=predicted)
  #Below is to calculate the true positive i.e same class and same cluster
  TP<-0
  for(i in 1:length(level))
  {
    truepositive<-classify[classify$predicted==level[i],1]
    truePositive_vector <- as.matrix(table(truepositive))
    #get the combination for each class in the same component, and add them together ...the second value is always 2
    TP<-TP+sum(apply(truePositive_vector,1,function(x) combinations(x)))
    
  }
  #The number of pairs in the same cluster, we need to know the total pair of data in the same cluster
  same_cluster <- as.matrix(table(classify$predicted))
  same_cluster<-sum(apply(same_cluster,1,function(x) combinations(x)))
  FP<-same_cluster - TP
  
  #To get the false negative and true negative, we need to know the total pair of data in the same class
  same_class<- as.matrix(table(classify$actual))
  same_class<-sum(apply(same_class,1,function(x) combinations(x)))
  
  # same_class = TP + FN
  FN<- same_class - TP
  #This is the total number that is TP+FP+TN+FN
  Total_number = combinations(length(actual))
  
  #From total number calculate TN
  TN<- Total_number - (TP+FP+FN)
  #Below are the False Positive and False Negative Rates
  #False Positive Rate is FP/FP+TP
  FPR <- FP/(FP+TN)
  #False Negative Rate is FN/FN+TN
  FNR <- FN/(TP+FN)
  RI<-(TP+TN)/(TP+FP+TN+FN)
  
  P<- TP/(TP+FP)
  R<-TP/(TP+FN)

  

  F_measure<-((beta^2)+1)*((P*R)/(((beta^2)*P)+R))
  FM<- sqrt(P*R)
  
  return(list(FPR=FPR,FNR=FNR,RI=RI,F=F_measure,mallow=FM))
  
  
}

combinations<-function(n)
{
  perm<- (n*(n-1))/2
  return(perm)
}
