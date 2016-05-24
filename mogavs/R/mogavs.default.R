mogavs.default <-
function(x,y,maxGenerations=10*ncol(x),popSize=ncol(x),noOfOffspring=ncol(x),crossoverProbability=0.9,mutationProbability=1/ncol(x),kBest=1,plots=F,additionalPlots=F, ...){
  if(missing(x))
    stop("Arg x is missing")
  if(missing(y))
    stop("Arg y is missing")
  if(missing(maxGenerations))
    stop("Arg maxGenerations is missing")
  
  if(is.matrix(x)==FALSE) x<-as.matrix(x)

  if(is.matrix(y)==FALSE) y<-as.matrix(y)
  

#Initializing the parent poulation named popMembers. Moreover, calculating MSE and the number of variables(NumVar) of each string
temp<-.initializePop(popSize, x, y)
popMembers<-temp[[1]]
numOfVariables<-temp[[2]]
MSE<-temp[[3]]
rm(temp) 

#Finding the non-dominated set of popMembers
temp<-.nonDomination(numOfVariables,MSE,popMembers)
eliteMembers<-temp[[1]]
dominatedSet<-temp[[2]]
#Note: ignoring other output values
rm(temp)

numOfVariablesInitial <- numOfVariables
MSEInitial <- MSE

archiveSet<-popMembers
obj1ArchiveSet<-t(numOfVariables)
obj2ArchiveSet<-t(MSE)

#initialize empty values
obj1Augmented<-""
obj2Augmented<-""


for(genCount in 1:maxGenerations){
  #Producing offspring
  offspring<-.crossover(noOfOffspring,popMembers,eliteMembers,crossoverProbability)

  offspring<-.mutation(offspring,mutationProbability)

  #Evaluating offspring
  temp<-.evaluateOffspring(offspring,x,y)
  numOfVariablesOffspring<-temp[[1]]
  MSEOffspring<-temp[[2]]
  rm(temp)

  archiveSet<-rbind(archiveSet,offspring)
  obj1ArchiveSet<-rbind(obj1ArchiveSet,t(numOfVariablesOffspring))
  obj2ArchiveSet<-rbind(obj2ArchiveSet,t(MSEOffspring))
  
  #Adding the offspring to the parent population
  popMembersAugmented<-rbind(popMembers,offspring)
  obj1Augmented<-rbind(t(numOfVariables),t(numOfVariablesOffspring))
  obj2Augmented<-rbind(t(MSE),t(MSEOffspring))
  
  #Finding the non-dominated set from Augmented set

  temp<-.removeDuplicates(popMembersAugmented,obj1Augmented,obj2Augmented)
  popMembersAugmented<-temp[[1]]
  obj1Augmented<-temp[[2]]
  obj2Augmented<-temp[[3]]
  rm(temp)
  
  #popMembersAugmented is updated as [nonDominatedSet; dominatedSet]
  #Note: First part of obj1Augmented and obj2Augmented store information about non-dominated set
  #Note: Second part of obj1Augmented and obj2Augmented store information about dominated set

  temp<-.nonDomination(obj1Augmented,obj2Augmented,popMembersAugmented)

  nonDominatedSet<-temp[[1]]
  dominatedSet<-temp[[2]]
  popMembersAugmented<-temp[[3]]
  obj1Augmented<-temp[[4]]
  obj2Augmented<-temp[[5]]
  rm(temp)

  temp<-.sortByComplexity(nonDominatedSet,obj1Augmented[1:nrow(nonDominatedSet)],obj2Augmented[1:nrow(nonDominatedSet)])
  nonDominatedSet<-temp[[1]]
  obj1Augmented[1:nrow(nonDominatedSet)]<-temp[[2]]
  obj2Augmented[1:nrow(nonDominatedSet)]<-temp[[3]]
  rm(temp)
  popMembersAugmented[1:nrow(nonDominatedSet),]<-nonDominatedSet

  #Update population members for next generation
  temp<-.updatePopulationMembers(nonDominatedSet,dominatedSet,popMembersAugmented, obj1Augmented, obj2Augmented, popSize)
  popMembers<-temp[[1]]
  numOfVariables<-temp[[2]]
  MSE<-temp[[3]]
  eliteMembers<-temp[[4]]
  rm(temp)
  
  print(paste("Generation ",genCount," / ",maxGenerations,sep=""))
  
  if(plots==TRUE){
    plot(numOfVariables,MSE,pch=8,col="blue")
    legend("topright",paste("Population Members: Gen = ",toString(genCount)))
  }

  
}
#update MSE and numOfVariables to only dominated members
MSE<-MSE[1:nrow(nonDominatedSet)]
numOfVariables<-numOfVariables[1:nrow(nonDominatedSet)]
if(plots==TRUE){
  #plot the initial members from the archive set
  plot(numOfVariablesInitial,MSEInitial,pch=8,col="red")
  points(numOfVariables,MSE,pch=8,col="blue")
  legend("topright",c("Population Members","Initial Members"),col=c("blue","red"),pch=c(8,8))
}

#Remove duplicate members from the archive set
temp<-.removeDuplicates(archiveSet, obj1ArchiveSet, obj2ArchiveSet)
archiveSet<-temp[[1]]
obj1ArchiveSet<-temp[[2]]
obj2ArchiveSet<-temp[[3]]

temp<-list("nonDominatedSet"=nonDominatedSet,"numOfVariables"=numOfVariables,"MSE"=MSE,"archiveSet"=archiveSet,"kBest"=kBest,"maxGenerations"=maxGenerations,"crossoverProbability"=crossoverProbability,"noOfOffspring"=noOfOffspring,"popSize"=popSize, "n_obs"=nrow(x), "obj1ArchiveSet"=obj1ArchiveSet, "obj2ArchiveSet"=obj2ArchiveSet)

#set returned var to be of class "mogavs"
class(temp)<-"mogavs"
if(additionalPlots==TRUE){
  createAdditionalPlots(temp,kBest,method="kbest")
}
return(temp)
}
