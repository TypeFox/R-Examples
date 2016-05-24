############## Perform DPA Analysis ###############
dpa.analysis.performDPA<- function(){


#############  Function to Specify the model based on the relations ################
dpa.specifyModel<-function(){
row<-NULL
n=nrow(relations)
arrowUni <- "->"
arrowBi <- "<->"
paramNumber=0
parameters <- NULL
baseVariables <- NULL

# for each of the rows, create the right format
for(i in 1:n){ 
  numberOfLags=as.numeric(relations[i,5])-as.numeric(relations[i,4])+1

  # for each relation, for each lag we create a relationship in the right format
  for (j in 1:numberOfLags){
    fromVar<-relations[i,1]
    toVar<-relations[i,2]
    arrowtype<-NULL
    relationshipFixedToOne<-FALSE
    lagNumber=as.numeric(relations[i,4])+j-1

    #if we're discussing a lagged variable now, we have to override one of the variable names
    if(lagNumber>0){

      #if we're lagging the TO variable...
      if(relations[i,3]=="To"){
        # change the toVar-name
        toVar<-paste(relations[i,2],"_L",lagNumber,sep="")
      } else {
        # change the fromVar-name
        fromVar<-paste(relations[i,1],"_L",lagNumber,sep="")
      }
    }

    #select one of the types of arrows
    if(relations[i,6]=="UniDirectional"){
      #one way
      arrowtype<-arrowUni
    } else {
      #two ways
      arrowtype<-arrowBi
    } 

    #is the fromVar already in the set of variables?
    foundFrom=FALSE
    if(!is.null(baseVariables)){
      for(k in 1:length(baseVariables)){
        if(as.character(fromVar)==baseVariables[k]){
          foundFrom=TRUE
        }
      }
    }
    if(!foundFrom){
      #put it in then...
      baseVariables <- rbind(baseVariables,as.character(fromVar))
    }

    #is the toVar already in the set of variables?
    foundTo=FALSE
    if(!is.null(baseVariables)){
      for(k in 1:length(baseVariables)){
        if(as.character(toVar)==baseVariables[k]){
          foundTo=TRUE
        }
      }
    }
    if(!foundTo){
      #put it in then...
      baseVariables <- rbind(baseVariables,as.character(toVar))
      #In this case, we need this relationship set to one...
      relationshipFixedToOne=TRUE
    }

    if(relationshipFixedToOne){
       row <- rbind(row,cbind(paste(fromVar," ",arrowtype," ",toVar,sep=""),NA,1))
       #print("added one... is fixed")       
       #print(row)
    }else{
       paramNumber=paramNumber+1;
       thisparameter<-paste("param",paramNumber,sep="")
       row <- rbind(row,cbind(paste(fromVar," ",arrowtype," ",toVar,sep=""),thisparameter,NA))
       parameters<-rbind(parameters,thisparameter)
       #print("added one... not fixed")       
       #print(row)
    }
  }
}

#print("Done...")
#print(baseVariables)

################# To add variance relations #######################
for(i in 1:length(baseVariables)){
  paramNumber=paramNumber+1
  thisparameter<-paste("param",paramNumber,sep="")
  row <- rbind(row,cbind(paste(baseVariables[i]," ",arrowBi," ",baseVariables[i],sep=""),thisparameter,NA))
  parameters<-rbind(parameters,thisparameter)
}

assign("row",row,env=.GlobalEnv)
assign("variables",baseVariables,env=.GlobalEnv)
}



################## Start data selection and performing DPA ###################
dpa.analysis.doPerform <- function(thecurrenttick=NULL){
  sem.DPA<-NULL
  sem.standardized<-NULL
  assign("sem",sem,env=.GlobalEnv)
  assign("sem.DPA",sem.DPA,env=.GlobalEnv)
  assign("sem.standardized",sem.standardized,env=.GlobalEnv)

  relevantData <- NULL

if (rbVal=="time_irrespective"){
  print(paste("Analysing irrespective of time"))
  for(i in 1:nrow(variables)){
    if(i==1){
       relevantData <- rbind(relevantData,e[variables[i]])
    }else{
       relevantData <- cbind(relevantData,e[variables[i]])
    }
  }
}


if (rbVal=="every_timeStep"){
  #print("ANALYSIS FOR SINGLE TICK: ")
  print(paste("Analysing tick ",thecurrenttick))
  for( i in 1:nrow(variables)){
    if(i==1){
      relevantData <- rbind(relevantData,e[which(e$tick==thecurrenttick),][variables[i]])
    }else{
      relevantData <- cbind(relevantData,e[which(e$tick==thecurrenttick),][variables[i]])
    }
  }
}

if (rbVal=="time_interval"){
  print("ANALYSIS FOR Time interval: ")
  for( i in 1:nrow(variables)){
    if(i==1){
      relevantData <- rbind(relevantData,e[which(e$tick==NumTick),][variables[i]])
    }else{
      relevantData <- cbind(relevantData,e[which(e$tick==NumTick),][variables[i]])
    }
  }
}


assign("relevantData",relevantData,env=.GlobalEnv)
#print("Number of cases selected from data:")
#print(nrow(relevantData))


################  Create covariance matrix from the relevant data  #####################
covMatrix<-cov(relevantData)

##################  Perform sem with the model resulting from relations and covariance matrix  ######################
try(sem.DPA<-sem(row,covMatrix,nrow(relevantData),maxiter=10000),TRUE)

#check if it worked, then save the results
if(!is.null(sem.DPA)){
  #record standardized results/parameters
  sem.standardized<-std.coef(sem.DPA)
  sem.results.parameters<-t(t(sem.standardized[,3]))
  sem.results.statistics<-cbind(sem.results.statistics,rbind(
     summary(sem.DPA)$iterations,
     summary(sem.DPA)$df,
     summary(sem.DPA)$GFI,
     summary(sem.DPA)$AGFI,
     summary(sem.DPA)$RMSEA[1],
     summary(sem.DPA)$SRMR,
     summary(sem.DPA)$NFI,
     summary(sem.DPA)$NNFI
  ))
  sem.results.coefficients<-cbind(sem.results.coefficients,sem.standardized[,2])
  if(rbVal=="every_timeStep"){
    listOfTicks<-cbind(listOfTicks,thecurrenttick)
  }

  # assign all results globally
  assign("sem.DPA",sem.DPA,env=.GlobalEnv)
  assign("sem.results.coefficients",sem.results.coefficients,env=.GlobalEnv)
  assign("sem.results.parameters",sem.results.parameters,env=.GlobalEnv)
  assign("sem.standardized",sem.standardized,env=.GlobalEnv)
  assign("sem.results.statistics",sem.results.statistics,env=.GlobalEnv)
  assign("listOfTicks",listOfTicks,env=.GlobalEnv)
}else{
  print("No solution found for this tick")
}

}


### Start here
assign("listOfTicks",NULL,env=.GlobalEnv)
assign("sem.DPA",NULL,env=.GlobalEnv)
assign("sem.standardized",NULL,env=.GlobalEnv)
assign("sem.results.parameters",NULL,env=.GlobalEnv)
assign("sem.results.statistics",NULL,env=.GlobalEnv)
assign("sem.results.coefficients",NULL,env=.GlobalEnv)


dpa.specifyModel()

if(rbVal=="time_irrespective"){
  dpa.analysis.doPerform()
  if(!is.null(sem.DPA)){
    if(sem.DPA$convergence==1){
      #print("DPA was successful and converged")
      dpa.results.viewRelationsPlots()
    } else {
      #print("DPA was unsuccessful")
    }
  }
}


if (rbVal=="every_timeStep"){
for( i in min(e$tick):max(e$tick)){
dpa.analysis.doPerform(i)
  if(!is.null(sem.DPA)){
    if(sem.DPA$convergence==1){
      #print("DPA was successful and converged")
      dpa.results.viewRelationsPlots(i)
    } else {
      #print("DPA was unsuccessful")
    }
  }
}


}
}