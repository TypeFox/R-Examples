# the function body:
getStdRate= function(relativeRate, model, referencePopulation, scale=100000) {

# check for negatives
if(any(relativeRate < -0.0001,na.rm=TRUE ))
  warning("negative numbers in rates, make sure they're not on the log scale")

  
newpop <- formatCases(referencePopulation)

newpop <- newpop[newpop$POPULATION!=0,]

newpop$logpop <- log(newpop$POPULATION) + log(scale)



if(is.numeric(model)) {
    # model is a vector of rates
        # check breaks for groups, make sure they line up
        rateBreaks =getBreaks(names(model))
        names(model) = rateBreaks$newNames
        newpop = formatCases(referencePopulation, ageBreaks=rateBreaks)
        
        referenceRate =sum(newpop$POPULATION * model[paste(newpop$sex, newpop$ageNumeric, sep=".")]) * scale

        
    } else {
# model is a glm object

# predict.glmZeros

# find the ages used in the model and get rid of those that do not exist in model
modelAges = model$xlevels$age
newpop= newpop[newpop$age %in% modelAges,]
#remove main effects that not in the model
for(Dlevel in names(model$xlevels)) {
alllevels = levels(newpop[[Dlevel]])
if(!all(alllevels %in% model$xlevels[[Dlevel]])) {
	# remove the rows from the missing levels
	tokeep = newpop[[Dlevel]] %in% model$xlevels[[Dlevel]]
	newpop = newpop[tokeep,]
 }
}


#Take out interactions with NA, work with model without interactions
#temp<-row.names(summary(model)$coefficients)[-1]
interactNA<-names(model$coefficients)[is.na(model$coefficients)]
#if there is any coeff with NA, remove
if(length(interactNA)>0){
  interact<-grep(":",interactNA,value=TRUE)
  #keep the rows that the interaction is not NA,
  newpop$param<-paste(paste("age",newpop$age,sep=""),paste("sex",newpop$sex,sep=""),sep=":")
  newpop = newpop[!newpop$param %in% interact,]
  newpop$param<-NULL
}

#find expected


		
referenceRate <- sum(predict(model, newpop[,c("age","sex","logpop")], type = "response"))

}



	newpop= (relativeRate*referenceRate)
	
newpop
}



#### function of list format: 

#getstdRate.list()<-function(relativeRate, model, referencePopulation, scale=100000){
#lapply(referencePopulation,getstdRate,relativeRate=relativeRate, scale=scale, model=model)
#}


#### the most original idea of the function: 
#getStdRate= function(relativeRate, model, referencePopulation, scale=100000) {

#relativeRate * sum(predict(model, referencePopulation)) * scale

#}
