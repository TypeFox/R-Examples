ClassifyByDecisionBoundaries=function(Data,DecisionBoundaries,ClassLabels){
# Cls = ClassifyByDecisionBoundaries(Data,DecisionBoundaries)
# Classify Data according to decision Boundaries
#
# INPUT
# Data(1:n,1)                 vector of Data
# DecisionBoundaries(1:L)    decision boundaries 
# OPTIONAL
# ClassLabels(1:L+1)        numbered class labels that are assigned to the classes. default (1:L) 
# OUTPUT
# Cls(1:n,1:d)               classiffication of Data
# Author MT 04/2015


if(!is.vector(Data)){
  warning('Data converted to vector')
 Data=as.vector(Data) 
}
if(is.list(DecisionBoundaries)){
  DecisionBoundaries=as.vector(DecisionBoundaries$DecisionBoundaries)
  print('DecisionBoundaries was a list, assuming usage of BayesDecisionBoundaries()')
}
AnzBounds= length(DecisionBoundaries)
if(missing(ClassLabels)){
ClassLabels=seq(from=1,by=1,to=(AnzBounds+1))
}

Cls=rep(1,length(Data)) # default alles in Klasse 1

for(b in 1:AnzBounds){
  
  ind=Data>DecisionBoundaries[b]

    Cls[ind] = rep(ClassLabels[b+1],sum(ind))
} # for c

return(Cls)
}
