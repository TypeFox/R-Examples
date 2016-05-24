## GNU License: written by Moreno I. Coco, 05/2013 

## Recurrence is calculating using contingency tables
## ref: Dale, R., Warlaumont, A. S., & Richardson, D. C. (2011).
## Nominal cross recurrence as a generalized lag sequential analysis for
## behavioral streams. International Journal of Bifurcation and Chaos,
## 21, 1153-1161. (special issue on recurrence)

## pro: every lag we have full co-occurence matrices of objects 
## cons: we do not actually use this info: only look at diagonal recurrence 

## Shift = a numerical vector for the delays, e.g., seq(1,100, 1)
## lags = a numerical vector, equal to shift, with the actual ms
## time stamps, e.g., seq(-3000, 3000, 25);

 ##CT; slower but can track objects co-occurences
 ## needs categorical data
              

## Build a contingency table at every delay given a set
## of objects (categorical states) of the two series

.packageName <- 'crqa'

CT.build <- function(ts1, ts2, set.object, lag){
 
CT = matrix(0,ncol=length(set.object),nrow=length(set.object))
T= length(ts1)

for (t in 1:(T-lag)){
  
  ind.col = grep( paste( "^", ts1[t],"$",sep = "" ), set.object)
  ind.row = grep( paste ("^", ts2[t+lag],"$",sep = "" ), set.object)
  CT[ ind.row, ind.col ] = CT[ ind.row, ind.col ] + 1

}

return(CT)

}

## Extract the recurrence of a particular pair of objects
## from the contingency table.

RR.pair <- function(ts1,ts2,obj1,obj2,set.object,lag){

  T = length(ts1);
  CT = CT.build(ts1,ts2,set.object,lag)
     i = grep(paste(obj1,"$",sep=""), set.object);
     j = grep(paste(obj2,"$",sep=""), set.object);

     RR = sum(CT[i,j])/(T-lag)

if (length(RR) == 0){
  RR =0}
  
return(RR)

}

### Diagonal-wise recurrence: calculating the recurrence along the diagonal of common pairs sequences that have to be compared

RecMatch = function(ts1,ts2,set.object,lag){

Rec = 0

for (i in 1:length(set.object)){
  obj=set.object[i]
  Rec = Rec + RR.pair(ts1,ts2, obj, obj, set.object, lag)
  }

return(Rec)

}


## recombine the recurrences such that it unfolds over negative and positive
## lags: i.e., speaker vs listener leading interpretation
## make sure that lags are correctly assigned.

combineside <- function(rec1, rec2){
  
  t.len = length(rec2)
  effind = seq(t.len, 1, -1)
  rec1 = rec1[effind]
  rec2 = rec2[2:length(rec2)]
  
  recTot = c(rec1,rec2)
      
  return (recTot)

}


################## use the functions above #################


CTcrqa <- function(ts1, ts2, par){
  datatype = thrshd = lags = NULL
  for (v in 1:length(par)) assign(names(par)[v], par[[v]]) ## assign parameters

  tryCatch({
      
      res = checkts(ts1, ts2, datatype, thrshd) ## first check that sequences
                                                ## have the same length 
      if ( res[[2]] == TRUE ){
      
          tsnorm = res[[1]]
          ts1 = tsnorm[,1]; ts2 = tsnorm[,2]
          objects = unique( c(ts2, ts1) ) ##take set of unique objects
  
          rec1 = rec2 = vector() 
                                        #  here there might be an issue to save info: check out
          
          for (lg in lags){
              
              if (lg == length(lags)/2 ){
                  
                  print ("Half way through the lags") }
              
              rec1 = c(rec1, RecMatch(ts1, ts2, objects, lg) )          
              rec2 = c(rec2, RecMatch(ts2, ts1, objects, lg) )  
      
          }
    
          print( paste( "Mean recurrence of pair:",
                       mean( c( rec1,rec2 ), na.rm = T ) ) )
    
          res = combineside(rec1, rec2) ## order and stitch together recurrence with the convention recurrence at rec2 is observed at positive lags
          
          return(res)
    
      } else {
    
          res = vector()
          print ( paste( "Sequences had a difference of", res[[1]], "units") )

          return(res)
   
      }
  }, warning = function(war) { ## here exception can be handled
    
      print(paste("WARNING:  ", war))
   
  }, error = function(err) {
 
      print(paste("ERROR:  ", err))
  })
  
  
}

