# GNU License: written by Moreno I. Coco (moreno.cocoi@gmail.com)
# based on code of Rick Dale (rdale@ucmerced.edu)

# compute diagonal recurrence profile for two
# categorical time series, t1 t2 (x / y)
# the ws is the total nr. of lags over which the series are evaluated

.packageName <- 'crqa'

drpdfromts <- function(t1, t2, ws, datatype, radius){

drpd = vector();

if (datatype == "categorical"){ ## just to make sure that data is given as vector

    t1 = as.character( as.matrix( t1 ) )
    t2 = as.character( as.matrix( t2 ) )

}

if (datatype == "continuous"){
    
    t1 = as.numeric( as.matrix( t1 ) )
    t2 = as.numeric( as.matrix( t2 ) )

}


for (i in (-ws-1):-2){
  
    ix = abs(i);
    y = t2[ix:length(t2)];
    x = t1[1:length(y)];

    if (datatype == "categorical"){

        drpd = c(drpd, sum(y==x)/length(y))}
    
    if (datatype == "continuous"){
        
        dif = abs(x-y)
        recpoint = length( which(dif <= radius) )
        drpd = c(drpd, recpoint / length(y))

    }
    
    ## drpd = [drpd ; sum(y==x)/sum((y<9|x<9))]; if we have to conditionalize the probabilities
}

if (datatype == "categorical"){
    drpd = c(drpd,  sum(t1==t2)/length(t1)); ## this is lag 0    
}

if (datatype == "continuous"){
    dif = abs(t1-t2)
    recpoint = length( which(dif <= radius) )
    drpd = c(drpd, recpoint / length(t1))
}



for (i in 2:(ws+1)){
    
    x = t1[i:length(t1)]
    y = t2[1:length(x)]
    
    if (datatype == "categorical"){
        drpd = c(drpd,  sum(y==x)/length(y)) ## this is lag 0    
    }
    
    if (datatype == "continuous"){
        dif = abs(x-y)
        recpoint = length( which(dif <= radius) )
        drpd = c(drpd, recpoint / length(y))
    }

}

## extract  max recurrence and the lag at which it occurred

maxrec = max(drpd);
maxlag = which(drpd == maxrec);

return( list(profile = drpd, maxrec = maxrec, maxlag = maxlag) )

}
