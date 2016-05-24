lagdata <-
function(data,varlist, laglength) {
 
    	mydata=data

    length(varlist)
    for (iv in 1:length(varlist))
    { 
        varname=varlist[iv] #"meantemp"
        ncc=dim(mydata)[2]
        nrr=dim(mydata)[1]
        names(mydata)
        
        #if not available variable
        if(sum(names(mydata)==varname)!=1) {
            print (paste("Check the name of variable ", varname, sep=""))
            stop
        }#if
        
        ### single and mva lag
        prefix=c(paste(varname, "_s", 1:laglength, sep=""), #single
        paste(varname, "_m", 1:laglength, sep=""))
        
        var_lag<-array(,c(nrr,2*laglength),dimnames = list(1:nrr,prefix))
        #dim(var_lag); names(var_lag)
        
        lag_ftn<-function(v, k){
          n=length(v)
          if(n>k){
            v_lag<-v
            v_lag[1:k]<-NA
            v_lag[(k+1):n]<-v[1:(n-k)]
          }
          return(v_lag)
        }
        
        ### single lag
        for (i in 1:laglength)
        {
            var_lag[,i]=lag_ftn(mydata[, varname],k=i)
        }#i
        
        ### mva
        sum_lag<-array(,dim=nrr)
        for (i in (laglength+1):(2*laglength))
        {
            i2=i-laglength
            sum_lag=0
            for (k2 in 0:i2)
            {
                sum_lag=sum_lag + lag_ftn(mydata[, varname],k=k2)
            }#k2
            var_lag[,i]=sum_lag/(i2+1)
        }#i
        
        mydata=cbind(mydata, var_lag)
    } #iv
    
    return(mydata)

}
