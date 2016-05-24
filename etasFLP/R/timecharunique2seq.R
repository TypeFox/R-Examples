timecharunique2seq <-function(timestring){
        timesecv1.v2= as.double(difftime(timestring,as.POSIXct("1900-01-01",tz="GMT"),tz="GMT",units="secs"))
        t1  =timestring
        t3  =as.double(timesecv1.v2)
        t2  =trunc(timesecv1.v2/86400)
        return(list(char=t1,sec=t3,day=t2))
}
