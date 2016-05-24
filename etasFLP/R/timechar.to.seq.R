timechar.to.seq <-function(v1,v2,v3){
        i1  =as.character(v1)
        j1  =as.character(v2)
        norm    =paste(as.character(((as.numeric(substr(i1,1,2))<10)+19)),i1,sep="")
        norm    =paste(norm,j1,sep="")
        dt  =strptime(norm,format="%Y%m%d%H%M",tz="GMT")
        timesecv1.v2= as.double(difftime(dt,as.POSIXct("1900-01-01",tz="GMT"),tz="GMT",units="secs"))
        t1  =norm
        t3  =as.double(timesecv1.v2+v3)
        t2  =trunc(timesecv1.v2/86400)
        return(list(char=t1,sec=t3,day=t2))
}
