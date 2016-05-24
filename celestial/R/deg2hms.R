deg2hms <-
function(deg,type='mat',sep=':',digits=2){
if(any(deg< 0 | deg>360)){stop('All deg values should be 0<=d<=360')}
    deg[deg < 0] = deg[deg < 0] + 360
    HRS = floor(deg/15)
    MIN = floor((deg/15 - HRS) * 60)
    SEC = (deg/15 - HRS - MIN/60) * 3600
    SEC[SEC<0]=0; SEC[SEC>60]=60
    MIN[SEC == 60] = MIN[SEC == 60] + 1
    SEC[SEC == 60] = 0
    HRS[MIN == 60] = HRS[MIN == 60] + 1
    MIN[MIN == 60] = 0
    HRS = HRS%%24
    if(digits<0){width=9}
    if(digits==0){width=2}
    if(digits>0){width=digits+3}
    SEC = formatC(SEC, format = "f", width = width, flag = 0, digits = digits)
    MIN = formatC(MIN, format = "f", width = 2, flag = 0, digits = 0)
    HRS = formatC(HRS, format = "f", width = 2, flag = 0, digits = 0)
    if(type=='mat'){output = cbind(HRS, MIN, SEC)}
    if(type=='cat' & sep!='HMS' & sep!='hms'){output=apply(cbind(HRS, MIN, SEC),1,paste,collapse=sep)}
    if(type=='cat' & sep=='HMS'){output=paste(paste(paste(HRS,MIN,sep='H'),SEC,sep='M'),c('','',''),sep='S')}
    if(type=='cat' & sep=='hms'){output=paste(paste(paste(HRS,MIN,sep='h'),SEC,sep='m'),c('','',''),sep='s')}
return(output)
}
