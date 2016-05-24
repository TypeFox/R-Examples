deg2dms <-
function(deg,type='mat',sep=':',digits=2){
if(any(deg< -90 | deg>90)){stop('All deg values should be -90<=deg<=90')}
    temp = sign(deg)
    deg = abs(deg)
    DEG = floor(deg)
    MIN = floor((deg - DEG) * 60)
    SEC = (deg - DEG - MIN/60) * 3600
    SEC[SEC<0]=0; SEC[SEC>60]=60
    MIN[SEC == 60] = MIN[SEC == 60] + 1
    SEC[SEC == 60] = 0
    DEG[MIN == 60] = DEG[MIN == 60] + 1
    MIN[MIN == 60] = 0
    if(digits<0){width=9}
    if(digits==0){width=2}
    if(digits>0){width=digits+3}
    SEC = formatC(SEC, format = "f", width = width, flag = 0, digits = digits)
    MIN = formatC(MIN, format = "f", width = 2, flag = 0, digits = 0)
    DEG = formatC(DEG, format = "f", width = 2, flag = 0, digits = 0)
    DEG[temp == -1] = paste("-", DEG[temp == -1], sep = "")
    DEG[temp == 1] = paste("+", DEG[temp == 1], sep = "")
    DEG[temp == 0] = paste("+", DEG[temp == 0], sep = "")
    if(type=='mat'){output = cbind(DEG, MIN, SEC)}
    if(type=='cat' & sep!='DMS' & sep!='dms'){output=apply(cbind(DEG, MIN, SEC),1,paste,collapse=sep)}
    if(type=='cat' & sep=='DMS'){output=paste(paste(paste(DEG,MIN,sep='D'),SEC,sep='M'),c('','',''),sep='S')}
    if(type=='cat' & sep=='dms'){output=paste(paste(paste(DEG,MIN,sep='d'),SEC,sep='m'),c('','',''),sep='s')}
return(output)
}
