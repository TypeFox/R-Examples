IAUID <-
function(ra,dec,name='GAMA',epoch='J'){

HMSiau=function(deg){
deg[deg<0]=deg[deg<0]+360
HRS=floor(deg/15)
MIN=floor((deg/15-HRS)*60)
SEC=(deg/15-HRS-MIN/60)*3600
SEC=floor(SEC*100*1.00001)/100
MIN[SEC==60]=MIN[SEC==60]+1
SEC[SEC==60]=0
HRS[MIN==60]=HRS[MIN==60]+1
MIN[MIN==60]=0
HRS=HRS%%24

SEC=formatC(SEC, format="f", width=5, flag=0, digits=2)
MIN=formatC(MIN, format="f", width=2, flag=0, digits=0)
HRS=formatC(HRS, format="f", width=2, flag=0, digits=0)

output=cbind(HRS,MIN,SEC)
}

DMSiau=function(deg){
temp=sign(deg)
deg=abs(deg)
DEG=floor(deg)
MIN=floor((deg-DEG)*60)
SEC=(deg-DEG-MIN/60)*3600
SEC=floor(SEC*10*1.00001)/10
MIN[SEC==60]=MIN[SEC==60]+1
SEC[SEC==60]=0
DEG[MIN==60]=DEG[MIN==60]+1
MIN[MIN==60]=0

SEC=formatC(SEC, format="f", width=4, flag=0, digits=1)
MIN=formatC(MIN, format="f", width=2, flag=0, digits=0)
DEG=formatC(DEG, format="f", width=2, flag=0, digits=0)
DEG[temp==-1]=paste('-',DEG[temp==-1],sep='')
DEG[temp==1]=paste('+',DEG[temp==1],sep='')
DEG[temp==0]=paste('+',DEG[temp==0],sep='')
output=cbind(DEG,MIN,SEC)
}

temphms=HMSiau(ra)
tempdms=DMSiau(dec)

IAUname=paste(name,epoch, temphms[,1], temphms[,2], temphms[,3], tempdms[,1], tempdms[,2], tempdms[,3],sep='')

return(IAUname)
}
