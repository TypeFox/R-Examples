swat_objective_function<-function(x,calib_range,calib_params,flowgage){

calib_params$current<-x
file.remove(list.files(pattern="output."))
alter_files(calib_params)

libarch = if (nzchar(version$arch)) paste("libs", version$arch, sep = "/") else "libs"
swatbin <- "rswat2012.exe"
system(shQuote(paste(path.package("SWATmodel"), libarch, swatbin, sep = "/")))

start_year=read.fortran(textConnection(readLines("file.cio")[9]),"f20")
test=readLines(file("output.rch"))
rchcolname=sub(" ","",(substr(test[9],50,61)))
flow=data.frame(as.numeric(as.character(substr(test[10:length(test)],50,61))))
colnames(flow)=rchcolname
reach=data.frame(as.numeric(as.character(substr(test[10:length(test)],8,10))))
rchcolname=sub(" ","",(substr(test[9],8,10)))
colnames(reach)=rchcolname
outdata=cbind(reach,flow)
test2=subset(outdata,outdata$RCH==3)
test2$mdate=as.Date(row(test2)[,1],origin=paste(start_year-1,"-12-31",sep=""))
test3=merge(flowgage$flowdata,test2[c(2,length(test2))],all=F)
test3$Qm3ps=test3$flow/3600/24
NS = NSeff(test3$Qm3ps[730:(length(test3$Qm3ps))], test3$FLOW_OUTcms[730:(length(test3$Qm3ps))])
plot(test3$mdate, test3$Qm3ps,type="l")
lines(test3$mdate, test3$FLOW_OUTcms,col="red")
print(NS)
return(abs(NS-1))
}


