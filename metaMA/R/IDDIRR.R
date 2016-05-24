`IDDIRR` <-
function(finalde,deindst)
{
DE=length(finalde)
gains=finalde[which(!(finalde %in% deindst))]
IDD=length(gains)
IDR=IDD/DE*100
perte=which(!(deindst %in% finalde))
Loss=length(perte)
IRR=Loss/length(deindst)*100
res=c(DE,IDD,Loss,round(IDR,2),round(IRR,2))
names(res)=c("DE","IDD","Loss","IDR","IRR")
res
}

