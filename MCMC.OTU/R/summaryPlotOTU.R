summaryPlotOTU <-
function(data,xgroup,facet=NA,type="bar",x.order=NA,whiskers="ci",otus=NA,log.base=10){
#	smry=smm0$summary;xgroup="bank";type="line";whiskers="ci";facet="species";otus=NA;x.order=NA;log.base=10
	smry=data
	if(!is.na(otus[1])){ smry=smry[smry$otu %in% otus,] }
	if (!is.na(x.order[1])) { smry[,xgroup]=factor(smry[,xgroup],levels=x.order) }
	if (dimnames(table(smry$mean>0))[[1]]=="FALSE" && table(smry$mean>0)[[1]]==length(smry$mean)) { 
		ylabb=paste("log",log.base,"(proportion)",sep="")
	} else {
		ylabb=paste("log",log.base,"(abundance)",sep="")
	}
	if (is.na(facet)) {
		if (type=="line" & whiskers=="ci") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",group="otu",colour="otu"))+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),lwd=0.4,
				width=0.7,position=pd)+
				geom_line(aes_string(group="otu"),position=pd)+
				geom_point(aes_string(group="otu"),position=pd,size=2.5)+
				theme_bw()+xlab(xgroup)+ylab(ylabb)
		}
		if (type=="bar" & whiskers=="ci") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="otu"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),lwd=0.4,width=0.5,
				colour="grey40",position=pd)+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
		if (type=="line" & whiskers=="sd") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",group="otu",colour="otu"))+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
				lwd=0.4,width=0.7,position=pd)+
				geom_line(aes_string(group="otu"),position=pd)+
				geom_point(aes_string(group="otu"),position=pd,size=2.5)+
				theme_bw()+xlab(xgroup)+ylab(ylabb)
		}
		if (type=="bar" & whiskers=="sd") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="otu"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
				lwd=0.4,width=0.5,colour="grey40",position=pd)+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
	} else {
		if (type=="line" & whiskers=="ci") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",colour="otu"))+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),
				lwd=0.4,width=0.7,position=pd)+
				geom_line(aes_string(group="otu"),position=pd)+
				geom_point(aes_string(group="otu"),position=pd,size=2.5)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(ylabb)
		}
		if (type=="bar" & whiskers=="ci") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="otu"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),
				lwd=0.4,width=0.5,colour="grey40",position=pd)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
		if (type=="line" & whiskers=="sd") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",colour="otu"))+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
				lwd=0.4,width=0.7,position=pd)+
				geom_line(aes_string(group="otu"),position=pd)+
				geom_point(aes_string(group="otu"),position=pd,size=2.5)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(ylabb)
		}
		if (type=="bar" & whiskers=="sd") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="otu"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),lwd=0.4,width=0.5,colour="grey40",position=pd)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
	}
	print(gpp)
	return(gpp)
}
