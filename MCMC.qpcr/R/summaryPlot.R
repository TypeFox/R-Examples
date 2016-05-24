summaryPlot <-
function(data,xgroup,facet=NA,type="bar",x.order=NA,whiskers="ci",genes=NA,log.base=2){
#	smry=su01;xgroup="condition";type="line";whiskers="ci";facet="timepoint";genes=NA;x.order=NA
	smry=data
	if(!is.na(genes[1])){ smry=smry[smry$gene %in% genes,] }
	if (!is.na(x.order[1])) { smry[,xgroup]=factor(smry[,xgroup],levels=x.order) }
	if (is.na(facet)) {
		if (type=="line" & whiskers=="ci") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",group="gene",colour="gene"))+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),lwd=0.4,
				width=0.7,position=pd)+
				geom_line(aes_string(group="gene"),position=pd)+
				geom_point(aes_string(group="gene"),position=pd,size=2.5)+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(abundance)",sep=""))
		}
		if (type=="bar" & whiskers=="ci") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="gene"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),lwd=0.4,width=0.5,
				colour="grey40",position=pd)+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
		if (type=="line" & whiskers=="sd") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",group="gene",colour="gene"))+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
				lwd=0.4,width=0.7,position=pd)+
				geom_line(aes_string(group="gene"),position=pd)+
				geom_point(aes_string(group="gene"),position=pd,size=2.5)+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(abundance)",sep=""))
		}
		if (type=="bar" & whiskers=="sd") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="gene"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
				lwd=0.4,width=0.5,colour="grey40",position=pd)+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
	} else {
		if (type=="line" & whiskers=="ci") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",colour="gene"))+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),
				lwd=0.4,width=0.7,position=pd)+
				geom_line(aes_string(group="gene"),position=pd)+
				geom_point(aes_string(group="gene"),position=pd,size=2.5)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(abundance)",sep=""))
		}
		if (type=="bar" & whiskers=="ci") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="gene"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes_string(ymin="lower",ymax="upper"),
				lwd=0.4,width=0.5,colour="grey40",position=pd)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
		if (type=="line" & whiskers=="sd") {
			pd=position_dodge(0.3)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",colour="gene"))+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
				lwd=0.4,width=0.7,position=pd)+
				geom_line(aes_string(group="gene"),position=pd)+
				geom_point(aes_string(group="gene"),position=pd,size=2.5)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(abundance)",sep=""))
		}
		if (type=="bar" & whiskers=="sd") {
			pd=position_dodge(0.85)
			gpp=ggplot(smry,aes_string(x=xgroup,y="mean",fill="gene"))+
				geom_bar(stat="identity",position=pd)+
				geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),lwd=0.4,width=0.5,colour="grey40",position=pd)+
				facet_wrap(as.formula(paste("~", facet)))+
				theme_bw()+xlab(xgroup)+ylab(paste("log",log.base,"(fold change)",sep=""))
		}
	}
	print(gpp)
	return(gpp)
}
