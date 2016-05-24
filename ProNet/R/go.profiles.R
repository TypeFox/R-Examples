##' @title GO enrichment comparison
##' 
##' @description GO enrichment comparison between networks.
##' 
##' @param g_ref Vector of the first set of molecules. 
##' @param g_comp	Vector of the second set of molecules. 
##' @param color Color represent of the two sets.
##' @param onto GO categories, three possible values are \code{MF} for GO function, \code{BP} for GO process, and \code{CC} for GO componet
##' @param mode	Mode of the GO vlues, either \code{Frequency} or \code{Number}.
##' @param demonstration.number	Number of the top GOs to display.
##' @param plot Logical value, whether to plot the profiling result (if \code{TRUE}) or not (if \code{FALSE}). 
##' @param ...	other arguments.
##' @return	A data frame of the profiling result and plots.
##' @export
##' @examples
##' entrez1<-c("11067","414157","196477","147339","642")
##' entrez2<-c("121549","51160","83878","11338","196477","9319","608","7015")
##' go.profiles(entrez1,onto="MF",main="Only Network 1")
##' go.profiles(entrez1,entrez2,onto="MF",main=c("Network 1 vs 2"))
go.profiles<-function(g_ref,g_comp,color=c("red","white"),onto=c("MF","BP","CC"),
                      mode=c("number","frequency"),demonstration.number=10,plot=TRUE,...)
{	
	if(missing(g_ref)){
		stop("Missing of g_ref")
	}
	mode=match.arg(mode)
	percentage<-ifelse(mode=="number",FALSE,TRUE)
	res<-NULL
	if(onto=="MF"){
    GO_function<-NULL
		data("GO_function", envir = environment())
		if(missing(g_comp))
		  res<-go.profile.table1(g_ref,GOset=GO_function[,c("GeneID","GO_ID","GO_term")])
		else
		  res<-go.profile.table2(g_ref,g_comp,GOset=GO_function[,c("GeneID","GO_ID","GO_term")])
	}else if(onto=="BP"){
    GO_process<-NULL
		data("GO_process", envir = environment())
		if(missing(g_comp))
		  res<-go.profile.table1(g_ref,GOset=GO_process[,c("GeneID","GO_ID","GO_term")])
		else
		  res<-go.profile.table2(g_ref,g_comp,GOset=GO_process[,c("GeneID","GO_ID","GO_term")])
	}else if(onto=="CC"){
    GO_component<-NULL
		data("GO_component", envir = environment())
		if(missing(g_comp))
		  res<-go.profile.table1(g_ref,GOset=GO_component[,c("GeneID","GO_ID","GO_term")])
		else
		  res<-go.profile.table2(g_ref,g_comp,GOset=GO_component[,c("GeneID","GO_ID","GO_term")])
	}
	if(demonstration.number<nrow(res))
    res<-res[order(res[,3],decreasing=TRUE)[1:demonstration.number],]
	else
    res<- res[order(res[,3],decreasing=TRUE),]
	rownames(res)<-NULL
	if(mode=="frequency"){
		if(ncol(res)==3)
			res[,3]<-round(res[,3]/(sum(res[,3])+0.00000001)*100,2)
		else if(ncol(res)==4){
			res[,3]<-round(res[,3]/(sum(res[,3])+0.00000001)*100,2)
			res[,4]<-round(res[,4]/(sum(res[,4])+0.00000001)*100,2)
		}
	}
	if(plot){
		if(!is.null(res)){
			if(ncol(res)==3){
				go.profile.plot1(res,mode=mode,...)
			}else if(ncol(res)==4){
				go.profile.plot2(res,mode=mode,legend.name=c(substitute(g_ref),substitute(g_comp)),...)
			}
		}
	}
	if(ncol(res)==3)
    colnames(res)<-c("GO_term","GOID",substitute(g_ref))
	else if(ncol(res)==4)
    colnames(res)<-c("GO_term","GOID",substitute(g_ref),substitute(g_comp))

  return(res)
}

## GO profiling for a single set of molecules.
go.profile.table1<-function(gene,GOset)
{
	gene2go<-GOset[!is.na(match(GOset$GeneID,gene)),]
	GOids<-table(gene2go$GO_ID)
	GOids<-data.frame(GOids)
	GOids[,1]<-as.character(GOids[,1])
	GOids<-cbind(GOset[match(GOids[,1],GOset$GO_ID),"GO_term"],GOids)
	colnames(GOids)<-c("GO_term","GOID","Freq")
	GOids$GO_term<-as.character(GOids$GO_term)
	return(GOids)
}

## GO profiling for two sets of molecules.
go.profile.table2<-function(gene1,gene2,GOset)
{
	gene2go1<-GOset[!is.na(match(GOset$GeneID,gene1)),]
	gene2go2<-GOset[!is.na(match(GOset$GeneID,gene2)),]
	GOids<-table(c(gene2go1$GO_ID,gene2go2$GO_ID))
	freq1<-rep(0,length(GOids))
	for(i in seq_along(freq1)){
		freq1[i]<-sum(gene2go1$GO_ID==names(GOids)[i])
	}
	GOids<-GOids-freq1
	GOids<-as.data.frame(cbind(freq1,GOids))
	GOids<-cbind(GOset[match(rownames(GOids),GOset$GO_ID),"GO_term"],rownames(GOids),GOids)
	colnames(GOids)<-c("GO_term","GOID","Freq1","Freq2")
	GOids$GO_term<-as.character(GOids$GO_term)
	GOids$GOID<-as.character(GOids$GOID)
	return(GOids)	
}

## Plot the first GO profiling result.
go.profile.plot1<-function(x,mode,...)
{
	chars<-paste(x$GOID,x$GO_term,sep=":")
	nchars<-max(unlist(lapply(chars,function(item){nchar(item)})))
	opt<-par(mar=c(4,31*0.38,4,4),xpd=TRUE, cex.axis=0.01)
	y<-x$Freq[order(x$Freq)]
	rm(x)
	names(y)<-charseg(chars)
	if(mode=="frequency") 
    xlim<-c(0,100) 
  else 
    xlim<-c(0,max(y)+min(3,max(y)*0.4))
	bp<-barplot(y,horiz=TRUE,beside=T,col=cm.colors(length(y)),
              xlim=xlim,names.arg="",xaxt="n",...)
	text(y,round(bp,1),y,pos=4,cex=0.8)
	axis(2,at=bp,labels=rev(names(y)),cex.axis=0.8,las=2)
	if(mode=="frequency")
		axis(1,cex.axis=0.8,labels=seq(0,100,by=20),at=seq(0,100,by=20))
	else
    axis(1,cex.axis=0.8,labels=seq(0,xlim[2],length=5),at=seq(0,xlim[2],length=5))
	par(opt)
}

## Plot the second GO profiling result.
go.profile.plot2<-function(x,mode,color=c("red","white"),legend.name,...)
{
	chars<- paste(x$GOID,x$GO_term,sep=":")
	nchars<-max(unlist(lapply(chars,function(item){nchar(item)})))
	y<-x[,c("Freq2","Freq1")]
	y<-y[order(y$Freq1),]
	rownames(y)<-charseg(chars)
	y<-as.matrix(t(y))
	opt<-par(mar=c(4,31*0.35,4,4),fig=c(0,0.9,0,1),cex.axis=0.01)
	if(mode=="frequency")
    xlim<-c(0,100) 
  else 
    xlim<-c(0,round(max(y)+min(3,max(y)*0.4)))
	bp<-barplot(y,horiz=TRUE,beside=T,col=color,xlim=xlim,
              names.arg="",axisnames=F,xaxt="n",...)
	text(y-xlim[2]*0.015,round(bp,1),y,pos=4,cex=0.7)
	axis(2,at=(bp[1,]+bp[2,])/2,labels=rev(colnames(y)),cex.axis=0.7,las=2)
	if(mode!="frequency")
    axis(1,cex.axis=0.8,labels=round(seq(xlim[1],xlim[2],length=5),0),at=round(seq(xlim[1],xlim[2],length=5),0))
	else
    axis(1,cex.axis=0.8,labels=seq(0,100,by=20),at=seq(0,100,by=20))
	legend("bottomright",legend=legend.name,fill=rev(color),cex=0.8,pt.cex=0.7)
	par(opt)
}

charseg<-function(chars)
{
  unlist(lapply(chars,function(name){
				len<-nchar(name)
				if(len>121)
					paste(substring(name,c(1,31,61,91,121),c(30,60,90,120,1000000L)),collapse="\n")
				else if(len>91)
					paste(substring(name,c(1,31,61,91),c(30,60,90,120)),collapse="\n")
				else if(len>61)
					paste(substring(name,c(1,31,61),c(30,60,90)),collapse="\n")
				else if(len>31)
					paste(substring(name,c(1,31),c(30,60)),collapse="\n")
				else
          name
	}))
}
