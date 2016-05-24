############## Generate results as graphs #############
####  set results directory first  ####
dpa.results.setGraphDir <- function(graphDir=NULL)
{
 if(is.null(graphDir)){
   graphDir <-tclvalue(tkchooseDirectory())
 }
 assign("graphDir",graphDir,env=.GlobalEnv)
}

################ Node plot ######################
dpa.results.viewNodePlots <-function()
{
plot <- tktoplevel()
Font1<-tkfont.create(family="times",size=10)
Font2<-tkfont.create(family="times",size=12)
tktitle(plot) <- "  Nodes Plot  "

tkgrid(tklabel(plot,text="x-axis parameter",font=Font2),row=0,column=0,padx=5,pady=5)
source <- names(e)
xColumn <- tkwidget(plot,"ComboBox",editable=FALSE,values=source)
tkconfigure(xColumn)
tkgrid(xColumn,row=1,column=0,padx=5,pady=5)

tkgrid(tklabel(plot,text="y-axis parameter",font=Font2),row=0,column=1,padx=5,pady=5)
source <- names(e)
yColumn <- tkwidget(plot,"ComboBox",editable=FALSE,values=source)
tkconfigure(yColumn)
tkgrid(yColumn,row=1,column=1,padx=5,pady=5)

tkgrid(tklabel(plot,text="tick number",font=Font2),row=0,column=4,padx=5,pady=5)
ticks <- c(1:max(e["tick"]))
tickNumber <- tkwidget(plot,"ComboBox",editable=TRUE,values=ticks)
tkconfigure(tickNumber)
tkgrid(tickNumber,row=1,column=4,padx=5,pady=5)

dpa.results.plotNodes<- function()
{ 
  dpa.incrementValue(i)
  x_column <- source[as.numeric(tclvalue(tcl(xColumn,"getvalue")))+1]
  #print(x_column)
  tkgrid(tklabel(plot,text=paste(x_column),font=Font1),row=i+1,column=0,padx=5,pady=5)
  y_column <- source[as.numeric(tclvalue(tcl(yColumn,"getvalue")))+1]
  #print(y_column)
  tkgrid(tklabel(plot,text=paste(y_column),font=Font1),row=i+1,column=1,padx=5,pady=5)
 
###########  scatter plot function  #################
 dpa.plot.scatter <- function(dataframe=NULL,xcolumn=NULL,ycolumn=NULL)
 {
   attach(dataframe)
   on.exit(detach(dataframe))
   colors=densCols(eval(parse(text=xcolumn)),eval(parse(text=ycolumn)))
   #graphDir <- "C:/Users/amit/Desktop/dynamic_path_approach/exercise_R/graphs"
   png(file=paste(graphDir,"/scatter_",xcolumn,"_",ycolumn,"_",".png",sep=""),res=300,width=3000,height=2000)
   smoothScatter(
	eval(parse(text=xcolumn)),
	eval(parse(text=ycolumn)),
	xlab=xcolumn,
      ylab=ycolumn,
	pch=20,
      col=colors,
	nrpoints=Inf
    )
  devs<-dev.off()
 }
dpa.plot.scatter(e,x_column,y_column)
}
nodePlot.but <-tkbutton(plot,text="   plot irrespective of time   ",command=dpa.results.plotNodes)
tkgrid(nodePlot.but,row=1,column=2,padx=5,pady=0) 

######################  plotting graphs for a specific tick entered by user  ###########################
dpa.results.plotNodesWithTime <-function()
{
  dpa.incrementValue(i)
  x_column <- source[as.numeric(tclvalue(tcl(xColumn,"getvalue")))+1]
  #print(x_column)
  tkgrid(tklabel(plot,text=paste(x_column),font=Font1),row=i+1,column=0,padx=5,pady=5)
  y_column <- source[as.numeric(tclvalue(tcl(yColumn,"getvalue")))+1]
  #print(y_column)
  tkgrid(tklabel(plot,text=paste(y_column),font=Font1),row=i+1,column=1,padx=5,pady=5)
 
  tick_number <- ticks[as.numeric(tclvalue(tcl(tickNumber,"getvalue")))+1]
  #print(tick_number)
  tkgrid(tklabel(plot,text=paste(tick_number),font=Font1),row=i+1,column=4,padx=5,pady=5)

  newData <-NULL
  newData <- rbind(newData,e[which(e["tick"]==tick_number),])

#scatter plot function 1
 dpa.plot.scatter <- function(dataframe=NULL,xcolumn=NULL,ycolumn=NULL,tick=NULL)
 {
   attach(dataframe)
   on.exit(detach(dataframe))
   colors=densCols(eval(parse(text=xcolumn)),eval(parse(text=ycolumn)))
   #graphDir <- "C:/Users/amit/Desktop/dynamic_path_approach/exercise_R/graphs"
   png(file=paste(graphDir,"/scatter_",xcolumn,"_",ycolumn,"_",tick,".png",sep=""),res=300,width=3000,height=2000)
   smoothScatter(
	eval(parse(text=xcolumn)),
	eval(parse(text=ycolumn)),
	xlab=xcolumn,
      ylab=ycolumn,
	pch=20,
      col=colors,
	nrpoints=Inf
    )
  devs<-dev.off()
 }
 dpa.plot.scatter(newData,x_column,y_column,tick_number)
}
nodePlotWithTime.but <-tkbutton(plot,text="   plot for tick number   ",command=dpa.results.plotNodesWithTime)
tkgrid(nodePlotWithTime.but,row=1,column=3,padx=5,pady=0) 
}

###############################  Relations plot  ################################
dpa.results.viewRelationsPlots <-function(tickNumber=NULL)
{

dpa.plot.graph <- function(dataframe=NULL,parameters,relations,tickNumber=NULL,extension="pdf")
{
 attach(dataframe)
 on.exit(detach(dataframe))
 if (rbVal=="time_irrespective"){
   if(extension=="png"){
     png(file=paste(graphDir,"/graph_","relation",".png",sep=""),res=92,width=1000,height=1000)
   }else{
     pdf(file=paste(graphDir,"/graph_","relation",".pdf",sep=""),width=6,height=6)
   }
 }else{
   if(extension=="png"){
     png(file=paste(graphDir,"/graph_","relation_",format(tickNumber,width=4),".png",sep=""),res=92,width=1000,height=1000)
   }else{
     pdf(file=paste(graphDir,"/graph_","relation_",tickNumber,".pdf",sep=""),width=6,height=6)
   }
 }

  strength<-NULL
  params=data.frame()
  for(i in 1:nrow(row)){
  
    #define from and to
    from<-unlist(strsplit(row[i,1]," "))[1]
    type<-unlist(strsplit(row[i,1]," "))[2]
    to<-unlist(strsplit(row[i,1]," "))[3]
    params[i,1]=from
    params[i,2]=to
  }

  relations<-params
  strengths<-NULL
  
  #override with standardized coefficients
  strengths<-sem.standardized[,2]
  parameters <- t(t(sort(parameters)))
  #print(parameters)
  #print(relations)
  graph1 <- graph.data.frame(relations, directed=TRUE, vertices=parameters)
  #lay1 <- layout.circle(graph1,parameters)
  #lay1 <- layout.sphere(graph1)
  #lay1 <- layout.fruchterman.reingold(graph1,niter=3000,verbose=igraph.par("verbose"))
  lay1 <- layout.reingold.tilford(graph1,parameters)
  #lay1 <- layout.lgl(graph1)
  title=""
  if(!is.null(tickNumber)){
     title<-paste("Tick: ",tickNumber)
  }

  #dashed lines for negative relations, solid for positive ones.
  typeofArrow<-NULL
  typeofArrow<-(strengths>0)
  typeofArrow<-match(typeofArrow,TRUE,FALSE)
  typeofArrow<-replace(typeofArrow,typeofArrow<1,2)
  typeofArrow<-replace(typeofArrow,strengths==1,0)

  maximum=max(abs(strengths))
  if(maximum<1){maximum=1}

  if(extension=="pdf"){
    pdffactor=.7
  } else {
    pdffactor=1.9
  }

  #replace strengths of exactly one to 0
  labels<-round(strengths,digits=3) 
  labels<-replace(labels,strengths==1,"")
  #print(labels)
  strengths<-replace(strengths,strengths==1,0.00000001)
  #print(strengths)
  #print("TypeofArrow:")
  #print(typeofArrow)

  plot.igraph(graph1,axes=FALSE,layout=lay1,
	#Titles and properties
	main=title,
	xlab="",ylab="",

	#Edge propeties:
	edge.label=labels,
	edge.label.cex=1.2*pdffactor, 
      edge.label.family="sans",
	edge.color="black",
	edge.label.color="black",
   	edge.arrow.size=1*pdffactor,
      edge.width=7*abs(strengths),
      edge.lty=typeofArrow,

	#Vertex propeties:
	vertex.label=parameters,
	vertex.label.dist=1, 
	vertex.label.cex=1.2*pdffactor, 
	vertex.size=3*pdffactor,
      vertex.color="black",
      vertex.label.family="sans",
	vertex.label.color="black"

  )
   devs<-dev.off()
}

parameters <- variables
relationsNew<-data.frame(rbind(relations[,1],relations[,2]),strength=c(10,5))
dpa.plot.graph(e,parameters,relations,tickNumber,'png')
dpa.plot.graph(e,parameters,relations,tickNumber,'pdf')
}


###############################  Coefficients plot  ################################
dpa.results.generateCoefficientsPlots <-function(filename=NULL,colors=NULL,indices=NULL,legend=NULL){

dpa.plot.graph<-function(extension=NULL,filename=NULL,colors=NULL,indices=NULL,legend=NULL){

nrofrows=0
for(i in 1:nrow(sem.results.coefficients)){
  if(min(sem.results.coefficients[i,]==1)&max(sem.results.coefficients)==1) {} 
  else{
    nrofrows=nrofrows+1
  }
}

sem.results.stripped.coefficients<-matrix(0,ncol=ncol(sem.results.coefficients),nrow=nrofrows)
sem.results.stripped.parameters<-matrix(0,ncol=1,nrow=nrofrows)
rows=0
for(i in 1:nrow(sem.results.coefficients)){
  if(min(sem.results.coefficients[i,]==1)&max(sem.results.coefficients)==1) {} 
  else{
    rows=rows+1
    sem.results.stripped.coefficients[rows,]<-sem.results.coefficients[i,]
    sem.results.stripped.parameters[rows,1]<- as.character(sem.results.parameters[rows])
  }
}

#make an additional selection
if(is.null(indices)){
  numbers=sem.results.stripped.coefficients
  labels=sem.results.stripped.parameters[,1]
}else{
  numbers=sem.results.stripped.coefficients[indices,]
  labels=sem.results.stripped.parameters[indices,1]
}

if(extension=='png'){
  fileinlcudingpath=paste(graphDir,"/",filename,".",extension,sep="")
  png(fileinlcudingpath,res=300,width=3000,height=2000)
}
if(extension=='pdf'){
  fileinlcudingpath=paste(graphDir,"/",filename,".",extension,sep="")
  pdf(fileinlcudingpath,width=5,height=5)
}

#plot.new()
matplot(
  type="l",
  as.numeric(listOfTicks),
  t(numbers),
  xlab="time",
  ylab="standardized coefficient",
  col=colors,
  lwd=2.5,
  lty=c(1:6),
  bty='n',
  ylim=c(-1,1),
  xlim=c(min(e$tick),max(e$tick))
)
legend(legend,t(labels),cex=1,col=colors,lwd=2.5,bty='n',lty=c(1:6));
devs<-dev.off()
}


  if(is.null(filename)){
    filename="coefficients"
  }
  if(is.null(colors)){
    colors=c('#FE9800','#FD97CA','#98CBFE','#BFBFBF','#FECB98','#FEFE98','#98CB00','#DEDEDE')
  }
  if(is.null(legend)){
    legend="bottomright"
  }
  extension='pdf'
  dpa.plot.graph(extension,filename,colors,indices,legend)
  extension='png'
  dpa.plot.graph(extension,filename,colors,indices,legend)


}

##############################Generate Fit Plots ###########################################
dpa.results.generateFitPlots <-function(filename=NULL,colors=NULL,indices=NULL,legend=NULL){

dpa.plot.graph<-function(extension=NULL,filename=NULL,colors=NULL,indices=NULL,legend=NULL){

labels=c("Iterations","Deg. freedom","GFI","Adj. GFI","RMSEA","SRMR","NFI","NNFI")

#make an additional selection
if(is.null(indices)){
  numbers=sem.results.statistics
}else{
  numbers=sem.results.statistics[indices,]
  labels=labels[indices]
}

if(extension=='png'){
  fileinlcudingpath=paste(graphDir,"/",filename,".",extension,sep="")
  png(fileinlcudingpath,res=300,width=3000,height=2000)
}
if(extension=='pdf'){
  fileinlcudingpath=paste(graphDir,"/",filename,".",extension,sep="")
  pdf(fileinlcudingpath,width=5,height=5)
}

#bound the y axis if appropriate
ylimits=c(min(numbers),max(numbers))
if(max(numbers)<=1){
  ylimits=c(ylimits[1],1)
}
if(min(numbers)>=1){
  ylimits=c(0,ylimits[2])
}

#plot.new()
matplot(
  type="l",
  as.numeric(listOfTicks),
  t(numbers),
  xlab="time",
  ylab="statistic",
  col=colors,
  lwd=2.5,
  lty=c(1:6),
  bty='n',
  ylim=ylimits,
  xlim=c(min(e$tick),max(e$tick))
)
legend(legend,t(labels),cex=1,col=colors,lwd=2.5,bty='n',lty=c(1:6));
devs<-dev.off()
}


  if(is.null(filename)){
    filename="statistics"
  }
  if(is.null(colors)){
    colors=c('#FE9800','#FD97CA','#98CBFE','#BFBFBF','#FECB98','#FEFE98','#98CB00','#DEDEDE')
  }
  if(is.null(legend)){
    legend="bottomright"
  }
  extension='pdf'
  dpa.plot.graph(extension,filename,colors,indices,legend)
  extension='png'
  dpa.plot.graph(extension,filename,colors,indices,legend)

}
