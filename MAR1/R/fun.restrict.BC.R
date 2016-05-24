fun.restrict.BC<-function(namesvar,namescovar){

taxa<-namesvar
#taxa<-c("krill","chaetognath","copepod","tunicate","phyto")
drivers<-namescovar
#drivers<-c("fish","fish2")

###################################################

res<-rep(0,length(taxa)^2+length(taxa)*length(drivers))

tt<-tktoplevel()
#tkwm.iconbitmap(tt,"mar1.ico")
tkwm.title(tt,"Specify Restrictions")
frame.overall<-tkframe(tt)
frame.var<-tkframe(frame.overall,relief="groove",borderwidth=3)
frame.covar<-tkframe(frame.overall,relief="groove",borderwidth=3)
Sp<-tkfont.create(size=1)
Spp<-tkfont.create(size=4)
Ph<-tkfont.create(size=12) # place holder for covariate buttons
f.tit<-tkfont.create(size=10,weight="bold") # title font
change<-tclVar(0)
changeVal<-as.numeric(tclvalue(change))

#### create buttons ####

for(c in 1:length(taxa)^2){

	func<-paste("function(){tclvalue(change)<-",c,"}",sep="")
	assign(paste("func",c,sep="."),eval(parse(text=func)))

assign( paste("var","b",c,sep="."),tkbutton(frame.var,text="  0.5  ",bg="gray80") )
tkconfigure( get(paste("var","b",c,sep=".")),command=get(paste("func",c,sep=".")) )#FUNCTION!!!
}

if(length(drivers)>0){
for(c in 1:(length(taxa)*length(drivers))){

	func<-paste("function(){tclvalue(change)<-",c+length(taxa)^2,"}",sep="")
	assign(paste("func",c,sep="."),eval(parse(text=func)))

assign( paste("var","b",c+length(taxa)^2,sep="."),tkbutton(frame.covar,text="  0.5  ",bg="gray80") )
tkconfigure( get(paste("var","b",c+length(taxa)^2,sep=".")),command=get(paste("func",c,sep=".")) )#FUNCTION!!!
}
}

done.b<-tkbutton(tt,text="  Done  ",command=function() tclvalue(change)<--999,font=f.tit)
tkbind(tt,"<Destroy>",function() tclvalue(change)<--999)

#### instruction lines ####

tkgrid(tklabel(tt,text=" ",font=Sp))
tkgrid(tklabel(tt,text="  Set any interaction restrictions that you would like incorporated into the model  ",font=f.tit))
tkgrid(tklabel(tt,text="(column variables affect row variables)"))
tkgrid(tklabel(tt,text=" "))

#### set column widths ####

colwidth<-paste(rep(" ",max(nchar(c(taxa,drivers)))),collapse=" ")

line1.var<-paste("tkgrid(",
	paste(
		rep("tklabel(frame.var,text=colwidth)",length(taxa)+1),
		collapse=","),
	")",sep="")

if(length(drivers)>0){
line1.covar<-paste("tkgrid(",
	paste(
		rep("tklabel(frame.covar,text=colwidth)",length(drivers)),
		collapse=","),
	")",sep="")
} else line1.covar<-"tkgrid(tklabel(frame.covar,text=colwidth))"
eval(parse(text=line1.var))
eval(parse(text=line1.covar))

tkwm.geometry(tt,"-30+30")

#### make row of taxa & driver labels ####

	taxa.labs<-"tklabel(frame.var,text=\" \")" #blank column
for(i in 1:length(taxa)){
	taxa.labs<-paste(taxa.labs,paste("tklabel(frame.var,text=\"",taxa[i],"\")",sep=""),sep=",")} #taxa names
	eval(parse( text=  paste("tkgrid(",taxa.labs,")")  ))

if(length(drivers)>0){
	driver.labs<-NULL
for(i in 1:length(drivers)){
	driver.labs<-paste(driver.labs,paste("tklabel(frame.covar,text=\"",drivers[i],"\")",sep=""),sep=",")}
	driver.labs<-substr(driver.labs,2,nchar(driver.labs))
	eval(parse( text=  paste("tkgrid(",driver.labs,")")  ))
} else tkgrid(tklabel(frame.covar,text="NA"))

#### arrange button rows ####

for(i in 1:length(taxa)){
	buttons<-NULL
for(j in 1:length(taxa)){
	buttons<-paste(buttons,paste(",","var.b.",j+((i-1)*length(taxa)),sep=""),sep="")
}
	line<-paste("tkgrid(tklabel(frame.var,text=\"",taxa[i],"\")",buttons,")",sep="")#;print(line)
	eval(parse(text=line))
}

if(length(drivers)>0){
for(i in 1:length(taxa)){
	buttons<-NULL
for(j in 1:length(drivers)){
	buttons<-paste(buttons,paste(",","var.b.",(length(taxa)^2+j)+((i-1)*length(drivers)),sep=""),sep="")
}
	buttons<-substr(buttons,2,nchar(buttons))
	line<-paste("tkgrid(",buttons,")",sep="")
	eval(parse(text=line))
}
} else {
for(i in 1:length(taxa)){tkgrid(tklabel(frame.covar,text="-",font=Ph))
		tkgrid.rowconfigure(frame.covar,i,pad=2)}
}

#### print frames ####

tkgrid(tklabel(frame.var,text=""))
tkgrid(tklabel(frame.covar,text=""))
tkgrid(tklabel(frame.overall,text=""),tklabel(frame.overall,text="Variates",font=f.tit),
	tklabel(frame.overall,text="Covariates",font=f.tit))
tkgrid(tklabel(frame.overall,text=""),frame.var,frame.covar,tklabel(frame.overall,text=""))
tkgrid(frame.overall)
tkgrid(tklabel(tt,text=" "))
tkgrid(done.b)
tkgrid(tklabel(tt,text=" "))
tkfocus(tt)

#### set up button interactions ####

while(changeVal>=0){
tkwait.variable(change)
changeVal<-as.numeric(tclvalue(change))
if(changeVal>=1) {
	res[changeVal]<-(res[changeVal]+1)
	res[changeVal]<-res[changeVal]-floor(res[changeVal]/3)*3
#	print(res)
#	print(changeVal)
if(res[changeVal]==0) {prob<-"  0.5  ";bgcol<-"gray80"}
if(res[changeVal]==1) {prob<-"   0   ";bgcol<-"lightcoral"}
if(res[changeVal]==2) {prob<-"   1   ";bgcol<-"darkolivegreen1"}
tkconfigure(get(paste("var.b",changeVal,sep=".")),text=prob,bg=bgcol)  }
} 

tkdestroy(tt)

#### make interaction restriction matrices ####

res[res==0]<-0.5;res[res==1]<-0;res[res==2]<-1
res.var<-res[1:length(taxa)^2]
if (length(drivers)>0) res.covar<-res[(length(taxa)^2+1):length(res)]
indexBGlobal<-matrix(res.var,nrow=length(taxa),ncol=length(taxa),byrow=T)
if(length(drivers)>0) {
indexCGlobal<-matrix(res.covar,nrow=length(taxa),ncol=length(drivers),byrow=T)} else {
indexCGlobal<-NULL}

cbind(indexBGlobal,indexCGlobal)->indexBCGlobal
indexBCGlobal

}

