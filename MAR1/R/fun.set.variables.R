fun.set.variables<-function(heads){

tt<-tktoplevel()
tkwm.geometry(tt,"-30+30")
tkwm.title(tt,"Select variates and covariates")
frame.overall<-tkframe(tt) # relief= raised, sunken, flat, ridge, solid, or groove
frame.b<-tkframe(frame.overall,relief="groove",borderwidth=3)
frame.sub<-tkframe(frame.overall,relief="flat",borderwidth=1)
frame.done<-tkframe(frame.overall,relief="flat",borderwidth=1)
Sp<-tkfont.create(size=1)
Spp<-tkfont.create(size=4)
f.tit<-tkfont.create(size=10,weight="bold") # title font

tkgrid(tklabel(frame.b,text="",font=Sp))

track<-rep(0,length(heads)-2)
change<-tclVar(0)
changeVal<-as.numeric(tclvalue(change))

func.1<-function(){tclvalue(change)<-1}
tog.b.1<-tkbutton(frame.b,text="   not included   ",
		bg="gray90",foreground="gray60")
tkconfigure( tog.b.1,command=func.1 )

for(i in 2:length(track)){
	func<-paste("function(){tclvalue(change)<-",i,"}",sep="")
	assign(paste("func",i,sep="."),eval(parse(text=func)))

	assign( paste("tog","b",i,sep="."),tkbutton(frame.b,text="   not included   ",
		bg="gray90",foreground="gray60") )
	tkconfigure( get(paste("tog","b",i,sep=".")),command=get(paste("func",i,sep=".")) )
}

for(i in 1:length(heads)){
head<-tklabel(frame.b,text=heads[i])
colnum<-tklabel(frame.b,text=paste("column",i,":",sep=" "))
if(i==1) {tkgrid(tklabel(frame.b,text="    "),colnum,tklabel(frame.b,text="   "),head,tklabel(frame.b,text="    "),
	tklabel(frame.b,text="time blocks *",foreground="gray40"),tklabel(frame.b,text="       "))}
if(i==2) tkgrid(tklabel(frame.b,text=" "),colnum,tklabel(frame.b,text=" "),head,tklabel(frame.b,text=" "),tklabel(frame.b,text="dates **",foreground="gray40"))
if(i>2) tkgrid(tklabel(frame.b,text=" "),colnum,tklabel(frame.b,text=" "),head,tklabel(frame.b,text=" "),get(paste("tog","b",i-2,sep=".")))
#tkgrid.configure(head,sticky="w")
tkgrid.configure(colnum,sticky="e")
}

tkgrid.columnconfigure(frame.b,tog.b.1,minsize=100)

tkgrid(tklabel(frame.b,text=" ",font=Spp))

tkgrid(tklabel(frame.overall,text="",font=Sp),columnspan=3)
tkgrid(tklabel(frame.overall,
	text="Specify the variates and covariates that\nyou would like to include in the model",font=f.tit),
	columnspan=3)
tkgrid(tklabel(frame.overall,text="",font=Spp),columnspan=3)

tkgrid(tklabel(frame.overall,text="  "),frame.b,tklabel(frame.overall,text="  "))

tkgrid(tklabel(frame.overall,text=" ",font=Sp))
done.b<-tkbutton(frame.done,text="  Done  ",font=f.tit,command=function() tclvalue(change)<--999)
cancel.b<-tkbutton(frame.done,text=" Cancel ",font=f.tit,command=function() tclvalue(change)<--9999)
tkgrid(done.b,tklabel(frame.done,text="    "),cancel.b)
tkgrid(frame.done,columnspan=3)
tkgrid(tklabel(frame.overall,text=" "))

tkgrid(tklabel(frame.sub,text="\n*  The first column should be a dummy variable coding ",foreground="gray40"),sticky="w")
tkgrid(tklabel(frame.sub,text="             for blocks of data to be considered continuous.",foreground="gray40"),sticky="w")
tkgrid(tklabel(frame.sub,text=" ",font=Sp))
tkgrid(tklabel(frame.sub,text="** The second column should be dates ",foreground="gray40"),sticky="w")
tkgrid(tklabel(frame.sub,text="             in yyyy-mm-dd format.",foreground="gray40"),sticky="w")
tkgrid(tklabel(frame.sub,text=" "))
tkgrid(frame.sub,columnspan=3)

tkgrid(frame.overall)
tkbind(tt,"<Destroy>",function() tclvalue(change)<--9999)

while(changeVal>=0){
tkwait.variable(change)
changeVal<-as.numeric(tclvalue(change)) #;print(changeVal)
if(changeVal>0){
	track[changeVal]<-track[changeVal]+1
	if(track[changeVal]==3) track[changeVal]<-0
if(track[changeVal]==0) {b.lab<-"   not included   ";textcol<-"gray60";bgcol<-"gray90"}
if(track[changeVal]==1) {b.lab<-"        variate         ";textcol<-"blue3";bgcol<-"gray75"}
if(track[changeVal]==2) {b.lab<-"      covariate       ";textcol<-"firebrick";bgcol<-"gray75"}
tkconfigure(get(paste("tog.b",changeVal,sep=".")),text=b.lab,foreground=textcol,bg=bgcol) } #;print(track)
}

tkdestroy(tt)

track<-c(0,0,track)
if(changeVal<(-999)) track<-NULL
track

}

