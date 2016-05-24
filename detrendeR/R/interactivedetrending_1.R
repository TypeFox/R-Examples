InteractiveDetrending <-
function (rwl, Detrend=NULL, method="Mean",  n=32, nPerc=0.67, p=0.5, ...){


.curve <- DETRENDING_INTERACTIVE_OUTPUT <-Detrend
.assign("DETRENDING_INTERACTIVE_OUTPUT", .curve)
.seriesNames<-colnames(rwl)
.x = as.numeric(rownames(rwl))
done <- tclVar(0)

    FLAG=FALSE
   # DETRENDING_INTERACTIVE_FLAG<<-FALSE 
    .assign("DETRENDING_INTERACTIVE_FLAG", FALSE)
    DETRENDING_METHOD = GetDetrendMethod(method=method,  n=n, nPerc=nPerc, p=p)
    DETRENDING_INTERACTIVE_OUTPUT_CHANGE <- matrix(NA, ncol=3, nrow=ncol(rwl))
    rownames(DETRENDING_INTERACTIVE_OUTPUT_CHANGE) <- (colnames(rwl))
    DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,1]<- 1:ncol(rwl)
    DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,2]<- (colnames(rwl))
    DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,3]<- DETRENDING_METHOD
    
	#DETRENDING_INTERACTIVE_OUTPUT_CHANGE <<- DETRENDING_INTERACTIVE_OUTPUT_CHANGE    
    .assign("DETRENDING_INTERACTIVE_OUTPUT_CHANGE", DETRENDING_INTERACTIVE_OUTPUT_CHANGE) 
if (is.null(Detrend)) {

    DETRENDING_INTERACTIVE_OUTPUT <- try(apply(rwl, 2, RemoveTrend, method=method, BandwidthPerc=nPerc, Bandwidth = n, P=p),silent=TRUE)
    } 


METHOD.value<-tclVar(method)
BANDWIDTH.value<-tclVar(n)
BANDWIDTH.P.value<-tclVar(nPerc)
P.value<-tclVar(p)
cbValue <- tclVar("0")

draw <-function(...) {
show.all.series="l"
    currentMETHOD<-method
    seriesIndex <- as.integer(tkcurselection(TL))+1

    .assign("currentBANDWIDTH.P", nPerc) 	#currentBANDWIDTH.P <<- nPerc
    .assign("currentBANDWIDTH" , n ) 		#currentBANDWIDTH <<- n
    .assign("currentP", p )			#currentP <<- p
    FLAG<-as.logic(tclvalue(cbValue))
    .curve<-DETRENDING_INTERACTIVE_OUTPUT

 if(!FLAG) show.all.series ="n"
 
    matplot(.x, rwl, lty=1, col="grey", type=show.all.series, main=.seriesNames[seriesIndex], ylab="", xlab="Year", las=1)
    lines(.x, rwl[,seriesIndex], type="l", col="blue")
  
   METHOD_CHANGED <- GetDetrendMethod(method=currentMETHOD,  n=currentBANDWIDTH, nPerc=currentBANDWIDTH.P, p=currentP)
    temp = rwl[,seriesIndex]
    
    #lines(.x, rwl[,panel$Series], col="blue")
    mtext(paste(METHOD_CHANGED, sep=""), line=0.5, side =3, adj =1, cex=0.90, col="blue",font=1)
    mtext(DETRENDING_INTERACTIVE_OUTPUT_CHANGE[seriesIndex,3], line=0.5, side =3, adj =0, cex=0.90, col="black",font=1)
    
   lines(.x, .curve[ ,seriesIndex], col="black")
   lines(.x, temp, col=2)
   cv<-try(RemoveTrend (temp, method = currentMETHOD, BandwidthPerc=currentBANDWIDTH.P, Bandwidth=currentBANDWIDTH, P=currentP) , silent=TRUE)
   lines(.x, cv , col=4)
    } 
    
re.draw <-function(...) {
show.all.series="l"
    currentMETHOD<-as.character(tclvalue(METHOD.value))
    seriesIndex <- as.integer(tkcurselection(TL))+1
    #currentBANDWIDTH.P <<- as.numeric(tclvalue("BANDWIDTH.P.value"))
    #currentBANDWIDTH <<- as.numeric(tclvalue("BANDWIDTH.value"))
    #currentP <<- as.numeric(tclvalue("P.value"))
    
    .assign("currentBANDWIDTH.P", as.numeric(tclvalue("BANDWIDTH.P.value")))
    .assign("currentBANDWIDTH", as.numeric(tclvalue("BANDWIDTH.value")))
    .assign("currentP", as.numeric(tclvalue("P.value")))
    
    FLAG<-as.logic(tclvalue(cbValue))
    .curve<-.get("DETRENDING_INTERACTIVE_OUTPUT")

 if(!FLAG) show.all.series ="n"
 
    matplot(.x, rwl, lty=1, col="grey", type=show.all.series, main=.seriesNames[seriesIndex], ylab="", xlab="Year", las=1)
    lines(.x, rwl[,seriesIndex], type="l", col="blue")
  
   METHOD_CHANGED <- GetDetrendMethod(method=currentMETHOD,  n=currentBANDWIDTH, nPerc=currentBANDWIDTH.P, p=currentP)
    temp = rwl[,seriesIndex]
    
    #lines(.x, rwl[,panel$Series], col="blue")
    mtext(paste(METHOD_CHANGED, sep=""), line=0.5, side =3, adj =1, cex=0.90, col="blue",font=1)
    mtext(DETRENDING_INTERACTIVE_OUTPUT_CHANGE[seriesIndex,3], line=0.5, side =3, adj =0, cex=0.90, col="black",font=1)
   
   lines(.x, temp, col=2)
   lines(.x, .curve[ ,seriesIndex], col="black")
   
   cv<-try(RemoveTrend (temp, method = currentMETHOD, BandwidthPerc=currentBANDWIDTH.P, Bandwidth=currentBANDWIDTH, P=currentP) , silent=TRUE)
   lines(.x, cv , col=4)
    } 
    

rwl.draw.change = function(...){
        currentMETHOD<-as.character(tclvalue(METHOD.value))
        seriesIndex <- as.integer(tkcurselection(TL))+1
        currentBANDWIDTH.P <- as.numeric(tclvalue("BANDWIDTH.P.value"))
    currentBANDWIDTH <- as.numeric(tclvalue("BANDWIDTH.value"))
    currentP <- as.numeric(tclvalue("P.value"))
     DETRENDING_INTERACTIVE_OUTPUT_CHANGE[seriesIndex,3]<<- GetDetrendMethod(method=currentMETHOD,  n=currentBANDWIDTH, nPerc=currentBANDWIDTH.P, p=currentP)
      temp = rwl[,seriesIndex]
         cv<-try(RemoveTrend (temp, method = currentMETHOD, BandwidthPerc=currentBANDWIDTH.P, Bandwidth=currentBANDWIDTH, P=currentP) , silent=TRUE)
        DETRENDING_INTERACTIVE_OUTPUT_CHANGE[seriesIndex,3] <<- GetDetrendMethod(method=currentMETHOD,  n=currentBANDWIDTH, nPerc=currentBANDWIDTH.P, p=currentP)
      DETRENDING_INTERACTIVE_OUTPUT[,seriesIndex]<<-cv
           .assign("DETRENDING_INTERACTIVE_OUTPUT", DETRENDING_INTERACTIVE_OUTPUT)
       
         re.draw()
          
      }

  
 RadioButton = function (FRAME, variable= NULL, BUTTON=c("b.r1", "b.r2"), VALUE=c(1,2), command=invisible,...){
 BUTTON<-as.vector(BUTTON)
 for (i in 1:length(BUTTON)){
 tkpack(tkradiobutton(FRAME, text = BUTTON[i], value = VALUE[i], variable = variable, command=command), anchor = "w")
 }
 }

detrend.types = c("Neg Exp", "Spline", "Spline%", "Mean")
detrend.values = c("ModNegExp", "Spline", "Spline%", "Mean")


ttt <- tktoplevel() #bg="white") 
tkwm.title(ttt, "Interactive detrending")
tkwm.deiconify(ttt)
tkgrab.set(ttt)
  tkwm.resizable(ttt,1,0)

tt<-tkframe(ttt, relief = "groove", borderwidth = 2)
   frame2 <- tkframe(ttt, relief = "groove", borderwidth = 2)
   frame2.1 <- tkwidget(frame2, "labelframe", text=" Series: ", relief = "groove", borderwidth = 2)
   frame2.2 <- tkwidget(frame2, "labelframe", text=" Detrend method: ",relief = "groove", borderwidth = 2)
   frame2.3 <- tkwidget(frame2, "labelframe", text=" Spline options: ",relief = "groove", borderwidth = 2)
   frame2.4 <- tkwidget(frame2, "labelframe", text="",relief = "groove", borderwidth = 2)
   frame2.4.1 <- tkframe(frame2.4)
################################################################################
SRC <- tkscrollbar(frame2.1, repeatinterval=25, command=function(...)tkyview(TL,...))
TL<-tklistbox(frame2.1,height=10,width= 20, selectmode="single", yscrollcommand=function(...)tkset(SRC,...),background="white")

#tkgrid(tklabel(frame2.1,text="Series:", foreground = "blue"))
tkgrid(TL,SRC)
tkgrid.configure(SRC,rowspan=10,sticky="nsw")


try(series <- colnames(rwl), silent=TRUE)->flag
if (class(flag)!="try-error"){
.assign("series", colnames(rwl))
for (i in (1:length(series))) { tkinsert(TL,"end",series[i]) }
  }
tkselection.set(TL,0)

tkbind(TL, "<ButtonRelease-1>", function(...) re.draw())


tkbind(TL, "<Up>", function(...) {    
i<-as.numeric(tkcurselection(TL))
if(i==0) i=1
tkinsert(TL, i+1, .seriesNames[i+1])
tkdelete(TL,i)
tkselection.set(TL,i-1) 
re.draw()})
        
tkbind(TL, "<Down>", function(...) {
i<-as.numeric(tkcurselection(TL)) 
tkselection.clear(TL, i)
if(i+1==length(.seriesNames))  i= i-1
tkselection.set(TL,i+1) 
re.draw()})


################################################################################ 
#tkpack(tklabel(frame2.2,text="Detrend method:", foreground = "blue"))
RadioButton(frame2.2, variable=METHOD.value, BUTTON = detrend.types, VALUE = detrend.values, command=re.draw )

################################################################################
#tkpack(tklabel(frame2.3, text ="Bandwidth", foreground = "blue"))
#tkpack(tkscale(frame2.3, label="Bandwidth (years)", command=re.draw, from=10, to=1000, variable="BANDWIDTH.value", showvalue=TRUE, resolution=5, orient="hor"), fill="x")
bandwidth.scale <- tkscale(frame2.3, label="Bandwidth (years)", from=n, to=n, variable="BANDWIDTH.value", showvalue=TRUE, resolution=5, orient="hor")#, command=re.draw)

tkpack(bandwidth.scale, fill="x")
#tkbind(bandwidth.scale, "<ButtonRelease-1>",function (...) re.draw())

#tkpack(tklabel(frame2.3, text =   "Bandwidth percentage", foreground = "blue"))
#tkpack(tkscale(frame2.3, label="Bandwidth (ratio)", command=re.draw, from=0.05, to=2.00, variable="BANDWIDTH.P.value", showvalue=TRUE, resolution=0.05, orient="hor"), fill="x")
bandwidth.p.scale <- tkscale(frame2.3, label="Bandwidth (ratio)",  from=nPerc, to=nPerc, variable="BANDWIDTH.P.value", showvalue=TRUE, resolution=0.05, orient="hor")#, command=re.draw)
tkpack(bandwidth.p.scale, fill="x")
#tkbind(bandwidth.p.scale, "<ButtonRelease-1>"=function(...) { re.draw() })

#tkpack(tklabel(frame2.3, text = "P", foreground = "blue"), side="left")
#tkpack(tkscale(frame2.3,label =   "P", command=re.draw, from=0, to=1.00, variable="P.value", showvalue=TRUE, resolution=0.1, orient="hor"), fill="x")
p.scale <- tkscale(frame2.3, label="P",  from=p, to=p, variable="P.value", showvalue=TRUE, resolution=0.05, orient="hor")#, command=re.draw)
tkpack(p.scale, fill="x")
#tkbind(p.scale, "<ButtonRelease-1>" =function(...) { re.draw() })


#   Change = function() {}
#next.series<-(tkbutton(frame2.4,text=" Next series",command=NextSeries)) #, width=10)


save.function = function(){
if(length(unique(DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,3]))>1){
     DETRENDING_INTERACTIVE_OUTPUT_CHANGE<<-rbind(c("Seq", "Series    ", "Detrending            "),DETRENDING_INTERACTIVE_OUTPUT_CHANGE)
   cat(rep("=", 38), "\n",sep="")
   DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,2] <-(format(DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,2], justify="left") )
   DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,3] <-(format(DETRENDING_INTERACTIVE_OUTPUT_CHANGE[,3], justify="left") )
   WriteMatrix(DETRENDING_INTERACTIVE_OUTPUT_CHANGE, col.width=3, ID=FALSE,row.names=F, col.names=F, sep="|")
   cat(rep("=", 38), "\n",sep="")
   }
 #  DETRENDING_INTERACTIVE_FLAG<<-TRUE
   tclvalue(done) <-1
}

  cancel.function = function(){   
  # DETRENDING_INTERACTIVE_FLAG<<-FALSE
     tclvalue(done) <-2
   }



change.bt<-(tkbutton(frame2.4,text=  "Change",command=rwl.draw.change))
close.bt<-(tkbutton(frame2.4,text=  "Close without saving",command=cancel.function))
save.bt<-(tkbutton(frame2.4,text=  "Save changes",command=save.function))
#tkpack(next.series)


cb <- tkcheckbutton(frame2.4.1,command=re.draw)

tkconfigure(cb,variable=cbValue, text="Show all series")


#tkpack(cb, tklabel(frame2.4,text="Show all series"), fill="x")
 
tkpack(cb, side="left", fill="x")
tkpack(frame2.4.1, fill="both")
tkpack(change.bt, fill="both")
tkpack(close.bt, fill="both")
tkpack(save.bt, fill="both")
            
tkpack(frame2.1, frame2.2, frame2.3, frame2.4,fill="x")
tkpack(frame2, fill="x")
tkpack(tt, fill="x")
 
#re.draw()
tkconfigure(bandwidth.scale, from=5, to=1000)
tkconfigure(bandwidth.p.scale, from=0.05, to=2)
tkconfigure(p.scale, from=0.0, to=1)

#bandwidth.p.scale
     tkbind(bandwidth.p.scale, "<ButtonRelease-1>", function(...) re.draw())
     tkbind(bandwidth.p.scale, "<ButtonRelease-2>", function(...) re.draw())
     tkbind(bandwidth.p.scale, "<Left>", function(...) re.draw())
     tkbind(bandwidth.p.scale, "<Right>", function(...) re.draw()) 
     tkbind(bandwidth.scale, "<ButtonRelease-1>", function(...) re.draw())
     tkbind(bandwidth.scale, "<ButtonRelease-2>", function(...) re.draw())
     tkbind(bandwidth.scale, "<Left>", function(...) re.draw())
     tkbind(bandwidth.scale, "<Right>", function(...) re.draw())
     tkbind(p.scale, "<ButtonRelease-1>", function(...) re.draw())
     tkbind(p.scale, "<ButtonRelease-2>", function(...) re.draw())
     tkbind(p.scale, "<Left>", function(...) re.draw())
     tkbind(p.scale, "<Right>", function(...) re.draw())
    draw()
   #  



tkbind(ttt, "<Destroy>", function() tclvalue(done) <- 2)
    #tkbind(ttt, "<KeyPress-Return>", function() OnOk())
    tkbind(ttt, "<KeyPress-Escape>", function() tclvalue(done) <- 2)
    #tkfocus(top_detrending)
    tkwait.variable(done)
    tkgrab.release(ttt)
    
    if (tclvalue(done) == "2") {
    tkdestroy(ttt)
    #dev.off(.Internal(dev.cur()))->a
   # try(dev.off(.Internal(dev.cur())),silent=T)
   try(dev.off(as.vector(dev.cur())),silent=T)
    .assign("DETRENDING_INTERACTIVE_FLAG", FALSE)
      }
     if (tclvalue(done) == "1") {
    tkdestroy(ttt)
    #dev.off(.Internal(dev.cur()))->a
    #try(dev.off(.Internal(dev.cur())),silent=T) 
      try(dev.off(as.vector(dev.cur())),silent=T)
    .assign("DETRENDING_INTERACTIVE_FLAG", TRUE)
    .assign("DETRENDING_INTERACTIVE_OUTPUT_CHANGE",DETRENDING_INTERACTIVE_OUTPUT_CHANGE)
    #.assign("DETRENDING_INTERACTIVE_OUTPUT", DETRENDING_INTERACTIVE_OUTPUT)
  }
}

