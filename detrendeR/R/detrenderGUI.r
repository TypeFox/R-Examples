detrenderGUI = function() {
#require(tcltk2)
#require(tcltk) || stop("tcltk support is absent") 

tt <- tktoplevel()
tkwm.resizable(tt,1,0)
done <- tclVar(0)
tkwm.title(tt, "DetrendeR")

#size=c(440,415,0,0)
#geo <- paste(size[1], "x", size[2], "+", size[3],"+", size[4], sep = "")
#    tkwm.geometry(tt, geo)
    tkwm.deiconify(tt)
    tkgrab.set(tt)
    #tkwm.resizable(tt)
    #tkfocus(tt)
    
#################################
# TREE MASK
#################################

frame1.parent <- tkframe(tt, relief = "groove", borderwidth = 2)
frame1 <- tkframe(frame1.parent, relief = "groove", borderwidth = 2)
frame1.1 <- tkframe( frame1, relief = "groove", borderwidth = 0)
frame1.2 <- tkframe( frame1, relief = "groove", borderwidth = 0)
frame1.3 <- tkframe( frame1, relief = "groove", borderwidth = 0)
  tkpack(tklabel(frame1, text = "    - Tree Mask -    ", foreground = "blue"))

site.value=tclVar(stc[1])
tree.value=tclVar(stc[2])
core.value=tclVar(stc[3])

site.entry <- tkentry(frame1.1, textvariable = site.value, width = 3)
site.lab<-tklabel(frame1.1, text = "Site: ")
tkpack(site.lab,site.entry, side ="left")
 
tree.entry <- tkentry(frame1.2, textvariable = tree.value, width = 3)
tree.lab<-tklabel(frame1.2, text = "Tree:")
tkpack(tree.lab,tree.entry, side ="left")

core.entry <- tkentry(frame1.3, textvariable = core.value, width = 3)
core.lab<-tklabel(frame1.3, text = "Core:")
tkpack(core.lab,core.entry, side ="left")

tkpack(frame1.1,frame1.2, frame1.3, side="top")
tkpack(frame1, side="left", fill="both")
tkpack(frame1.parent, fill="both")

############################
#  SELECT SERIES
############################
 frame.1 <- tkframe(frame1.parent, relief = "groove", borderwidth = 2)
   
    remove.shorter.series.value <- tclVar(remove.shorter.series)
    delete.shorter.series.value <- tclVar(delete.shorter.series)
 
    tkgrid(tklabel(frame.1, text = "          - Select Series - ", foreground = "blue"))
    select.series.cbut <- tkcheckbutton(frame.1, text = "Delete series shorter than:", variable = remove.shorter.series.value)
    select.series.entry <- tkentry(frame.1, textvariable = delete.shorter.series.value, width = 4)
 
    #select.series.corr.cbut <- tkcheckbutton(frame.1, text = "Delete series corr with master lower:", variable = remove.shorter.series.value)
    tkgrid(select.series.cbut, select.series.entry)
    tkpack(frame.1, fill = "x")

############################
#  MAKE SEGMENT PLOT
############################
 frame.1.1 <- tkframe(frame1.parent, relief = "groove", borderwidth = 2)
 makeSegPlot.value <- tclVar(makeSegPlot)
 tkgrid(tklabel(frame.1.1, text = " - Plot Series - ", foreground = "blue"))

 makeSegPlot.cbut <- tkcheckbutton(frame.1.1, text = "Make graph with length series", variable = makeSegPlot.value)
   tkgrid(makeSegPlot.cbut)
   tkpack(frame.1.1, fill = "x")

##################
# DETRENDING #####
##################
    frame2 <- tkframe(tt, relief = "groove", borderwidth = 2)
 tkpack(tklabel(frame2, text = " - Detrending - ", foreground = "blue"))
 
  interactive.detrend.value <- tclVar(interactive.detrend)
interactive.detrend.cbut <- tkcheckbutton(frame2, text = "Interactive detrending", variable = interactive.detrend.value)
tkpack(interactive.detrend.cbut, side="bottom", anchor = "w")

    #Det1frame <- tkframe(frame2, relief = "groove", borderwidth = 2)
Det1frame <- tkwidget( frame2, "labelframe", foreground="blue",text= "First detrend: ",relief = "groove", borderwidth = 2)
   
    #Det2frame <- tkframe(frame2, relief = "groove", borderwidth = 2)
Det2frame <- tkwidget( frame2, "labelframe", foreground="blue",text= "Second detrend: ",relief = "groove", borderwidth = 2)
  
Det1.1frame <- tkframe( Det1frame, relief = "groove", borderwidth = 0)
#Det1.1frame <- tkwidget( Det1frame, "labelframe", foreground="blue",text= "First detrend: ",relief = "groove", borderwidth = 0)
Det1.2frame <- tkwidget( Det1frame, "labelframe", foreground="blue",text= "Spline options: ", relief = "groove", borderwidth = 2)
#Det1.2frame <- tkframe( Det1frame, relief = "groove", borderwidth = 2)
Det1.2.1frame <- tkframe( Det1.2frame, relief = "groove", borderwidth = 0)
Det1.2.2frame <- tkframe( Det1.2frame, relief = "groove", borderwidth = 0)
Det1.2.3frame <- tkframe( Det1.2frame, relief = "groove", borderwidth = 0)
#    tkpack(tklabel(Det1.1frame, text = " -  First  Detrend - ", foreground = "blue"), anchor = "w")
#    tkpack(tklabel(Det1.2frame, text = " -  Spline options - ", foreground = "blue"), anchor = "w")
  
Det2.1frame <- tkframe( Det2frame, relief = "groove", borderwidth = 0)
#Det2.2frame <- tkframe( Det2frame, relief = "groove", borderwidth = 0)
Det2.2frame <- tkwidget( Det2frame, "labelframe", foreground = "blue", text= "Spline options: ", relief = "groove", borderwidth = 2)
#Det2.2frame <- tkframe( Det2frame, relief = "groove", borderwidth = 2)
Det2.2.1frame <- tkframe( Det2.2frame, relief = "groove", borderwidth = 0)
Det2.2.2frame <- tkframe( Det2.2frame, relief = "groove", borderwidth = 0)
Det2.2.3frame <- tkframe( Det2.2frame, relief = "groove", borderwidth = 0)

#    tkpack(tklabel(Det2.1frame, text = " -  Second  Detrend - ", foreground = "blue"), anchor = "w")
#    tkpack(tklabel(Det2.2frame, text = " -  Spline options - ", foreground = "blue"), anchor = "w")

#RadioButton = function (FRAME, variable= NULL, BUTTON=c("b.r1", "b.r2"), VALUE=c(1,2)){
#BUTTON<-as.vector(BUTTON)
#for (i in 1:length(BUTTON)){
# tkpack(tkradiobutton(FRAME, text = BUTTON[i], value = VALUE[i], variable = variable), anchor = "w")
#}
#}

method1.value=tclVar(method1)
method2.value=tclVar(method2)
n1.value = tclVar(n1)
#if (is.numeric(nPerc1)) {nPerc1.value = tclVar(round(nPerc1,3))} else {nPerc1.value = tclVar(nPerc1)} 
#nPerc1.value = tclVar(round(nPerc1,3))
nPerc1.value = tclVar(nPerc1)
p1.value = tclVar(p1)
n2.value = tclVar(n2)
nPerc2.value = tclVar(nPerc2)
#nPerc2.value = tclVar(round(nPerc2,3))
p2.value = tclVar(p2)

detrend.types = c("Neg Exp", "Spline", "Spline%", "Mean", "No Detrending")
detrend.values = c("ModNegExp", "Spline", "Spline%", "Mean", "No Detrending")
RadioButton(Det1.1frame, variable=method1.value, BUTTON = detrend.types, VALUE = detrend.values )
RadioButton(Det2.1frame, variable=method2.value, BUTTON = detrend.types, VALUE = detrend.values )

n1.entry <- tkentry(Det1.2.1frame, textvariable = n1.value, width = 5)
Det1.2.1lab<-tklabel(Det1.2.1frame, text = "Spline length:")
tkpack(Det1.2.1lab,n1.entry, side ="left")
 
nPerc1.entry<-tkentry(Det1.2.2frame, textvariable = nPerc1.value, width = 5)
Det1.2.2lab<-tklabel(Det1.2.2frame, text = "Spline ratio:  ")

tkpack(Det1.2.2lab,nPerc1.entry, side ="left", anchor="w")

p1.entry<-tkentry(Det1.2.3frame, textvariable = p1.value, width = 5)
Det1.2.3lab<-tklabel(Det1.2.3frame, text = "Value of p:    ")
tkpack(Det1.2.3lab,p1.entry, side ="left", anchor="w")

tkpack(Det1.2.1frame,Det1.2.2frame, Det1.2.3frame,side="top")

n2.entry <- tkentry(Det2.2.1frame, textvariable = n2.value, width = 5)
Det2.2.1lab<-tklabel(Det2.2.1frame, text = "Spline length:")
tkpack(Det2.2.1lab,n2.entry, side ="left")
 
nPerc2.entry<-tkentry(Det2.2.2frame, textvariable = nPerc2.value, width = 5)
Det2.2.2lab<-tklabel(Det2.2.2frame, text = "Spline ratio:  ")

tkpack(Det2.2.2lab,nPerc2.entry, side ="left", anchor="w")

p2.entry<-tkentry(Det2.2.3frame, textvariable = p2.value, width = 5)
Det2.2.3lab<-tklabel(Det2.2.3frame, text = "Value of p:    ")
tkpack(Det2.2.3lab,p2.entry, side ="left", anchor="w")

tkpack(Det2.2.1frame,Det2.2.2frame, Det2.2.3frame,side="top")

tkpack(Det1.1frame, Det1.2frame, side = "left", expand = 1, fill="x")
 
 tkpack(Det2.1frame, Det2.2frame, side = "left", expand = 1, fill="x")
 tkpack(Det1frame, Det2frame, side = "left", expand = 1, fill="x")

    tkpack(frame2, fill="x")

###############################
#AUTORREGRESIVE MODEL
###############################
    frame3 <- tkframe(tt, relief = "groove", borderwidth = 2)
    #frame3.1 <- tkframe(frame3, relief = "groove", borderwidth = 2)
   
    makeAr.value <- tclVar(makeAr)
    arMAX.value <- tclVar(arMAX)

    #tkgrid(tklabel(frame3.1, text = " - Autorregresive model - ", foreground = "blue"), columnspan = 2)
frame3.1 <- tkwidget( frame3, "labelframe", foreground="blue",text= "Autorregresive model: ",relief = "groove", borderwidth = 2)
    makeAr.cbut <- tkcheckbutton(frame3.1, text = "Autorregresive model of order:", variable = makeAr.value)
    arMAX.entry <- tkentry(frame3.1, textvariable = arMAX.value, width = 2)
   
    tkgrid( makeAr.cbut,  arMAX.entry)
   # tkpack(frame3.1, fill = "x")
   # tkpack(frame3, fill = "x")

# frame3.2 <- tkframe(frame3, relief = "groove", borderwidth = 2)
#frame4.0 <- tkframe(frame1.parent, relief = "groove") 
frame3.2 <- tkwidget( frame3, "labelframe", foreground="blue",text= "Mean: ",relief = "groove", borderwidth = 2)

          #tkgrid(tklabel(frame3.2,text="Mean:", foreground = "blue"))
          #tkpack(frame4.0, fill = "x") 
          #frame4 <- tkframe(frame1.parent, relief = "groove") 
          rb1 <- tkradiobutton(frame3.2)
          rb2 <- tkradiobutton(frame3.2)
          rbValue <- tclVar(biweightMean)
          tkconfigure(rb1,variable=rbValue,value=TRUE)
          tkconfigure(rb2,variable=rbValue,value=FALSE)
          tkgrid(tklabel(frame3.2,text="Robust     "),rb1)
          tkgrid(tklabel(frame3.2,text="Arithmetic "),rb2)
          #tkpack(frame3.1, fill = "x")
 tkpack(frame3.1, frame3.2, side = "left", expand = 1, fill="both")

          tkpack(frame3, fill = "x") 
          


##############################
#EPS
################################

frame4 <- tkframe(tt, relief = "groove", borderwidth = 2)
frame4.1 <- tkframe(frame4, relief = "groove", borderwidth = 2)
frame4.1.1 <- tkframe( frame4.1, relief = "groove", borderwidth = 0)
frame4.1.2 <- tkframe( frame4.1, relief = "groove", borderwidth = 0)
frame4.1.3 <- tkframe( frame4.1, relief = "groove", borderwidth = 0)

frame4.2 <- tkframe(frame4, relief = "groove", borderwidth = 2)
frame4.2.1 <- tkframe( frame4.2, relief = "groove", borderwidth = 0)
frame4.2.2 <- tkframe( frame4.2, relief = "groove", borderwidth = 0)
frame4.2.3 <- tkframe( frame4.2, relief = "groove", borderwidth = 0)

run.win.analysis.value=tclVar(run.win.analysis)
winLength.value=tclVar(winLength)
stepWin.value=tclVar(stepWin)
make.common.EPS.value=tclVar(make.common.EPS)
make.select.period.EPS.value=tclVar(make.select.period.EPS)
first.year.common.value=tclVar(first.year.common)
last.year.common.value=tclVar(last.year.common)
 
frame4.label<-tklabel(frame4, text = " - EPS analysis - ", foreground = "blue")

run.win.analysis.cbut <- tkcheckbutton(frame4.1.1, text = "EPS window analysis", variable = run.win.analysis.value)
tkpack(frame4.label, fill = "x")
tkpack(run.win.analysis.cbut, fill = "x", side ="left")
tkpack(frame4.1.1, fill = "x")

winLength.entry <- tkentry(frame4.1.2, textvariable = winLength.value, width = 3)
winLength.lab<-tklabel(frame4.1.2, text = "Window length: ")
tkpack(winLength.lab,winLength.entry, side ="left") 
tkpack(frame4.1.2, fill = "x")

stepWin.entry <- tkentry(frame4.1.3, textvariable = stepWin.value, width = 3)
stepWin.lab<-tklabel(frame4.1.3, text = "                  Lag: ")
tkpack(stepWin.lab,stepWin.entry, side ="left") 
tkpack(frame4.1.3, fill = "x")

make.common.EPS.cbut <- tkcheckbutton(frame4.2.1, text = "Common interval", variable = make.common.EPS.value)
#tkpack(frame4.2.label, fill = "x")
tkpack(make.common.EPS.cbut, fill = "x", side ="left")
#tkpack(frame4.2.1, fill = "x")

   make.select.period.EPS.cbut <- tkcheckbutton(frame4.2.2, text = "Period        First year:", variable = make.select.period.EPS.value)
   first.year.common.entry <- tkentry(frame4.2.2, textvariable =  first.year.common.value, width = 4)
   tkpack(make.select.period.EPS.cbut,first.year.common.entry, fill = "x", side ="left")
#tkpack(frame4.2.2, fill = "x")

last.year.common.entry <- tkentry(frame4.2.3, textvariable =  last.year.common.value, width = 4)
last.year.common.label<-tklabel(frame4.2.3, text = "                         Last year:")

tkpack(last.year.common.label, last.year.common.entry, side="left")

tkpack(frame4.2.1,frame4.2.2,frame4.2.3, fill = "x", side="top")
tkpack(frame4.1, fill = "both", expand=1, side="left")
tkpack( frame4.2,fill = "both", expand=1, side="right")
tkpack(frame4, fill = "x")

##########################################################
# OnOK
##########################################################
#tclvalue=function(x) tclvalue(x)

OnOk = function(){

#remove.shorter.series<<-as.logic(tclvalue(remove.shorter.series.value))
#delete.shorter.series<<-toNumber(tclvalue(delete.shorter.series.value))
.assign("remove.shorter.series",as.logic(tclvalue(remove.shorter.series.value)))
.assign("delete.shorter.series",toNumber(tclvalue(delete.shorter.series.value)))

#makeSegPlot<<-as.logic(tclvalue(makeSegPlot.value))
.assign("makeSegPlot", as.logic(tclvalue(makeSegPlot.value)))
#interactive.detrend<<-as.logic(tclvalue(interactive.detrend.value))
.assign("interactive.detrend",as.logic(tclvalue(interactive.detrend.value)))

#makeFirstDetrending <<- TRUE
#method1<<-tclvalue(method1.value)
#n1<<-toNumber(tclvalue(n1.value))
#nPerc1<<-toNumber(tclvalue(nPerc1.value))
#p1<<-toNumber(tclvalue(p1.value))
#if(method1=="No Detrending")  makeFirstDetrending <<- FALSE
#first.detrending.method <<- GetDetrendMethod(method1,n1 ,nPerc1, p1)

.assign("makeFirstDetrending", TRUE)
.assign("method1",tclvalue(method1.value))
.assign("n1",toNumber(tclvalue(n1.value)))
.assign("nPerc1",toNumber(tclvalue(nPerc1.value)))
.assign("p1",toNumber(tclvalue(p1.value)))
if(method1=="No Detrending")  .assign("makeFirstDetrending", FALSE)
.assign("first.detrending.method", GetDetrendMethod(method1,n1 ,nPerc1, p1))

#makeSecondDetrending <<- TRUE
#method2<<-tclvalue(method2.value)
#n2<<-toNumber(tclvalue(n2.value))
#nPerc2<<-toNumber(tclvalue(nPerc2.value))
#p2<<-toNumber(tclvalue(p2.value))

.assign("makeSecondDetrending", TRUE)
.assign("method2",tclvalue(method2.value))
.assign("n2",toNumber(tclvalue(n2.value)))
.assign("nPerc2",toNumber(tclvalue(nPerc2.value)))
.assign("p2",toNumber(tclvalue(p2.value)))

if(method2=="No Detrending")  .assign("makeSecondDetrending", FALSE) #makeSecondDetrending <<- FALSE
#second.detrending.method <<- GetDetrendMethod(method2,n2 ,nPerc2, p2)
.assign("second.detrending.method", GetDetrendMethod(method2,n2 ,nPerc2, p2))

.assign("makeAr",as.logic(tclvalue(makeAr.value))) #makeAr<<-as.logic(tclvalue(makeAr.value))
.assign("arMAX",toNumber(tclvalue(arMAX.value))) #arMAX<<-toNumber(tclvalue(arMAX.value))
.assign("biweightMean", as.logic(tclvalue(rbValue))) 

#run.win.analysis<<-as.logic(tclvalue(run.win.analysis.value))
#winLength<<-toNumber(tclvalue(winLength.value))
#stepWin<<-toNumber(tclvalue(stepWin.value))
.assign("run.win.analysis", as.logic(tclvalue(run.win.analysis.value)))
.assign("winLength", toNumber(tclvalue(winLength.value)))
.assign("stepWin", toNumber(tclvalue(stepWin.value)))

#make.common.EPS<<-as.logic(tclvalue(make.common.EPS.value))
#make.select.period.EPS<<-as.logic(tclvalue(make.select.period.EPS.value))
#first.year.common<<-toNumber(tclvalue(first.year.common.value))
#last.year.common<<-toNumber(tclvalue(last.year.common.value))
.assign("make.common.EPS",as.logic(tclvalue(make.common.EPS.value)))
.assign("make.select.period.EPS",as.logic(tclvalue(make.select.period.EPS.value)))
.assign("first.year.common",toNumber(tclvalue(first.year.common.value)))
.assign("last.year.common",toNumber(tclvalue(last.year.common.value)))

if(make.select.period.EPS){
if (any(is.na(first.year.common), is.na(last.year.common))) {    
#make.select.period.EPS<<-FALSE
#make.common.EPS<<-TRUE}
.assign("make.select.period.EPS", FALSE)
.assign("make.common.EPS",TRUE)
}
}

stc.<-c(toNumber(tclvalue(site.value)),toNumber(tclvalue(tree.value)),toNumber(tclvalue(core.value)))
if (sum(stc.)!=8) {
tk_messageBox(type="ok", "Please correct the tree mask!")
} else {
#stc[1]<<-stc.[1]
#stc[2]<<-stc.[2]
#stc[3]<<-stc.[3]
.assign("stc", stc.)
tkdestroy(tt)
#detrender()
tclvalue(done) <- 1
}
}

###########################################################
#
###########################################################
 frame5 <- tkframe(tt, relief = "groove", borderwidth = 2)
# Save.but <- tkbutton(frame5, text = "Save settings", command = function() tkdestroy(tt))
 Cancel.but <- tkbutton(frame5, text = "Cancel", command = function() tkdestroy(tt))
 Ok.but <- tkbutton(frame5, text = "Ok", command = OnOk) 
# tkpack(Save.but, Cancel.but , Ok.but, side = "left", expand = "TRUE",   fill = "x")
 tkpack( Cancel.but , Ok.but, side = "left", expand = "TRUE",   fill = "x")
 tkpack(frame5, fill="x")
###########################################################
#tk.block(tt)
    tkfocus(tt)
    tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(tt, "<KeyPress-Return>", function() OnOk())
    tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
    tkwait.variable(done)
    if (tclvalue(done) == "2") {
    tkdestroy(tt)
    return(FALSE)  }
     if (tclvalue(done) == "1") {
    tkdestroy(tt)
    return(TRUE)  }
}
# detrenderGUI()
