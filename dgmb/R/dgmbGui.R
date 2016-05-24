dgmbGui <-
function()
{
cat (" --------------------------------------------------\n")
cat ("| Simulating Data for PLS Mode B Structural Models |\n")
cat ("| Version: 1.2                                     |\n")                                 
cat ("| Depends: abind, tcltk, MASS                      |\n")
cat ("| Licence: GPL >= 2                                |\n")
cat ("| System Requirements: tktable, BWidget            |\n")
cat (" --------------------------------------------------\n")
#cat ("Linux users must to download the bwidget and Tktable2.10 packages\n from http://sourceforge.net/projects/tktable/files/tktable/\n")

#packages 
#library (tcltk)
tclRequire("BWidget")
tclRequire("Tktable")
                                
#Main Frame
tt <-tktoplevel()
tkwm.geometry(tt, "780x360") 
tkwm.resizable(tt, FALSE, FALSE)
tktitle(tt) <- "dgmb: Simulating Data for PLS Mode B Structural Models"
topMenu <- tkmenu(tt)
tkconfigure(tt, menu=topMenu)

#Main Menu
fileMenu <- tkmenu(topMenu,tearoff=FALSE)
AcercadeMenu <- tkmenu(topMenu,tearoff=TRUE)
ayudaMenu <- tkmenu(topMenu,tearoff=TRUE)
tkadd(topMenu,"cascade",label="File",menu=fileMenu)
tkadd(topMenu,"cascade",label="About dgmbGui",menu=AcercadeMenu)
tkadd(topMenu,"cascade",label="Help",menu=ayudaMenu)

#File Menu
tkadd(fileMenu,"command",label="Exit",command=function() tkdestroy(tt))

#About Menu
tkadd(AcercadeMenu,"command",label="Version",command = function () tkmessageBox(title="Version",message="dgmbGui (Beta 1.2)", icon="info", type="ok"))
tkadd(AcercadeMenu,"command",label="Developed by",command = function () tkmessageBox(title="Developed by",message="Dr. Alba Ester Martinez Ruiz, amartine@ucsc.cl \nMg. Claudia Loreto Martinez Araneda, cmartinez@ucsc.cl"))


#Help Menu
tkadd(ayudaMenu,"command",label="dgmb-package",command= function () ViewHelp())
tkadd(ayudaMenu,"command",label="R-package",command= function() help.start())

#Number of data set
Ntxt <- tclVar("500")#default value

#Sample size
ntxt <- tclVar("250")#default value
entry.Ntxt <-tkentry(tt,width="4",textvariable=Ntxt)
entry.ntxt <-tkentry(tt,width="4",textvariable=ntxt)

tkgrid(tklabel(tt,text="I N P U T   P A R A M E T E R S", foreground="RED"),row=0, columnspan=5)
tkgrid(tklabel(tt,text="Number of data sets",  foreground="BLACK", padx=0),row=1, column=0, ipadx=0, sticky="w",columnspan=1 )
tkgrid(entry.Ntxt, row=1, column= 1, sticky="w", columnspan=1)
tkgrid(tklabel(tt,text="Sample size", foreground="BLACK", padx=0),row=2, column=0, ipadx=0, sticky="w", columnspan=1)
tkgrid(entry.ntxt, row=2, column= 1, sticky="w")

#Number of indicators for each block of variables
indicators <- c(2,4,6,8)

comboBox <- tkwidget(tt,"ComboBox",editable=FALSE,values=indicators)
tkgrid(tklabel(tt,text="Indicators for each block of variables"), row=3, column=0, sticky="w")
tkgrid(comboBox)

#Distribution
distribution <- c("Normal")  
tl1<-tklistbox(tt,height=2,selectmode="single",background="white")
tkgrid(tklabel(tt,text="Distribution (read only)"), row=3, column=2, sticky="w")
tkgrid(tl1, sticky="e", row=4, column=2, sticky="w")
tkinsert(tl1,"end",distribution[1])
tkselection.set(tl1,1)  #Default selection
tkconfigure(tl1, state = "disabled") #read only

#Outer weights
rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)
rbValue <- tclVar(TRUE)
tkconfigure(rb1,variable=rbValue, value=TRUE, text="Equal")
tkconfigure(rb2,variable=rbValue, value=FALSE, text="Different")
tkgrid(tklabel(tt,text="Outer weights"), row=5, column=0, sticky="w")
tkgrid(rb1, sticky="w")
tkgrid(rb2, sticky="w")

#Matrix of structural relationships (r.s)
tkgrid(tklabel(tt,text="Matrix of linear effects (read only)", foreground="BLACK"),row=1, column=3 )

rs <- matrix(c(0,0,0,1,
                0,0,0,1,
                0,0,0,1,
                1,1,1,0),4,4,byrow=TRUE)

rss <- matrix(c(".","Y1","Y2","Y3","Y4",
                "Y1",0,0,0,1,
                "Y2",0,0,0,1,
                "Y3", 0,0,0,1,
                "Y4",1,1,1,0),5,5,byrow=TRUE)
                              
tclrs <- tclArray()
for (i in (1:5))
      for (j in (1:5)) {tclrs[[i-1,j-1]] <- rss[i,j]}

table1 <- tkwidget(tt,"table",variable=tclrs,rows=5,cols=5,titlerows=1,titlecols=1,selectmode="extended",colwidth=5,foreground="white")
tkgrid(table1, row=2, column=3)
tkconfigure(table1,variable=tclrs,background="grey",selectmode="extended", state = "disabled", resizeborders="none")
 

#Matrix of no linear effect indicators (r.ie)
tkgrid(tklabel(tt,text="Matrix of nonlinear and interaction effects (read only)",foreground="BLACK"),row=3, column=3)
rie <- matrix(c(0,0,1,
                0,1,0,
	        0,0,0),3,3,byrow=TRUE)

riee <- matrix(c(".", "Y1","Y2","Y3",
                      "Y1",0,0,1,
	              "Y2",0,1,0,
		      "Y3",0,0,0),4,4,byrow=TRUE)
		             
tclrie <- tclArray()
for (i in (1:4))
     for (j in (1:4)){tclrie[[i-1,j-1]] <- riee[i,j]}        

table2 <- tkwidget(tt,"table",variable=tclrie,rows=4,cols=4,titlerows=1,titlecols=1,selectmode="extended",colwidth=5,foreground="white")
tkgrid(table2, row=4, column=3)
tkconfigure(table2,variable=tclrie,background="gray",selectmode="extended", state = "disabled", resizeborders="none")

#trigger ViewParameters function
ViewParameters.but <- tkbutton(tt,text="View fixed parameters",command= function() ViewParameters(comboBox,rbValue))
tkgrid(ViewParameters.but, row=8, column=2)
tkgrid(tklabel(tt,text="     "))
tkfocus(tt)

#trigger to Process
OK.but <- tkbutton(tt,text="Data Generation", command=function() ToProcess(as.numeric(tclvalue(Ntxt)), as.numeric(tclvalue(ntxt)), "F", indicators[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1], as.integer(tclvalue(rbValue)) )  )   
tkgrid(OK.but, row=8, column=3)
tkpack.propagate(tt, FALSE)
tkfocus(tt)
    
}#End dgmbGui

