inpbox_computational_options <-
function(inp){
	table<-gtkTable(rows=12,columns=4,homogeneous=FALSE)
	diag_label<-gtkVBox()
	labelgp<-list()
	labelgp$t1<-gtkLabel('Initial Temperature (C)                       ')
	labelgp$t2<-gtkLabel('Heat Transfer Coeff(cal/k min)                ')
	labelgp$t3<-gtkLabel('Simulation End Time (min)                     ')
	labelgp$t4<-gtkLabel('Conversion Limit (%)                          ')
	labelgp$t5<-gtkLabel('Induction Time (min)                          ')
	labelgp$t6<-gtkLabel('Initial pH                                                 ')
	labelgp$t7<-gtkLabel('Initial Cationic Conc. (mol/L)                     ')
	labelgp$t8<-gtkLabel('Initial Pressure (bar)                                     ')
	labelgp$t9<-gtkLabel('Inflow Temperature (C)                           ')
	labelgp$t10<-gtkLabel('Jacket Temperature (C)                          ')
	labelgp$t1$SetWidthChars(30)
	labelgp$t2$SetWidthChars(30)
	labelgp$t3$SetWidthChars(30)
	labelgp$t4$SetWidthChars(30)
	labelgp$t5$SetWidthChars(30)
	labelgp$t6$SetWidthChars(30)
	labelgp$t7$SetWidthChars(30)
	labelgp$t8$SetWidthChars(30)
	labelgp$t9$SetWidthChars(30)
	labelgp$t10$SetWidthChars(30)
	sapply(labelgp,diag_label$packStart,expand=TRUE,fill=TRUE)
	diag_entry<-gtkVBox()
	entrygp<-list()
	entrygp$e1<-gtkEntryNew()
	entrygp$e1$SetText(inp[1])
	entrygp$e2<-gtkEntryNew()
	entrygp$e2$SetText(inp[2])
	entrygp$e3<-gtkEntryNew()
	entrygp$e3$SetText(inp[3])
	entrygp$e4<-gtkEntryNew()
	entrygp$e4$SetText(inp[4])
	entrygp$e5<-gtkEntryNew()
	entrygp$e5$SetText(inp[5])
	entrygp$e6<-gtkEntryNew()
	entrygp$e6$SetText(inp[6])
	entrygp$e7<-gtkEntryNew()
	entrygp$e7$SetText(inp[7])
	entrygp$e8<-gtkEntryNew()
	entrygp$e8$SetText(inp[8])
	entrygp$e9<-gtkEntryNew()
	entrygp$e9$SetText(inp[9])
	entrygp$e10<-gtkEntryNew()
	entrygp$e10$SetText(inp[10])
	sapply(entrygp,diag_entry$packStart)
	diag_check<-gtkVBox()
	checkgp<-list()
	checkgp$t1<-gtkCheckButtonNewWithLabel(label='Diffusion Controlled Propagation')
	checkgp$t1$SetActive(inp[11])
	checkgp$t2<-gtkCheckButtonNewWithLabel(label='Segmental Diffusion Termination ')
	checkgp$t2$SetActive(inp[12])
	checkgp$t3<-gtkCheckButtonNewWithLabel(label='Diffusion Controlled Termination')
	checkgp$t3$SetActive(inp[13])
	checkgp$t4<-gtkCheckButtonNewWithLabel(label='Reaction Diffusion Termination  ')
	checkgp$t4$SetActive(inp[14])
	checkgp$t5<-gtkCheckButtonNewWithLabel(label='Variable Initiator Efficiency   ')
	checkgp$t5$SetActive(inp[15])
	checkgp$t6<-gtkCheckButtonNewWithLabel(label='Thermal Initiation              ')
	checkgp$t6$SetActive(inp[16])
	checkgp$t6$SetInconsistent(TRUE)
	checkgp$t7<-gtkCheckButtonNewWithLabel(label='pH dependent reaction           ')
	checkgp$t7$SetActive(inp[17])
	checkgp$t8<-gtkCheckButtonNewWithLabel(label='Gas Phase Calculation           ')
	checkgp$t8$SetActive(inp[18])
	checkgp$t8$SetInconsistent(TRUE)
	checkgp$t9<-gtkCheckButtonNewWithLabel(label='Liquid Viscosity Prediction     ')
	checkgp$t9$SetInconsistent(TRUE)
	checkgp$t10<-gtkCheckButtonNewWithLabel(label='Isothermal System              ')
	checkgp$t10$SetActive(inp[20])
	sapply(checkgp,diag_check$packStart)
	submit_hbox<-gtkHBox()
	button_ok<-gtkButton("OK")
	button_canc<-gtkButton("Cancel")
	submit_hbox$packStart(button_ok,expand=TRUE)
	submit_hbox$packEnd(button_canc,expand=TRUE)
	gSignalConnect(button_ok,"clicked",f=function(button_ok){
	button_ok$setData('ans',list(
		c(
			as.numeric(entrygp$e1$getText()),
			as.numeric(entrygp$e2$getText()),
			as.numeric(entrygp$e3$getText()),
			as.numeric(entrygp$e4$getText()),
			as.numeric(entrygp$e5$getText()),
			as.numeric(entrygp$e6$getText()),
			as.numeric(entrygp$e7$getText()),
			as.numeric(entrygp$e8$getText()),
			as.numeric(entrygp$e9$getText()),
			as.numeric(entrygp$e10$getText()),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t1))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t2))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t3))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t4))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t5))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t6))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t7))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t8))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t9))),
			as.numeric(as.logical(gtkToggleButtonGetActive(checkgp$t10)))
		)))
		gtkMainQuit()
		}
	)
	gSignalConnect(button_canc,"clicked",f=function(button_canc){
	button_ok$setData('ans',NULL)
	gtkMainQuit()
	}
	)
	diag_align<-gtkAlignment(xalign=0)
	diag_align$add(diag_label)
	diag_align1<-gtkAlignment(xalign=0)
	diag_align1$add(diag_entry)
	diag_align2<-gtkAlignment(xalign=0)
	diag_align2$add(diag_check)
	table$attach(diag_align,left.attach=1,2,top.attach=1,2,xoptions=c('expand','fill'),yoptions=c('expand','fill'))
	table$attach(diag_align1,left.attach=2,3,top.attach=1,2,xoptions=c('expand','fill'),yoptions=c('expand','fill'))
	table$attach(diag_align2,left.attach=3,4,top.attach=1,2,xoptions=c('expand','fill'),yoptions=c('expand','fill'))
	table$attach(submit_hbox,left.attach=1,4,top.attach=11,12,xoptions=c('expand','fill'),yoptions=c('expand','fill'))
	table$setColSpacing(0,5)
	window<-gtkWindow(show=FALSE)
	window['border-width']<-14
	window$setTitle('Computational Options')
	window$SetDeletable(FALSE)
	window$SetResizable(FALSE)
	window$add(table)
	window$showAll()
	gtkMain()
	ans<-button_ok$GetData('ans')
	window$Destroy()
	return(ans)
}
