###########################################################
########### demo.Gamlss Demonstration Routine #############
###########################################################
# author:   Konstantinos Pateras

gamlss.demo<-function()
{

  tmaster <- tktoplevel() # Core Static Window with type of Distributions.
   tkwm.resizable(tmaster, 0, 0)  # Cancels Maximization
   type<-tclVar(0)
   tkwm.title(tmaster,"Demonstration of Gamlss.demo")
   tkfocus(tmaster)
   frame.base<-tkframe(tmaster,relief="ridge",borderwidth=8)
   frame.2nd<-tkframe(frame.base,relief="ridge",borderwidth=8,cursor="hand1")
   tkgrid(tklabel(tmaster,text="               Choose one of the following:              ",font = tkfont.create(family="times",size=12,weight="bold")) ,columnspan=2)
   for (pa in c("Demos for gamlss.family distributions", "Demos for smoothing techniques", "Demos for local polynomial smoothing")) {
     master.radio.box.options<-tkradiobutton(frame.2nd,text=pa,value=pa,variable=type)
     tkpack(master.radio.box.options,anchor="w")
   }
   tclvalue(type)<-"Demos for gamlss.family distributions"
   tkgrid(frame.2nd,columnspan=2)
   tkgrid(frame.base,columnspan=2)
   tkgrid(tklabel(tmaster,text=" "))
   proceedtmaster<-function(){
    if(tclvalue(type)=="Demos for gamlss.family distributions"){demoDist()}
    if(tclvalue(type)=="Demos for smoothing techniques"){demoPsplines()}
    if(tclvalue(type)=="Demos for local polynomial smoothing"){demoLpolyS()}
    tkdestroy(tmaster) 
   }
   button.proceed<-tkbutton(tmaster,text="Proceed",cursor="hand2",font=tkfont.create(family="times",size=12,weight="bold",slant="italic"),command=proceedtmaster)
   button.exit<-tkbutton(tmaster,text="EXIT",cursor="X_cursor",font=tkfont.create(family="times",size=12),command=function()tkdestroy(tmaster))
   tkgrid(button.exit,button.proceed)
   tkgrid.configure(button.exit,sticky="ews",column=0,row =3)
   tkgrid.configure(button.proceed,sticky="wes",column=1,row =3)

} # End of Gamlss.demo

###########################################################
########### demo.Gamlss Lpoly Routine #####################
###########################################################

 demoLpolyS<-function()
 {
  tpoly  <-  tktoplevel() # Core Static Window with type of Psplines.
  tkwm.resizable(tpoly, 0, 0)  # Cancels Maximization                                                                        
  type <- tclVar(0)
  tkwm.title(tpoly,"Demos for local polynomial smoothing")
  tkfocus(tpoly)
  frame.base <- tkframe(tpoly,relief="ridge",borderwidth=8)
   frame.2nd <- tkframe(frame.base,relief="ridge",borderwidth=8,cursor="hand1")
  tkgrid(tklabel(tpoly,text="               Choose one of the following:              ",font = tkfont.create(family="times",size=12,weight="bold")) ,columnspan=2)
  for (pa in c("Local mean", "Weighted Local mean", "Local polynomial", "Weighted Local polynomial", "Demo for Local Fits"))
   {
    master.radio.box.options<-tkradiobutton(frame.2nd,text=pa,value=pa,variable=type)
    tkpack(master.radio.box.options,anchor="w")
   }
  tclvalue(type) <- "Local mean"
  tkgrid(frame.2nd, columnspan=2)
  tkgrid(frame.base,columnspan=2)
  tmpoly<-function(){
   if(tclvalue(type)=="Local mean"){demo.Locmean()}
   if(tclvalue(type)=="Local polynomial"){demo.Locpoly()}
   if(tclvalue(type)=="Weighted Local mean"){demo.WLocmean()} 
   if(tclvalue(type)=="Weighted Local polynomial"){demo.WLocpoly()} 
   if(tclvalue(type)=="Demo for Local Fits"){demo.LocalRegression()} 
  } 
  Dpoly<-function(){
   tkdestroy(tpoly)
   gamlss.demo()
  }                            
  tkgrid(tklabel(tpoly,text=" ")) # Add extra spaces between the buttons and the box of choices.
  button.proceed<-tkbutton(tpoly,text="Proceed",cursor="hand2", font=tkfont.create(family="times",size=12,weight="bold",slant="italic"),command=tmpoly)
  button.exit<-tkbutton(tpoly,text="Back",cursor="sb_left_arrow",font=tkfont.create(family="times",size=12,slant="italic"),command=Dpoly)
  tkgrid(button.exit,button.proceed)
  tkgrid.configure(button.exit,sticky="ews",column=0,row =3)
  tkgrid.configure(button.proceed,sticky="wes",column=1,row =3) 
 } # End of Poly Function
