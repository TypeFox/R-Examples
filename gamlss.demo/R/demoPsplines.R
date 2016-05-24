###########################################################
########### demo.Gamlss Splines Routine ###################
###########################################################
# author :   Konstantinos Pateras                                                           
#library(gamlss.demo)
demoPsplines<-function()
{
# require(gamlss.demo) || stop("gamlss support is absent")
# require(tcltk) || stop("tcltk support is absent")
if(exists("tsplines",envir=.GlobalEnv)) {tkdestroy(tsplines)}
#if(exists("tlist",envir=.GlobalEnv)) {tkdestroy(tlist)}

   tsplines  <-  tktoplevel() # Core Static Window with type of Psplines.
   tkwm.resizable(tsplines, 0, 0)  # Cancels Maximization                                                                        
   type <- tclVar(0)
   tkwm.title(tsplines,"Demo Psplines")
   tkfocus(tsplines)
  frame.base <- tkframe(tsplines,relief="ridge",borderwidth=8)
   frame.2nd <- tkframe(frame.base,relief="ridge",borderwidth=8,cursor="hand1")
   tkgrid(tklabel(tsplines,text="          Choose one of the following:         ",font = tkfont.create(family="times",size=12,weight="bold")) ,columnspan=2)
   for (pa in c("demo.BSplines", "demo.RandomWalk", "demo.histSmo", "demo.interpolateSmo", "demo.PSplines"))
     {
      master.radio.box.options<-tkradiobutton(frame.2nd,text=pa,value=pa,variable=type)
      tkpack(master.radio.box.options,anchor="w")
     }
 tclvalue(type) <- "demo.BSplines"
 tkgrid(frame.2nd, columnspan=2)
 tkgrid(frame.base,columnspan=2)
                     tmsplines<-function(){
  if(tclvalue(type)=="demo.BSplines"){demo.BSplines()}
  if(tclvalue(type)=="demo.RandomWalk"){demo.RandomWalk()}
  if(tclvalue(type)=="demo.histSmo"){demo.histSmo()}
  if(tclvalue(type)=="demo.interpolateSmo"){demo.interpolateSmo()}       
  if(tclvalue(type)=="demo.PSplines"){demo.PSplines()}              
 } 
  Dsplines<-function(){
   tkdestroy(tsplines)
 }                            
 tkgrid(tklabel(tsplines,text=" ")) # Add extra spaces between the buttons and the box of choices.
 button.proceed<-tkbutton(tsplines,text="Proceed",cursor="hand2", font=tkfont.create(family="times",size=12,weight="bold",slant="italic"),command=tmsplines)
 button.exit<-tkbutton(tsplines,text="EXIT",cursor="X_cursor",font=tkfont.create(family="times",size=12),command=Dsplines)
 tkgrid(button.exit,button.proceed)
 tkgrid.configure(button.exit,sticky="ews",column=0,row =3)
 tkgrid.configure(button.proceed,sticky="wes",column=1,row =3)  
}
# demoPsplines()
