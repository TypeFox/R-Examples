###########################################################
########### demo.Gamlss Demonstration Routine #############
###########################################################
# author:   Konstantinos Pateras
demoDist<-function()
{
# require(gamlss.demo) || stop("gamlss support is absent")
# require(tcltk) || stop("tcltk support is absent")
 if(exists("tmain",envir=.GlobalEnv)) {tkdestroy(tmain)}
 if(exists("tcore",envir=.GlobalEnv)) {tkdestroy(tcore)}
 dists<-c("Normal","Logistic","Normal Vs Logistic","Gumbel","Reverse Gumbel","Exponential Gaussian","Generalised t","Power Exponential","PE & Normal","t family","t Vs Normal","Exponential Generalized Beta type 2","Generalised t","Jonhson's SU","Jonhson's SU original","Normal-Exponential-t","Sinh-arcsinh","Skew exponential power 1","Skew exponential power 2","Skew exponential power 3","Skew exponential power 4","Skew t 1","Skew t 2","Skew t 3","Skew t 4","Skew t 5","Exponential","Gamma","Log Normal","Normal Vs Log Normal","Inverse Gaussian","Weibull","Weibull 2","Weibull 3","Zero adjusted gamma","Zero adjusted inverse Gaussian","Box-Cox t","Box-Cox power exponential","Box Cox Cole and Green","Generalized gamma","Generalized Inverse Gaussian","Generalized Beta type 1","Generalised Beta type 2","Beta","Normal Truncated","Gamma Truncated","Binomial","Beta binomial","Delaporte","Logarithmic","Negative binomial type I","Negative binomial type II","Poisson","Sichel","Poisson Inverse Gaussian","Zero inflated binomial","Zero altered binomial","Beta Inflated","Beta Inflated at zero","Beta Inflated at one","Zero altered Logarithmic","Zero altered Poisson","Zero inflated Poisson","Zero Inflated Poisson 2","Zero altered beta binomial","Zero inflated beta binomial","Zero altered negative binomial","Zero inflated negative binomial","Zero inflated PIG", "Skew Normal type 1", "Skew Normal type 2", "t family type 2", "SST", "Logit Normal", "Log Nornal type 2", "Yule", "Waring","Geometric", "Inverse Gamma", "Parelo type 2", "Parelo type 2 original","Sinh-arcsinh original", "Sinh-arcsinh original type 2", "Double Poisson" )
 
 dists.func<-c("demo.NO","demo.LO","demo.NO.LO","demo.GU","demo.RG","demo.exGAUS","demo.GT","demo.PE","demo.PE.NO","demo.TF","demo.TF.NO","demo.EGB2","demo.GT","demo.JSU","demo.JSUo","demo.NET","demo.SHASH","demo.SEP1","demo.SEP2","demo.SEP3","demo.SEP4","demo.ST1","demo.ST2","demo.ST3","demo.ST4","demo.ST5","demo.EXP","demo.GA","demo.LOGNO","demo.NO.LOGNO","demo.IG","demo.WEI","demo.WEI2","demo.WEI3","demo.ZAGA","demo.ZAIG","demo.BCT","demo.BCPE","demo.BCCG","demo.GG","demo.GIG","demo.GB1","demo.GB2","demo.BE","demo.NOtr","demo.GAtr","demo.BI","demo.BB","demo.DEL","demo.LG","demo.NBI","demo.NBII","demo.PO","demo.SICHEL","demo.PIG","demo.ZIBI","demo.ZABI","demo.BEINF","demo.BEINF0","demo.BEINF1","demo.ZALG","demo.ZAP","demo.ZIP","demo.ZIP2","demo.ZABB","demo.ZIBB","demo.ZANBI","demo.ZINBI","demo.ZIPIG", "demo.SN1", "demo.SN2", "demo.TF2", "demo.SST", "demo.LOGITNO", "demo.LOGNO2", "demo.YULE", "demo.WARING", "demo.GEOM", "demo.IGAMMA", "demo.PARETO2", "demo.PARETO2o", "demo.SHASHo", "demo.SHASHo2", "demo.DPO")
 dists.type<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1, 1, 0,0,0,0,0,1 )
#  , , , , , demo.SHASHo, demo.SHASHo2,  , 
 dists.type<-as.numeric(dists.type)
 dists.all<-as.data.frame(cbind(dists,dists.func,dists.type)) 

  tcore <<- tcore <-tktoplevel() # Core Static Window with type of Distributions.
   tkwm.resizable(tcore, 0, 0)  # Cancels Maximization
   type<-tclVar(0)
   tkwm.title(tcore,"Demonstration of Gamlss.demo")
   tkfocus(tcore)
   frame.base<-tkframe(tcore,relief="ridge",borderwidth=8)
   frame.2nd<-tkframe(frame.base,relief="ridge",borderwidth=8,cursor="hand1")
   tkgrid(tklabel(tcore,text="               Choose one of the following:              ",font = tkfont.create(family="times",size=12,weight="bold")) ,columnspan=2)
   for (pa in c("Continuous Distributions", "  Discrete Distributions", "   Mixed Distributions")) {
     master.radio.box.options<-tkradiobutton(frame.2nd,text=pa,value=pa,variable=type)
     tkpack(master.radio.box.options,anchor="w")
   }
   tclvalue(type)<-"Continuous Distributions"
   tkgrid(frame.2nd,columnspan=2)
   tkgrid(frame.base,columnspan=2)
   tkgrid(tklabel(tcore,text=" "))
   proceedtcore<-function(){
     tkdestroy(tcore)
     tmainwin()
   }
   button.proceed<-tkbutton(tcore,text="Proceed",cursor="hand2",font=tkfont.create(family="times",size=12,weight="bold",slant="italic"),command=proceedtcore)
   tcl("wm", "protocol", tcore, "WM_DELETE_WINDOW")
   button.exit<-tkbutton(tcore,text="EXIT",cursor="X_cursor",font=tkfont.create(family="times",size=12),command=function()tkdestroy(tcore))
   tkgrid(button.exit,button.proceed)
   tkgrid.configure(button.exit,sticky="ews",column=0,row =3)
   tkgrid.configure(button.proceed,sticky="wes",column=1,row =3)
   tmain <-0
   tmainwin<-function(){
    tmain  <<- tmain <- tktoplevel()   # Dynamic Window changes according to selected type.
     tkwm.resizable(tmain, 0, 0)  # Cancels Maximization
     tkwm.title(tmain,"Demonstration of Gamlss.demo")
     object.list.scroll <- tkscrollbar(tmain,bg = "light blue",command=function(...)tkyview(object.list,...))
            object.list <- tklistbox(tmain,bg="light blue",cursor="hand1", font= tkfont.create(family="times",size=13), yscrollcommand=function(...)tkset(object.list.scroll,...))
     m=NULL
     if(tclvalue(type)=="Continuous Distributions"){m=0}
     if(tclvalue(type)=="  Discrete Distributions"){m=1}
     if(tclvalue(type)=="   Mixed Distributions"){m=2} 
     temp1<-as.character(dists.all[,1][dists.all[,3]==m])
     temp2<-as.character(dists.all[,2][dists.all[,3]==m])
     temp3<-dists.all[,3][dists.all[,3]==m]
     temp.all<-as.data.frame(cbind(temp1,temp2,temp3))
     for( i in temp1)
     {
       tkinsert(object.list,"end",i)
     }
     call.dist<-function(){
       tempora.dist=NULL                                            
       tempora.dist<-as.numeric(tkcurselection(object.list))+1
       eval(call(paste(temp.all[tempora.dist,2],sep="")))
     }
     ddist<-function(){
       tkdestroy(tmain)                                            
       demoDist()
     }
     tkselection.set(object.list,0)
     tkgrid(tklabel(tmain,text=paste("           ",tclvalue(type)," :         "),font= tkfont.create(family="times",size=12,weight="bold")),column=1)
     tkgrid(object.list,sticky="wnse", column=1)
     tkgrid.configure(object.list.scroll,sticky="", column=2, row=1)  
     button.ok <- tkbutton(tmain,text="Select!",cursor="hand2",state="normal", font=tkfont.create(family="times", weight="bold",size=13,slant="italic"),command=call.dist)
     tkconfigure(button.ok,cursor="hand2")
     tkconfigure(tmain,cursor="")
      button.cancel <- tkbutton(tmain,text="Back",cursor="sb_left_arrow", font=tkfont.create(family="times", weight="bold",size=13,slant="italic"), command=ddist)
     tkgrid(button.cancel,button.ok)
     tkgrid.configure(button.cancel,sticky="snw",row = 2, column=0, columnspan=2)
     tkgrid.configure(button.ok,sticky="sne", row = 2, column=1, columnspan=2)
     tkgrid.configure(object.list.scroll,sticky="wns")
     
  }
}
#demoDist()
