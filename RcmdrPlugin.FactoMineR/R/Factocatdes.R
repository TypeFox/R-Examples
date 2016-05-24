Factocatdes <-
function(){
top<-tktoplevel(borderwidth=10)
tkwm.title(top,.Facto_gettext("Description of categories"))


####### FRAMES -------------------------------------------------------------------------------------------------------------------------------
    descopFrame<-tkframe(top)
    descFrame <- tkframe(descopFrame)
    listdesc<-tklistbox(descFrame,selectmode="extended",exportselection=FALSE,yscrollcommand=function(...) tkset(scr,...))
    scr <- tkscrollbar(descFrame,repeatinterval=5,command=function(...)tkyview(listdesc,...))
   
    factFrame <- tkframe(descopFrame)
    listfact<-tklistbox(factFrame,selectmode="single",height=7,width=20,yscrollcommand=function(...) tkset(scrfact,...))
   
    scrfact <- tkscrollbar(factFrame,repeatinterval=5,command=function(...)tkyview(listfact,...))
    tkgrid(listfact, scrfact,sticky = "nw")
    tkgrid.configure(scrfact, sticky = "wns")
    tkgrid.configure(listfact,sticky = "ew")

    descopFrame2<-tkframe(top)
    descFrame2 <- tkframe(descopFrame2)
    listdesc2<-tklistbox(descFrame2,selectmode="extended",exportselection=FALSE,yscrollcommand=function(...) tkset(scr,...))
    scr2 <- tkscrollbar(descFrame2,repeatinterval=5,command=function(...)tkyview(listdesc,...))

    donnee<-get(getRcmdr(".activeDataSet"))
    nomdonnee<-activeDataSet()
    vars<-colnames(donnee)
    vars.fact = NULL
    vars.desc = NULL
    vars.desc2 = NULL
    for (i in (1:ncol(donnee))){
      if (is.numeric(donnee[,i])){
        tkinsert(listdesc,"end",vars[i])
        vars.desc = c(vars.desc,vars[i])
      }
      else {
        vars.fact = c(vars.fact,vars[i])
        tkinsert(listfact,"end",vars[i])
        vars.desc2 = c(vars.desc2,vars[i])
        tkinsert(listdesc2,"end",vars[i])
      }
    }
    tkgrid(listdesc, scr,sticky = "nw")
    tkgrid.configure(scr, sticky = "wns")
    tkgrid.configure(listdesc,sticky = "ew")
    
    tkgrid(listdesc2, scr2,sticky = "nw")
    tkgrid.configure(scr2, sticky = "wns")
    tkgrid.configure(listdesc2,sticky = "ew")
  optionsFrame <- tkframe(descopFrame,borderwidth=2,relief="ridge")
  
    probab.val<-tclVar("0.05")
    probab<-tkentry(optionsFrame,width=10,textvariable=probab.val)
    probab.lab<-tklabel(optionsFrame,text=.Facto_gettext("P-value:"))
    tkgrid(probab.lab,probab,sticky="w")
   
    resu.val<-tclVar("results")
    resu<-tkentry(optionsFrame,width=10,textvariable=resu.val)
    resu.lab<-tklabel(optionsFrame,text=.Facto_gettext("Name of the result object:"))
    tkgrid(resu.lab,resu,sticky="w")


#################### Fonctions liées au bouton 'Sorties' --------------------------------------------------------------------------------------------------

env<-environment()
if (length(vars.fact)>1){
  Gquali<-TRUE
  Gcateg<-TRUE
}
else{
  Gquali<-FALSE
  Gcateg<-FALSE
}
if (length(vars.desc)>0) Gquanti<-TRUE
else  Gquanti<-FALSE

onSortie<-function(){
    sortiestop<-tktoplevel(borderwidth=10)
    tkwm.title(sortiestop,"Outputs")

    onOKsortie<-function(){
        if(tclvalue(qualiValue)=="1") assign("Gquali", TRUE, envir=env)
        else assign("Gquali", FALSE, envir=env)
        if(tclvalue(categValue)=="1") assign("Gcateg", TRUE, envir=env)
        else assign("Gcateg", FALSE, envir=env)
        if(tclvalue(quantiValue)=="1") assign("Gquanti", TRUE, envir=env)
        else assign("Gquanti", FALSE, envir=env)
        tkdestroy(sortiestop)
    }

    quali.check <- tkcheckbutton(sortiestop)
    if (Gquali) qualiValue <- tclVar("1")
    else qualiValue <- tclVar("0")
    quali.lab<-tklabel(sortiestop,text=.Facto_gettext("Description by qualitative variables"))
    tkconfigure(quali.check,variable=qualiValue)

    categ.check <- tkcheckbutton(sortiestop)
    if (Gcateg) categValue <- tclVar("1")
    else categValue <- tclVar("0")
    tkconfigure(categ.check,variable=categValue)
    categ.lab<-tklabel(sortiestop,text=.Facto_gettext("Description by categories"))

    quanti.check <- tkcheckbutton(sortiestop)
    if (Gquanti) quantiValue <- tclVar("1")
    else quantiValue <- tclVar("0")
    tkconfigure(quanti.check,variable=quantiValue)
    quanti.lab<-tklabel(sortiestop,text=.Facto_gettext("Description by quantitative variables"))

      tkgrid(tklabel(sortiestop,text=.Facto_gettext("Options for the outputs"),fg="red"))
      tkgrid(quali.lab,quali.check,sticky="w")
      tkgrid(categ.lab,categ.check,sticky="w")
      tkgrid(quanti.lab,quanti.check,sticky="w")

    tkgrid(tklabel(sortiestop,text=""))
    tkgrid(tkbutton(sortiestop,text="OK",width=9,command=onOKsortie))
    tkfocus(sortiestop)
}

App<-function(){
    nbitemlist<-c(tclvalue(tkcurselection(listdesc)))
    nbitem<-unlist(strsplit(nbitemlist,"\\ "))

nbitemlist2<-c(tclvalue(tkcurselection(listdesc2)))
nbitem2<-unlist(strsplit(nbitemlist2,"\\ "))

    variable <- vars.desc[as.numeric(tkcurselection(listdesc))+1]
    variable2 <- vars.desc2[as.numeric(tkcurselection(listdesc2))+1]
    resultat<-tclvalue(resu.val)
    probabltat<-tclvalue(probab.val)
    done = 0  
    if (length(vars.fact[as.numeric(tkcurselection(listfact))+1])==0) tkmessageBox(message=.Facto_gettext("No variable selected for the variable to describe"),icon="warning",type="ok")         
    else {
      done=1
      if (length(nbitem2)>0)  {      
        if(any(ind<-grep(vars.fact[as.numeric(tkcurselection(listfact))+1],variable2))) variable2=variable2[-ind]   
        if(length(variable2!=0)) {     
          if (length(nbitem)>0) {  command3=paste(resultat,'=catdes(',nomdonnee,'[,c("', paste(vars.fact[as.numeric(tkcurselection(listfact))+1], collapse='", "'),'", "',   paste(variable, collapse='", "'),'", "', paste(variable2, collapse='", "'),'")] ,num.var=1,proba=',as.numeric(tclvalue(probab.val)),')',sep='')  }
          else{
            if(length(vars.desc>0)){ command3=paste(resultat,'=catdes(',nomdonnee,'[,c("', paste(vars.fact[as.numeric(tkcurselection(listfact))+1], collapse='", "'),'", "',paste(vars.desc, collapse='", "'),'", "', paste(variable2, collapse='", "'),'")] ,num.var=1,proba=',as.numeric(tclvalue(probab.val)),')',sep='')  }
            else{command3=paste(resultat,'=catdes(',nomdonnee,'[,c("', paste(vars.fact[as.numeric(tkcurselection(listfact))+1], collapse='", "'),'", "', paste(variable2, collapse='", "'),'")] ,num.var=1,proba=',as.numeric(tclvalue(probab.val)),')',sep='') }
          }
          justDoIt(command3)
          logger(command3)
        }
        else{
          if (length(nbitem)>0) command3=paste(resultat,'=catdes(',nomdonnee,'[,c("', paste(vars.fact[as.numeric(tkcurselection(listfact))+1], collapse='", "'),'", "',   paste(variable, collapse='", "'),'")] ,num.var=1,proba=',as.numeric(tclvalue(probab.val)),')',sep='')  
          else{
            i=grep(vars.fact[as.numeric(tkcurselection(listfact))+1],vars)
            command3=paste(resultat,'=catdes(',nomdonnee,',num.var=',i,',proba=',as.numeric(tclvalue(probab.val)),')',sep='') 
          }
          justDoIt(command3)
          logger(command3)
        }
      }
      else{
        if(length(vars.desc2)==1){ 
          if (length(nbitem)>0) { command3=paste(resultat,'=catdes(',nomdonnee,'[,c("', paste(vars.fact[as.numeric(tkcurselection(listfact))+1], collapse='", "'),'", "', paste(variable, collapse='", "'),'")] ,num.var=1,proba=',as.numeric(tclvalue(probab.val)),')',sep='')  }
          else{
            i=grep(vars.fact[as.numeric(tkcurselection(listfact))+1],vars)
            command3=paste(resultat,'=catdes(',nomdonnee,',num.var=',i,',proba=',as.numeric(tclvalue(probab.val)),')',sep='')
          }
          justDoIt(command3)
          logger(command3)
        }
        else{          
          if(any(ind<-grep(vars.fact[as.numeric(tkcurselection(listfact))+1],vars.desc2))) vars.desc2=vars.desc2[-ind]  
          if (length(nbitem)>0) command3=paste(resultat,'=catdes(',nomdonnee,'[,c("', paste(vars.fact[as.numeric(tkcurselection(listfact))+1], collapse='", "'),'", "', paste(variable, collapse='", "'),'", "', paste(vars.desc2, collapse='", "'),'")] ,num.var=1,proba=',as.numeric(tclvalue(probab.val)),')',sep='')  
          else{
            i=grep(vars.fact[as.numeric(tkcurselection(listfact))+1],vars)
            command3=paste(resultat,'=catdes(',nomdonnee,',num.var=',i,',proba=',as.numeric(tclvalue(probab.val)),')',sep='')
          }  
          justDoIt(command3)
          logger(command3) 
        }
      }
               
      if (Gquali==1){                
        if ((length(nbitem2)>0) & (length(variable2)>0 ))  {           
          doItAndPrint(paste(resultat,'$test.chi',sep=''))
        }
        if((length(nbitem2)==0) & (length(vars.desc2)>=1))  {       
          doItAndPrint(paste(resultat,'$test.chi',sep=''))
        }
      }
      if ((Gquanti==1)& (length(vars.desc)>0)){
        doItAndPrint(paste(resultat,'$quanti',sep=''))
      }
      if ((Gcateg==1)){
        doItAndPrint(paste(resultat,'$category',sep=''))
      }
    }
        
    return(done)
  }


onOK <- function(){
  done = App()
  if (done >0) tkdestroy(top)
}

sorties<-tkbutton(optionsFrame,text=.Facto_gettext("Outputs"),borderwidth=3,width=12,fg="darkred",command=onSortie)
appel<-tkbutton(top,text=.Facto_gettext("Apply"),borderwidth=3,width=12,fg="blue",command=App)
OKCancelHelp(helpSubject="catdes")

tkgrid(tklabel(top,text=""))
tkgrid(tklabel(top, text = .Facto_gettext("Choose the variable to describe:"), fg = "blue"), columnspan = 2, sticky = "w")
tkgrid(descopFrame,columnspan = 2, sticky="w")
tkgrid(factFrame,tklabel(descopFrame,text="    "),optionsFrame,sticky="w")
tkgrid(sorties,sticky="e")
tkgrid.configure(optionsFrame,sticky="s")
 
tkgrid(tklabel(descopFrame,text=""))
if (length(vars.desc)>0){ 
  tkgrid(tklabel(descopFrame, text = .Facto_gettext("Choose the quantitative variables (by default all)"), fg = "blue"), columnspan = 2, sticky = "w")

  tkgrid(descopFrame,columnspan = 2, sticky="w")
  tkgrid(descFrame,tklabel(descopFrame,text="    "),optionsFrame,sticky="w")
}

if (length(vars.fact)>1){
  tkgrid(tklabel(top,text=""))
  tkgrid(tklabel(top, text = .Facto_gettext("Choose the qualitative variables (by default all)"), fg = "blue"), columnspan = 2, sticky = "w")
  tkgrid(descopFrame2,columnspan = 2, sticky="w")
  tkgrid(descFrame2,tklabel(descopFrame2,text="    "),sticky="w")
}

tkgrid(tklabel(top,text=""))
tkgrid(buttonsFrame,appel)
tkgrid.configure(buttonsFrame, columnspan=2)
tkgrid.configure(appel, column=2)
tkfocus(top)
}
