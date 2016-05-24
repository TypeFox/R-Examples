"sliderv" <-
function(refresh.code,names,minima,maxima,resolutions,starts,
   title="control",no=0, set.no.value=0) {

  if(no!=0)
    return(as.numeric(tcltk::tclvalue(get(paste("slider",no,sep=""),envir=slider.env))))

  if(set.no.value[1]!=0){
    try(eval(parse(text=paste("tcltk::tclvalue(slider",set.no.value[1],")<-",
                     set.no.value[2],sep="")),envir=slider.env));
    return(set.no.value[2])
   }

  if(!exists("slider.env")) slider.env<<-new.env()
  #library(tcltk);
  nt<-tcltk::tktoplevel(); tcltk::tkwm.title(nt,title); tcltk::tkwm.geometry(nt,"+0+0")
  for(i in seq(names))
    eval(parse(text=paste("assign(\"slider",i,"\",tcltk::tclVar(starts[i]),envir=slider.env)",sep="")))
  for(i in seq(names)){
    tcltk::tkpack(fr<-tcltk::tkframe(nt),side='left');
    lab<-tcltk::tklabel(fr, text=names[i], width="1")
    sc<-tcltk::tkscale(fr, command=refresh.code, from=minima[i],
                to=maxima[i], showvalue=T, resolution=resolutions[i])
    assign("sc",sc,envir=slider.env); tcltk::tkpack(lab,sc,side="top")
    eval(parse(text=paste("tcltk::tkconfigure(sc,variable=slider",i,")",sep="")),
         envir=slider.env)
  }

  tcltk::tkpack(fr<-tcltk::tkframe(nt),fill="x")
  tcltk::tkpack(tcltk::tkbutton(fr, text="Exit", command=function()tcltk::tkdestroy(nt)),
         side="right")
  tcltk::tkpack(tcltk::tkbutton(fr, text="Reset",
                  command=function(){ for(i in seq(starts))
                                        eval(parse(text=paste("tcltk::tclvalue(slider",i,")<-",starts[i],sep="")),envir=slider.env)
                                      refresh.code()
                                    } ),side="left")
}

