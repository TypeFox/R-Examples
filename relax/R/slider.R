slider<-
function (sl.functions, sl.names, sl.mins, sl.maxs, sl.deltas, 
    sl.defaults, but.functions, but.names, no, set.no.value, 
    obj.name, obj.value, reset.function, title, prompt=FALSE,
    sliders.frame.vertical=TRUE) 
{ # pwolf 080310 080603 100917
  slider.env<-"1"; rm("slider.env")
  if (!exists("slider.env")) slider<-slider.env<<-new.env(parent=.GlobalEnv)
  if (!missing(no)) 
          return(as.numeric(tclvalue(get(paste("slider", no, sep = ""), 
              envir = slider.env))))
  if (!missing(set.no.value)) {
          try(eval(parse(text = paste("tclvalue(slider", set.no.value[1], 
              ")<-", set.no.value[2], sep = "")), envir = slider.env))
          return(set.no.value[2])
  }
  if (!missing(obj.name)) {
          if (!missing(obj.value)) 
              assign(obj.name, obj.value, envir = slider.env)
          else obj.value <- get(obj.name, envir = slider.env)
          return(obj.value)
  }
  if (missing(title))        title<-"slider control widget"
  if (missing(sl.names))     {sl.defaults <- sl.names <- NULL}
  if (missing(sl.functions)) sl.functions <- function(...) {  }
   

  # require(tcltk) # 140306
  nt <- tktoplevel()
  tkwm.title(nt, title)
  tkwm.geometry(nt, "+0+15")
  assign("tktop.slider",nt,envir=slider.env)

  "relax"  
  tkpack(f.slider<-tkframe(nt)) ##vertical
  for (i in seq(sl.names)) {
     "relax"
     eval(parse(text=paste("assign('slider", i, 
                           "',tclVar(sl.defaults[i]),envir=slider.env)",sep="")
     ))
  }
  for (i in seq(along=sl.names)) {
     tkpack(fr <- tkframe(f.slider), 
            side = if (sliders.frame.vertical) "top" else "right") ##vertical
     lab <- tklabel(fr, text = sl.names[i], width = "25")
     sc <- tkscale(fr, from = sl.mins[i], to = sl.maxs[i], 
     showvalue = TRUE, resolution = sl.deltas[i], orient = "horiz")
     tkpack(lab, sc, side=if(sliders.frame.vertical) "right" else "top") ##vertical
     assign("sc", sc, envir = slider.env)
     eval(parse(text=paste("tkconfigure(sc,variable=slider", 
                           i, ")", sep = "")), envir = slider.env)
     sl.fun<-if(length(sl.functions)>1)sl.functions[[i]] else sl.functions
     if(!is.function(sl.fun)) 
       sl.fun <- eval(parse(text = paste("function(...){", sl.fun, "}")))
     if(prompt) tkconfigure(sc,command=sl.fun) else tkbind(sc,"<ButtonRelease>",sl.fun)
  }

  assign("slider.values.old", sl.defaults, envir = slider.env)
  tkpack(f.but <- tkframe(nt), fill = "x")
  tkpack(tkbutton(f.but, text = "Exit", 
         command = function() tkdestroy(nt)), side = "right")
  if(!missing(reset.function)){
    # reset.function <- function(...) print("relax")
    if(!is.function(reset.function)) 
      reset.function <- eval(parse(text = 
          paste("function(...){",reset.function, "}")))
    tkpack(tkbutton(f.but,text="Reset", command = function(){
         for (i in seq(sl.names)) 
           eval(parse(text = 
             paste("tclvalue(slider", 
                   i, ")<-", sl.defaults[i], 
                   sep = "")), envir = slider.env)
         reset.function()
    }), side = "right")
  }

  if (missing(but.names)) but.names <- NULL
  for(i in seq(but.names)) {
    but.fun<-if(length(but.functions)>1)but.functions[[i]] else but.functions
    if(!is.function(but.fun)) 
      but.fun<-eval(parse(text = c("function(...){",but.fun,"}")))
    tkpack(tkbutton(f.but, text = but.names[i], command = but.fun), 
                    side = "left")
  }

  invisible(nt)
}

gslider<-
function (sl.functions, sl.names, sl.mins, sl.maxs, sl.deltas, 
    sl.defaults, but.functions, but.names, no, set.no.value, 
    obj.name, obj.value, reset.function, title, prompt=FALSE,
    sliders.frame.vertical=TRUE, hscale=1, vscale=1,
    pos.of.panel = c("bottom","top","left","right")[1]) 
{ # pwolf 100915 / 121206
  require(tkrplot) 
  slider.env<-"1"; rm("slider.env")
  if (!exists("slider.env")) slider<-slider.env<<-new.env(parent=.GlobalEnv)
  if (!missing(no)) 
          return(as.numeric(tclvalue(get(paste("slider", no, sep = ""), 
              envir = slider.env))))
  if (!missing(set.no.value)) {
          try(eval(parse(text = paste("tclvalue(slider", set.no.value[1], 
              ")<-", set.no.value[2], sep = "")), envir = slider.env))
          return(set.no.value[2])
  }
  if (!missing(obj.name)) {
          if (!missing(obj.value)) 
              assign(obj.name, obj.value, envir = slider.env)
          else obj.value <- get(obj.name, envir = slider.env)
          return(obj.value)
  }
  if (missing(title))        title<-"slider control widget"
  if (missing(sl.names))     {sl.defaults <- sl.names <- NULL}
  if (missing(sl.functions)) sl.functions <- function(...) {  }
   

  # require(tcltk) # 140306
  nt <- tktoplevel()
  tkwm.title(nt, title)
  tkwm.geometry(nt, "+0+15")
  assign("tktop.slider",nt,envir=slider.env)

  "relax"  
  nt.bak <- nt
  sl.frame <- tkframe(nt); gr.frame <- tkframe(nt)
  if( !any(pos.of.panel == c("bottom","top","left","right"))) pos.of.panel <- "bottom"
  tkpack(sl.frame,gr.frame,side=pos.of.panel)
  # gslider start:  
  newpl<-function(...){ plot(0:2,0:2,type="n",xlab="",ylab=""); text(1,1,"dummy plot") }
  img <- tkrplot::tkrplot(gr.frame, newpl, vscale=vscale, hscale=hscale ); tkpack(img,side="top") 
  assign("img",img,envir=slider.env)



  # :gslider end
  ## sliders.frame.vertical ##vertical
  tkpack(f.slider<-tkframe(sl.frame)) ##vertical
  for (i in seq(along=sl.names)) {
     "relax"
     eval(parse(text=paste("assign('slider", i, 
                           "',tclVar(sl.defaults[i]),envir=slider.env)",sep="")
     ))
  }
  # gslider start:
  parent.env<-sys.frame(sys.nframe()-1)
  # :gslider end
  for (i in seq(along=sl.names)) {
     tkpack(fr <- tkframe(f.slider), # 
            side = if( pos.of.panel %in% c("left","right")) "top" else {
                     if (sliders.frame.vertical) "top" else "right" ##vertical
                   })
     lab <- tklabel(fr, text = sl.names[i], width = "25")
     sc <- tkscale(fr, from = sl.mins[i], to = sl.maxs[i], 
                   showvalue = TRUE, resolution = sl.deltas[i], orient = "horiz")
     tkpack(lab, sc, 
            side=  if( pos.of.panel %in% c("left","right")) "top" else {
                     if(sliders.frame.vertical) "right" else "top" ##vertical
                   })
     assign("sc", sc, envir = slider.env)
     eval(parse(text=paste("tkconfigure(sc,variable=slider", 
                           i, ")", sep = "")), envir = slider.env)
     sl.fun<-if(is.list(sl.functions))sl.functions[[min(i,length(sl.functions))]]
             else                     sl.functions
     if(!is.function(sl.fun)) 
       sl.fun <- eval(parse(text = paste("function(...){", sl.fun, "}")))
    # gslider start:
      fname<-paste("tkrrsl.fun",i,sep="")
      eval(parse(text=c(paste(fname, " <-")," function(...){", 
                  "tkrreplot(get('img',envir=slider.env),fun=function()",
                  deparse(sl.fun)[-1],")", "}" )))       
      eval(parse(text=paste("environment(",fname,")<-parent.env")))     
      if (prompt) tkconfigure(sc, command = get(fname))
          else tkbind(sc, "<ButtonRelease>", get(fname))
  }
  if(exists("tkrrsl.fun1")) get("tkrrsl.fun1")() ## gslider only
  # :gslider end




  assign("slider.values.old", sl.defaults, envir = slider.env)
  tkpack(f.but <- tkframe(sl.frame), fill = "x")
  tkpack(tkbutton(f.but, text = "Exit", 
         command = function() tkdestroy(nt)), side = "right")
  if(!missing(reset.function)){
    # reset.function <- function(...) print("relax")
    if(!is.function(reset.function)) 
      reset.function <- eval(parse(text = 
          paste("function(...){",reset.function, "}")))
    fname<-"reset.function"
    idx<-seq(along=sl.names)
    hhh <- paste(sep = "","sl",idx,"<-get('slider",idx,"',envir=slider.env);",
                          "tclvalue(sl",idx,")<-",sl.defaults[idx],"\n")
    eval(parse(text = c(paste(fname, " <-"), " function(...){", hhh,
              "tkrreplot(get('img',envir=slider.env),fun=function()", 
              deparse(reset.function)[-1], ")", "}")))
    eval(parse(text = paste("environment(", fname, ")<-parent.env")))
    tkpack(tkbutton(f.but, text = "Reset", command = get(fname)), side = "right")
  }

  if (missing(but.names)) but.names <- NULL
  for(i in seq(along=but.names)) {
    but.fun<-if(is.list(but.functions)) but.functions[[min(i,length(but.functions))]]
             else                       but.functions

    if(!is.function(but.fun)) but.fun<-eval(parse(text = c("function(...){",but.fun,"}")))
    # gslider start:
    fname<-paste("tkrr.fun",i,sep="")
    eval(parse(text=c(paste(fname, " <-")," function(...){", 
                  "tkrreplot(get('img',envir=slider.env),fun=function()",
                  deparse(but.fun)[-1],")", "}" )))       
    eval(parse(text=paste("environment(",fname,")<-parent.env")))     
    
    tkpack(tkbutton(f.but, text = but.names[i], command = get(fname)), side = "left")
  }
  if(exists("tkrr.fun1")) get("tkrr.fun1")() ## gslider only
  # :gslider end


  invisible(img)
}

