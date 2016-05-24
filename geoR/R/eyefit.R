plot.eyefit <- function(x, ...){
  xv <- seq(0, x$max.dist, l=201)
  yv <- max(apply(summary(x)[,c("sigmasq", "tausq")],1,sum))
  plot(xv, yv, type="n", ...)
  lines(x)
  return(invisible())
}

lines.eyefit <- function(x, ...){
  lapply(x, lines.variomodel, ...)
  return(invisible())
}

summary.eyefit <- function(object, ...){
  if(length(object) == 0) return(NULL)
  else{
    res <- data.frame(t(rep(0,7)))
    names(res) <-c("cov.model","sigmasq","phi","tausq", "kappa", "kappa2", "practicalRange")
    for(i in 1:length(object)){
      tmp <- unclass(unlist(object[[i]]))
      if (tmp["cov.model"] != "gneiting.matern") tmp["kappa2"] <- NA
      if (all(tmp["cov.model"] != c("powered.exponential","cauchy","gencauchy","matern")))
        tmp["kappa"] <- NA
      res[i,] <- tmp[c("cov.model","cov.pars1","cov.pars2","nugget", "kappa", "kappa2", "practicalRange")]
    }
  }
  return(res)
}

print.summary.eyefit <- function(x, ...){
  print(x, ...)
  return(invisible())
}

print.eyefit <- function(x, ...){
  print(summary(x), ...)
  return(invisible())
}

eyefit <- function(vario, silent=FALSE){  
  if(!requireNamespace("tcltk", quietly=TRUE)) stop("package tcltk is required to run eyefit()")

  geterrmessage()  
  done <- tclVar(0)
  
  eyefit.env <- new.env()
  assign("eyefit.tmp", list(), envir=eyefit.env)
  ##  .eyefit.tmp <<- list()
  
  dmax <- max(vario$u)
 
  kappa1 <- tclVar("0.5")
  kappa2 <- tclVar("1.0")    
  kernel<- tclVar("exponential")
  mdist <- tclVar(max(vario$u))
  nugget <- tclVar(0.1*(max(vario$v)))
  sill <- tclVar(0.8*(max(vario$v)))	
  range <- tclVar(dmax/3)
  
  replot <- function(...) {        
    k <- as.character(tclObj(kernel))
    kp1 <- as.numeric(tclObj(kappa1))
    kp2 <- as.numeric(tclObj(kappa2))
    maxdist <- as.numeric(tclObj(mdist))
    sigmasq <- as.numeric(tclObj(sill))	
    phi <- as.numeric(tclObj(range))
    tausq <- as.numeric(tclObj(nugget))			
    
    eval(substitute(plot(vario)))    			

    fit <- get("eyefit.tmp", envir=eyefit.env)
    lapply(fit,
           function(x)
           lines.variomodel(seq(0, maxdist, l=100), cov.model=x$cov.model,
                            kappa=x$kappa, cov.pars=x$cov.pars,
                            nug=x$nug, max.dist=x$max.dist))

    if (k == "gneiting.matern" |k == "gencauchy" )
      {
        lines.variomodel(x=seq(0,maxdist,l=100), cov.model=k, kappa=c(kp1,kp2),
                         cov.pars = c(sigmasq,phi),nug=tausq,max.dist=maxdist)
      }
    else if (k == "powered.exponential" || k =="cauchy" || k == "matern")
      {	        
        lines.variomodel(x=seq(0,maxdist,l=100), cov.model=k, kappa=kp1,
                         cov.pars = c(sigmasq,phi),nug=tausq,max.dist=maxdist)
      }
    else
      lines.variomodel(x=seq(0,maxdist,l=100), cov.model=k, cov.pars=c(sigmasq,phi),
                       nug=tausq,max.dist=maxdist)

  }
  
  redraw <- function(...)
    {
      var <- as.character(tclObj(kernel))	
      if (var == "gneiting.matern" | var == "gencauchy") 
	{
          tkconfigure(entry.kappa1, state="normal")
    	  tkconfigure(ts5, state="normal")
          tkfocus(entry.kappa1)
          tkconfigure(entry.kappa2, state="normal")
    	  tkconfigure(ts6, state="normal")
	}
      else if (var == "powered.exponential" || var =="cauchy" || var == "matern") 
	{
          tkconfigure(entry.kappa1,state="normal")
  	  tkconfigure(ts5,state="normal")
          tkfocus(entry.kappa1)
          tkconfigure(entry.kappa2,state="disabled")
    	  tkconfigure(ts6,state="disabled")		
	}
      else
	{
	  tkconfigure(ts5,state="disabled")
  	  tkconfigure(ts6,state="disabled")
	  tkconfigure(entry.kappa1,state="disabled")
          tkconfigure(entry.kappa2,state="disabled")
	}
       replot()    
     }
  
  base <- tktoplevel()
  tkwm.title(base, "Eyefit 1.0")
  
  spec.frm <- tkframe(base,borderwidth=2)
  left.frm <- tkframe(spec.frm)
  right.frm <- tkframe(spec.frm)
  
  frame1 <- tkframe(left.frm, relief="groove", borderwidth=2)
  tkpack(tklabel(frame1,text="Parameters"),fill="both",side="top")    
  
  entry.mdist <-tkentry(frame1,width="4",textvariable=mdist)    
  tkpack(ts1<-tkscale(frame1, label="Max. Distance",command=replot, from=0.00, to=dmax,
                      showvalue=0, variable=mdist,
                      resolution=0.01, orient="horiz",relief="groove"),
         fill="both",expand=1,padx=3,pady=2,ipadx=3,ipady=2,side="left")
  tkpack(entry.mdist,fill="none",side="right")
  
  frame3 <- tkframe(left.frm, relief="groove", borderwidth=2)
  tkpack(tklabel(frame3,text="Cov. Parameters"),fill="both",side="top")    		   
  
  entry.sill <-tkentry(frame3,width="4",textvariable=sill)    
  tkpack(ts2<-tkscale(frame3, label="Sill (sigmasq):",command=replot, from=0.00, to=2*max(vario$v),
                      showvalue=0, variable=sill,
                      resolution=0.01, orient="horiz",relief="groove"),
         fill="none",expand=1,padx=3,pady=2,ipadx=3,ipady=2,side="left")
  tkpack(entry.sill,side="right")
  
  frame4 <- tkframe(left.frm, relief="groove", borderwidth=2)
  entry.range <-tkentry(frame4,width="4",textvariable=range)        
  tkpack(ts3<-tkscale(frame4, label="Range (phi):",command=replot, from=0.00, to=2*dmax,
                      showvalue=1, variable=range,
                      resolution=0.01, orient="horiz",relief="groove"),
         fill="x",expand=1,padx=3,pady=2,ipadx=3,ipady=2,side="left")	    
  tkpack(entry.range,side="right")		    
  
  frame5 <- tkframe(left.frm, relief="groove", borderwidth=2)
  tkpack(tklabel(frame5,text="Nugget"),fill="both",side="top")    		       
  entry.nugget <-tkentry(frame5,width="4",textvariable=nugget)            
  tkpack(ts4<-tkscale(frame5, label="Nugget (tausq):",command=replot,
                      from=0.00, to=2*max(vario$v),
                      showvalue=1, variable=nugget,
                      resolution=0.01, orient="horiz",relief="groove"),
         fill="x",expand=1,padx=3,pady=2,ipadx=3,ipady=2,side="left")
  tkpack(entry.nugget,side="right")
  
  frame6 <- tkframe(left.frm, relief="groove", borderwidth=2)
  tkpack(tklabel(frame6,text="Kappa"),fill="both",side="top")        
  entry.kappa1<-tkentry(frame6,width="4",textvariable=kappa1,state="disabled")    
  tkpack(ts5<-tkscale(frame6, label="Kappa 1:",command=replot,
                      from=0.00, to=10.00,
                      showvalue=1, variable=kappa1,state="disabled",
                      resolution=0.01, orient="horiz",relief="groove"),
         fill="x",expand=1,padx=3,pady=2,ipadx=3,ipady=2,side="left")
  tkpack(entry.kappa1,side="right",fill="none")
  
  
  frame7 <- tkframe(left.frm, relief="groove", borderwidth=2)        
  entry.kappa2 <-tkentry(frame7,width="4",textvariable=kappa2,state="disabled")
  tkpack(ts6<-tkscale(frame7, label="Kappa 2:",command=replot,
                      from=0.00, to=10.00,
                      showvalue=1, variable=kappa2,state="disabled",
                      resolution=0.01, orient="horiz",relief="groove"),
         fill="x",expand=1,padx=3,pady=2,ipadx=3,ipady=2,side="left")
  tkpack(entry.kappa2,side="right",fill="none")
  
  frame2 <- tkframe(right.frm, relief="groove", borderwidth=2)
  tkpack(tklabel(frame2, text="Function"))
  for ( i in c("cauchy","gencauchy","circular","cubic","exponential","gaussian", 
               "gneiting","gneiting.matern","linear","matern","power",
               "powered.exponential","pure.nugget","spherical","wave") ) {        
    
    tmp <- tkradiobutton(frame2, command=redraw,
                         text=i, value=i, variable=kernel)
    tkpack(tmp, anchor="w")        
  }
  
  OnOK <- function(){
    replot()	
  }    
  
  OnQuit <- function(){			
    tclvalue(done)<-2
  }
  
  OnClear <- function(aux=vario){
    assign("eyefit.tmp", list(), envir=eyefit.env)	
    plot(aux)
  }
  
  OnSave <- function(){
    k <- as.character(tclObj(kernel))
    kp1 <- as.numeric(tclObj(kappa1))
    if(k == "gneiting.matern")
      kp2 <-  as.numeric(tclObj(kappa2))
    else kp2 <- NULL
    maxdist <- as.numeric(tclObj(mdist))
    sigmasq <- as.numeric(tclObj(sill))	
    phi <- as.numeric(tclObj(range))
    tausq <- as.numeric(tclObj(nugget))
    aux <- list(cov.model=k, cov.pars=c(sigmasq, phi), 
                nugget=tausq, kappa=c(kp1,kp2), lambda = vario$lambda,
                trend = vario$trend,
                practicalRange = practicalRange(cov.model=k,
                  phi = phi, kappa = kp1),
                max.dist=maxdist)      
    oldClass(aux) <- "variomodel"
    assign("eyefit.tmp", c(get("eyefit.tmp", envir=eyefit.env),list(aux)), envir=eyefit.env)
    
    replot()
  }
  
  tkpack(frame1, frame3, frame4, frame5, frame6, frame7, fill="x")
  tkpack(frame2, fill="x")
  tkpack(left.frm, right.frm,side="left", anchor="n")
  
  c.but <- tkbutton(base,text="Clear",
                    command=function(){OnClear(vario)})

  q.but <- tkbutton(base, text="Quit", command=OnQuit)
  
  save.but <- tkbutton(base, text="Save", command=OnSave)		      
  tkpack(spec.frm)
  tkpack(q.but, side="right")
  tkpack(c.but, side="left")        
  tkpack(save.but, side="right")
  
  replot()
  
  tkbind(entry.kappa1, "<Return>", function(){replot()})
  tkbind(entry.kappa2, "<Return>", function(){replot()})
  tkbind(base, "<Destroy>", function() tclvalue(done)<-2)
  
  tkwait.variable(done)
    
  tkdestroy(base)
  
  if (!silent){
    fit <- get("eyefit.tmp", envir=eyefit.env)
    oldClass(fit) <- "eyefit"
    return(fit)
  }
  else
    return(invisible())
}
