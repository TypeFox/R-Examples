see.lma.tck<-function (){

local({
have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
if (have_ttk) {
    tkbutton <- ttkbutton
    tkcheckbutton <- ttkcheckbutton
    tkentry <- ttkentry
    tkframe <- ttkframe
    tklabel <- ttklabel
    tkradiobutton <- ttkradiobutton
}
tclServiceMode(FALSE) 
dialog.sd <- function() {
  tt <- tktoplevel()
  tkwm.title(tt,"Linear models")
  y.entry <- tkentry(tt, textvariable=Y, width =45)
  x0.entry <- tkentry(tt, textvariable=X0, width =45) 
  x1.entry <- tkentry(tt, textvariable=X1, width =45) 
  x2.entry <- tkentry(tt, textvariable=X2, width =45) 
  x3.entry <- tkentry(tt, textvariable=X3, width =45) 
 
  nvars <- tkentry(tt, textvariable=NV, width =5) 
showXY<-tclVar(1)  
  done <- tclVar(0)
 
reset<-function(){
Y<-"c(11,17,16,14,15,12,10,15,19,11,23,20,18,17,27,33,22,26,28)"
X0<-"c(rep(1,5),rep(14,1))"
X1<-"c(rep(0,5),rep(1,5),rep(0,9))"
X2<-"c(rep(0,10),rep(1,4),rep(0,5))"
X3<-"c(rep(0,14),rep(1,5))"

NV<-"4"
  tclvalue(showXY)<-"1"
}

reset.but <- tkbutton(tt, text = "Reset", command = reset)
submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)



build <- function() {
  showXY <- as.logical(tclObj(showXY))
  Y <-parse(text=tclvalue(Y))[[1]]
  NV<-parse(text=tclvalue(NV))[[1]]
  
if(NV == "2"){
X0 <- parse(text=tclvalue(X0))[[1]]
X1 <- parse(text=tclvalue(X1))[[1]]
X<-substitute(cbind(as.numeric(X0),as.numeric(X1)))}  
if(NV == "3"){
X0 <- parse(text=tclvalue(X0))[[1]]
X1 <- parse(text=tclvalue(X1))[[1]]
X2 <- parse(text=tclvalue(X2))[[1]]
X<-substitute(cbind(as.numeric(X0),as.numeric(X1),as.numeric(X2)))}
if(NV =="4"){
X0 <- parse(text=tclvalue(X0))[[1]]
X1 <- parse(text=tclvalue(X1))[[1]]
X2 <- parse(text=tclvalue(X2))[[1]]
X3 <- parse(text=tclvalue(X3))[[1]]
X<-substitute(cbind(as.numeric(X0),as.numeric(X1),as.numeric(X2),as.numeric(X3)))}

substitute(pm1(as.numeric(Y),as.numeric(X),as.numeric(.8),showXY=showXY))
}                

nc.cbut <- tkcheckbutton(tt, text="Show XY", variable=showXY)
  tkgrid(tklabel(tt, text = "Linear model (ANOVA)"), 
      columnspan = 2)
  tkgrid(tklabel(tt, text = ""))
  tkgrid(tklabel(tt, text = 'Y'), y.entry)
  tkgrid(tklabel(tt, text = 'X1'), x0.entry)
  tkgrid(tklabel(tt, text = 'X2'), x1.entry)
  tkgrid(tklabel(tt, text = "X3"), x2.entry)
  tkgrid(tklabel(tt, text = "X4"), x3.entry)
  
  tkgrid(tklabel(tt, text = "Number of factor levels"), nvars)
  tkgrid(tklabel(tt, text = ""))
    tkgrid(nc.cbut)  
  tkgrid(tklabel(tt, text = ""))
  tkgrid(submit.but, reset.but, sticky ="e")
 
  tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
  tkwait.variable(done)
  if (tclvalue(done) == "2") 
      stop("aborted")
  tkdestroy(tt)
  cmd <- build()
  eval.parent(cmd)
}
Y<-tclVar("c(11,17,16,14,15,12,10,15,19,11,23,20,18,17,27,33,22,26,28)")
X0<-tclVar("c(rep(1,5),rep(0,14))")
X1<-tclVar("c(rep(0,5),rep(1,5),rep(0,9))")
X2<-tclVar("c(rep(0,10),rep(1,4),rep(0,5))")
X3<-tclVar("c(rep(0,14),rep(1,5))")

NV<-tclVar("4")
dialog.sd()
})
}    
                                                            