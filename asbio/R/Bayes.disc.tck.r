Bayes.disc.tck<-function (){

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
  tkwm.title(tt,"Bayesian analysis -- discrete data and priors")
  prior.entry <- tkentry(tt, textvariable=Prior, width =15)
  data.entry <- tkentry(tt, textvariable=Data1, width =15) 
  data.name.entry <- tkentry(tt, textvariable=Data.name, width = 15) 
  c.data.entry <- tkentry(tt, textvariable=C.data, width = 15) 
  
  done <- tclVar(0)
  show.plot<-tclVar(1)

reset<-function(){
  tclvalue(Prior) <- ""
  tclvalue(Data1) <- ""
  tclvalue(Data.name) <- "Data"
  tclvalue(C.data) <- "seq(1,length(Prior))"
  tclvalue(show.plot)<-"1"
}

reset.but <- tkbutton(tt, text = "Reset", command = reset)
submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)


build <- function() {
  Prior <-parse(text=tclvalue(Prior))[[1]]
  Data <- parse(text=tclvalue(Data1))[[1]]
  data.name <- tclvalue(Data.name)
  show.plot <- as.logical(tclObj(show.plot))
  c.data <-parse(text=tclvalue(C.data))[[1]]

substitute(Bayes.disc(Likelihood=Data,Prior=Prior, 
data.name=data.name,c.data=c.data,plot=show.plot))
}

nc.cbut <- tkcheckbutton(tt, text="Show plot", variable=show.plot)
  tkgrid(tklabel(tt, text = " Bayesian analysis of discrete data "), 
      columnspan = 2)
  tkgrid(tklabel(tt, text = ""))
  tkgrid(tklabel(tt, text = 'Data|\u03b8 '), data.entry)
  tkgrid(tklabel(tt, text = 'Priors; \u03b8 '), prior.entry)
  tkgrid(tklabel(tt, text = "Data name"), data.name.entry)
  tkgrid(tklabel(tt, text = "Prior names"), c.data.entry)
  tkgrid(tklabel(tt, text = ""))
  tkgrid(nc.cbut)
  tkgrid(submit.but, reset.but, sticky ="e")
  
  tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
  tkwait.variable(done)
  if (tclvalue(done) == "2") 
      stop("aborted")
  tkdestroy(tt)
  cmd <- build()
  eval.parent(cmd)
invisible(tclServiceMode(TRUE))
}
Prior <- tclVar("c(1/3, 1/3, 1/3)")
Data1<-tclVar("c(1, 1/2, 0)")
Data.name <- tclVar("data")
C.data <- tclVar("c(1, 2, 3)")
dialog.sd()
})
}
