InitErgmm.rsender<-function(model, var=1, var.df=3){
  if (!is.directed(model[["Yg"]]))
    stop("Sender effects are not allowed with an undirected network; use 'sociality'", call.=FALSE)
  model[["sender"]]<-TRUE
  model[["prior"]][["sender.var"]]<-var
  model[["prior"]][["sender.var.df"]]<-var.df
  model
}
InitErgmm.rreceiver<-function(model, var=1, var.df=3){
  if (!is.directed(model[["Yg"]]))
    stop("receiver effects are not allowed with an undirected network; use 'sociality'", call.=FALSE)
  model[["receiver"]]<-TRUE
  model[["prior"]][["receiver.var"]]<-var
  model[["prior"]][["receiver.var.df"]]<-var.df
  model
}
InitErgmm.rsociality<-function(model, var=1, var.df=3){
  model[["sociality"]]<-TRUE
  model[["prior"]][["sociality.var"]]<-var
  model[["prior"]][["sociality.var.df"]]<-var.df
  model
}
