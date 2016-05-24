pirls.mean.pv <-
  function(pvlabel="ASRREA", by, data, export=FALSE, name= "output", folder=getwd()) {
    
    intsvy.mean.pv(pvnames = paste(pvlabel, "0", 1:5, sep=""), 
                   by=by, data=data, export=export,
                   name=name, folder=folder, config=pirls_conf)
    
}
