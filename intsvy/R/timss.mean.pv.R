timss.mean.pv <-
  function(pvlabel="BSMMAT", by, data, export=FALSE, name= "output", folder=getwd()) {
    
  intsvy.mean.pv(pvnames = paste(pvlabel, "0", 1:5, sep=""), 
                 by=by, data=data, export=export,
                 name=name, folder=folder, config=timss4_conf)

}
