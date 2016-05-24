pisa.mean.pv <-  
function(pvlabel, by, data, export=FALSE, name= "output", folder=getwd()) {
  
  intsvy.mean.pv(pvnames = paste("PV", 1:5, pvlabel, sep=""), 
                 by=by, data=data, export=export,
                 name=name, folder=folder, config=pisa_conf)
  
}
