timss.log.pv <- 
  function(pvlabel="BSMMAT", x, by, cutoff, data, 
           export=FALSE, name= "output", folder=getwd()) {
    
  intsvy.log.pv(x=x, pvlabel=pvlabel, cutoff=cutoff, by=by, data=data, export=export, 
                name= name, folder=folder, config=timss8_conf)
  
}