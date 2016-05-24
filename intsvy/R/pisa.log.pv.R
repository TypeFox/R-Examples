pisa.log.pv <- 
  function(pvlabel="READ", x, by, cutoff, data, export=FALSE, 
           name= "output", folder=getwd()) {
    
  intsvy.log.pv(x=x, pvlabel=pvlabel, cutoff=cutoff, by=by, data=data, export=export, 
                name= name, folder=folder, config=pisa_conf)
  
}