pisa.reg.pv <- 
  function(x, pvlabel="READ", by, data, export=FALSE, name= "output", folder=getwd(), std=FALSE) {
    
    intsvy.reg.pv(x=x, pvlabel=pvlabel, by=by, data=data, std=std, export=export, 
                  name= name, folder=folder, config=pisa_conf) 
}