piaac.reg <-
function(y, x, by, data, export=FALSE, name= "output", folder=getwd()) { 
  
  intsvy.reg(x=x, y=y, by=by, data=data, export=export,
             name=name, folder=folder, config=piaac_conf)
  
}
