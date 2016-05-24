piaac.mean <- 
function(variable, by, data, export=FALSE, name= "output", folder=getwd()) {
  
  intsvy.mean(variable=variable, by=by, data=data, 
              export=export, name=name, folder=folder,
              config=piaac_conf)
}

