pirls.rho <-
function(variables, by, data, export=FALSE, name= "output", folder=getwd()) {
  
  intsvy.rho(variables=variables, by=by, data=data, export=export, name=name,
             folder=folder, config = pirls_conf)
}
