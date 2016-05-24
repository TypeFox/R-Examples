timss.rho.pv <-
function(variable, pvlabels, by, data, export=FALSE, name= "output", folder=getwd()) {

  intsvy.rho.pv(variable=variable, pvlabels=pvlabels, 
                by=by, data=data, export=export, 
                name= name, folder=folder, config = pirls_conf)
}
