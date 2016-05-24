piaac.mean.pv <- 
  function(pvlabel, by, data, export=FALSE, name= "output", folder=getwd()) {

    intsvy.mean.pv(pvnames = paste("PV", pvlabel, 1:10, sep=""), 
                   by=by, data=data, export=export,
                   name=name, folder=folder, config=piaac_conf)

}

