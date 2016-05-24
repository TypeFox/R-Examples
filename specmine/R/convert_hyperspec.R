## functions to convert hyperspec objects into our own

convert_from_hyperspec = function(hsobj, type = "undefined", description = "") {
  datamatrix = t(hsobj$spc)
  x.values = hyperSpec::wl(hsobj)
  if (!is.null(hsobj$..) && ncol(hsobj$..) >= 1)
    metadata = hsobj$..
  else metadata = NULL
  labels = hyperSpec::labels(hsobj)
  dataset = create_dataset(datamatrix, type = type, metadata = metadata, description = description, 
                            x.axis.values = x.values, label.x = labels$.wavelength, label.values = labels$spc)
  dataset
}

## functions to convert our datasets into hyperspec objects

convert_to_hyperspec = function(dataset) {
  if (is.null(dataset$metadata))
    hyper.object = new("hyperSpec", spc = t(dataset$data), 
                       wavelength = as.numeric(rownames(dataset$data)) )  
  else
    hyper.object = new("hyperSpec", spc = t(dataset$data), 
                       wavelength = as.numeric(rownames(dataset$data)), data= dataset$metadata)
  if (!is.null(dataset$labels)) {
    hyper.object@label$.wavelength = dataset$labels$x
    hyper.object@label$spc = dataset$labels$val
  }
  hyper.object
}