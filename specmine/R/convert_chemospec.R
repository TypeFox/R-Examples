## functions to convert ChemoSpec objects into our own

convert_from_chemospec = function(csobj, type = "undefined", description = "") {
  datamatrix = t(csobj$data)
  x.values = csobj$freq
  metadata = data.frame(csobj$groups)
  samplenames = csobj$names
  label.x = csobj$unit[1]
  label.val = csobj$unit[2]
  dataset = create_dataset(datamatrix, type = type, metadata = metadata, description = description, 
                           x.axis.values = x.values, sample.names = samplenames, 
                           label.x = label.x, label.values = label.val)
  dataset
}

