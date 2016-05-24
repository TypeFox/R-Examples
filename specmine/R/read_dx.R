# functions to read JDX spectra files 

"read_dataset_dx" = function(folder.data, filename.meta= NULL, type = "undefined", description = "", 
                              label.x = NULL, label.values = NULL,  
                              header.col.meta = TRUE, header.row.meta = TRUE, sep.meta = ",")
{
  if (!is.null(filename.meta))
    metadata = read_metadata(filename.meta, header.col = header.col.meta, header.row = header.row.meta, sep = sep.meta)
  else metadata = NULL
  
  data.dx = read_data_dx(folder.data)
  
  freqs = data.dx[[1]][,1] # get frequencies from first spectrum
  datamat = matrix(data = NA, nrow = length(freqs), ncol = length(data.dx))

  for (i in 1:length(data.dx)) datamat[,i] = data.dx[[i]][,2]
  
  rownames(datamat) = as.character(freqs)
  colnames(datamat) = names(data.dx)
  
  dataset = create_dataset(datamat, type = type, metadata = metadata, description = description, 
                           label.x = label.x, label.values = label.values)
  dataset
}

"read_data_dx" = function(foldername, debug = F)
{
  filenames = dir(foldername, pattern=".[Dd][Xx]$", full.names=TRUE)
  sampleList = list()
  sampleNames = c()
  snames <- gsub("\\.[^.]*$", "", basename(filenames));
  for (i in 1:length(filenames)) {
    print(paste("Reading sample ", filenames[i]))
    sampleList[[i]] = ChemoSpec::readJDX(filenames[i], debug = debug)
  }
  sampleNames = snames
  names(sampleList) = sampleNames
  sampleList
}

get_samples_names_dx = function(foldername){
  files = list.files(foldername,pattern=".[Dd][Xx]$", recursive = TRUE, full.names= TRUE)
  samples.names = gsub("\\.[^.]*$", "",basename(files))
  samples.names
}