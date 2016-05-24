# Functions that allow to filter the dataset by some criteria

# SUBSET functions - allow to define criteria for the information to keep

# returns dataset with selected set of samples
# samples - vector with indexes or names of the samples to select
"subset_samples" = function(dataset, samples, rebuild.factors = T) {
  
  dataset$metadata = dataset$metadata[samples,,drop= F]
  if (rebuild.factors) dataset$metadata = rebuild_factors_df(dataset$metadata)
  
  dataset$data = dataset$data[,samples, drop=F]
  dataset
}

# selects set of samples by the value of a metadata variable
"subset_samples_by_metadata_values" = function(dataset, metadata.varname, values)
{
  indexes = which(dataset$metadata[,metadata.varname] %in% values)
  subset_samples(dataset, indexes)
}

# selects a random subset of nsamples from the dataset
# returns a new dataset with the selected samples
"subset_random_samples" = function(dataset, nsamples)
{
  indexes = sample(num_samples(dataset), nsamples)
  subset_samples(dataset, indexes)
}

"subset_x_values" = function(dataset, variables, by.index = FALSE) {
  if (!by.index) {
    variables = as.character(variables)
    indexes = which(rownames(dataset$data) %in% variables)
  }
  else indexes = variables
  dataset$data = dataset$data[indexes,,drop=F]
  dataset
}

"subset_x_values_by_interval" = function(dataset, min.value, max.value)
{
  x.values = get_x_values_as_num(dataset)
  indexes = which(x.values >= min.value & x.values <= max.value)
  subset_x_values(dataset, indexes, by.index = T)
}

"subset_by_samples_and_xvalues" = function(dataset, samples, variables = NULL, by.index = F, 
                                           variable.bounds = NULL, rebuild.factors = T)
{
  if (!by.index) {
    if (is.null(variables)) {
      if (is.null(variable.bounds)) 
        stop("One of variables or variable.bounds parameters needs to be defined")
      else {
        x.values = get_x_values_as_num(dataset)
        x.indexes = which(x.values >= variable.bounds[1] & x.values <= variable.bounds[2])
      }
    }
    else {
      variables = as.character(variables)
      x.indexes = which(rownames(dataset$data) %in% variables)
    }
  }
  else x.indexes = variables
  
  dataset$metadata = dataset$metadata[samples,,drop= F]
  if (rebuild.factors) dataset$metadata = rebuild_factors_df(dataset$metadata)
  
  dataset$data = dataset$data[x.indexes,samples, drop=F]
  dataset
}

"subset_metadata" = function(dataset, variables)
{
  dataset$metadata = dataset$metadata[,variables, drop = F]
  dataset
}

  
# REMOVE functions
# Used to remove SAMPLES, DATA VARIABLES or METADATA

# removes selected samples, data points or metadata variables from a dataset
# type: "sample", "data" or "metadata"
# data.to.remove - defines names or indexes of data to be removed
# by.index - specifies if data.to.remove has indexes or names
# rebuild.factors - if T, in metadata categorical variables are recalculated so that levels that
# disapear are removed (this does not maintain order of the levels)
# xaxis.num - indicates if data.to.remove for data variables is given as numeric vector (if T) or text

"remove_data" = function(dataset, data.to.remove, type = "sample", by.index = F, rebuild.factors = T) {
  if (type == "sample")
    dataset = remove_samples(dataset, data.to.remove, rebuild.factors)
  else if(type == "data")
    dataset = remove_data_variables(dataset, data.to.remove, by.index)
  else if(type == "metadata")
    dataset = remove_metadata_variables(dataset, data.to.remove)
  else stop("Type of data to remove is undefined")
  dataset
}   

"remove_samples" = function(dataset, samples.to.remove, rebuild.factors = T) {
  if (is.numeric(samples.to.remove))
    res = subset_samples(dataset, -samples.to.remove, rebuild.factors = rebuild.factors)
  else {
    indexes.to.remove = which(colnames(dataset$data) %in% samples.to.remove)
    res = subset_samples(dataset, -indexes.to.remove, rebuild.factors = rebuild.factors)
  }
  res
}

"remove_data_variables" = function(dataset, variables.to.remove, by.index = FALSE) {
  if (length(variables.to.remove) == 0) {
    warning("No variables to remove")
    return (dataset)
  }
    
  if (!by.index) {
    variables.to.remove = as.character(variables.to.remove)
    indexes.to.remove = which(rownames(dataset$data) %in% variables.to.remove)
  }
  else indexes.to.remove = variables.to.remove
  subset_x_values(dataset, -indexes.to.remove, by.index = T)
}

"remove_x_values_by_interval" = function(dataset, min.value, max.value)
{
  x.values = get_x_values_as_num(dataset)
  indexes.to.remove = which(x.values >= min.value & x.values <= max.value)
  subset_x_values(dataset, -indexes.to.remove, by.index = T)
}

"remove_metadata_variables" = function(dataset, variables.to.remove)
{
  if (!is.numeric(variables.to.remove))
    indexes.to.remove = which(colnames(dataset$metadata) %in% variables.to.remove)
  else indexes.to.remove = variables.to.remove
  
  if (!is.null(indexes.to.remove) & length(indexes.to.remove) > 0)
    dataset$metadata = dataset$metadata[,-indexes.to.remove, drop = F]
  else warning("No metadata variables removed since no fields matched the criteria")
  dataset
}

# functions to remove samples / variables with NAs

"remove_samples_by_nas" = function(dataset, max.nas = 0, by.percent = F)
{
  if (by.percent== T) max.nas = 100 * max.nas / num_x_values(dataset)
  res = apply(dataset$data, 2, function(x) sum(is.na(x)))
  to.remove = which(res > max.nas)
  remove_samples(dataset, to.remove)
}

"remove_samples_by_na_metadata" = function(dataset, metadata.var)
{
  to.remove = which(is.na(dataset$metadata[,metadata.var]))
  remove_samples(dataset, to.remove)
}

"remove_variables_by_nas" = function(dataset, max.nas = 0, by.percent = F)
{
  if (by.percent== T) max.nas = 100 * max.nas / num_samples(dataset)
  res = apply(dataset$data, 1, function(x) { sum(is.na(x)) } )
  to.remove = which(res > max.nas)
  remove_data_variables(dataset, to.remove, by.index = T)
}

# aggregate samples 
# can be used to merge replicates
# indexes - vector where each index will indicate different groups
# c(1,1,2,2,3,3) would mean that the 6 samples would be aggregated into 3 groups
# with the first with samples 1 and 2, the second with samples 3 and 4, and so on
# aggreg.fn - function used to aggregate data values: "mean", "median", "sum", "max", "min", ...
# metadata variables are handled as data variables if numeric; if factor, the most common value is taken
# meta.to.remove - metadata variables to be removed 

"aggregate_samples" = function(dataset, indexes, aggreg.fn = "mean", meta.to.remove = c()) {
  groups = unique(indexes)
  newdata = matrix(NA, nrow(dataset$data), length(groups))
  rownames(newdata) = rownames(dataset$data)
  colnames(newdata) = vector(mode = "character", length=length(groups) )
  
  newmeta = data.frame()
  for (i in 1:length(dataset$metadata)) {
    if (is.numeric(dataset$metadata[[i]]))
      newmeta[[i]] = vector(mode= "numeric", length=0 )
    else if (is.factor(dataset$metadata[[i]]))
      newmeta[[i]] = factor(levels=levels(dataset$metadata[[i]]) ) 
                            #labels = labels(dataset$metadata[[i]]) )
    else newmeta[[i]] = vector(mode="character", length=length(groups))
  }
  names(newmeta) = names(dataset$metadata)
  
  for (g in 1:length(groups)) {
    to.merge = which(indexes == groups[g])
    newdata[,g] = apply(dataset$data[,to.merge,drop=F], 1, aggreg.fn)
    for (i in 1:length(dataset$metadata)) {
      if (is.numeric(dataset$metadata[[i]])) {
        func= match.fun(aggreg.fn)
        newmeta[g,i] = func(dataset$metadata[to.merge,i])
      }
#        newmeta[g,i] = eval(parse(text=paste(aggreg.fn,"(",dataset$metadata[to.merge,i],")",sep="")))
      else if (is.factor(dataset$metadata[[i]]))
        newmeta[g,i] = levels(newmeta[[i]])[which.max(table(dataset$metadata[to.merge,i]))]
      else 
        newmeta[g,i] = dataset$metadata[to.merge[1],i]
    }
    first.index = to.merge[1]
    colnames(newdata)[g] = colnames(dataset$data)[first.index]
    rownames(newmeta)[g] = colnames(dataset$data)[first.index]
  }

  newdataset = list()
  newdataset$data = newdata
  newdataset$metadata = newmeta
  newdataset$labels = dataset$labels
  newdataset$type = dataset$type
  newdataset$description = dataset$description
  newdataset = remove_metadata_variables(newdataset, meta.to.remove)
  newdataset
}

# Merge data and metadata into a data frame (rows are samples; columns are x.values and metadata)
# samples - if defined, allows to specify which samples will be kept; if undefined, all samples
# x.values - if defined, allows to specify which x.values to keep; if undefined, all will be kept 
# metadata.vars - if defined, allows to specify which metadata to keep; if undefined, all will be kept
"merge_data_metadata" = function(dataset, samples = NULL, metadata.vars = NULL, x.values = NULL, 
                                 by.index = F)
{
  if (!is.null(samples) )
    dataset = subset_samples(dataset, samples, rebuild.factors = T)
  if (!is.null(x.values) )
    dataset = subset_x_values(dataset, x.values, by.index = by.index)
  if (!is.null(metadata.vars))
    dataset = subset_metadata(dataset, metadata.vars)
  
  df = as.data.frame(t(dataset$data))
  df = cbind(df, dataset$metadata)
  df
}
