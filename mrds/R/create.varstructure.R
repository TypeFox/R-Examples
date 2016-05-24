#' Creates structures needed to compute abundance and variance
#'
#' Creates samples and obs dataframes used to compute abundance and its
#' variance based on a structure of geographic regions and samples within each
#' region.  The intent is to generalize this routine to work with other
#' sampling structures.
#'
#' The function performs the following tasks: 1)tests to make sure that region
#' labels are unique, 2) merges sample and region tables into a samples table
#' and issue a warning if not all samples were used, 3) if some regions have no
#' samples or if some values of Area were not valid areas given then issue
#' error and stop, then an error is given and the code stops, 4) creates a
#' unique region/sample label in samples and in obs, 5) merges observations
#' with sample and issues a warning if not all observations were used, 6) sorts
#' regions by its label and merges the values with the predictions from the
#' fitted model based on the object number and limits it to the data that is
#' appropriate for the fitted detection function.
#'
#' @param model fitted ddf object
#' @param region region table
#' @param sample sample table
#' @param obs table of object #'s and links to sample and region table
#' @return List with 2 elements: \item{samples }{merged dataframe containing
#'   region and sample info - one record per sample} \item{obs}{merged
#'   observation data and links to region and samples}
#' @note Internal function called by \code{\link{dht}}
#' @author Jeff Laake
#' @keywords utility
create.varstructure <- function(model,region,sample,obs){

  # Test to make sure that region labels are unique
  if(length(unique(region$Region.Label))!=length(region$Region.Label)){
    stop("Region labels must be unique")
  }

  # Merge sample and region tables into samples; warn if not all samples used
  samples <- merge(region,sample,by.x="Region.Label",all.x=TRUE,all.y=TRUE)
  if(any(is.na(samples$Area))){
     warning("Some samples not included in the analysis")
  }

  # Test to make sure that sample labels are unique within region
  if(dim(samples)[1] != dim(unique(data.frame(region=samples$Region.Label,
                                            sample=samples$Sample.Label)))[1]){
    stop("Sample labels must be unique within region")
  }

  # Merge again but don't use all.y which ignores samples not used
  samples <- merge(region,sample,by.x="Region.Label",all.x=TRUE)
  samples <- samples[order(samples$Region.Label,samples$Sample.Label),]

  # If some regions have no samples then issue error and stop; also
  # if invalid areas given then issue error and stop
  if(any(is.na(samples$Sample.Label))){
    stop(paste("Following regions have no samples - ",
        paste(samples$Region.Label[is.na(samples$Sample.Label)],collapse=",")))
  }

  if(any(is.na(region$Area)) | any(!is.numeric(region$Area))){
    stop("Invalid or missing Area values for regions\n")
  }

  # Create a unique region/sample label in samples and in obs
  samples$Label <- paste(samples$Region.Label,samples$Sample.Label,sep="")
  obs$Label <- paste(obs$Region.Label,obs$Sample.Label,sep="")

  # we only want the following columns from obs:
  #  Label, object, Region.Label, Sample.Label
  # so remove everything else so we don't end up with .x and .ys in the 
  # merges that follow...
  obs <- obs[,c("Label", "object", "Region.Label", "Sample.Label")]

  # Merge observations with sample; warn if not all observations used
  data <- merge(obs,samples,by.x="Label",by.y="Label",all.x=TRUE,sort=FALSE)

  if(any(is.na(data$Region.Label.y))){
    warning("Some observations not included in the analysis")
  }

  data <- data[!is.na(data$Region.Label.y),]
  data$Region.Label.y <- NULL
  data$Sample.Label.y <- NULL
  data$Sample.Label <- data$Sample.Label.x
  data$Region.Label <- data$Region.Label.x
  data$Region.Label.x <- NULL
  data$Sample.Label.x <- NULL

  # Sort regions by label
  region <- region[order(region$Region.Label),]

  # Merge with data from model and limit to data appropriate for method
  data <- merge(data,model$data,by.x="object",by.y="object",sort=FALSE)

  # observer =1 to avoid problems with merge; this forces abundance 
  #  to always be estimated using observer 1 as the primary
  obs <- 1
  if(!model$method %in% c("io","io.fi","rem","rem.fi")){
     if(model$method!="ds"){
       data <- data[data$observer==obs,]
       data <- data[data$detected==1,]
     }
  }else{
    data <- data[data$observer==obs,]
  }

  # Return vectors and lists for computation
  return(list(samples=samples,
              obs=data,
              region=region))
}
