#' @title Function to return full download records using \code{site}s, \code{dataset}s, or dataset IDs.
#' @description Using the dataset ID, site object or dataset object, return all records associated with the data as a \code{download_list}.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr content GET
#' @importFrom stats na.omit
#' @param x A single numeric dataset ID or a vector of numeric dataset IDs as returned by \code{get_datasets}, or a \code{site}, \code{dataset}, or \code{dataset_list}.
#' @param verbose logical; should messages on API call be printed?
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return This command returns either object of class \code{"try-error"}' (see \code{\link{try}}) definined by the error returned from the Neotoma API call, or an object of class \code{download_list}, containing a set of \code{download} objects, each with relevant assemblage information and metadata:
#' The \code{download} object is a list of lists and data frames that describe an assemblage, the constituent taxa, the chronology, site and PIs who contributed the data. The following are important components:
#'
#'  \item{ \code{dataset} }{A table describing the collection, including dataset information, PI data compatable with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
#'  \item{ \code{sample.meta} }{Dataset information for the core, primarily the age-depth model and chronology.  In cases where multiple age models exist for a single record the most recent chronology is provided here.}
#'  \item{ \code{taxon.list} }{The list of taxa contained within the dataset, unordered, including information that can be used in \code{\link{get_taxa}}}
#'  \item{ \code{counts} }{The assemblage data for the dataset, arranged with each successive depth in rows and the taxa as columns.  All taxa are described in \code{taxon.list}, the chronology is in \code{sample.data}}
#'  \item{ \code{lab.data} }{A data frame of laboratory data, such as exotic pollen spike, amount of sample counted, charcoal counts, etc.}
#'  \item{ \code{chronologies} }{A list of existing chronologies.  If only a single chronology exists for a record then this is the same as the age-model in \code{sample.meta}.}
#'
#' @section Note:
#' The function returns a warning in cases where single taxa are defined by multiple taphonomic characteristics, for example grains that are identified seperately as crumpled and torn in the same sample and sums these values within a sample.
#' In the case that a geochronology dataset is passed to \code{get_download} the function returns a message and a NULL object (that is later excized).  Use \code{get_geochron} for these objects.
#' The chronologies can be augmented using the function \code{get_chroncontrol}, where the individual chronology objects in \code{chronologies} will consist of a table equivalent to \code{sample.meta} and a \code{chroncontrol} object.
#'
#' @examples \dontrun{
#' #  Search for sites with "Pseudotsuga" pollen that are older than 8kyr BP and
#' #  that are roughly within western British Columbia:
#' t8kyr.datasets <- get_dataset(taxonname='*Picea*', loc=c(-90, 41, -89, 44), 
#'                               ageold = 20000, ageyoung=10000)
#'
#' #  Returns 20 records (as of 04/04/2013), get the dataset for all records:
#' pollen.records <- get_download(t8kyr.datasets)
#'
#' #  Standardize the taxonomies for the different records using the WS64 taxonomy.
#' compiled.sites <- compile_taxa(pollen.records, list.name='WS64')
#'
#' #  Extract the Pseudotsuga curves for the sites:
#' get.curve <- function(x, taxa) {
#'                if (taxa %in% colnames(x$counts)) {
#'                  count <- x$counts[,taxa]/rowSums(x$counts, na.rm=TRUE)
#'                } else {
#'                  count <- rep(0, nrow(x$count))
#'                }
#'                data.frame(site = x$dataset$site.data$site.name,
#'                age = x$sample.meta$age,
#'                count = count)
#'              }
#'
#' curves <- do.call(rbind.data.frame,
#'                   lapply(compiled.sites, get.curve, taxa = 'Larix/Pseudotsuga'))
#'
#' #  For illustration, remove the sites with no Pseudotsuga occurance:
#' curves <- curves[curves$count > 0, ]
#'
#' smooth.curve <- predict(loess(sqrt(count)~age, data=curves),
#'                         data.frame(age=seq(20000, 0, by = -100)))
#'                         
#' plot(sqrt(count) ~ age, data = curves,
#'      ylab = '% Pseudotsuga/Larix', xlab='Calibrated Years BP', pch=19,
#'      col=rgb(0.1, 0.1, 0.1, 0.1), xlim=c(0, 20000))
#' lines(seq(20000, 0, by = -100), smooth.curve, lwd=2, lty=2, col=2)
#'
#' #  This figure shows us an apparent peak in Larix/Pseudotsuga pollen in the early-Holocene that
#' #  lends support to a warmer, drier early-Holocene in western North America.
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords IO connection
#' @export
get_download <- function(x, verbose = TRUE) {
  UseMethod('get_download')
}

#' @title Function to return full download records using \code{numeric} dataset IDs.
#' @description Using the dataset ID, return all records associated with the data as a \code{download_list}.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr content GET
#' @param x A single numeric dataset ID or a vector of numeric dataset IDs as returned by \code{get_datasets}.
#' @param verbose logical; should messages on API call be printed?
#' @export
get_download.default <- function(x, verbose = TRUE) {

  # Updated the processing here. There is no need to be fiddling with
  # call. Use missing() to check for presence of argument
  # and then process as per usual

  if (missing(x)) {
      stop(paste0("Either a ",sQuote("dataset id"), " (x) or a dataset object must be provided."))
  }
  if (!missing(x) & !is.numeric(x)) {
          stop('The dataset id (x) must be numeric.')
  }

  get.sample <- function(x) {
    # query Neotoma for data set
    base.uri <- 'http://api.neotomadb.org/v1/data/downloads'
    
    neotoma_content <- httr::content(httr::GET(paste0(base.uri, '/', x)), as = "text")
    if (identical(neotoma_content, "")) stop("")
    aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
    
    # Might as well check here for error and bail
    if (inherits(aa, "try-error"))
        return(aa)

    # if no error continue processing
    if (isTRUE(all.equal(aa[[1]], 0))) {
        stop(paste('Server returned an error message:\n', aa[[2]]),
             call. = FALSE)
    }

    if (isTRUE(all.equal(aa[[1]], 1))) {

      if (aa[[1]] == 1 & length(aa[[2]]) == 0) {
        # Is this the best way to deal with it?
        message(paste0("Dataset ID ", x, " has no associated record in Neotoma."))
        return(NULL)
      }
        aa <- aa[[2]]
        
        rep_NULL <- function(x) { 
          # small function to recursively fill all NULL values with NAs.
          if (is.null(x)) {NA}
          else{
            if (class(x) == 'list') {
              lapply(x, rep_NULL)
            } else {
              return(x)
            }
          }
        }
        
        aa <- lapply(aa, function(x)rep_NULL(x))
        
        if (verbose) {
            message(strwrap(paste0("API call was successful. ",
                                   "Returned record for ",
                                   aa[[1]]$Site$SiteName)))
        }

        # Here the goal is to reduce this list of lists to as
        # simple a set of matrices as possible.
        nams <- names(aa[[1]])
        aa1 <- aa[[1]]

        # If there are actual stratigraphic samples with data
        # in the dataset returned.

        if ('Samples' %in% nams  & length(aa1$Samples) > 0) {

          # Build the metadata for the dataset.
          dataset <- list(
            site.data = data.frame(site.id = aa1$Site$SiteID,
                                   site.name = aa1$Site$SiteName,
                                   long = mean(unlist(aa1$Site[c('LongitudeWest', 'LongitudeEast')]),
                                               na.rm = TRUE),
                                   lat = mean(unlist(aa1$Site[c('LatitudeNorth', 'LatitudeSouth')]),
                                              na.rm = TRUE),
                                   elev = aa1$Site$Altitude,
                                   description = aa1$Site$SiteDescription,
                                   long.acc = abs(aa1$Site$LongitudeWest - aa1$Site$LongitudeEast),
                                   lat.acc = abs(aa1$Site$LatitudeNorth - aa1$Site$LatitudeSouth),
                                   row.names = aa1$Site$SiteName,
                                   stringsAsFactors = FALSE),
            dataset.meta = data.frame(dataset.id = aa1$DatasetID,
                                      dataset.name = aa1$DatasetName,
                                      collection.type = aa1$CollUnitType,
                                      collection.handle = aa1$CollUnitHandle,
                                      dataset.type =  aa1$DatasetType,
                                      stringsAsFactors = FALSE),
            pi.data = do.call(rbind.data.frame,
                                aa1$DatasetPIs),
            submission = data.frame(submission.date = strptime(aa1$NeotomaLastSub,
                                                               '%m/%d/%Y'),
                                    submission.type = 'Last submission to Neotoma',
                                    stringsAsFactors = FALSE),
            access.date = Sys.time())

          # Assign classes
          
          class(dataset) <- c('dataset', 'list')
          class(dataset$site) <- c('site', 'data.frame')
          
          #  Geochronological datasets behave differently than any other dataset.
          #  It's frustrating.  This is the only way we can figure it out in a
          #  general way.
          if (dataset$dataset.meta$dataset.type == 'geochronologic') {
            
            message(paste0('The dataset ID ', dataset$dataset.meta$dataset.id,
                           ' is associated with a geochronology object, not count data.'))
            return(NULL)

          } else {
            
            # copy to make indexing below easier?
            samples <- aa1$Samples

            # Build the metadata for each sample in the dataset.
            sample.meta <- do.call(rbind.data.frame,
                                   lapply(samples, `[`,
                                          c("AnalysisUnitDepth",
                                            "AnalysisUnitThickness",
                                            "SampleID", "AnalysisUnitName"
                                            )))

            # Sample age data
            # Not all depths have the same number of chronologies,
            # which is a bit annoying, and actually more
            # complicated than I originally thought.

            # There may be multiple chronologies associated with each record,
            # and not all chronologies will cover the entire core,
            # which makes things frustrating and difficult.

            # First, get all unique chronology names.
            # Some cores don't have age models, so we use a try function.
            
            # There is an issue in about 31 records in the US database (where I tested)
            # where there are samples without named sample ages, but not a second chronology.
            # The issue is that there are undated samples within a record that has dated records.
            # So some samples are NAPD1 and then the undated samples, deeper in the core are just 'NA'
            # This means that NA is just a continuous flow from NAPD1.
            
            # In records with two chronologies we get an array with two rows.
            # If there are unidentified samples within the chronology they aren't named.
            # This becomes problematic when we have `n` chronologies, plus samples without
            # named chronologies.
            
            chron_list <- lapply(samples, '[[', 'SampleAges')
            
            chron_names <- unique(unlist(sapply(chron_list, function(x)unique(sapply(x, '[[', 'ChronologyName')))))
            chron_lengths <- sapply(chron_list, length)
            # When there is an NA aged sample in a record it fucks everything up
            # because it's not inherently a part of a chronology.
            
            if (any(is.na(chron_names)) & !all(diff(chron_lengths) == 0)) {
              # This implies that many records have multiple chronology
              # coverage, but that some don't & have NA coverage.
              # In that case we want to duplicate the record, and then
              # assign names to the chronologies.
              reassign <- function(x, chron_names) {
                if (length(x) == length(stats::na.omit(chron_names))) {
                  return(x)
                } else {
                  x <- rep(x, max(chron_lengths))
                  return(lapply(1:max(chron_lengths), 
                          function(y) {x[[y]]$ChronologyName <- chron_names[y]; x[[y]]}))
                }
              }
            } else if (all(is.na(chron_names))) {
              chron_names <- 'No chronology'
              chron_list <- lapply(chron_list, function(x) {x[[1]]$ChronologyName <- chron_names; x})
            } else if (any(is.na(chron_names)) & all(diff(chron_lengths) == 0)) {
                # This implies that many records have multiple chronology
                # coverage, but that some don't & have NA coverage.
                # In that case we want to duplicate the record, and then
                # assign names to the chronologies.
                reassign <- function(x, chron_names) {
                  if (length(x) == length(na.omit(chron_names)) & !is.na(x[[1]]$ChronologyName)) {
                    return(x)
                  } else {
                    x <- rep(x, max(chron_lengths))
                    return(lapply(1:max(chron_lengths), 
                                  function(y) {x[[y]]$ChronologyName <- chron_names[y]; x[[y]]}))
                  }
                }
                chron_list <- lapply(chron_list, reassign, chron_names = chron_names)   
            }
            
            chron_vectors <- sapply(chron_list, function(x)unique(sapply(x, '[[', 'ChronologyName')))

            if (is.list(chron_vectors)) {
              #  If it's a list then there are levels without common chronologies
              #  in practice this is more likely for faunal/macrofossil remains.
              unique_vectors <- stats::na.omit(unique(unlist(chron_vectors)))
              
              # Now, we inject an empty age for the missing chronology:
              for (i in unique_vectors) {
                for (j in 1:length(chron_list)) {
                  if (!i %in% sapply(chron_list[[j]], '[[', 'ChronologyName')) {
                    if (any(is.na(sapply(chron_list[[j]], '[[', 'ChronologyName')))) {
                      # Clear any empty records:
                      chron_list[[j]][[which(is.na(sapply(chron_list[[j]], '[[', 'ChronologyName')))]] <- NULL
                    }
                    add_length <- length(chron_list[[j]]) + 1
                    
                    chron_list[[j]][[add_length]] <- list(AgeOlder = NA,
                                                          Age = NA,
                                                          AgeYounger = NA,
                                                          ChronologyName = i,
                                                          AgeType = NA,
                                                          ChronologyID = NA)
                    
                  }
                }
              }
              
              # Now we just need to re-sort the lists:
              chron_list <- lapply(chron_list, function(x){
                lapply(match(sapply(x, '[[', 'ChronologyName'), unique_vectors), function(y)x[[y]])
              })
              
              chron_vectors <- sapply(chron_list, function(x)unique(sapply(x, '[[', 'ChronologyName')))
              
            }
            # This fills in the end of a set of sample ages if there are NAs in a series,
            # but with the fix above it shouldn't be necessary.
            if (!is.null(nrow(chron_vectors))) {
              # There's a problem here that the things aren't in the right order.
              
              chron_vectors <- t(apply(chron_vectors, 1, function(x) {
                                        if (any(is.na(x)) & !all(is.na(x))) {
                                          x[is.na(x)] <- unique(x[!is.na(x)])
                                        }
                                    x } ))
            } else {
              chron_vectors[is.na(chron_vectors)] <- unique(unlist(chron_vectors[!is.na(chron_vectors)]))
            }
            
            chrons <- try(unique(as.vector(unlist(chron_vectors))), silent = TRUE)

            base.frame <- data.frame(age.older = rep(NA, nrow(sample.meta)),
                                     age = rep(NA, nrow(sample.meta)),
                                     age.younger = rep(NA, nrow(sample.meta)),
                                     chronology.name = rep(NA, nrow(sample.meta)),
                                     age.type = rep(NA, nrow(sample.meta)), 
                                     chronology.id = rep(NA, nrow(sample.meta)), 
                                     dataset.id = rep(NA, nrow(sample.meta)))
            
            colnames(base.frame) <- c('age.older', 'age',
                                      'age.younger', 'chronology.name',
                                      'age.type', 'chronology.id', 'dataset.id')

            if (!class(chrons) == 'try-error') {
              # Now we create the chronologies, so long as samples have assigned "SampleAges"
              # If they don't, then we stick in the empty `base.frame` and assign it a name "1"
              # Create the list:
              chron.list <- lapply(1:length(chrons), function(x) base.frame)
              
              if (!(is.null(chrons) | any(is.na(chrons))) & length(chrons) > 0) {
                # If there is a chronology and there are no NAs in it:
                names(chron.list) <- chrons

                for (i in 1:length(samples)) {
                  
                  if (length(chron_list[[i]]) < length(chrons) & length(chron_list[[i]]) == 1) {
                    #  If there's a sample with nothing in it then it gets only a single
                    # chron_list object, which makes the `j` loop choke.
                    #  This current implementation only accounts for a completely
                    # undated sample within the core, but not if one model spans a different length
                    # of the core.
                    
                    for (j in 1:length(chrons)) {
                      # Fix the `j` placeholder at 1
                      chron.list[[j]][i, ] <- data.frame(chron_list[[i]][[1]],
                                                         stringsAsFactors = FALSE)
                      chron.list[[j]]$dataset.id <- dataset$dataset.meta$dataset.id
                    }
                    
                  } else if (!is.null(chrons) & any(is.na(chrons)) & length(chrons) > 0) {
                    # There is a chronology, but one of levels has an NA:
                    for (i in 1:length(chron_list)) {
                      chron_list
                    }
                    
                  } else {
                    for (j in 1:length(chrons)) {
                      # Some of the new datasets are passing data without any chronology information.
                      # Here we're filling in the dataset metadata
                      
                      chron.list[[j]][i, ] <- data.frame(chron_list[[i]][[j]],
                                                         stringsAsFactors = FALSE)
                      chron.list[[j]]$dataset.id <- dataset$dataset.meta$dataset.id
                    }
                  }
                }
              } else {
                chron.list[[ 1 ]][1, ] <-
                  data.frame(samples[[1]]$SampleAges[[1]],
                             stringsAsFactors = FALSE)
                chron.list[[1]]$dataset.id <- dataset$dataset.meta$dataset.id
              }
              
              default_chron <- which(sapply(chron.list, function(x)x$chronology.id[1]) == aa1$DefChronologyID)
              if (length(default_chron) > 1) { default_chron <- default_chron[1] }
              
              if (length(default_chron) == 0) {
                warning(paste0('This dataset has no defined default chronology.  Please use caution.\n',
                        'Using the first chronology, ', names(chron.list)[1],' as the default.'),
                        immediate. = TRUE)
                default_chron <- 1
              }
              
            } else {
              chron.list <- list(base.frame)
              default_chron <- 1
            }

            # sample names - can be NULL hence replace with NA if so
            tmp <- sapply(sample.names <-
                          lapply(samples, `[[`, "SampleUnitName"), is.null)
            sample.names[tmp] <- NA

            # stick all that together, setting names, & reordering cols
            # the most recent age model is provided as the default.
            # 
            sample.meta <- cbind.data.frame(sample.meta,
                                            chron.list[[default_chron]],
                                            unlist(sample.names))
            names(sample.meta) <- c("depth", "thickness",
                                    "sample.id", "unit.name",
                                    colnames(chron.list[[length(chron.list)]]),
                                    "sample.name")

            #  re-ordering the columns so they make sense.
            sample.meta <- sample.meta[, c(1:2, 5:10, 3, 11, 4)]

            # sample data/counts
            # 1) extract each SampleData component & then rbind. Gives a
            # list of data frames
            sample.data <- lapply(lapply(samples, `[[`, "SampleData"),
                                  function(x) do.call(rbind.data.frame, x))
            
            if (length(unlist(sample.data)) == 0) {
              warning('This record contains no count data.  Returning a NULL record.')
              return(NULL)
            }
            
            # 2) How many counts/species in each data frame?
            nsamp <- sapply(sample.data, nrow)
            
            # 3) bind each data frame - result is a data frame in long format
            sample.data <- do.call(rbind, sample.data)
            
            # 4) add a Sample column that is the ID from sample.meta
            sample.data$sample.id <- rep(sample.meta$sample.id, times = nsamp)

            # We're going to isolate the count data and clean it up by
            # excluding lab data and charcoal.  The charcoal exclusion
            # needs some further consideration.
            
            colnames(sample.data) <- c('taxon.name', 'variable.units',
                                       'variable.element', 'variable.context',
                                       'taxon.group', 'value',
                                       'ecological.group', 'sample.id')
            
            # get the table:
            cast_table <- reshape2::dcast(sample.data, 
                                taxon.name + variable.units + variable.element +
                                  variable.context + taxon.group +
                                  ecological.group ~ sample.id, value.var = "value", fun.aggregate = sum)
            
            taxon.list <- cast_table[ ,1:6]
            
            # Now we check for any duplicated names and give them an alias
            # that binds the name with the variable units.
            if (any(duplicated(cast_table$taxon.name))) {
              
              which_dup <- as.character(taxon.list$taxon.name[duplicated(taxon.list$taxon.name)])
              
              # A variable may be duplicated because the same variable was mapped in
              # different units, or because it comes from different contexts.
              # So we want to replace with the appropriate name:
              
              taxon.list$alias <- as.character(taxon.list$taxon.name)
              
              for (i in which_dup) {
                # Choose the column to resolve the problem.
                dup_rows <- data.frame(taxon.list[taxon.list$taxon.name %in% i,],
                                       stringsAsFactors = FALSE)
                
                # The pick order should be:
                # 1. variable.context
                # 2. variable.element - it's clear
                # 3. variable.unit - not sure how often the same element is expressed in different units.
                # 4. context|element
                dup_rows <- dup_rows[,c('variable.context','variable.element',
                                        'variable.units')]
                
                # Remove any columns where all elements are NA:
                dup_rows <- dup_rows[, apply(dup_rows, 2, function(x) { !all(is.na(x)) } )]
                
                #  To resolve the duplication issue we want to find the shortest combination of columns for
                #  which the sum of duplicates is zero:
                if (is.null(nrow(dup_rows))) {
                  # In some cases there are lots of NAs. . . 
                  taxon.list$alias[taxon.list$taxon.name %in% i] <-  paste0(taxon.list$alias[taxon.list$taxon.name %in% i], '|',
                                                                            as.character(dup_rows))
                } else {
                  taxon.list$alias[taxon.list$taxon.name %in% i] <- sapply(1:nrow(dup_rows),
                                                                           function(x) {
                                                                             paste0(taxon.list$alias[taxon.list$taxon.name %in% i][x], '|',
                                                                                             paste0(as.character(t(dup_rows)[,x]), collapse = '|')) })
                }
              }
              
              message <- paste0('\nThere were multiple entries for ',
                                unique(which_dup),
                                sapply(unique(which_dup), function(x)ifelse(length(grep("\\.$", 'abcd', perl = TRUE)) == 1, '.', '')),
                                ' \nget_download has mapped aliases for the taxa in the taxon.list.')
              warning(immediate. = TRUE, message, call. = FALSE)
            }  

            # This pulls out charcoal & lab data to stick in a different data.frame.
            take <- !(taxon.list$taxon.group == "Laboratory analyses" |
                        taxon.list$taxon.group == "Charcoal")

            count.data <- t(cast_table[take, 7:ncol(cast_table)])
            if ('alias' %in% colnames(taxon.list)) {
              colnames(count.data) <- taxon.list$alias[take]
            } else {
              colnames(count.data) <- taxon.list$taxon.name[take]
            }

            # Pull out the lab data and treat it in
            # the same way as the previous:
            
            if (sum(take) > 0) {
              lab.data <- t(cast_table[!take, 7:ncol(cast_table)])
              colnames(lab.data) <- taxon.list$alias[!take]
            } else {
                lab.data <- NULL
            }

            # stick all this together
            aa <- list(dataset = dataset,
                       sample.meta = sample.meta,
                       taxon.list = taxon.list,
                       counts = count.data,
                       lab.data = lab.data,
                       chronologies = chron.list)
          
            class(aa) <- c('download', 'list')
          }
        }
        if ((!'Samples' %in% nams) | length(aa1$Samples) == 0) {
          message('Dataset contains no sample data.')
          return(NULL)
        }
    
    }
    
    aa
    
  }

  aa <- lapply(x, get.sample)

  drop.any <- sapply(aa,is.null)

  if (any(drop.any > 0)) {
    if (all(drop.any > 0)) {
      stop('All datasets return non-download objects\n')
    } else {
      warning('Some datasetids returned empty downloads, be aware that length(datasetid) is longer than the download_list.\n')
      aa <- aa[-(which(sapply(aa,is.null),arr.ind=TRUE))]
    }
  }
  
  names(aa) <- sapply(lapply(lapply(aa, '[[', 'dataset'), '[[', 'dataset.meta'), '[[', 'dataset.id')
  
  class(aa) <- c('download_list', 'list')

  aa
  
}

#' @title Function to return full download records using a \code{dataset}.
#' @description Using a \code{dataset}, return all records associated with the data as a \code{download_list}.
#'
#' @param x An object of class \code{dataset}.
#' @param verbose logical; should messages on API call be printed?
#' @export
get_download.dataset <- function(x, verbose = TRUE) {

  # Updated the processing here. There is no need to be fiddling with
  # call. Use missing() to check for presence of argument
  # and then process as per usual
  
  datasetid <- x$dataset.meta$dataset.id
  
  if (!x$dataset.meta$dataset.type %in% 'geochronologic') {
    aa <- get_download(datasetid, verbose = verbose)
  } else {
    cat('Dataset is a geochronological data object.  Defaulting to get_geochron.\n')
    geochron <- get_geochron(datasetid)
    aa <- list(dataset = x,
               geochronology = geochron)
  }

  aa
}

#' @title Function to return full download records using a \code{dataset_list}.
#' @description Using a \code{dataset_list}, return all records associated with the data as a \code{download_list}.
#'
#' @param x An object of class \code{dataset_list}.
#' @param verbose logical; should messages on API call be printed?
#' @export
get_download.dataset_list <- function(x, verbose = TRUE) {
  
  # Updated the processing here. There is no need to be fiddling with
  # call. Use missing() to check for presence of argument
  # and then process as per usual
  
  datasetid <- sapply(x, FUN=function(x)x$dataset$dataset.id)
  
  aa <- get_download(datasetid, verbose = verbose)
  
  aa
}

#' @title Function to return full download records using a \code{site}.
#' @description Using a \code{site}, return all records associated with the data as a \code{download_list}.
#'
#' @param x An object of class \code{site}.
#' @param verbose logical; should messages on API call be printed?
#' @export
get_download.site <- function(x, verbose = TRUE) {
  
  message('Fetching datasets for the site(s)')
  dataset <- get_dataset(x)
  
  datasetid <- unlist(lapply(dataset, FUN=function(x)x$dataset$dataset.id))
  
  message('Getting downloads:')
  aa <- get_download(datasetid, verbose = verbose)
  
  aa
}
