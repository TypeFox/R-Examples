##' Read proxy data from a Tilia TLX format file.
##'
##' @importFrom xml2 xml_attr read_xml as_list xml_text xml_find_one xml_find_all
##' @title Read proxy data from Tilia TLX files
##'
##' @param file a string representing a Tilia TLX format file.
##' @return Return a `download` object.
##'
##' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
##'
##' @export
##' @rdname read.tilia
##'
##' @examples
##' \dontrun{
##' marion <- read.tilia('crystal.tlx')
##'
##' western.cnt <- counts(western.dl)
##' sapply(western.cnt, dim)
##' marion.cnt<- counts(western.dl[[1]])
##' dim(marion.cnt)
##' }
##' 
`read.tilia` <- function(file) {
    tilia_xml <- xml2::read_xml('inst/crystal.tlx')
    tilia_list <- xml2::as_list(tilia_xml)
    
    find_NA <- function(x,y) {
      wrap <- try(xml2::xml_text(xml2::xml_find_one(x,y)), silent=TRUE)
      if (class(wrap) == 'try-error') wrap <- NA
      return(wrap)
    }
    
    # Making the site data:
    tilia_site <- xml2::xml_find_all(tilia_xml, '//Site')
    
    site <- data.frame(site.id     = NA,
                       site.name   = find_NA(tilia_site,'.//SiteName'),
                       long        = mean(as.numeric(find_NA(tilia_site,'.//LongEast')),
                                          as.numeric(find_NA(tilia_site,'.//LongWest')), na.rm=TRUE),
                       lat         = mean(as.numeric(find_NA(tilia_site,'.//LatNorth')), 
                                          as.numeric(find_NA(tilia_site,'.//LatSouth')), na.rm=TRUE),
                       elev        = as.numeric(find_NA(tilia_site,'.//Altitude')),
                       description = as.character(find_NA(tilia_site,'.//Notes')),
                       long.acc    = abs(as.numeric(find_NA(tilia_site,'.//LongEast')) - 
                                              as.numeric(find_NA(tilia_site,'.//LongWest'))),
                       lat.acc     = abs(as.numeric(find_NA(tilia_site,'.//LatNorth')) - 
                                           as.numeric(find_NA(tilia_site,'.//LatSouth'))),
                       row.names   = find_NA(tilia_site,'.//SiteName'),
                       stringsAsFactors = FALSE)
    
    class(site) <- c('site', 'data.frame')
    
    #######################################################################
    # Pull the contact objects:
    contact_nodes <- xml2::xml_find_all(tilia_xml, '//Contact')
    
    contacts <- do.call(rbind.data.frame,lapply(contact_nodes, function(x) {
      data.frame(contact.name     = find_NA(x,'.//FullContactName'),
                 contact.status   = find_NA(x,'.//Status'),
                 family.name      = find_NA(x,'.//FamilyName'),
                 leading.initials = find_NA(x,'.//LeadingInitials'),
                 given.names      = find_NA(x,'.//GivenNames'),
                 suffix           = find_NA(x,'.//Suffix'),
                 title            = find_NA(x,'.//Title'),
                 phone            = find_NA(x,'.//Phone'),
                 fax              = find_NA(x,'.//Fax'),
                 email            = find_NA(x,'.//Email'),
                 url              = find_NA(x,'.//URL'),
                 address          = find_NA(x,'.//Address'),
                 notes            = find_NA(x,'.//Notes'),
                 contact.id       = find_NA(x,'.//NeotomaContactID'),
                 alias.id         = find_NA(x,'.//NeotomaAliasID'))
      }))
  
    contacts <- contacts[!rowSums(is.na(contacts)) == ncol(contacts),]
    class(contacts) <- c('contact', 'data.frame')
    
    #####################################################################
    
    dataset <- list(
      site.data = site,
      dataset.meta = data.frame(dataset.id = NA,
                                dataset.name      = find_NA(tilia_xml, 
                                                            './/Datasets//Dataset//Name'),
                                collection.type   = find_NA(tilia_xml, 
                                                            './/CollectionUnit//CollectionType'),
                                collection.handle = find_NA(tilia_xml,
                                                            './/CollectionUnit//Handle'),
                                dataset.type      = find_NA(tilia_xml, 
                                                            './/Datasets//Dataset//DatasetType'),
                                stringsAsFactors = FALSE),
      pi.data      = data.frame(ContactID = contacts$contact.id,
                                ContactName = contacts$contact.name),
      submission   = data.frame(submission.date = NA,
                                submission.type = NA,
                                stringsAsFactors=FALSE),
      access.date = NA)
    
    class(dataset) <- c('dataset', 'list')
    
    ###############################################
    #
    spreadsheet <- xml2::xml_find_all(tilia_xml,'//SpreadSheet')
    
    empty.frame <- rep(NA, max(sapply(1:length(xml2::xml_children(spreadsheet)), 
                                      function(x)length(xml2::xml_children(xml2::xml_children(spreadsheet)[x]))), na.rm=TRUE))
    
    # Find the "Data" table:
    # x comes in as (xml_children(spreadsheet)[1])
    
    sample_pull <- function(x) {
      cells <- as.numeric(xml2::xml_attr(xml2::xml_children(x), 'row'))
      
      empty.frame[cells] <- xml2::xml_text(xml2::xml_children(x))
      
      empty.frame
    }
    
    all_sample <- do.call(rbind.data.frame, 
                    lapply(1:length(xml2::xml_children(spreadsheet)),
                      function(x)as.character(sample_pull(xml2::xml_children(spreadsheet)[x]))))
    
    all_sample <- apply(all_sample, 2, as.character)
    colnames(all_sample) <- NULL
    
    ##############################################################
    #
    #  Build the sample meta:
    
    chrons <- grep("Chron", all_sample[1,])
    chron_nos <- regexpr('[0-9]',all_sample[1,chrons], perl = TRUE)
    unique_chrons <- substr(all_sample[1,chrons], chron_nos, chron_nos)
    
    # Now we have unique chronology numbers.
    
    sample.meta <- list()
    
    for(i in unique_chrons) {
      
      chron_set <- all_sample[,chrons[unique_chrons == i]]
      
      # Clear the 
      depths <- !is.na(all_sample[,1])
      
      if (!is.null(ncol(chron_set))) {
        # Get age elements if there are multiple elements for the chronology:  
        if (length(grep('old', chron_set[1,], ignore.case = TRUE)) > 0) {
          age_older <- suppressWarnings(as.numeric(gsub('\n', '', 
                                                    chron_set[depths,grep('old', chron_set[1,], 
                                                                      ignore.case = TRUE)])))
        } else { age_older <- rep(NA, sum(depths)) }
        
        if (length(grep('young', chron_set[1,], ignore.case = TRUE)) > 0) {
          age_younger <- suppressWarnings(as.numeric(gsub('\n', '', 
                                                    chron_set[depths,grep('young', chron_set[1,], 
                                                                      ignore.case = TRUE)])))
        } else { age_younger <- rep(NA, sum(depths)) }
        
        if (length(grep(paste0('^\n#Chron',i,'\n$'), chron_set[1,])) > 0) {
          age <- suppressWarnings(as.numeric(gsub('\n', '', 
                                                  chron_set[depths,grep(paste0('^\n#Chron',i,'\n$'), 
                                                                    chron_set[1,], 
                                                                    ignore.case = TRUE)])))  
        } else {
          age <- rep(NA, sum(depths))
        }
          
      } else {
        # There's only one vector of ages, check if they're the age, 
        # older or younger
        if (length(grep(paste0("Chron", i, ".Old"), chron_set[1])) == 0) {
          age_older <- rep(NA, sum(depths))
        } else {
          age_older <- suppressWarnings(as.numeric(chron_set[depths]))
        }
        if (length(grep(paste0("Chron", i, ".Young"), chron_set[1])) == 0) {
          age_younger <- rep(NA, sum(depths))
        } else {
          age_younger <- suppressWarnings(as.numeric(chron_set[depths]))
        }
        if (length(grep(paste0("^\n#Chron", i, "\n$"), chron_set[1])) == 0) {
          age <- rep(NA, sum(depths))
        } else {
          age <- suppressWarnings(as.numeric(chron_set[depths]))
        }
      }
      
      if (is.null(ncol(chron_set))) {
        chron_name <- gsub('\n', '', chron_set[2])
        age_type   <- gsub('\n', '', chron_set[4])
      } else {
        chron_name <- gsub('\n', '', chron_set[2, which.min(nchar(chron_set[2,]))])
        age_type   <- gsub('\n', '', chron_set[4, which.min(nchar(chron_set[2,]))])
      }
      
      sample.meta[[i]] <- data.frame(depth = as.numeric(gsub('\n', '', all_sample[depths,1])),
                                thick = NA,
                                age.older = age_older,
                                age = age,
                                age.younger = age_younger,
                                chronology.name = chron_name,
                                age.type = age_type,
                                chronology.id = NA,
                                dataset.id = NA)
      
    }
    
    # There can be only one sample.meta though.  The default will be (from now on) the
    chronologies <- sample.meta
    sample.meta <- sample.meta[[which.max(names(sample.meta))]]
    # chronology with the highest number in the chron index.
    
    # Here we want to push the chroncontrols for each model into a list:
    
    models <- xml2::xml_find_all(tilia_xml, ".//AgeModel")
    
    get_controls <- function(x) {
      
      #default = which(sapply(tilia_list$AgeModels, function(x)x$Default == "True"))
      controls <- xml2::xml_find_one(x, ".//ChronControls")
      
      controls <- data.frame(age.older = as.numeric(xml2::xml_text(xml2::xml_find_all(controls, ".//AgeLimitOlder"))),
                             age = as.numeric(xml2::xml_text(xml2::xml_find_all(controls, ".//Age"))),
                             age.younger =  as.numeric(xml2::xml_text(xml2::xml_find_all(controls, ".//Age"))),
                             chronology.name = (xml2::xml_text(xml2::xml_find_all(x, ".//ChronologyName"))),
                             age.type = (xml2::xml_text(xml2::xml_find_all(x, ".//AgeUnits"))),
                             chronology.id = NA,
                             dataset.id = NA)
      controls
    }
    
    chron_controls <- lapply(models, get_controls)
    
    #  Everything in the spreadsheet table,
    #  but everything (chrons &cetera) are all together.
    #  First we need to pull out the taxonomy data (it doesn't start with a hash mark)
    
    count_cols <- regexpr('#', all_sample[1,]) < 0 & !is.na(all_sample[1,])
    
    taxon_sub <- all_sample[1:5, count_cols]
    
    taxon_list = data.frame(taxon.name = gsub('\n', '', taxon_sub[2,]),
                            variable.units = gsub('\n', '', taxon_sub[4,]),
                            variable.element = gsub('\n', '', taxon_sub[3,]),
                            variable.context = NA,
                            taxon.group = NA,
                            ecological.group = gsub('\n', '', taxon_sub[5,]))
    
    count_data <- matrix(as.numeric(gsub('\n', '', all_sample[6:nrow(all_sample), count_cols])),
                         ncol = sum(count_cols))
    
    lab_data   <- count_data[, taxon_list$ecological.group == 'LABO']
    count_data <- count_data[, !taxon_list$ecological.group == 'LABO']
    
    aa <- list(dataset        = dataset,
               sample.meta    = sample.meta,
               taxon_list     = taxon_list,
               counts         = count_data,
               lab.data       = lab_data,
               chron_controls = chron_controls,
               chronologies   = chronologies)
    
    class(aa) <- 'download'
    
    aa
}

