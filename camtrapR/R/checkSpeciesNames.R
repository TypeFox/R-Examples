checkSpeciesNames <- function(speciesNames,
                              searchtype,
                              accepted = TRUE,
                              ask = TRUE
){

  # check input
  if(searchtype  %in% c("scientific", "common") == FALSE) stop ("'searchtype' must be 'scientific' or 'common'")
  stopifnot(is.logical(accepted))
  stopifnot(is.character(speciesNames) | is.factor(speciesNames))
  speciesNames <- unique(as.character(speciesNames))

  file.sep <- .Platform$file.sep
  
  # query ITIS TSN (taxnonomic serial number)
  tsns <- get_tsn(searchterm = speciesNames,
                  searchtype = searchtype,
                  verbose    = FALSE,
                  accepted   = accepted,
                  ask        = ask)

  # warning if a name was not found
  if(any(is.na(tsns))){
    not.matched <- which(is.na(tsns))
    warning(paste("found no matches for", length(not.matched), "name(s):\n",  paste(speciesNames[not.matched], collapse = ", ")), immediate. = TRUE, call. = FALSE)
    tsns_worked <- as.tsn(tsns[-not.matched], check = FALSE)
  } else {
    tsns_worked <- tsns
  }

  if(any(!is.na(tsns))){

  scientific <- common <- author <- rankname  <- taxon_status <- data.frame()

  for(i in 1:length(tsns_worked)){
    scientific   <- rbind(scientific,   getscientificnamefromtsn (tsns_worked[i]))       # get scientific names
    common       <- rbind(common,       getcommonnamesfromtsn (tsns_worked[i]))          # get common names
    author       <- rbind(author,       gettaxonauthorshipfromtsn (tsns_worked[i]))      # get author
  	rankname     <- rbind(rankname,     gettaxonomicranknamefromtsn (tsns_worked[i]))    # get rank name
	  if(accepted == FALSE) taxon_status <- rbind(taxon_status, getcoremetadatafromtsn (tsns_worked[i]))         # get taxonomic status
  }

  # if more than 1 common name, condense
  if(any(table(common$tsn) > 1)) {
    common2 <- tapply(common$comname, INDEX = common$tsn, FUN = paste, collapse = file.sep)
    common <- data.frame(comname = common2,
                         tsn     = rownames(common2))
  }

  # make outtable
  dat.out <- data.frame(user_name       = speciesNames,
                        tsn             = as.numeric(tsns),
                        scientific_name = NA,
                        taxon_author    = NA,
                        common_name     = NA,
                        rankname        = NA,
                        taxon_status    = NA,
                        itis_url        = NA)

  dat.out$scientific_name[match(scientific$tsn, dat.out$tsn)] <- scientific$combinedname
  dat.out$taxon_author[match(author$tsn, dat.out$tsn)]        <- author$authorship
  dat.out$common_name[match(common$tsn, dat.out$tsn)]         <- common$comname
  dat.out$rankname[match(rankname$tsn, dat.out$tsn)]          <- rankname$rankname
  dat.out$itis_url[match(tsns_worked, dat.out$tsn)]           <- attributes(tsns_worked)$uri

  if(accepted == FALSE){
    dat.out$taxon_status[match(taxon_status$tsn, dat.out$tsn)]  <- taxon_status$taxonusagerating
  } else {
    dat.out$taxon_status[!is.na(dat.out$tsn)]  <- "accepted"
  }

  return(dat.out)

  } else {stop("found no TSNs for speciesNames")}   # error if no tsns found
}