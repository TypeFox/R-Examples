
# As dBase is an old DOS format, Umlaute  are  stored  using  a  different  code  table
#    (namely ASCII) than most modern unices (namely ANSI).
taxname.abbr <- function(x, hybrid = c('ignore', 'TRUE', 'preserve', 'FALSE', 'substitute'), species = FALSE, cf = FALSE, ...) {
    hybrid <- as.character(hybrid)
    hybrid <- match.arg(hybrid, c('ignore', 'TRUE', 'preserve', 'FALSE', 'substitute'))
  #  loc <- Sys.getlocale(category='LC_CTYPE')
#  Sys.setlocale("LC_ALL","C")
#  print('Executing taxname.abbr ...')
    x <- sub('\ ag[.]', ' agg.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ aggr[.]', ' agg.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ species group[\ ]', ' agg.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ ssp[.]', ' subsp.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ v[.]\ ', ' var. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sv[.]\ ', ' subvar. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ Sec[.]\ ', ' sect. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sect[.][(]bot[.][)]\ ', ' sect. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ Ser[.]\ ', ' ser. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ Subs[.]\ ', ' subsect. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]l[.]', ' s. l.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]str[.]', ' s. str.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]\ str[.]', ' sensustricto', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]\ l[.]', ' sensulato', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]lat[.]', ' sensulato', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.] lat[.]', ' sensulato', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]\ ', ' subsp. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sensustricto', ' s. str.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sensulato', ' s. l.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ f[.]\ ', ' fo. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ species', ' spec.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ nothosubsp[.]' , '\ nothossp.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ nothossp[.]' , '\ nssp.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ nothovar[.]' , '\ nvar.', x, perl=TRUE, useBytes=TRUE)
    if(hybrid %in% c('ignore', 'TRUE')) {
      x <- sub('\ x.' , '\ ', x, perl=TRUE, useBytes=TRUE)
      x <- sub('\ nssp[.]' , '\ ssp.', x, perl=TRUE, useBytes=TRUE)
      x <- sub('\ nvar[.]' , '\ var.', x, perl=TRUE, useBytes=TRUE)
    }
    if(hybrid %in% c('preserve', 'FALSE')) {
    }
    if(hybrid == 'substitute') {
      x <- sub('\ x.' , ' \u00d7', x, perl=TRUE, useBytes=TRUE)
    }

    if(cf) x <- sub('^cf.\ ', '', x, ignore.case=TRUE)
		if(species)  {
      x <- sub('\ sp[.]', ' species', x, perl=TRUE, useBytes=TRUE)
      x <- sub('\ sp$', ' species', x, perl=TRUE, useBytes=TRUE)
      x <- sub('\ spec[.]', ' species', x, perl=TRUE, useBytes=TRUE) 
			} else {
			x <- sub('\ sp[.]' , '', x, perl=TRUE, useBytes=TRUE)
			x <- sub('\ sp$', '', x, perl=TRUE, useBytes=TRUE)
			x <- sub('\ spec[.]' , '', x, perl=TRUE, useBytes=TRUE)
			x <- sub('\ species$' , '', x, perl=TRUE, useBytes=TRUE)
			}
    x <- sub("\\s+$", "", x) # trim trailing leading spaces
    #  Sys.setlocale(category='LC_CTYPE', locale=loc)
   return(x)  
}


taxname.simplify <- function(x, genus=TRUE, epithet=TRUE, hybrid=TRUE, ...) {
#    x <- 'S\U00EBlixae calcarae subsp. holdae'
    x <- gsub('\U00EB', 'e', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('\U00CF', 'i', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('ii', 'i', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('nn', 'n', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('fa', 'pha', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('ph', 'p', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('rh', 'h', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('th', 't', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('tt', 't', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('y', 'i', x, perl=TRUE, useBytes=TRUE)
if(epithet) {
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ae\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('arum\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ea\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ei\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('eos\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ia\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ium\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ius\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('orum\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')

    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('a\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('e\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ens\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('es\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('i\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('is\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('on\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('um\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('us\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ae', 'e', substr(x, regexpr('\ ', x), nchar(x))), sep='')    
}
if(genus) {
    x <-	paste(sub('a$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('as$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('e$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('es$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('eus$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('is$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('on$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('u$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('um$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('us$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
}
if(hybrid) {
	x <- gsub(' x ', ' ', x)
	x <- gsub('\U00D7', '', x)
}
return(x)
}
# taxname.simplify(x, Gattungsendung=TRUE, Artendung=TRUE)


TCS.replace <- function(x) {
  x <- toupper(x)
## Turboveg & ## Florkart Germany (BfN lists)
  x <- replace(x, x %in% c('SPECIES_NR', 'TAXNR', 'NAMNR', 'NAMEID', 'TAXONUSAGEID'), 'TaxonUsageID')
  x <- replace(x, x %in% c('ABBREVIAT','TAXONNAME','TAXON','TAXNAME'), 'TaxonName')
  x <- replace(x, x %in% c('VALID_NR', 'SIPNR', 'SYNNAMEID', 'TAXONCONCEPTID'), 'TaxonConceptID')
  x <- replace(x, x %in% c('VALID_NAME', 'VALIDNAME', 'TAXONCONCEPT'), 'TaxonConcept')
  x <- replace(x, x %in% c('AGG', 'AGGNR', 'NAMEPARENTID', 'ISCHILDTAXONOFID'), 'IsChildTaxonOfID')
  x <- replace(x, x %in% c('AGG_NAME', 'ISCHILDTAXONOF'), 'IsChildTaxonOf')
  x <- replace(x, x %in% c('SECUNDUM', 'ACCORDINGTO'), 'AccordingTo')
  x <- replace(x, x %in% c('NATIVENAME', "COMMONNAME", 'VERNACULARNAME'), 'VernacularName')
  x <- replace(x, x %in% c('CLASSIFICA', 'CLASSIFICATION'), 'Classification')
  x <- replace(x, x %in% c('RANG', 'RANK', "TAXONOMICRANK", 'TAXONRANK'), 'TaxonRank')
  x <- replace(x, x %in% c('author', 'AUTHOR', "Author"), 'NameAuthor')

## ESveg
  x <- replace(x, x %in% c("TAXONCODE"), 'TaxonUsageID')
  x <- replace(x, x %in% c("OBSERVATIONCODE"), "RELEVE_NR")
  x <- replace(x, x %in% c("OBSERVATIONID"), "RELEVE_NR")
  x <- replace(x, x %in% c("STRATUMCODE"), "LAYER")
  x <- replace(x, x %in% c("STRATUM"), "LAYER")
  x <- replace(x, x %in% c("PERCENTAGE_MEAN"), "COVER_PERC")
  x <- replace(x, x %in% c("COVERPERCENT"), "COVER_PERC")

## CDM
  x <- replace(x, x %in% c('TAXON'), 'TaxonName')
  x <- replace(x, x %in% c('SYNUUID'), 'TaxonUsageID')
  x <- replace(x, x %in% c('ACCTAXONID'), 'TaxonConceptID')
  return(x)
}


TV.replace <- function(x) {
  x <- toupper(x)
  ## Turboveg & ## Florkart Germany (BfN lists)
  x <- replace(x, x %in% c('TaxonUsageID', 'TAXNR', 'NAMNR', 'NAMEID', 'TAXONUSAGEID'), 'SPECIES_NR')
  x <- replace(x, x %in% c('TaxonName','TAXONNAME','TAXON','TAXNAME'), 'ABBREVIAT')
  x <- replace(x, x %in% c('TaxonConceptID', 'SIPNR', 'SYNNAMEID', 'TAXONCONCEPTID'), 'VALID_NR')
  x <- replace(x, x %in% c('TaxonConcept', 'VALIDNAME', 'TAXONCONCEPT'), 'VALID_NAME')
  x <- replace(x, x %in% c('IsChildTaxonOfID', 'AGGNR', 'NAMEPARENTID', 'ISCHILDTAXONOFID'), 'AGG')
  x <- replace(x, x %in% c('IsChildTaxonOf', 'ISCHILDTAXONOF'), 'AGG_NAME')
  x <- replace(x, x %in% c('AccordingTo', 'ACCORDINGTO'), 'SECUNDUM')
  x <- replace(x, x %in% c('VernacularName', "COMMONNAME", 'VERNACULARNAME'), 'NATIVENAME')
  x <- replace(x, x %in% c('Classification', 'CLASSIFICATION'), 'CLASSIFICA')
  x <- replace(x, x %in% c('TaxonRank', 'RANK', "TAXONOMICRANK", 'TAXONRANK'), 'RANG')
  x <- replace(x, x %in% c('NameAuthor', 'AUTHOR', "Author"), 'author')
    
  ## CDM
  x <- replace(x, x %in% c('TAXON', 'TaxonName'), 'ABBREVIAT')
  x <- replace(x, x %in% c('SYNUUID', 'TaxonUsageID'), 'SPECIES_NR')
  x <- replace(x, x %in% c('ACCTAXONID', 'TaxonConceptID'), 'VALID_NR')
  return(x)
}

