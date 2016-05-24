#' Taxonomic class methods
#'
#' @name taxon_classes
#' @details
#' The taxonomic classes:
#'
#' \itemize{
#'  \item \strong{binomial} - A binomial, via \code{\link{binomial}}
#'  \item \strong{taxonref} - A single taxonref, via \code{\link{taxonref}}
#'  \item \strong{taxonrefs} - A list of taxonrefs, via \code{\link{taxonrefs}}
#'  \item \strong{grouping} - A grouping (classification) object, via \code{\link{grouping}}
#'  \item \strong{taxon} - A taxon, via \code{\link{taxon}}
#'  \item \strong{taxa} - List of taxon objects, via \code{\link{taxa}}
#' }
NULL

#' A class to represent a taxonomic binomial
#'
#' @export
#' @param genus A genus name
#' @param epithet A specific epithet name
#' @param canonical Canonical name
#' @param species Species, genus plus epithet
#' @param authority Authority name
#' @examples
#' binomial("Poa")
#' binomial("Poa", "annua", authority="L.")
binomial <- function(genus=NULL, epithet=NULL, canonical=NULL, species=NULL, authority=NULL){
  res <- tc(list(genus=genus, epithet=epithet, canonical=canonical, species=species, authority=authority))
  structure(res, class="binomial")
}

#' @export
print.binomial <- function(x, ...){
  cat("<binomial>", sep = "\n")
  cat(paste0("  genus: ", x$genus), sep = "\n")
  cat(paste0("  epithet: ", x$epithet), sep = "\n")
  cat(paste0("  canonical: ", x$canonical), sep = "\n")
  cat(paste0("  species: ", x$species), sep = "\n")
  cat(paste0("  authority: ", x$authority), sep = "\n")
}

#' A class to represent a taxonomic reference
#'
#' @export
#' @param rank (character) Taxonomic rank
#' @param name (character) A name
#' @param id (character,numeric) Identifier
#' @param uri (character) Source of name
#' @examples
#' taxonref("genus", "Poa", 56, "http://scottchamberlain.info/")
#'
#' # many names input
#' splist <- c('Litsea bindoniana', 'Rubus ghanakantae', 'Desmanthus palmeri',
#'  'Leptinella longipes', 'Asarum sakawanum', 'Cistanche compacta',
#'  'Ormosia nanningensis', 'Claoxylon physocarpum', 'Hedycarya arborea',
#'  'Hypnum gracile')
#' lapply(splist, function(x) taxonref("species", name = x))
taxonref <- function(rank="none", name="none", id="none", uri="none"){
  res <- list(rank = rank, name = name, id = as.character(id), uri = uri)
  check_type(res, "character")
  structure(res, class = "taxonref")
}

#' @export
print.taxonref <- function(x, ...){
  cat("<taxonref>", sep = "\n")
  cat(paste0("  rank: ", x$rank), sep = "\n")
  cat(paste0("  name: ", x$name), sep = "\n")
  cat(paste0("  id: ", x$id), sep = "\n")
  cat(paste0("  uri: ", x$uri), sep = "\n")
}

#' A class to represent a list of taxonomic references
#'
#' @export
#' @param ... One or more taxonref objects
#' @examples
#' a <- taxonref("genus", "Poa", 56, "http://scottchamberlain.info/")
#' b <- taxonref("genus", "Quercus", 32343, "http://scottchamberlain.info/")
#' taxonrefs(a, b)
taxonrefs <- function(...){
  res <- list(...)
  check_type(res, "taxonref")
  structure(res, class = "taxonrefs")
}

check_type <- function(x, type){
  resis <- sapply(x, class)
  if (any(resis != type)) {
    stop(sprintf("One or more inputs was not of class %s", type),
         call. = FALSE)
  }
}

check_one_type <- function(x, type){
  if (!is(x, type)) {
    stop(sprintf("One or more inputs was not of class %s", type),
         call. = FALSE)
  }
}

#' A class to represent a taxonomic classification
#'
#' @export
#' @param kingdom A kingdom name
#' @param subkingdom A subkingdom name
#' @param infrakingdom A infrakingdom name
#' @param division A division name
#' @param phylum A phylum name
#' @param subdivision A subdivision name
#' @param infradavision A infradavision name
#' @param superclass A superclass name
#' @param clazz A clazz name
#' @param subclass A subclass name
#' @param infraclass A infraclass name
#' @param superorder A superorder name
#' @param order A order name
#' @param suborder A suborder name
#' @param infraorder A infraorder name
#' @param superfamily A superfamily name
#' @param family A family name
#' @param subfamily A subfamily name
#' @param tribe A tribe name
#' @param subtribe A subtribe name
#' @param genus A genus name
#' @param subgenus A subgenus name
#' @param section A section name
#' @param subsection A subsection name
#' @param species A species name
#' @param subspecies A subspecies name
#' @param variety A variety name
#' @param race A race name
#' @param subvariety A subvariety name
#' @param stirp A stirp name
#' @param morph A morph name
#' @param form A form name
#' @param aberration A aberration name
#' @param subform A subform name
#' @param unspecified A unspecified name
#' @examples
#' grouping(kingdom=taxonref("kingdom", "Animalia"),
#'                species=taxonref("species", "Homo sapiens"))
grouping <- function(kingdom=NULL,subkingdom=NULL,infrakingdom=NULL,division=NULL,
  phylum=NULL,subdivision=NULL,infradavision=NULL,superclass=NULL,clazz=NULL,subclass=NULL,
  infraclass=NULL,superorder=NULL,order=NULL,suborder=NULL,infraorder=NULL,superfamily=NULL,
  family=NULL,subfamily=NULL,tribe=NULL,subtribe=NULL,genus=NULL,subgenus=NULL,section=NULL,
  subsection=NULL,species=NULL,subspecies=NULL,variety=NULL,race=NULL,subvariety=NULL,
  stirp=NULL,morph=NULL,form=NULL,aberration=NULL,subform=NULL,unspecified=NULL){

  res <- tc(list(kingdom=kingdom,subkingdom=subkingdom,infrakingdom=infrakingdom,
    division=division,phylum=phylum,subdivision=subdivision,
    infradavision=infradavision,superclass=superclass,clazz=clazz,
    subclass=subclass,infraclass=infraclass,superorder=superorder,
    order=order,suborder=suborder,infraorder=infraorder,
    superfamily=superfamily,family=family,subfamily=subfamily,
    tribe=tribe,subtribe=subtribe,genus=genus,subgenus=subgenus,
    section=section,subsection=subsection,species=species,subspecies=subspecies,
    variety=variety,race=race,subvariety=subvariety,stirp=stirp,
    morph=morph,form=form,aberration=aberration,subform=subform,unspecified=unspecified))
  check_type(res, "taxonref")
  structure(res, class = "grouping")
}

#' @export
print.grouping <- function(x, ...){
  cat("<grouping>", sep = "\n")
  clss <- x[x != "none"]
  for (i in seq_along(clss)) {
    cat(sprintf("  %s: %s", clss[[i]]$rank, clss[[i]]$name), sep = "\n")
  }
}

#' A class to represent a single taxon
#'
#' @export
#' @param binomial A binomial name
#' @param grouping A grouping object
#' @examples
#' bin <- binomial("Poa", "annua", authority="L.")
#' class <- grouping(kingdom=taxonref("kingdom", "Plantae"),
#'    species=taxonref("family", "Poaceae"))
#' taxon(bin, class)
#'
#' # many names input
#' splist <- c('Litsea bindoniana', 'Rubus ghanakantae', 'Desmanthus palmeri',
#'  'Leptinella longipes', 'Asarum sakawanum', 'Cistanche compacta',
#'  'Ormosia nanningensis', 'Claoxylon physocarpum', 'Hedycarya arborea',
#'  'Hypnum gracile')
#' lapply(splist, binomial)
taxon <- function(binomial, grouping){
  check_one_type(binomial, "binomial")
  check_one_type(grouping, "grouping")
  res <- list(binomial = binomial, grouping = grouping)
  structure(res, class = "taxon")
}

#' @export
print.taxon <- function(x, ...){
  cat("<taxon>", sep = "\n")
  cat(sprintf("  binomial: %s %s", x$binomial$genus, x$binomial$epithet), sep = "\n")
  cat("  grouping: ", sep = "\n")
  clss <- x$grouping[x$grouping != "none"]
  for (i in seq_along(clss)) {
    cat(sprintf("    %s: %s", clss[[i]]$rank, clss[[i]]$name), sep = "\n")
  }
}

#' A class to represent a list of taxa
#'
#' @export
#' @param ... An object of class taxon
#' @examples
#' bin <- binomial("Poa", "annua", authority="L.")
#' class <- grouping(kingdom=taxonref("kingdom", "Plantae"),
#'    species=taxonref("family", "Poaceae"))
#' taxa(taxon(bin, class), taxon(bin, class))
taxa <- function(...){
  res <- list(...)
  for (i in seq_along(res)) {
    if (is(res[[i]], "list")) {
      restmp <- res[i]
      res[[i]] <- NULL
      res <- c(res, unlist(restmp, recursive = FALSE))
    } else {
      res[[i]] <- res[[i]]
    }
  }
  check_type(res, "taxon")
  structure(res, class = "taxa")
}
