# setClass("binomial", slots = c(genus="character", epithet="character", canonical="character", species="character", authority="character"))

library("R6")

poa <- binomial("Poa")
poa$tojson()

binomial <- function(genus, epithet=NULL, canonical=NULL, species=NULL, authority=NULL){
  Binomial$new(genus = genus,
               epithet=epithet,
               canonical=canonical,
               species=species,
               authority=authority)
}

Binomial <- R6::R6Class("Binomial",
  public = list(
    genus=NULL,
    epithet=NULL,
    canonical=NULL,
    species=NULL,
    authority=NULL,

    initialize = function(genus = NULL, epithet = NULL, canonical = NULL, species = NULL, authority = NULL){
      self$genus <- genus
      self$epithet <- epithet
      self$canonical <- canonical
      self$species <- species
      self$authority <- authority
    },

    tojson = function(...){
      jsonlite::toJSON(list(genus=self$genus,
                            epithet=self$epithet,
                            canonical=self$canonical,
                            species=self$species,
                            authority=self$authority), auto_unbox = TRUE, ...)
    }
  )
)

Taxonref$new("genus", "Poa", 54, "http://scottchamberlain.info/")

Taxonref <- R6::R6Class("Taxonref",
    public = list(
      rank=NULL,
      name=NULL,
      id=NULL,
      uri=NULL,

      initialize = function(rank = NULL, name = NULL, id = NULL, uri = NULL){
        self$rank <- rank
        self$name <- name
        self$id <- id
        self$uri <- uri
      }
    )
)

Classification <- R6::R6Class("Classification",
  public = list(
    kingdom=NULL,
    subkingdom=NULL,
    infrakingdom=NULL,
    division=NULL,
    phylum=NULL,
    subdivision=NULL,
    infradavision=NULL,
    superclass=NULL,
    clazz=NULL,
    subclass=NULL,
    infraclass=NULL,
    superorder=NULL,
    order=NULL,
    suborder=NULL,
    infraorder=NULL,
    superfamily=NULL,
    family=NULL,
    subfamily=NULL,
    tribe=NULL,
    subtribe=NULL,
    genus=NULL,
    subgenus=NULL,
    section=NULL,
    subsection=NULL,
    species=NULL,
    subspecies=NULL,
    variety=NULL,
    race=NULL,
    subvariety=NULL,
    stirp=NULL,
    morph=NULL,
    form=NULL,
    aberration=NULL,
    subform=NULL,
    unspecified=NULL,

    initialize = function(kingdom=NULL, subkingdom=NULL, infrakingdom=NULL, division=NULL,
        phylum=NULL, subdivision=NULL, infradavision=NULL, superclass=NULL, clazz=NULL,
        subclass=NULL, infraclass=NULL, superorder=NULL, order=NULL, suborder=NULL,
        infraorder=NULL, superfamily=NULL, family=NULL, subfamily=NULL, tribe=NULL,
        subtribe=NULL, genus=NULL, subgenus=NULL, section=NULL, subsection=NULL,
        species=NULL, subspecies=NULL, variety=NULL, race=NULL, subvariety=NULL,
        stirp=NULL, morph=NULL, form=NULL, aberration=NULL, subform=NULL, unspecified=NULL){
      self$kingdom <- kingdom
      self$subkingdom <- subkingdom
      self$infrakingdom <- infrakingdom
      self$division <- division
      self$phylum <- phylum
      self$subdivision <- subdivision
      self$infradavision <- infradavision
      self$superclass <- superclass
      self$clazz <- clazz
      self$subclass <- subclass
      self$infraclass <- infraclass
      self$superorder <- superorder
      self$order <- order
      self$suborder <- suborder
      self$infraorder <- infraorder
      self$superfamily <- superfamily
      self$family <- family
      self$subfamily <- subfamily
      self$tribe <- tribe
      self$subtribe <- subtribe
      self$genus <- genus
      self$subgenus <- subgenus
      self$section <- section
      self$subsection <- subsection
      self$species <- species
      self$subspecies <- subspecies
      self$variety <- variety
      self$race <- race
      self$subvariety <- subvariety
      self$stirp <- stirp
      self$morph <- morph
      self$form <- form
      self$aberration <- aberration
      self$subform <- subform
      self$unspecified <- unspecified
    }
  )
)

Taxon <- R6::R6Class("Taxon",
  public = list(
    binomial = NULL,
    classification = NULL,
    x = NULL,

    initialize = function(binomial = NULL, classification = NULL){
      self$binomial <- binomial
      self$classification <- classification
    },

    hier = function(){
      tmp <- self$classification
      names(tmp)
#       nn <- slotNames(tmp)[-1]
#       vals <- vapply(nn, function(g) slot(tmp, g)@name, "", USE.NAMES = FALSE)
#       data.frame(rank=nn, value=vals, stringsAsFactors = FALSE)
    }
#     ,
#     `[[` = function(x) {
#       self$classification[[x]]
#     }
  )
)

maketax <- function(binomial, classification){
  res <- Taxon$new(binomial = binomial, classification = classification)
  structure(res, class=c("Taxon","stuff"))
}

foo <- function(x) UseMethod("foo")
foo.stuff <- function(x) x

# `[[`.Taxon <- function(x, y) {
#   x$classification[y]
# }
