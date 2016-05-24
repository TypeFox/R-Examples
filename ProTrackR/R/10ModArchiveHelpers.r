genre.table <- read.table(text = "
genre,                 genre.id
                          unset,                    unset
                          Alternative,                 48
                          Gothic,                      38
                          Grunge,                     103
                          Metal - Extreme,             37
                          Metal (general),             36
                          Punk,                        35
                          Chiptune,                    54
                          Demo Style,                  55
                          One Hour Compo,              53
                          Chillout,                   106
                          Electronic - Ambient,         2
                          Electronic - Breakbeat,       9
                          Electronic - Dance,           3
                          Electronic - Drum and Bass,   6
                          Electronic - Gabber,         40
                          Electronic - Hardcore,       39
                          Electronic - House,          10
                          Electronic - IDM,            99
                          Electronic - Industrial,     34
                          Electronic - Jungle,         60
                          Electronic - Minimal,       101
                          Electronic - Other,         100
                          Electronic - Progressive,    11
                          Electronic - Rave,           65
                          Electronic - Techno,          7
                          Electronic (general),         1
                          Trance - Acid,               63
                          Trance - Dream,              67
                          Trance - Goa,                66
                          Trance - Hard,               64
                          Trance - Progressive,        85
                          Trance - Tribal,             70
                          Trance (general),            71
                          Big Band,                    74
                          Blues,                       19
                          Jazz - Acid,                 30
                          Jazz - Modern,               31
                          Jazz (general),              29
                          Swing,                       75
                          Bluegrass,                  105
                          Classical,                   20
                          Comedy,                      45
                          Country,                     18
                          Experimental,                46
                          Fantasy,                     52
                          Folk,                        21
                          Fusion,                     102
                          Medieval,                    28
                          New Ages,                    44
                          Orchestral,                  50
                          Other,                       41
                          Piano,                       59
                          Religious,                   49
                          Soundtrack,                  43
                          Spiritual,                   47
                          Video Game,                   8
                          Vocal Montage,               76
                          World,                       42
                          Ballad,                      56
                          Disco,                       58
                          Easy Listening,             107
                          Funk,                        32
                          Pop - Soft,                  62
                          Pop - Synth,                 61
                          Pop (general),               12
                          Rock - Hard,                 14
                          Rock - Soft,                 15
                          Rock (general),              13
                          Christmas,                   72
                          Halloween,                   82
                          Hip-Hop,                     22
                          R and B,                     26
                          Reggae,                      27
                          Ska,                         24
                          Soul,                        25", header = T, sep = ",", stringsAsFactors = F, strip.white = T)

empty.line <- data.frame(mod.id        = integer(0),
                         title         = character(0),
                         filename      = character(0),
                         downloads     = integer(0),
                         favorited     = integer(0),
                         MD5           = character(0),
                         format        = character(0),
                         size          = character(0),
                         genre         = character(0),
                         member.rating = character(0),
                         revwr.rating  = character(0))

#html escape codes

htmlcodes <- read.table(text = "
                        char, code
                        8364, 'euro'
                        32, 'nbsp'
                        34, 'quot'
                        38, 'amp'
                        60, 'lt'
                        62, 'gt'
                        160, 'nbsp'
                        161, 'iexcl'
                        162, 'cent'
                        163, 'pound'
                        164, 'curren'
                        165, 'yen'
                        166, 'brvbar'
                        167, 'sect'
                        168, 'uml'
                        169, 'copy'
                        170, 'ordf'
                        172, 'not'
                        173, 'shy'
                        174, 'reg'
                        175, 'macr'
                        176, 'deg'
                        177, 'plusmn'
                        178, 'sup2'
                        179, 'sup3'
                        180, 'acute'
                        181, 'micro'
                        182, 'para'
                        183, 'middot'
                        184, 'cedil'
                        185, 'sup1'
                        186, 'ordm'
                        187, 'raquo'
                        188, 'frac14'
                        189, 'frac12'
                        190, 'frac34'
                        191, 'iquest'
                        192, 'Agrave'
                        193, 'Aacute'
                        194, 'Acirc'
                        195, 'Atilde'
                        196, 'Auml'
                        197, 'Aring'
                        198, 'AElig'
                        199, 'Ccedil'
                        200, 'Egrave'
                        201, 'Eacute'
                        202, 'Ecirc'
                        203, 'Euml'
                        204, 'Igrave'
                        205, 'Iacute'
                        206, 'Icirc'
                        207, 'Iuml'
                        208, 'ETH'
                        209, 'Ntilde'
                        210, 'Ograve'
                        211, 'Oacute'
                        212, 'Ocirc'
                        213, 'Otilde'
                        214, 'Ouml'
                        215, 'times'
                        216, 'Oslash'
                        217, 'Ugrave'
                        218, 'Uacute'
                        219, 'Ucirc'
                        220, 'Uuml'
                        221, 'Yacute'
                        222, 'THORN'
                        223, 'szlig'
                        224, 'agrave'
                        225, 'aacute'
                        226, 'acirc'
                        227, 'atilde'
                        228, 'auml'
                        229, 'aring'
                        230, 'aelig'
                        231, 'ccedil'
                        232, 'egrave'
                        233, 'eacute'
                        234, 'ecirc'
                        235, 'euml'
                        236, 'igrave'
                        237, 'iacute'
                        238, 'icirc'
                        239, 'iuml'
                        240, 'eth'
                        241, 'ntilde'
                        242, 'ograve'
                        243, 'oacute'
                        244, 'ocirc'
                        245, 'otilde'
                        246, 'ouml'
                        247, 'divide'
                        248, 'oslash'
                        249, 'ugrave'
                        250, 'uacute'
                        251, 'ucirc'
                        252, 'uuml'
                        253, 'yacute'
                        254, 'thorn'
                        ", sep = ",", header = T, quote = "'", strip.white = T)
htmlcodes$char <- intToUtf8(htmlcodes$char, T)
htmlcodes$code <- paste0("&", htmlcodes$code, ";")

htmlcodes <- rbind(htmlcodes, data.frame(
  char = intToUtf8(32:383, T),
  code = sprintf("&#%03i;", 32:383)
))

# An oversimplistic function to unescape html escape codes
.htmlUnescape <- function(text) {
  for (i in 1:nrow(htmlcodes))
    text <- gsub(htmlcodes$code[i], htmlcodes$char[i], text, fixed = T)
  Encoding(text) <- "UTF-8"
  return(text)
}

#' ModArchive helper functions
#'
#' \url{http://ModArchive.org} is the largest online archive of module files. These functions
#' will assist in accessing this archive.
#'
#' The \code{modArchive.info} function will retrieve info on a specific module from the
#' ModArchive. The \code{modArchive.download} will download modules from the archive.
#' The \code{modArchive.search.mod} function will search the archive for modules. Note
#' that the ModArchive also contains file formats other that ProTracker's MOD format.
#' This package will only handle the MOD format.
#'
#' @param mod.id An \code{integer} code used as module identifier in the ModArchive database.
#' A \code{mod.id} can be obtained by performing a search with \code{modArchive.search.mod}.
#' When downloading a module, make sure that the identifier represents a MOD file, as
#' other types will result in an error.
#' @param search.text A \code{character} string to be used as terms to search
#' in the ModArchive.
#' @param search.where A \code{character} string indicating where in the module files
#' to search for the \code{search.text}. See usage section for the available options.
#' @param format.filter File format filter to be used in a search in the ModArchive.
#' See the usage section for all possible options. Default is "unset" (meaning that
#' it will search for any file format). Note that only the `MOD' format
#' is supported by this package.
#' @param size.filter File size filter to be used in a search in the ModArchive.
#' Needs to be a \code{character} string representation of a file size
#' category as specified on ModArchive.org.
#' See the usage section for all possible options. Default is "unset" (meaning that
#' it will search for any file size). Note that the maximum file size of a
#' module is approximately 4068 kilobytes, meaning that the largest file size
#' category is irrelevant for `MOD' files. Also note that the category names are
#' inconsistant, these are the literal catagories used by ModArchive
#' @param genre.filter Genre filter to be used in a search in the ModArchive.
#' Needs to be a \code{character} string representation of a genre
#' as specified on ModArchive.org.
#' See the usage section for all possible options. Default is "unset" (meaning that
#' it will search for any genre).
#' @return \code{modArchive.info} will return a \code{\link{list}} with module info.
#' \code{modArchive.download} will download a module and return it as a
#' \code{\link{PTModule}} object. \code{modArchive.download} will search on the modArchive
#' and return the first page of search results in the form of a \code{\link{data.frame}}.
#' @name modArchive
#' @aliases modArchive.info
#' @aliases modArchive.download
#' @aliases modArchive.search.mod
#' @rdname modArchive
#' @note The `\code{modArchive}' functions were created for your convenience.
#' However, users and developers should not rely too heavily on them for
#' the following two reasons:
#' \itemize{
#' \item{The ModArchive is developed and maintained by a third party. Changes
#' in the structure of that database or the dissemination of the data may
#' cause the functions in this package to fail. I don't have any control on
#' these kind of changes.}
#' \item{Some of the \code{modArchive} functions implemented here, rely on regular
#' expressions to lookup information in ModArchive html code. I do realise
#' that using regular expressions in html code is the root of all evil (as
#' html is not a regular language). However, this approach is used to keep this package
#' light (i.e., not to rely on too many other packages), without having to write
#' a custom html-parser. That being said, this approach should work in most
#' cases, but may fail in some.}
#' }
#' @examples
#' \dontrun{
#' ## Search for the module that is also used as
#' ## an example in this package:
#' search.results <- modArchive.search.mod("ProTrackR", "module_instruments", "MOD")
#'
#' ## apparently there are multiple modules in
#' ## database that use 'ProTrackR' as a
#' ## sample name. Select the actual module from the
#' ## list:
#' search.select <- subset(search.results, title == "intro")
#'
#' ## get all available details for this mod from
#' ## ModArchive.org:
#' modArchive.info(search.select$mod.id)
#'
#' ## download the selected module from ModArchive.org:
#' mod <- modArchive.download(search.select$mod.id)
#' }
#' @author Pepijn de Vries
#' @export
modArchive.info <- function(mod.id)
{
  mod.id <- as.integer(mod.id[[1]])
  result = list()
  con <- url(paste("http://modarchive.org/module.php?", mod.id, sep = ""), "rb")
  page.source <- readLines(con, warn = F)
  close(con)
  if (any(grepl("<h1 class='notification'>Error</h1>", page.source)))
  {
    stop(unlist(regmatches(page.source, gregexpr("(?<=<p class='notification'>).*?(?=</p>)", page.source, perl = TRUE))))
  }
  result$mod.id        <- mod.id
  result$title         <- unlist(regmatches(page.source, gregexpr("(?<=<h1>).*?(?= <span class=\"module-sub-header\">[(])", page.source, perl = TRUE)))
  result$title         <- .htmlUnescape(result$title)
  result$filename      <- unlist(regmatches(page.source, gregexpr("(?<=<span class=\"module-sub-header\">[(]).*?(?=[)]</span></h1>)", page.source, perl = TRUE)))
  result$downloads     <- unlist(regmatches(page.source, gregexpr("(?<=<li class=\"stats\">Downloads: ).*?(?=</li>)", page.source, perl = TRUE)))
  result$favorited     <- unlist(regmatches(page.source, gregexpr("(?<=<li class=\"stats\">Favourited: ).*?(?= times</li>)", page.source, perl = TRUE)))
  result$MD5           <- unlist(regmatches(page.source, gregexpr("(?<=<li class=\"stats\">MD5: ).*?(?=</li>)", page.source, perl = TRUE)))
  result$format        <- unlist(regmatches(page.source, gregexpr("(?<=<li class=\"stats\">Format: ).*?(?=</li>)", page.source, perl = TRUE)))
  result$size          <- unlist(regmatches(page.source, gregexpr("(?<=<li class=\"stats\">Uncompressed Size: ).*?(?=</li>)", page.source, perl = TRUE)))
  result$genre         <- unlist(regmatches(page.source, gregexpr("(?<=<li class=\"stats\">Genre: ).*?(?=</li>)", page.source, perl = TRUE)))
  result$url           <- unlist(regmatches(page.source, gregexpr("(?<=&nbsp;<a href=\").*?(?=\" class=\"standard-link\">Download</a>)", page.source, perl = TRUE)))
  result$member.rating <- unlist(regmatches(page.source, gregexpr("(?<=Member Rating: ).*?(?=\\))", page.source, perl = TRUE)))
  result$member.rating <- gsub("<.*?>", "", result$member.rating, perl = T)
  result$member.rating <- gsub(" (", "", result$member.rating, fixed = T)
  result$revwr.rating  <- unlist(regmatches(page.source, gregexpr("(?<=Reviewer Rating: ).*?(?=\\))", page.source, perl = TRUE)))
  result$revwr.rating  <- gsub("<.*?>", "", result$revwr.rating, perl = T)
  result$revwr.rating  <- gsub(" (", "", result$revwr.rating, fixed = T)
  result               <- format.modarchive.table(result)
  return(result)
}

#' @rdname modArchive
#' @export
modArchive.download <- function(mod.id)
{
  mod.id <- as.integer(mod.id[[1]])
  con <- url(paste("http://api.modarchive.org/downloads.php?moduleid=", mod.id, sep = ""), "rb")
  mod <- read.module(con)
  close(con)
  return (mod)
}

#' @rdname modArchive
#' @export
modArchive.search.mod <- function(search.text,
                                  search.where  = c("filename_or_songtitle", "filename_and_songtitle", "filename", "songtitle", "module_instruments", "module_comments"),
                                  format.filter = c("unset", "669", "AHX", "DMF", "HVL", "IT", "MED", "MO3", "MOD", "MTM", "OCT", "OKT", "S3M", "STM", "XM"),
                                  size.filter   = c("unset", "0-99", "100-299", "300-599", "600-1025", "1025-2999", "3072-6999", "7168-100000"),
                                  genre.filter  = c("unset", "Alternative", "Gothic", "Grunge", "Metal - Extreme", "Metal (general)", "Punk", "Chiptune", "Demo Style",
                                                    "One Hour Compo", "Chillout", "Electronic - Ambient", "Electronic - Breakbeat", "Electronic - Dance",
                                                    "Electronic - Drum and Bass", "Electronic - Gabber", "Electronic - Hardcore", "Electronic - House", "Electronic - IDM",
                                                    "Electronic - Industrial", "Electronic - Jungle", "Electronic - Minimal", "Electronic - Other",
                                                    "Electronic - Progressive", "Electronic - Rave", "Electronic - Techno", "Electronic (general)", "Trance - Acid",
                                                    "Trance - Dream", "Trance - Goa", "Trance - Hard", "Trance - Progressive", "Trance - Tribal", "Trance (general)",
                                                    "Big Band", "Blues", "Jazz - Acid", "Jazz - Modern", "Jazz (general)", "Swing", "Bluegrass", "Classical", "Comedy",
                                                    "Country", "Experimental", "Fantasy", "Folk", "Fusion", "Medieval", "New Ages", "Orchestral", "Other", "Piano",
                                                    "Religious", "Soundtrack", "Spiritual", "Video Game", "Vocal Montage", "World", "Ballad", "Disco", "Easy Listening",
                                                    "Funk", "Pop - Soft", "Pop - Synth", "Pop (general)", "Rock - Hard", "Rock - Soft", "Rock (general)", "Christmas",
                                                    "Halloween", "Hip-Hop", "R and B", "Reggae", "Ska", "Soul"))
{
  search.where  <- match.arg(search.where)
  format.filter <- match.arg(format.filter)
  size.filter   <- match.arg(size.filter)
  genre.filter  <- match.arg(genre.filter)
  genre.filter  <- genre.table$genre.id[genre.table$genre == genre.filter]
  ## Set filters:
  con <- url(paste("http://modarchive.org/filter.php?genreid=", genre.filter, sep = ""))
  page.source   <- readLines(con, warn = F)
  close(con)
  ## Set filters:
  con <- url(paste("http://modarchive.org/filter.php?format=", format.filter, sep = ""))
  page.source   <- readLines(con, warn = F)
  close(con)
  ## Set filters:
  con <- url(paste("http://modarchive.org/filter.php?size=", size.filter, sep = ""))
  page.source   <- readLines(con, warn = F)
  close(con)
  ## make sure that details are returned:
  con <- url("http://modarchive.org/index.php?detail=1")
  page.source   <- readLines(con, warn = F)
  close(con)
  ## perform search:
  con <- url(paste("http://modarchive.org/index.php?request=search&query=",
                   utils::URLencode(search.text),
                   "&submit=Find&search_type=",
                   search.where, sep = ""))
  page.source   <- readLines(con, warn = F)
  close(con)
  ## analyse search results:
  result        <- lines.between(page.source, "Search Results</h1>", "Search</h1>")
  if (is.null(result)) return (empty.line)
  result        <- paste(result, collapse = "\n")
  ## split per table row:
  result        <- unlist(strsplit(result, "</tr>", fixed = T))
  result        <- as.list(result[c(-1, -length(result))])
  result        <- lapply(result, extract.search.results)
  result        <- do.call(rbind, result)
  result        <- format.modarchive.table(result)
  result$filename <- gsub("^.*?[#$]", "", result$url)
  return(result)
}

extract.search.results <- function(page.source)
{
  result <- list()
  result$mod.id        <- unlist(regmatches(page.source, gregexpr("(?<=<span class=\"module-detail\">Module ID:</span> <span class='module-detail'>).*?(?=</span>)", page.source, perl = TRUE)))
  result$title         <- unlist(regmatches(page.source, gregexpr("(?<=listing\">\n).*?(?=\n</span>)", page.source, perl = TRUE)))
  result$title         <- .htmlUnescape(result$title)
  result$filename      <- NA
  result$downloads     <- unlist(regmatches(page.source, gregexpr("(?<=<span class=\"module-detail\">Downloads:</span> <span class='module-detail'>).*?(?=</span>)", page.source, perl = TRUE)))
  result$favorited     <- NA
  result$MD5           <- NA
  result$format        <- unlist(regmatches(page.source, gregexpr("(?<=<span class=\"format-icon\">).*?(?=</span>)", page.source, perl = TRUE)))
  result$size          <- unlist(regmatches(page.source, gregexpr("(?<=<span class=\"module-detail\">Size:</span> <span class='module-detail'>).*?(?=</span>)", page.source, perl = TRUE)))
  result$genre         <- unlist(regmatches(page.source, gregexpr("(?<=<span class=\"module-detail\">Genre:</span> <span class='module-detail'>).*?(?=</span>)", page.source, perl = TRUE)))
  result$url           <- unlist(regmatches(page.source, gregexpr("(?<=<a href=\").*?(?=\" title)", page.source, perl = TRUE)))

  result$member.rating <- unlist(regmatches(page.source, gregexpr("(?<=<span class='module-listing'>).*?(?=</span>)", page.source, perl = TRUE)))
  result$member.rating <- gsub("<.*?>", "", result$member.rating, perl = T)
  result$member.rating <- gsub("Rated ", "", result$member.rating, fixed = T)
  result$revwr.rating  <- NA
  return (unlist(result))
}

lines.between <- function(source.text, start, end)
{
  start <- which(grepl(start, source.text))[1]
  end   <- which(grepl(end, source.text))
  end   <- end[length(end)]
  if (length(start) == 0 || length(end) == 0 ||
      is.na(start) || is.na(end)) return (NULL)
  if (start > end) return (NULL)
  return(source.text[start:end])
}

format.modarchive.table <- function(x)
{
  if (is.null(dim(x))) x <- t(as.matrix(x))
  x[,"mod.id"]    <- as.integer(as.character(x[,"mod.id"]))
  x[,"downloads"] <- as.integer(as.character(x[,"downloads"]))
  x[,"favorited"] <- as.integer(as.character(x[,"favorited"]))
  x               <- as.data.frame(x, stringsAsFactors = F)
}
