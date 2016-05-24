#' Downloads data of Brazilian government bonds directly from the website
#'
#' This function looks into the tesouro direto website
#' (http://www.tesouro.fazenda.gov.br/tesouro-direto-balanco-e-estatisticas) and
#' downloads all of the files containing prices and yields of government bonds.
#' You can use input asset.codes to restrict the downloads to specific bonds
#'
#' @param asset.codes Strings that identify the assets (1 or more assets) in the
#'   names of the excel files. E.g. asset.codes = 'LTN'. When set to NULL, it
#'   will download all available assets
#' @param dl.folder Name of folder to save excel files from tesouro direto (will
#'   create if it does not exists)
#' @param do.clean.up Clean up folder before downloading? (TRUE or FALSE)
#' @param do.overwrite Overwrite excel files? (TRUE or FALSE). If FALSE, will
#'   only download the new data for the current year
#'
#' @return TRUE if successful
#' @export
#'
#' @examples
#' # only download file where string LTN_2015 is found in its name
#' # (only 1 file for simplicity)
#' download.TD.data(asset.codes = 'LTN_2015')
#'
#' # The excel file shoulbe available in folder 'TD Files' (default name)
#'
download.TD.data <- function(asset.codes = 'LTN',
                             dl.folder = 'TD Files',
                             do.clean.up = T,
                             do.overwrite = F) {
  # check folder

  if (!dir.exists(dl.folder)) {
    warning(paste('Folder ', dl.folder, 'was not found. Creating a it..'))
    dir.create(dl.folder)
  }

  # clean up folders
  if (do.clean.up) {
    list.f <- dir(dl.folder, pattern = '*.xls', full.names = T)
    file.remove(list.f)
  }

  # check if user has internet
  test.internet <- curl::has_internet()

  if (!test.internet){
    stop('No internet connection found...')
  }

  base.url <-
    "http://www.tesouro.fazenda.gov.br/tesouro-direto-balanco-e-estatisticas"

  # read html
  html.code <- paste(readLines(base.url, warn = F), collapse = "\n")

  # fixing links strings

  my.links <-
    stringr::str_extract_all(html.code,pattern = '<a href=\"(.*?).xls\"')[[1]][c(-1,-2)]
  my.links <-
    stringr::str_replace_all(my.links,'<a href=\"',replacement = '')
  my.links <-
    stringr::str_replace_all(my.links,'\"',replacement = '')

  fixLinks <- function(linkIn) {
    if (!stringr::str_detect(linkIn, 'www3.tesouro.gov.br')) {
      linkOut <- paste('http://www.tesouro.fazenda.gov.br',linkIn,sep = '')

    } else {
      linkOut <- linkIn

    }

    return(linkOut)

  }

  my.links <- lapply(my.links ,FUN = fixLinks)
  my.links <- paste(my.links)

  # find asset code in links

  if (!is.null(asset.codes)) {
    idx <- logical(length = length(my.links))
    for (i.asset in asset.codes) {
      temp.idx <- stringr::str_detect(string = my.links,pattern = i.asset)
      idx <- idx | temp.idx
    }

    my.links <- my.links[idx]

  }


  n.links <- length(my.links)

  my.c <- 1
  for (i.link in my.links) {
    splitted.str <- stringr::str_split(i.link,'/')[[1]]
    out.file <-
      paste0(dl.folder,'/',splitted.str[length(splitted.str)])

    cat(paste0('\nDownloading file ', out.file, ' (',my.c, '-', n.links, ')'))

    # check if file exists and if it does not contain the current year
    # in its name (thats how tesouro direto stores new data)

    test.current.year <- stringr::str_detect(out.file,format(Sys.Date(),'%Y'))

    if (file.exists(out.file)&(!test.current.year)&(!do.overwrite)){

      cat(' Found file in folder, skipping it.')
      my.c <- my.c + 1
      next()
    }

    cat(' Downloading...')
    utils::download.file(
      url = i.link,
      method = 'internal',
      mode = 'wb',
      destfile = out.file,
      quiet = T )

    my.c <- my.c + 1

  }

  return(T)

}
