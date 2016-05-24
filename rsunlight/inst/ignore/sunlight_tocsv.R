#' Write data from any rsunlight function output to a csv file on your machine.
#'
#' @keywords internal
#' @param x Output from any of the rsunlight functions.
#' @param file File name, with path.
#' @param ... Further arguments passed to read.csv
#' @details This function attemps to coerce the raw output from each rsunlight function to a
#' data.frame to write to csv, but it may fail in some cases. You can always make your own
#' data.frame.
#' @examples \dontrun{
#' out <- cg_getcommittees(id = 'JSPR')
#' sunlight_tocsv(out, "~/myfile.csv")
#'
#' out <- cg_getcommitteesallleg(bioguide_id = 'S000148')
#' sunlight_tocsv(out, "~/myfile.csv")
#'
#' out <- cg_getcommitteeslist(chamber = 'Joint')
#' sunlight_tocsv(out, "~/myfile.csv")
#'
#' out <- cg_getdistrictlatlong(latitude = 35.778788, longitude = -78.787805)
#' sunlight_tocsv(out, "~/myfile.csv")
#'
#' out <- cg_getdistrictzip(zip = 27511)
#' sunlight_tocsv(out, "~/myfile.csv")
#'
#' out <- cg_getlegislator(last_name = 'Pelosi')
#' sunlight_tocsv(out, "~/myfile.csv")
#' out <- cg_getlegislator(party = 'D')
#' sunlight_tocsv(out, "~/myfile.csv")
#' }
sunlight_tocsv <- function(x, file="~/", ...){
  UseMethod("sunlight_tocsv")
}

#' @method sunlight_tocsv cg_getcommittees
#' @rdname sunlight_tocsv
#' @keywords internal
sunlight_tocsv.cg_getcommittees <- function(x, file="~/", ...){
  stopifnot(is(x, "cg_getcommittees"))
  notmembers <- x$committee[!names(x$committee) %in% "members"]
  notmembers <- replacemissing(notmembers)
  members <- do.call(rbind.fill, lapply(x$committee$members, data.frame, stringsAsFactors = FALSE))
  dat <- cbind(notmembers, members)
  write.csv(dat, file=file, ..., row.names=FALSE)
}

#' @method sunlight_tocsv cg_getcommitteesallleg
#' @rdname sunlight_tocsv
#' @keywords internal
sunlight_tocsv.cg_getcommitteesallleg <- function(x, file="~/", ...){
  stopifnot(is(x, "cg_getcommitteesallleg"))
  dat <- lapply(x$committees, function(b){
    df <- do.call(rbind.fill, lapply(b$committee$subcommittees, data.frame, stringsAsFactors = FALSE))
    tmp <- b$committee[!names(b$committee) %in% "subcommittees"]
    tmp <- replacemissing(tmp)
    tmp <- data.frame(tmp, stringsAsFactors = FALSE)
    names(tmp) <- c('committee_chamber','committee_id','committee_name')
    if(!is.null(df)){
      names(df) <- c('subcommittee_chamber','subcommittee_id','subcommittee_name')
      cbind(tmp, df)
    } else{ tmp }
  })
  dat2 <- do.call(rbind.fill, dat)
  write.csv(dat2, file=file, ..., row.names=FALSE)
}

#' @method sunlight_tocsv cg_getcommitteelist
#' @rdname sunlight_tocsv
#' @keywords internal
sunlight_tocsv.cg_getcommitteelist <- function(x, file="~/", ...){
  stopifnot(is(x, "cg_getcommitteelist"))
  iter <- x$committees
  iter <- replacemissing(iter)
  df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
  write.csv(df, file=file, ..., row.names=FALSE)
}

#' @method sunlight_tocsv cg_getdistrictlatlong
#' @rdname sunlight_tocsv
#' @keywords internal
sunlight_tocsv.cg_getdistrictlatlong <- function(x, file="~/", ...){
  stopifnot(is(x, "cg_getdistrictlatlong"))
  iter <- x$districts
  df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
  write.csv(df, file=file, ..., row.names=FALSE)
}

#' @method sunlight_tocsv cg_getdistrictzip
#' @rdname sunlight_tocsv
#' @keywords internal
sunlight_tocsv.cg_getdistrictzip <- function(x, file="~/", ...){
  stopifnot(is(x, "cg_getdistrictzip"))
  iter <- x$districts
  df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
  write.csv(df, file=file, ..., row.names=FALSE)
}

#' @method sunlight_tocsv cg_getlegislator
#' @rdname sunlight_tocsv
#' @keywords internal
sunlight_tocsv.cg_getlegislator <- function(x, file="~/", ...){
  stopifnot(is(x, "cg_getlegislator"))
  if(x$found == 1){
    iter <- x$data
    df <- data.frame(iter, stringsAsFactors = FALSE)
  } else {
    iter <- x$data
    iter <- lapply(iter, replacemissing)
    df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
  }
  write.csv(df, file=file, ..., row.names=FALSE)
}

replacemissing <- function(x){
  lapply(x, function(y){
    ifelse(length(y)==0, "", y)
  })
}
