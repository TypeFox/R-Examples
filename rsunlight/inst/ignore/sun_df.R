#' Write data from any rsunlight function output to a csv file on your machine.
#' 
#' @import assertthat
#' @export
#' @param x Output from any of the rsunlight functions.
#' @param ... Further args passed to read.csv
#' @details This function attemps to coerce the raw output from each rsunlight function to a 
#' data.frame to write to csv, but it may fail in some cases. You can always make your own 
#' data.frame. 
#' @examples \dontrun{
#' out <- cg_legislators(latitude = 35.778788, longitude = -78.787805)
#' sun_df(out)
#' 
#' out <- cg_committees(member_ids='L000551')
#' sun_df(out)
#' }
sun_df <- function(x, ...){
  UseMethod("sun_df")
}

#' @export
sun_df.cg_legislators <- function(x){
  assert_that(class(x) == "cg_legislators")
  if(x$page$count == 1){
    iter <- x$results[[1]]
    iter <- replacemissing(iter)
    df <- data.frame(iter, stringsAsFactors = FALSE)
  } else {
    iter <- x$results
    iter <- lapply(iter, replacemissing)
    iter <- 
      lapply(iter, function(v){
        v[sapply(v, class) %in% "list"] <- 
          as.list(unlist(v[sapply(v, class) %in% "list"]))
        v
      })
    df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
  }
  return( df )
}

#' @export
sun_df.cg_committees <- function(x){
  assert_that(is(x, "cg_committees"))
  if(any(names(x$results[[1]]) %in% 'subcommittees')){
    iter <- x$results[[1]]$subcommittees
    iter <- lapply(iter, replacemissing)
    df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
  } else if(any(names(x$results[[1]]) %in% 'members')){
    iter <- x$results[[1]]$members
    iter <- lapply(iter, function(w) as.list(unlist(w)))
    iter <- lapply(iter, replacemissing)
    df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
    names(df) <- gsub("\\.", "_", names(df))
  } else {
    #   notmembers <- x$committee[!names(x$committee) %in% "members"]
    #   notmembers <- replacemissing(notmembers)
    #   members <- do.call(rbind.fill, lapply(x$committee$members, data.frame, stringsAsFactors = FALSE))
    #   cbind(notmembers, members)
    if(x$count == 1){
      iter <- x$results[[1]]
      iter <- replacemissing(iter)
      df <- data.frame(iter, stringsAsFactors = FALSE)
    } else {
      iter <- x$results
      iter <- lapply(iter, replacemissing)
      df <- do.call(rbind.fill, lapply(iter, data.frame, stringsAsFactors = FALSE))
    }
  }
  return( df )
}

replacemissing <- function(z){
  lapply(z, function(y){
    ifelse(length(y)==0, "", y)
  })
}