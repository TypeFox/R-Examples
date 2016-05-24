# Access to Banku Danych Lokalnych
# through the API in https://mojepanstwo.pl/api/bdl/

getAllPages <- function(url0, debug=0) {
  page = 1
  result <- list()
  repeat {
    if (debug>0) cat("page ",page,"\n")

    url <- paste0(url0, page)
    document <- rjson::fromJSON(file=url, method='C')

    if (length(document$Dataobject) == 0) break()
    result <- c(result, document$Dataobject)
    page <- page + 1
  }
  result
}

getMPgminy <- function(debug = 0) {
  url0 <- 'https://api-v3.mojepanstwo.pl/dane/gminy?limit=500&page='
  result <- getAllPages(url0, debug=debug)

  tmp <- lapply(result, function(d) {
    c(id = d$id,
      nazwa = d$data$gminy.nazwa,
      teryt = d$data$gminy.teryt,
      powiat.id = d$data$powiaty.id,
      wojewodztwo = d$data$wojewodztwa.nazwa,
      powierzchnia = d$data$gminy.powierzchnia,
      wydatki_roczne = d$data$gminy.wydatki_roczne,
      zadluzenie_roczne = d$data$gminy.zadluzenie_roczne,
      typ_nazwa = d$data$gminy.typ_nazwa,
      liczba_ludnosci = d$data$gminy.liczba_ludnosci,
      adres = d$data$gminy.adres)
  })
  data.frame(do.call(rbind, tmp), stringsAsFactors = FALSE)
}

getMPpowiaty <- function(debug = 0) {
  url0 <- 'https://api-v3.mojepanstwo.pl/dane/powiaty?limit=500&page='
  result <- getAllPages(url0, debug=debug)

  tmp <- lapply(result, function(d) {
    c(id = d$id,
      nazwa = d$data$powiaty.nazwa,
      wojewodztwo = d$data$wojewodztwa.nazwa)
  })
  data.frame(do.call(rbind, tmp), stringsAsFactors = FALSE)
}

getMPwojewodztwa <- function(debug = 0) {
  url0 <- 'https://api-v3.mojepanstwo.pl/dane/wojewodztwa?limit=500&page='
  result <- getAllPages(url0, debug=debug)

  tmp <- lapply(result, function(d) {
    c(id = d$id,
      wojewodztwo = d$data$wojewodztwa.nazwa)
  })
  data.frame(do.call(rbind, tmp), stringsAsFactors = FALSE)
}

getBDLtree <- function(raw = FALSE, debug = 0) {
  url0 <- 'https://api-v3.mojepanstwo.pl/dane/bdl_wskazniki?limit=500&page='
  result <- getAllPages(url0, debug=debug)

  if (raw) return(result)

  tmp <- lapply(result, function(d) {
    c(id = d$id,
      slug = d$slug,
      opis = d$data$bdl_wskazniki.grupa_tytul,
      data_aktualizacji = d$data$bdl_wskazniki.data_aktualizacji)
  })
  data.frame(do.call(rbind, tmp), stringsAsFactors = FALSE)
}

getBDLsearch <- function(query = "", debug = 0, raw = FALSE) {
  url <- paste0('https://api-v3.mojepanstwo.pl/bdl/search?q=', htmlEscape(query))
  if (raw) {
    document <- jsonlite::fromJSON(txt = url,simplifyVector=FALSE)
    return(document)
  }
  else {
    document <- jsonlite::fromJSON(txt = url,simplifyDataFrame=TRUE)
    return(document)
  }
}

getBDLseries <- function (metric_id = "", slice = NULL, time_range = NULL, wojewodztwo_id = NULL,
          powiat_id = NULL, gmina_id = NULL, meta = NULL, debug = 0,
          raw = FALSE)
{
  url <- paste0("https://api-v3.mojepanstwo.pl/bdl/series?metric_id=",
                metric_id)
  if (!is.null(slice))
    url <- paste0(url, "&slice=", htmlEscape(slice))
  if (!is.null(time_range))
    url <- paste0(url, "&time_range=", htmlEscape(time_range))
  if (!is.null(wojewodztwo_id))
    url <- paste0(url, "&wojewodztwo_id=", htmlEscape(wojewodztwo_id))
  if (!is.null(powiat_id))
    url <- paste0(url, "&powiat_id=", htmlEscape(powiat_id))
  if (!is.null(gmina_id))
    url <- paste0(url, "&gmina_id=", htmlEscape(gmina_id))
  if (!is.null(meta))
    url <- paste0(url, "&meta=", htmlEscape(meta))

  if (raw) {
    document <- jsonlite::fromJSON(txt = url,
                                   simplifyDataFrame=FALSE)
    return(document)
  }
  else {
    document <- jsonlite::fromJSON(txt = url,
                                   simplifyDataFrame=FALSE)
    met <- t(sapply(document$slices, function(s) s$slice))
    if (nrow(met) == 1) met <- t(met)
    fullmet <- do.call(what = rbind, lapply(1:ncol(met), function(d) {
      tmp <- do.call(what = rbind, document$meta$dimensions[[d]]$options)
      rownames(tmp) <- tmp[,1]
      met[,d] <<- unlist(tmp[met[,d],"value"])
      data.frame(dim=d, tmp)
    }) )

    dgs <- lapply(seq_along(document$slices), function(sn) {
      s <- document$slices[[sn]]
      tmp <- data.frame(do.call(what = rbind, args = s$series), units = s$units)
      for (i in ncol(met):1)
        tmp <- data.frame(dimension = met[sn,i], tmp, row.names = 1:nrow(tmp))
      for (i in 1:ncol(tmp))
        if (class(tmp[,i]) == "list")
          tmp[,i] <- sapply(tmp[,i], '[', 1)
        tmp
    })
    document <- do.call(what = rbind, dgs)
    return(document)
  }
}
