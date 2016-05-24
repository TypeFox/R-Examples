### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2013-03-28 19:02:06 umusus>
### romanization.R
### Convert from Hiragana or Katakana to Romanji without kakasi

.syllabicate.hira <- function(y){
  ye <- strsplit(y, NULL)
  res <- lapply(ye, function(z){
    i <- grep('[\u3083\u3085\u3087]', z)
    if(any(i)){
      z[i-1] <- sapply(i, function(j){
        paste(z[j-1], z[j], sep = '')
      })
      return(z[-i])
    }else{
      return(z)
    }
  })
  return(res)
}

kata2hira <- function(x){
  stopifnot(is.character(x))
  tbl <- c(jpn.syllabary$Hiragana, jpn.syllabary.add$Hiragana)
  names(tbl) <- c(jpn.syllabary$Katakana, jpn.syllabary.add$Katakana)
  s <- strsplit(x, NULL)
  sapply(s,  function(y){
    i <- is.na(r <- tbl[y])
    r[i] <- y[i]
    paste(r, collapse="")
  })
}

hira2kata <- function(x){
  stopifnot(is.character(x))
  tbl <- c(jpn.syllabary$Katakana, jpn.syllabary.add$Katakana)
  names(tbl) <- c(jpn.syllabary$Hiragana, jpn.syllabary.add$Hiragana)
  s <- strsplit(x, NULL)
  sapply(s,  function(y){
    i <- is.na(r <- tbl[y])
    r[i] <- y[i]
    paste(r, collapse="")
  })
}

kana2roma <- function(x, type = c("Hepburn", "Nippon.shiki", "Kunrei.shiki"),
                        cap = FALSE, ascii.only = TRUE){
  stopifnot(is.character(x))
  type <- match.arg(type)
  i <- switch(type,
              Hepburn = 3,
              Nippon.shiki = 4,
              Kunrei.shiki = 5)
  tbl <- c(jpn.syllabary[, i], jpn.syllabary.add$Romaji)
  names(tbl) <- c(jpn.syllabary$Hiragana, jpn.syllabary.add$Hiragana)
##  print(tbl)
  y <- .syllabicate.hira(kata2hira(x))
  s <- sapply(y, function(z){
    j <- is.na(r <- tbl[z])
    r[j] <- z[j]
    paste(r, collapse="")
  })
  s <- gsub('ltu(.)', '\\1\\1', s)      # geminate consonant
  if(type == "Hepburn" && ascii.only){
    s <- gsub('cc', 'tc', s)           # exception for geminate consonant
    s <- gsub('n([bmp])', 'm\\1', s)    # exception for the (Japanese) syllabic nasal
    s <- gsub('o[uo]', 'o', s)    # long sound 
    s <- gsub('u[uo]', 'o', s)    # long sound 
  }
  if(cap){
    s <- paste(toupper(substring(s, 1,1)), substring(s, 2),
                 sep="", collapse=" ")
  }
  return(s)
}

hiragana <- function() {
  tbl <- rbind(matrix(c(rep(227, 63), rep(129, 63), 128 + 1:63), nrow = 63), 
    matrix(c(rep(227, 20), rep(130, 20), 127 + 1:20), nrow = 20))
  hira <- apply(tbl, 1, function(x) {
    rawToChar(as.raw(x))
  })
  return(hira)
}

katakana <- function() {
  tbl <- rbind(matrix(c(rep(227, 31), rep(130, 31), 160 + 1:31), nrow = 31), 
    matrix(c(rep(227, 55), rep(131, 55), 127 + 1:55), nrow = 55))
  kata <- apply(tbl, 1, function(x) {
    rawToChar(as.raw(x))
  })
  return(kata)
}

ya.hira2kata <- function(x) {
  paste2 <- function(x, ...) {
    paste(x, ..., sep = "", collapse = "")
  }
  hira <- hiragana()
  kata <- katakana()
  kata <- kata[1:length(hira)]
  return(chartr(paste2(hira), paste2(kata), x))
}

ya.kata2hira <- function(x) {
  paste2 <- function(x, ...) {
    paste(x, ..., sep = "", collapse = "")
  }
  hira <- hiragana()
  kata <- katakana()
  kata <- kata[1:length(hira)]
  return(chartr(paste2(kata), paste2(hira), x))
}
