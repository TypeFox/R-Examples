### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2015-06-25 14:30:47 umusus>
### kakasi function fimaly

kakasi <- function(x, kakasi.option="-Ha -Ka -Ja -Ea -ka",
                   ITAIJIDICTPATH=Sys.getenv("ITAIJIDICTPATH", unset = NA),
                   KANWADICTPATH=Sys.getenv("KANWADICTPATH", unset = NA)){
  if(is.na(ITAIJIDICTPATH)){
    Sys.setenv(ITAIJIDICTPATH=.set.dict("itaijidict"))
  }else{
    stopifnot(file.exists(ITAIJIDICTPATH))
    Sys.setenv(ITAIJIDICTPATH=ITAIJIDICTPATH)
  }
  if(is.na(KANWADICTPATH)){
    Sys.setenv(KANWADICTPATH=.set.dict("kanwadict"))
  }else{
    stopifnot(file.exists(KANWADICTPATH))
    Sys.setenv(KANWADICTPATH=KANWADICTPATH)
  }
  stopifnot(is.character(x))
  ops <- strsplit(kakasi.option, " ")[[1]]
  ## This is a trick to work correctly
  ops <- c(" ", ops)
  stopifnot(length(ops) != 0)
  lc.ctype <- Sys.getlocale("LC_CTYPE")
  if (lc.ctype != "ja_JP.UTF-8" || lc.ctype != ""){
    if (Sys.getlocale("LC_CTYPE") == "Japanese_Japan.932") x <- sjis2utf8(x)
  }else{
    warning("kakasi() assumes \"ja_JP.UTF-8\" for LC_TYPE")
  }
  u <- sapply(x, function(i){
    .Call("rkakasi", x = i, k = ops, PACKAGE = "Nippon")
  })
  return(u)
}

.set.dict <- function(dict){
  lib <- .libPaths()[sapply(.libPaths(), function(x){file.exists(file.path(x, "Nippon"))})]
  d.path <- file.path(lib, "Nippon", "share", dict)
  stopifnot(file.exists(d.path))
  return(d.path)
}
