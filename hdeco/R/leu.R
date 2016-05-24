"leu" <-
function (UJ="Y2",BE=names(.QND)) {
  ZEK <- UJ == BE

  if(sum(ZEK) == 0) {
    return(BE)
  }
  else {
    return(BE[!ZEK])
  }
}

