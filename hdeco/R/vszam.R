"vszam" <-
function (VARIAB="Y",WHICH=c(2,3)) {
  MIT <- paste(sep="",VARIAB,WHICH)
  HONNAN <- names(.QND)
  as.vector(charmatch(MIT,HONNAN))
}

