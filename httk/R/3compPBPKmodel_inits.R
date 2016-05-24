initparms3comp <- function(newParms = NULL){
  parms <- c(
    BW = 70,
    CLmetabolismc = 0.203,
    kgutabs = 1,
    Qcardiacc = 0,
    Qgfrc = 0.108,
    Qgutf = 0.205,
    Qliverf = 0.0536,
    Vgut = 0,
    Vliver = 0,
    Vrest = 0,
    Fraction_unbound_plasma = 0.0682,
    CLmetabolism = 0.0,
    Qcardiac = 0,
    Qgfr = 0.0,
    Qgut = 0.0,
    Qliver = 0.0,
    Kliver2plasma = 0,
    Krest2plasma = 0,
    Kgut2plasma = 0,
    Ratioblood2plasma = 0
  )
  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
  }
  if (!is.null(newParms)) parms[names(newParms)] <- newParms
  out <- .C("getParms_3comp",
   as.double(parms),
  out=double(length(parms)),
  as.integer(length(parms)))$out
  names(out) <- names(parms)
  out
}

Outputs3comp <- c(
    "Cgut",
    "Cliver",
    "Crest",
    "Cserum"
)


initState3comp <- function(parms, newState = NULL) {
  Y <- c(
    Agutlumen = 0.0,
    Agut = 0.0,
    Aliver = 0.0,
    Arest = 0.0,
    Ametabolized = 0.0,
    Atubules = 0.0,
    AUC = 0.0
  )
  Y <- with(as.list(parms), {  Y
  })

  if (!is.null(newState)) {
    if (!all(names(newState) %in% c(names(Y)))) {
      stop("illegal state variable name in newState")
    }
    Y[names(newState)] <- newState
  }
  Y
}
