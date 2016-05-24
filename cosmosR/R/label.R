#' Etiquetter un objet
#' 
#' Attribuer ou recuperer l'etiquette d'un objet
#' 
#' Methode par defaut pour acceder et modifier l'etiquette d'un objet ou d'un element d'un objet.
#' Il est possible de supprimer l'etiquette en passant NULL a la fonction.
#' @rdname label
#' @encoding UTF-8
#' @param objet L'objet a nommer
#' @param value Le texte de l'etiquette
#' @return Renvoie une chaine de caracteres contenant l'etiquette de l'objet
#' @examples
#' a <- c(18,25,23,32)
#' 
#' label(a) <- "Age"
#' label(a)
#' # Renvoie "Age"
#' 
#' label(a) <- NULL # Supprime le label
#' 
#' df <- data.frame(a=25, b="H")
#' label(df) <- c(a="Age",b="Sexe")
#' @export
label <- function(objet) UseMethod("label")

#' @rdname label
#' @export
label.default <- function(objet)
{
  if (missing(objet))
    NULL
  else if (is.null(attr(objet,"label")))
    deparse(substitute(objet))
  else
    attr(objet,"label")
}

#' @rdname label
#' @export
"label<-" <- function(objet,value) UseMethod("label<-")

#' @rdname label
#' @export
"label<-.default" <- function(objet,value)
{
  attr(objet,"label")<-value[[1]]
  if (length(value)>1) warning("label : vecteur de longueur 1")
  objet
}

#' @rdname label
#' @export
"label<-.data.frame" <- function(objet,value)
{
  if (length(value)==1 && is.null(names(value)))
  {
    objet<-unclass(objet)
    label(objet)<-value
    class(objet)="data.frame"
  }
  else
  {
    if (!is.null(names(value)))
    {
      for (nom in names(value))
      {
        if (is.null(objet[[nom]]))
          warning(paste(nom,"n'est pas un element"))
        else
          label(objet[[nom]])<-value[[nom]]
      }
    }
  }
  objet
}