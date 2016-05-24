.HTML.file <- NULL
.tabs <- 0
.ligne <- F

#' Initialiser un fichier HTML
#' 
#' Cree et remplit les headers pour un fichier hTML
#' 
#' Si aucun nom de fichier n'est fourni, cree un fichier temporaire dans le repertoire temporaire
#' Le nom du fichier actuel est stocke dans .HTML.file
#' @encoding UTF-8
#' @param file Nom du fichier HTML a creer, par defaut un fichier temporaire
#' @param title Titre de la page
#' @param CSSfile Fichier CSS a utiliser
#' @examples
#' \dontrun{
#' HTMLInit(file="sortie.html", title="Titre de la page", CSSfile="desc.css")
#' }
#' @export
HTMLInit <- function(file=tempfile(pattern="report", fileext=".html"), title="", CSSfile="")
{
  file.copy(from=file.path(path.package("cosmosR"),"desc.css"), to=file.path(dirname(file),"desc.css"),overwrite=T)
  file.copy(from=file.path(path.package("cosmosR"),"diag.css"), to=file.path(dirname(file),"diag.css"),overwrite=T)
  file.copy(from=file.path(path.package("cosmosR"),"cosmosR.js"), to=file.path(dirname(file),"cosmosR.js"),overwrite=T)
  unlockBinding(".HTML.file", env=asNamespace("cosmosR"))
  unlockBinding(".tabs", env=asNamespace("cosmosR"))
  unlockBinding(".ligne", env=asNamespace("cosmosR"))
  .HTML.file <<- file
  .tabs <<- 0
  .ligne <<- F
  
  HTML("<!DOCTYPE html>", append=F)
  HTML("<HTML>")
  inc()
    HTML("<HEAD>")
    inc()
      HTML("<meta charset='", localeToCharset()[1], "' />")
      HTML("<title>", title, "</title>")
      HTML("<link rel='stylesheet' href='", CSSfile, "' />")
      HTML("<script src='cosmosR.js'></script>")
    dec()
    HTML("</HEAD>")
    HTML("<BODY>")
    inc()
}

#' Termine et clos le fichier HTML
#' 
#' Ecrit le footer du fichier hTML
#' 
#' Ecrit le footer du fichier initialise par HTMLInit, ouvre le fichier dans le navigateur et supprime l'acces.
#' @encoding UTF-8
#' @export
HTMLEnd <- function()
{
  if (is.null(.HTML.file)) return
  
  unlockBinding(".HTML.file", env=asNamespace("cosmosR"))
  
    dec()
    HTML("</BODY>")
  dec()
  HTML("</HTML>")
  
  browseURL(.HTML.file)
  .HTML.file <<- NULL
}

#' Ecrit dans le fichier HTML
#' 
#' Ecrit dans le fichier HTML cree par HTMLInit
#' 
#' Ecrit dans le fichier initialise par HTMLInit dont le nom est contenu dans .HTML.file
#' @encoding UTF-8
#' @param x Contenu a ecrire
#' @param ... Contenu concatene sans espace a x
#' @param append Decide si x... doit etre ajoute a un fichier existant
#' @param sep Separateur de fin de ligne, modifier pour ecrire sur la meme ligne du fichier
#' @export
HTML <- function(x, ..., append=T,sep="\n")
{
  if (is.null(.HTML.file)) return
  
  unlockBinding(".ligne", env=asNamespace("cosmosR"))
  
  x <- paste0(x,...,collapse="")
  
  if (.ligne)
    tabs <- ""
  else
    tabs <- paste0(rep("\t", .tabs),collapse="")
  
  if (sep == "\n")
    .ligne <<- F
  else
    .ligne <<-T
  
  x <- paste0(tabs, x)
  
  cat(x,file=.HTML.file,sep=sep, append=append)
}

inc <- function()
{
  unlockBinding(".tabs", env=asNamespace("cosmosR"))
  .tabs <<- .tabs+1
}

dec <- function()
{
  unlockBinding(".tabs", env=asNamespace("cosmosR"))
  if (.tabs>0) .tabs <<- .tabs-1
}