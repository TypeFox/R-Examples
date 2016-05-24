#' Charger un fichier de donnees dans l'environnement
#' 
#' Charge un fichier texte ou excel contenant les donnees ainsi que les formats et labels
#' 
#' La fonction charge les donnees contenues dans le fichier de donnees et renvoie une table de valeurs.
#' Si des labels/formats sont definis ils seront appliques (fichiers labels.sas, formats.sas et attribformats.sas)
#' Les variables pour lesquelles un format est defini seront considerees comme des variables qualitatives.
#' La table de valeurs prend automatiquement le nom du fichier (suffixe par le numero de feuille).
#' Le fichier de donnees est charge depuis ../../data/, les formats depuis le repertoire courant.
#' @encoding UTF-8
#' @param fichier Fichier de donnees a charger
#' @param feuille Feuille a utiliser si fichier excel et en cas de feuilles mutiples (et qu'on veut acceder a une feuille au-dela de la premiere)
#' @return La data frame avec les labels et formats
#' @examples
#' \dontrun{Ma_table <- charger("donnees.xls", feuille=2)}
#' @export
charger <- function(fichier,feuille=1)
{
  # Verification de la presence du fichier
  fichier<-paste("../../data/",fichier,sep="")
  if (!file.exists(fichier))
  {
    warning("Le fichier source n'existe pas !")
    return -1
  }
  
  # Lecture du fichier selon l'extension
  if (grepl("\\.csv$",fichier) || grepl("\\.txt$",fichier))
    x<-read.csv2(fichier, na.strings="",encoding="native.enc")
  else if (grepl("\\.xlsx?$",fichier))
    x<-read.xlsx(fichier,feuille,encoding="native.enc")
  
  # Lecture des formats a partir du fichier SAS
  if (file.exists("formats.sas"))
  {
    formats=list(0)
    
    con=file("formats.sas","r",encoding="native.enc")
      formatsfile=readLines(con)
    close(con)
    
    # Nettoyage du fichier : espaces de debut et fin, et autour des '=', conservation uniquement des lignes significatives (value xxx et 'code'='level')
    formatsfile=sub("^[[:space:]]*","",formatsfile)
    formatsfile=sub("[[:space:]]*$","",formatsfile)
    formatsfile=formatsfile[grepl("^value *\\w",formatsfile) | grepl("\\d+ *= *\\'.*?\\'",formatsfile)]
    formatsfile=sub("[[:space:]]*=[[:space:]]","=",formatsfile)
    formatsfile=sub("^(\\d+)","\\'\\1\\'",formatsfile)
    
    # Creation de la liste de formats
    for (format in formatsfile)
    {
      if (grepl("^value *\\w",format))
        formats[[strsplit(format," ")[[1]][2]]]<-character(0)
      else
        eval(parse(text=paste0("formats[[length(formats)]]<-c(formats[[length(formats)]],",format,")")))
    }
    #formats<-lapply(formats,"attr<-",which="ordered",value=F)
    #attr(formats[[strsplit(formatsfile[grepl("^value *\\w* *order$",formatsfile)]," ")[[1]][2]]],"ordered")<-T
    formats[[1]]<-NULL
    
    # Lecture des attributions de format depuis le fichier SAS
    if (file.exists("attribformats.sas"))
    {
      con=file("attribformats.sas","r",encoding="native.enc")
        attribfile=readLines(con)
      close(con)
      
      # Nettoyage du fichier : espaces de debut et fin, espaces multiples, autour des '=', '.' final, conservation uniquement des lignes significatives (finissant en format=xxx)
      attribfile=sub("^[[:space:]]*","",attribfile)
      attribfile=sub("[[:space:]]*$","",attribfile)
      attribfile=sub("\\.?$","",attribfile)
      attribfile=attribfile[grepl("(\\w[[:space:]]+)+format *= *\\w",attribfile)]
      attribfile=gsub("[[:space:]]+"," ",attribfile)
      attribfile=sub("[[:space:]]*=[[:space:]]","=",attribfile)
      
      # Attribution des formats aux variables concernees
      for (attrib in attribfile)
      {
        attrib=strsplit(attrib," ")[[1]]
        format=strsplit(attrib[length(attrib)],"=")[[1]][2]
        for (var in attrib[-length(attrib)])
        {
          if (!is.null(x[[var]]))
              x[[var]]<-factor(x[[var]],levels=names(formats[[format]]),labels=formats[[format]])#,ordered=attr(formats[[format]],"ordered"))
        }
      }
    }
  }
  
  # Lecture des labels a partir du fichier SAS
  if (file.exists("labels.sas"))
  {
    con=file("labels.sas","r",encoding="native.enc")
    labelsfile=readLines(con)
    close(con)
    labels=labelsfile[grepl('^ *\\w*? *= *\\".*?\\" *$',labelsfile)]
    labels=paste(labels,collapse=",")
    label_exe = paste("label(x)<-c(",labels,")")
    eval(parse(text=label_exe))
  }
  
  # Remplacement des "" par NA dans les facteurs
  for (var in names(x))
  {
    if (is.factor(x[[var]]))
      levels(x[[var]])[levels(x[[var]])==""]<-NA
  }
  
  # Enregistrement de la table dans l'environnement
  x
}