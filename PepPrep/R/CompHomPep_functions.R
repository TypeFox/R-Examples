# MPL03 Funktionen

# Allgemeiner Ablauf
#  1.Lese Datenbanken ein
#  2.Finde alle Paare mit Hilfe des Entry_Names
#  3.Trypsiniere die Sequenzen der Paare
#  4.Finde identische Peptide
#  5.Füge diese zu einer sequenz zusammen
#  6.Speicher das Paar mit neuer Sequenz ab


# Hauptfunktion
# Eingabe:
#  path_o1(char) Pfad zu den Proteinen des 1.Organismus
#  path_o2(char) Pfad zu den Proteinen des 2.Organismus
#  path(char) Pfad zum Ziel für das Ergebnis
#  width(number) Formatierung der Breite der FASTA-Ausgabe
#  intermediate(binary) Ausgabe der Zwischenschritte in einer Liste zurückgeben
#  target(char) Aminosäure hinter der das Verdauungsenzym schneidet
#  exception(char) wenn target diese Aminosäure folgt wird nicht geschnitten
# Ausgabe:


compHomToPepFasta <- function(path_o1, path_o2, path, width = 60, intermediate = FALSE, target = "K|R", exception = "P")
{
  # fange fehlende Eingaben ab
  if(missing(path_o1))
  {
    message <- "compHomToPepFasta: Error, missing path for organism one! (path_o1)"
    print(message)
    return(message)
  }
  if(missing(path_o2))
  {
    message <- "compHomToPepFasta: Error, missing path for organism two! (path_o2)"
    print(message)
    return(message)
  }
  if(missing(path))
  {
    message <- "compHomToPepFasta: Error, missing path for the result! (path)"
    print(message)
    return(message)
  }

  # fange flasche Eingabe  von target und exception ab
  targettest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,target)
  exceptiontest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,exception)
  if(attr(targettest,"match.length") != nchar(target))
  {
    message <- "compHomToPepFatsa: Error target is wrong"
    print(message)
    return(message)
  }
  if(attr(exceptiontest,"match.length") != nchar(exception))
  {
    message <- "compHomToPepFasta: Error exception is wrong"
    print(message)
    return(message)
  }
  
  
  # Hauptroutine
  print("start findEntryPairs")
  tbl <- findEntryPairs(path_o1, path_o2)

  print("start addHeader")
  tbl <- addHeader(tbl, path_o1, path_o2)
   
  print("start addAllPeptides")
  tbl <- addAllPeptides(tbl, path_o1, path_o2, target, exception)
  
  print("start myWriteFasta")
  fasta <- myWriteFasta2(tbl, path, width)
  
  result<-tbl
  if(intermediate)
  {
    result <- list(tbl, fasta)
    return(result)
  }
  else
  {
    return(paste0("computation done: ", path))
  }
}



# Diese Funktion berechnet alle möglichen Paarungen für gleiche Gennamen und gibt eine
# Liste mit den HeaderPositionen in beiden Organismen zum Gennamen zurück
# Eingabe: Pfade zu den Fasta Files
# Ausgabe: data.frame (Entryname, org_one header position, org_two header position)
findEntryPairs <- function(path_o1, path_o2)
{
   # lese Quellen für Sequenzen
  org1fasta <- readLines(path_o1)
  org2fasta <- readLines(path_o2)
  
  # berechne die Positionen der header
  hpos_orgone <- grep(">", org1fasta)
  hpos_orgtwo <- grep(">", org2fasta)
  
  # Beispiel
  # >sp|A0M8Q6|LAC7_HUMAN Ig lambda-7 chain C region OS=Homo sapiens GN=IGLC7 PE=1 SV=2
  # split >sp  A0M8Q6  LAC7_HUMAN Ig lambda-7 chain C region OS=Homo sapiens GN=IGLC7 PE=1 SV=2
  # pos[1,1]=5
  # str_sub  LAC7
  
  # organismus 1
  # teile Header
  ool <- str_split(org1fasta[hpos_orgone],"\\|")
  # oo sammelt die Namen der Proteine
  oo <- c()
  # lese jetzt den generellen Teil aus der EntryID
  for (i in 1:length(ool))
  {
    # matrix mit start und stop des musters
    pos <- str_locate(ool[[i]][3],"_")
    oo <- c(oo,str_sub(ool[[i]][3],1,pos[1,1]-1))    
  }
  
  # organismus 2
  # teile Header
  otl <- strsplit(org2fasta[hpos_orgtwo],"\\|")
  ot <- c()
  # lese jetzt den generellen Teil aus der EntryID
  for (i in 1:length(otl))
  {
    # matrix mit start und stop des musters
    pos <- str_locate(otl[[i]][3],"_")
    ot <- c(ot,str_sub(otl[[i]][3],1,pos[1,1]-1))    
  }
  
  # Schnittmenge zwischen oo und ot enthält alle Proteine die ein Pärchen bilden
  pairs <- intersect(oo, ot)
  # an welchen Positionen stehen meine Namen
  # so nicht: !!!!!!!
  # moo<-match(pairs,oo)
  # mot<-match(pairs,ot)
  # wichtig ist jetzt alle Paarungen zwischen den Isoformen zu berechnen
  oohits <- (oo %in% pairs)
  othits <- (ot %in% pairs)
  ootbl <- data.frame(name=oo[oohits], hposorg1=hpos_orgone[oohits], stringsAsFactors=FALSE)
  ottbl <- data.frame(name=ot[othits], hposorg2=hpos_orgtwo[othits], stringsAsFactors=FALSE)
  result <- merge(ootbl, ottbl, by="name")
  
  # wo stehen die Header in der Eingabe?
  # result<-data.frame(name=pairs,hposorg1=hpos_orgone[moo],hposorg2=hpos_orgtwo[mot],stringsAsFactors=FALSE)
    
  # return(list(ool,oo,otl,ot,result))
  return(result)
}

# -------------------------------------------------------------
# Funktion, die zur Berechnung der Peptide führt
# 
# Teilschritte: 60er Abschnitte zusammen führen, Trypsinieren, Vergleichen und zusammenhängen
# 
# Eingabe: data.frame(name, posorgone,posorgtwo), Pfad zu Quelle org1 und Pfad zu Quelle org2
# Ausgabe: data.frame(name, posorgone,posorgtwo,peptideseq) Um die gesuchte Peptidsequenz erweiterte Tabelle

addAllPeptides<-function(tbl, path_o1, path_o2,target, exception)
{
  # lese Quellen für Sequenzen
  org1fasta <- readLines(path_o1)
  org2fasta <- readLines(path_o2)
  # Positionen aller Header im Quellfile
  allheadpos_org1 <- grep(">", org1fasta)
  allheadpos_org2 <- grep(">", org2fasta)
  lallhp_o1 <- length(allheadpos_org1)
  lallhp_o2 <- length(allheadpos_org2)
  
  # speichert die neuen Sequenzen
  newseq <- rep("seq", length(tbl$name))
  
  # gehe die Paare durch und suche Anfangs und Endstelle für die Sequenz
  for (i in seq(along=tbl$name))
  {
    # für Organismus one
    oohs <- match(tbl$hposorg1[i], allheadpos_org1)
    oostart <- tbl$hposorg1[i]+1
    if((oohs+1) <= lallhp_o1)
    {
      oohe <- allheadpos_org1[oohs+1]
      oostop <- oohe-1
    }
    else
    {      
      oostop <- length(org1fasta) 
    }
    
#     print(paste("i: ",i))
#     print(paste("oostart: ",oostart))    
#     print(paste("oohs: ",oohs))
#     print(paste("oohe: ",oohe))
#     print(paste("oostop: ",oostop))
    
    ooseq <- paste(org1fasta[oostart:oostop], sep="", collapse="")
        
    # für Organismus Two
    oths <- match(tbl$hposorg2[i],allheadpos_org2)
    otstart <- tbl$hposorg2[i]+1
    if((oths+1) <= lallhp_o2)
    {
      othe <- allheadpos_org2[oths+1]
      otstop <- othe-1
    }
    else
    {      
      otstop<-length(org2fasta) 
    }
    
#     print(paste("i: ",i))
#     print(paste("otstart: ",otstart))    
#     print(paste("oths: ",oths))
#     print(paste("othe: ",othe))
#     print(paste("otstop: ",otstop))
    
    otseq <- paste(org2fasta[otstart:otstop],sep="",collapse="")
    
#     if (tbl$name[i]=="UKD")
# {
#       print()
# }
    
    # Trypsinierung
    trypooseq <- trypsinateSeq(ooseq,target, exception)
    trypotseq <- trypsinateSeq(otseq,target, exception)
    
    # Vergleich 
    hits <- intersect(trypooseq,trypotseq)
    newseq[i] <- paste(hits, sep="", collapse="")
    
    
  }
  
  result <- cbind(tbl, data.frame(cpepseq=newseq, stringsAsFactors=FALSE))
  # entferne Zeilen mit leeren Sequenzen  
  result <- result[result$cpepseq != "", ]  

  return(result)
}

#------------
# Trypsinierung
# Eingabe: seq String/Sequenz, target char, exception char
# default für Trypsin
# Ausgabe: Vektor mit Sequenzstücken

trypsinateSeq<-function(seq, target= "K|R", exception= "P")
{
  # ueberpruefe ob 
  # Am Abkürzungen:ARNDCQEGHILKMFPSTWYV
  targettest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,target)
  exceptiontest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,exception)
  if(attr(targettest,"match.length") != nchar(target))
  {
    message <- "trypsinateSeq: Error target is wrong"
    print(message)
    return(message)
  }
  if(attr(exceptiontest,"match.length") != nchar(exception))
  {
    message <- "trypsinateSeq: Error exception is wrong"
    print(message)
    return(message)
  }
  
  
  # berechne zunächst die Abschnitte in denen geschnitten wird
  vrk <- str_locate_all(seq, target)[[1]]
  vp <- str_locate_all(seq, exception)[[1]]
  pminus <- vp-1
  # tryppos enthält jetzt alle Positionen von K und R denen kein P folgt
  tryppos <- vrk[!(vrk[ ,1] %in% pminus[ ,1]), ]
  
#   #für den Fall, dass keine Schnittpunkte vorliegen
#   if (class(tryppos)=="integer")
#   {
#     result<-seq    
#     return(result)
#   }
  #-----START-----------neuer Teil evtl. nicht brauchbar
  # entweder kein Schnittpunkt oder
  # nur einer
  if (class(tryppos) == "integer")  
  {
    # print("mutPepSeq: eine Sequenz")
    # für den Fall, dass keine Schnittpunkte vorliegen
    # ist die ganze Sequenz ein Block
    if (length(tryppos) == 0)
    {
      tryppos <- matrix( c(1,nchar(seq)), nrow=1, ncol=2, dimnames=list(c(), c("start","end")))
    }
    # tryppos wird zum Vektor mit Element 1 Wert vrk[1,1]
    else
    {
      tryppos <- matrix( c(tryppos["start"], tryppos["end"]), nrow=1, ncol=2, dimnames=list(c(), c("start","end")))      
    }
    
    
  }
  # tryppos wird zur leeren Matrix
  if(class(tryppos) == "matrix")
  {
    # für den Fall, dass keine Schnittpunkte vorliegen
    # ist die ganze Sequenz ein Block
    if(dim(tryppos)[1] == 0)
    {
      tryppos <- matrix( c(1,nchar(seq)), nrow=1, ncol=2, dimnames=list(c(), c("start","end")))
    }
    
  }
  #--------ENDE--------neuer Teil evtl. nicht brauchbar 
  
  # berechne Peptidbruchstücke
  # ACHTUNG hier kommen manchmal einige stücke zu klein an aber einfache Lösung, die Felder bleiben leer
  result <- rep("", length(tryppos[,1])+1)
  # result<-rep("bla01",length(tryppos)+1)
  
  # prints()
  # tryppos enthält die letzte Position eines Teilstrings
  first <- 1
  for (i in seq(along=tryppos[,1]))
  {
    last <- tryppos[i,1]
    pep <- substr(seq,first,last)
    result[i] <- pep
    first <- last+1    
  }
  
  # füge das letzte Element an
  last <- nchar(seq)
  if (first <= last)
  {
    result[length(result)] <- substr(seq,first,last)
  }  
    
  return(result)  
}


#------
# Header zur Tabelle hinzufügen
# Eingabe: Tabelle mit Name und Positionen (data.frame), Pfade zu den Quelldaten
# Ausgabe: um header erweiterte Tabelle

addHeader <- function(tbl,path_o1,path_o2)
{
  # lese Quellen für Sequenzen
  org1fasta <- readLines(path_o1)
  org2fasta <- readLines(path_o2)
  
  myheader <- rep("head",length(tbl$name))
  
  for (i in seq(along=tbl$name))
  {
  secheader <- gsub("^>"," org2:",org2fasta[tbl$hposorg2[i]])
  # myheader[i]<-paste(org1fasta[tbl$hposorg1[i]],org2fasta[tbl$hposorg2[i]],sep="")
  myheader[i] <- paste(org1fasta[tbl$hposorg1[i]], secheader, sep="")  
  }
  # print(paste("Länge tbl",length(tbl$name)))
  # print(paste("Länge header",length(myheader)))  
  result <- cbind(tbl, data.frame(header=myheader, stringsAsFactors=FALSE))
  
  return(result)
}


# -----------
# Funktion zum Schreiben der Fasta
# Eingabe: data.frame mit name,pos1,pos2,cpepseq; Pfad zum Speichern der Datei; Breite der Sequenz
# Ausgabe: Vektor (header,seq,header,seq usw...) und Datei mit Fastasequenz

myWriteFasta2 <- function(tbl, path, width)
{
    
  result <- c()
  for (i in seq(along=tbl$name))
  { 
    # print(i)
    # header
    result <- c(result, tbl$header[i])
    
    # Sequenz
    # teile die Sequenz
    len <- nchar(tbl$cpepseq[i])
    strt <- 1
    stp <- width
    test <- len-width
    # solange wir nur Mittelstücke betrachten
    while(strt <= test)
    {
      result <- c(result, substr(tbl$cpepseq[i], strt, stp))
      strt <- strt+width
      stp <- stp+width
    }
    # Endstück, Startposition ist schon richtig, nur Ende muss korrigiert werden
    stp <- len
    result <- c(result, substr(tbl$cpepseq[i], strt, stp))      
    
  }
  
  write(result, file=path) 
  
  return(result)
}