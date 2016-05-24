# Funktionen SNVtoFASTA

# Funktionen um AM-Sequenzen aus den ANNOVAR Annotierungen zu berechnen
# 1.Gene sortieren
# 2.Spalte AAchange auslesen
# 3.nach Transcript-ID anhand der mRNA-RefSeq forschen
# 4.Zur TranscriptID die AM Sequenz aus der Fasta Datei lesen
# 5.alle Änderungen an der fasta-Sequenz durchführen
# 6.neuen Fasta Header erzeugen und mit geänderter Sequenz in neue fasta-Datei schreiben

# Hauptfunktion
# Eingabe: 
# tbl   Annotationstabelle (data.frame)
# glst  Genliste (data.frame),
# spath   Pfad zur Proteindatenbank von ENSEMBL "BasisDaten/Homo_sapiens.GRCh37.70.pep.all.fa"
# tpath  Zielpfad fertige Fasta-Sequenz/Dateiname (character) 
# mymart  welches Mart soll es sein (character)
# myarchive   soll biomaRt auf ein Archive zugreifen? 
# width Breite der Sequenzen in der neuen FASTA-Datei (numeric)
# intermediate(binary) Ausgabe der Zwischenschritte in einer Liste zurückgeben
# target(char) Aminosäure hinter der das Verdauungsenzym schneidet
# exception(char) wenn target diese Aminosäure folgt wird nicht geschnitten

# Ausgabe: 
# list(aachanges=aachanges,transcripts=transcripts,mutfasta=mutfasta) 
# list[[1]] Annotationstabelle mit Angaben zu den Mutationen an der Aminosäuresequenz (data.frame) 
# list[[2]] Zuordnung zwischen mRNA-Refseq (NM_ID) und EnsembleTranskriptID
# list[[3]] Die Fasta als Vector
# list[[4]] Logbuch zur Mutation der Transkripte

# Ausgabe: FASTA-Transkripte (Textdatei an der Stelle Zielpfad)

snvToPepFasta <- function(tbl, glst, mymart = "ensembl", myarchive = FALSE, spath, tpath, width = 60, intermediate = FALSE, target = "K|R", exception = "P")
{
  if(missing(tbl))
  {
    message <- "snvToPepFasta: Error, missing annotationtable! (tbl)"
    print(message)
    return(message)    
  }
  if(missing(glst))
  {
    message <- "snvToPepFasta: Error, missing genlist! (glst)"
    print(message)
    return(message)
  }
  if(missing(spath))
  {
    message <- "snvToPepFasta: Error, missing path to ENSEMBL protein database! (spath)"
    print(message)
    return(message)
  }
  if(missing(tpath))
  {
    message <- "snvToPepFasta: Error, targetpath! (tpath)"
    print(message)
    return(message)
  }
  
  # fange flasche Eingabe  von target und exception ab
  targettest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,target)
  exceptiontest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,exception)
  if(attr(targettest,"match.length") != nchar(target))
  {
    message <- "snvToPepFasta: Error target is wrong"
    print(paste("snvToPepFasta: Error target is wrong",target))
    return(message)
  }
  if(attr(exceptiontest,"match.length") != nchar(exception))
  {
    message <- "snvToPepFasta: Error exception is wrong"
    print("snvToPepFasta: Error exception is wrong")
    return(message)
  }
  
  
  glst <- glst[order(glst$Genes),]
  # print(class(glst))
  # zum Testen
  # Füge SNP-Key hinzu
  print("start addKey")
  tbl <- addKey(tbl)
  
  # prozessiere die AM-Änderungen
  print("start getAMchanges")
  aachanges <- getAMchanges(tbl, glst)
  
  # bestimme die Ensembletranskriptids
  print("start getTranscriptID")
  transcripts <- getTranscriptIDv2(aachanges$nmid, mymart, myarchive)
  
  # mutiere jetzt die fasta Sequenzen aus Homo_sapiens.GRCh37.70.pep.all.fa und erstelle eine neue Fasta-Ausgabe
  # print(head(transcripts))
  print("start mutateProtToPep")
  mutfastalist <- mutateProtToPep(aachanges, transcripts, spath, target, exception)
  mutfasta <- mutfastalist[[1]]
  mutlog <- mutfastalist[[2]]
  
  # mutfasta<-mutateAM(aachanges,transcripts)
  # print(head(mutfasta))
  print("start myWriteFasta")
  # schreibe die fertig bearbeitete AM-Sequenz in eine Datei
  writefasta <- myWriteFasta(mutfasta, tpath, width)
  
  if(intermediate)
  {
    result <- list(aachanges=aachanges, transcripts=transcripts, mutfasta=mutfasta, mutlog=mutlog)
    return(result)  
  }
  else
  {
    return(paste0("computation done: ",tpath))    
  }
  
}


# SNP-Key hinzufügen
# Eingabe: data.frame Tabelle ohne Schlüssel
# Ausgabe: data.frame mit Schlüsselspalte snpkey
addKey <- function(tbl)
{
  result <- tbl
  result$snpkey <- paste(tbl$Chr, tbl$Start, tbl$Ref, tbl$Obs, sep="_")
  return(result)
}

# berechne die nötigen Änderungen für die Geneauswahl
# Eingabe: 
# tbl   Annotationstabelle (data.frame),
# glst  Genliste(character)
# Ausgabe:
# result  AM-Mutationsliste (data.frame) sprich die Annotationstabelle erweitert um die Änedrungsangaben aber 
#         beschränkt auf die Genlistenauswahl
getAMchanges <- function(tbl, glst)
{
  result <- data.frame()
  # welche SNVs passen zu unserer Genauswahl?
  snvlist <- tbl[tbl$Gene %in% glst, ]
  
  # splitte jetzt die AAChangespalte auf
  myindex <- c()
  refseq <- c()
  aachange <- c()
  # Prozessiere nun die Angabe zur Mutation
  for (i in seq(along=snvlist$Gene))
  {
    # print(i)
    details <- unlist(strsplit(snvlist[i, "AAChange"], ":"))
    # Hinweis beim check
    # details <- unlist(strsplit(snvlist[i, "AAChange"], "\:"))
    # Error: '\:' is an unrecognized escape in character string starting ""\:"
    # details<-unlist(strsplit(snvlist[i,"AAChange"],"\\:"))
    
    # print(details)
    if (sum(grepl("p\\.", details)) == 1)
    {
      refseq <- c(refseq, details[grep("NM", details)])
      aachange <- c(aachange, details[grep("p\\.", details)])
      myindex <- c(myindex, i)
    }  
    
  }
  
  # print(length(aachange))
  
  origAA <- c()
  posAA <- c()
  mutAA <- c()
  for (i in seq(along=aachange))
  {
    # entferne p.
    aachange[i] <- sub("p\\.", "", aachange[i])    
    # teile auf nach Original Position Mutation
    l <- nchar(aachange[i])
    origAA <- c(origAA, substr(aachange[i], 1, 1))
    posAA <- c(posAA, substr(aachange[i], 2, l-1))
    mutAA <- c(mutAA, substr(aachange[i], l, l))
  }
  #   print("refseq")  
  #   print(length(refseq))
  #   print(unique(refseq))
  #   print("origAA")  
  #   print(length(origAA))  
  #   print("posAA")
  #   print(length(posAA))
  #   print("mutAA")  
  #   print(length(mutAA))
  
  result <- cbind(snvlist[myindex,], data.frame(nmid=refseq, origAA=origAA, posAA=posAA, mutAA=mutAA)) 
  # result<-cbind(mytable,data.frame(nmid=refseq,origAA=origAA,posAA=posAA,mutAA=mutAA))
  #   print(refseq)
  #   print(aachange)
  #   print(origAA)
  #   print(posAA)
  #   print(mutAA)
  return(result)
}


# berechne nun die Ensemble-TranskriptIDs aus der refseq NM_ID
# diesmal vektorwertig
# Eingabe: 
#  nmlst  NM_ID liste (vector/data.frame?)
#  mymart  welches Mart soll es sein (character)
#  myarchive   soll biomaRt auf ein Archive zugreifen?
# Ausgabe: 
#  result Tabelle mit NM_ID und die passenden TranskriptIDs  (data.frame)

getTranscriptIDv2 <- function(nmlst, mymart, myarchive)
{
  # print(nmlst)
  # print(unique(nmlst))
  # print(order(unique(nmlst)))
  # auf einzigartige beschränken und sortieren
  nmlst <- unique(nmlst)
  nmlst <- nmlst[order(nmlst)]
  # print(nmlst)
  # biomart-Abfrage
  # library(biomaRt)
  
  ensembl <- useMart(mymart, archive=myarchive)
  ensembl <- useDataset("hsapiens_gene_ensembl",mart <- ensembl)  
  # Beispiel für Rückgabe
  # ensembl_transcript_id refseq_mrna                                                                           description
  # ENST00000327044   NM_015658 nucleolar complex associated 2 homolog (S. cerevisiae) [Source:HGNC Symbol;Acc:24517]
  myAttri <- c("ensembl_transcript_id", "refseq_mrna", "description")
  myFilter <- c("refseq_mrna")
  
  # result<-data.frame()  
  tbl <- getBM(attributes=myAttri, filters=myFilter, values=nmlst, mart=ensembl)   
  result <- data.frame(ensembl_transcript_id=tbl$ensembl_transcript_id, nmid=tbl$refseq_mrna, pname=tbl$description, stringsAsFactors=TRUE)
  
  return(result)
}


# Trypsinierung mit direkter Selektion mutierter Peptide
# Eingabe:AAsequenz, dataframe:nmid,origAA,posAA,mutAA,snpkey
# Ausgabe:peptidAAsequenz

mutPepSeq <- function(mutpos, seq, target= "K|R", exception= "P")
{
  
  # ueberpruefe ob 
  # Am Abkürzungen:ARNDCQEGHILKMFPSTWYV
  targettest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,target)
  exceptiontest <- regexpr("^[ARNDCQEGHILKMFPSTWYV]([|][ARNDCQEGHILKMFPSTWYV])*" ,exception)
  if(attr(targettest,"match.length") != nchar(target))
  {
    message <- "mutPepSeq: Error target is wrong"
    print(message)
    return(message)
  }
  if(attr(exceptiontest,"match.length") != nchar(exception))
  {
    message <- "mutPepSeq: Error exception is wrong"
    print(message)
    return(message)
  }
  
  
  result <- c("")
  # mutpos zur Sicherheit erneut aufsteigend sortieren
  mutpos <- mutpos[order(mutpos$posAA), ]
  # print(class(mutpos$posAA))
  # berechne zunächst die Abschnitte in denen geschnitten wird
  vrk <- str_locate_all(seq, target)[[1]]
  vp <- str_locate_all(seq, exception)[[1]]
  pminus <- vp-1
  tryppos <- vrk[ !(vrk[,1]%in%pminus[,1]), ]
  
#   print("tryppos")
#   print(class(tryppos))
#   print(tryppos)
  
  # entweder kein Schnittpunkt oder
  # nur einer
  if (class(tryppos) == "integer")  
  {
    # print("mutPepSeq: eine Sequenz")
    # für den Fall, dass keine Schnittpunkte vorliegen
    # ist die ganze Sequenz ein Block
    if (length(tryppos) == 0)
    {
      tryppos <- matrix(c(1, nchar(seq)), nrow=1, ncol=2, dimnames=list(c(), c("start", "end")))
    }
    # tryppos wird zum Vektor mit Element 1 Wert vrk[1,1]
    else
    {
      tryppos <- matrix( c(tryppos["start"], tryppos["end"]), nrow=1, ncol=2, dimnames=list(c(), c("start", "end")))      
    }
    
    
  }

  # tryppos wird zur leeren Matrix
  if(class(tryppos) == "matrix")
  {
    # für den Fall, dass keine Schnittpunkte vorliegen
    # ist die ganze Sequenz ein Block
    if(dim(tryppos)[1] == 0)
    {
      tryppos <- matrix(c(1, nchar(seq)), nrow=1, ncol=2, dimnames=list(c(), c("start", "end")))
    }
    
  }
  
  # überprüfe für jede Mutation zu welcher Position diese gehört
  # Start und Endposition vom Substring
  spos <- 1
  epos <- 0
  # merke den Einstiegspunkt für den Durchlauf
  my_i <- 1
  
#   print("vor der Schleife")
#   print("mutpos")
#   print(mutpos$posAA)
#   print("tryppos")
#   print(tryppos)
  
# Schleife für das dataframe
  for(k in seq(along=mutpos$posAA))
  {
    # print(c("k:",k))
    
    # Schleife für die Abschnittsmatrix    
    for(i in my_i:dim(tryppos)[1])
    {
      # print(c("i:",i))    
      
      # falls die Position vor einem Musterblock liegt
      if(as.numeric(levels(mutpos$posAA)[mutpos$posAA[k]]) <= tryppos[i,"start"]) 
      {
        # print("case1")
        epos <- tryppos[i, "start"]
        # str_sub(x,2,1)="" gibt einen leeren String zurück falls die Positionen falsch sind
        # daher funktioniert dies hier
        # result<-paste(result,str_sub(seq,spos,epos),"_start_",spos,"_",epos,"x",tryppos[i,"start"],"_",tryppos[i,"end"],"_",tryppos[i-1,"end"],"x",sep="")
        result <- paste(result, str_sub(seq, spos, epos), sep="")
        # print(result)
        spos <- tryppos[i, "start"]+1
        break
      }
      
             
      else 
      {
        # falls die Position in einem Musterblock liegt 
        if (i < dim(tryppos)[1])
        {
          #print("case2")
          spos <- tryppos[i, "start"]+1
          my_i = i+1
        }
        # letzten Musterblock überschritten und meine Mutationsposition ist größer als
        # das Ende vom letzten Musterblock
        else
        {
          # print("case3")
          spos <- tryppos[i, "start"]+1 
          epos <- nchar(seq)
          # result<-paste(result,str_sub(seq,spos,epos),"_Fin_",spos,"_",epos,"xxx",sep="")
          result <- paste(result, str_sub(seq, spos, epos), sep="")
          # print(result)
          # Achtung hier darf ich nur ein mal rein, sonst hängt er immer neue Sequenzen an
          return(result)
          
        }
        
      }
      
    }
    
    
  }
  
  return(result)
  
}


# Mutation 
# Eingabe: 
#  tbl     Um AA-Mutationen erweiterte Annotationstabelle (data.frame)
#  nm      Tabelle mit Zuordnung NM_ID zu EnsembleTranskriptID (data.frame)
#  spath   Datei-Pfad zur Quelle der Transcripte (character)
# Ausgabe (list): 
#  resfasta  Vector mit Mutierten fasta-Sequenzen 
#  log       Vector mit Fehlern im Mutationsteil

mutateProtToPep <- function(tbl, nm, spath, target, exception)
{
  # alle fertig bearbeiteten Fastas
  resfasta <- c()
  # Logbuch für die Fehler im Mutationsteil
  log <- c()
  # Quelle
  pepseq <- readLines(spath)
  # pepseq<-readLines("BasisDaten/Homo_sapiens.GRCh37.70.pep.all.fa")
  # fasse alle Positionen der Header zusammen
  allheaderpos <- grep(">", pepseq)
  # Vektor mit allen Headern
  allheader <- pepseq[allheaderpos]
  
  # print(head(allheaderpos))
  # gehe alle Transkripte durch
  for (i in seq(along=nm$ensembl_transcript_id ))
  {
    # header vom gesuchten Gen
    header <- allheaderpos[grep(nm$ensembl_transcript_id[i], allheader)]
    # print("i")
    # print(i)
    if (!length(header) == 0)
    {
      # print(nm$ensembl_transcript_id[i])
      # print(header)
      # print(class(header))
      # print(length(header))
      # print(allheaderpos)
      if (length(header) > 1) 
      {
        print("mehr als ein header")
        log <- c(log,"mehr als ein header")
        break
      }
      # falls es keinen weiteren header gibt
      if (match(header, allheaderpos) >= length(allheaderpos))
      {
        # bis zu letzten Zeile 
        nheader <- length(pepseq)+1
      }
      else
      {
        # oder bis zum nächsten header
        nheader <- allheaderpos[match(header,allheaderpos)+1]
      }
      #     #nächster Header
      #     nheader<-allheaderpos[match(header,allheaderpos)+1]
      #     print("nheader")
      #     print(nheader)
      #     print(class(nheader))
      #     print(length(nheader))
      #     #falls es keinen weiteren header gibt
      #     if (is.na(nheader))
      #     {
      #       #bis zu letzten Zeile 
      #       nheader<-length(pepseq)+1
      #     }
      
      #erstelle den neuen header
      workheader <- pepseq[header]
      en_ids <- enheaderToDetail(workheader)
      # falls die ENST bzw. die gewählte Accssesionnumber doppelt vor kommt muss diese um ein Suffix erweitert werden.
      if (length(grep(en_ids["ENST"], resfasta)) > 0)
      {
        mysuffix <- length(grep(en_ids["ENST"], resfasta))+1
        workheader <- paste(paste(">", paste(en_ids["ENST"], "x", mysuffix, sep=""), sep=""), "|", nm$pname[i], sep=" ")
        
      }
      else
      {
        workheader <- paste(paste(">", en_ids["ENST"], sep=""), "|", nm$pname[i], sep=" ")
        # workheader<-paste(paste(">",en_ids["ENST"],sep=""),"|",nm$pname[i],en_ids["ENSP"],en_ids["ENSG"],sep= " ")
        
      }
      
      
      # erstelle die neue AA-Sequenz
      workseq <- c()
      # print("header;nheader")
      # print(header)
      # print(nheader)
      for (j in (header+1):(nheader-1))
      {
        # print("j")
        # print(j)
        workseq <- paste(workseq, pepseq[j], sep="")        
      }
      # bearbeite die workseq mit den Mutationen
      # print(attributes(tbl))
      # mutpos<-subset(tbl, nmid == nm$nmid[i],select=c(nmid,origAA,posAA,mutAA,snpkey))
      # print(levels(nm$nmid))
     
      mutpos <- tbl[tbl$nmid == levels(nm$nmid)[nm$nmid[i]], c("nmid", "origAA", "posAA", "mutAA", "snpkey")]
      mutpos <- mutpos[order(as.numeric(levels(mutpos$posAA)[mutpos$posAA])), ]
      # print(mutpos)
      mutseq <- c()
      cutstart <- 1
      cutend <- 1
      for (k in seq(along=mutpos$nmid))
      {      
        # cutend<-as.numeric(mutpos$posAA[k])
        # cutend<-mutpos$posAA[k]
        cutend <- as.numeric(levels(mutpos$posAA[k]))[mutpos$posAA[k]]
        # vergleiche ob Orginial in der Sequenz stimmt
        # print(cutend)
        # print(class(cutend))
        # print("workseq")
        # print(class(workseq))
        # print(substr(workseq,cutend,cutend))
        if (substr(workseq, cutend, cutend) == mutpos$origAA[k])
        {
          # falls ja füge die Mutation ein
          mutseq <- paste(mutseq, substr(workseq, cutstart, cutend-1),mutpos$mutAA[k], sep="")
          
          #ergänze hier auch den Header um die Mutation der AA
          workheader <- paste(workheader, " ", mutpos$origAA[k], "->", mutpos$mutAA[k], "_", mutpos$posAA[k], sep="")
          
          cutstart <- cutend+1
        }
        else
        {
          # hänge einen Eintrag zu unkompatiblen Mutationen an den Header
          workheader <- paste(workheader, " wrongAA:", mutpos$origAA[k], "->", mutpos$mutAA[k], "_", mutpos$posAA[k], sep="")
          # falls nein, nehme die nächste Mutation
          # print(paste("wrongAA",mutpos$nmid[k]))
          log <- c(log, paste("wrongAA:", mutpos$nmid[k], mutpos$origAA[k], "->", mutpos$mutAA[k], "_", mutpos$posAA[k], "PeptidAA:", substr(workseq,cutend,cutend)))
          next
        }      
      }
      # füge jetzt den Sequenzrest an, da die obige Schleife nur bis zur letzten Mutation verlängert
      mutseq <- paste(mutseq, substr(workseq, cutstart, nchar(workseq)), sep="")
      
      # Trypsinierung
      # print("call mutPepSeq")
      # print(mutpos)
      # print(mutseq)
      mutseq <- mutPepSeq(mutpos, mutseq, target, exception)
      # print("done mutPepSeq")
      # print(mutseq)
      
      # Ausgabe
      resfasta <- c(resfasta, workheader, mutseq)
      # result<-c(result,workheader,workseq,mutseq)
      
    }
    else
    {
      # print(paste("TranskriptID nicht in Quelle" ,nm$ensembl_transcript_id[i]))
      log <- c(log, paste("TranskriptID nicht in Quelle", nm$ensembl_transcript_id[i]))
      next
    }
    
    
  }
  # fertigmutierte Sequenz
  # mutfasta
  result <- list(resfasta, log)
  return(result)
}





# Funktion zum Schreiben der Fasta-Sequenz in eine Datei
# Eingabe: 
# Vektor aus header und Sequenz, jeweils alternierend
#  tpath   Pfad zum Ziel für die FASTA-Datei
#  width   Breite der AM Ausgabe 
# Ausgabe: Aufgesplittete Sequenzliste

myWriteFasta <- function(lst, tpath, width)
{
  result <- c()
  for (i in seq(along=lst))
  {
    # falls ein header vorliegt    
    if (i %% 2 == 1)
    {
      result <- c(result, lst[i])
      next
    }
    # teile die Sequenz
    else 
    {
      len <- nchar(lst[i])
      strt <- 1
      stp <- width
      test <- len-width
      # solange wir nur Mittelstücke betrachten
      while(strt <= test)
      {
        result <- c(result, substr(lst[i], strt, stp))
        strt <- strt+width
        stp <- stp+width
      }
      # Endstück, Startposition ist schon richtig, nur Ende muss korrigiert werden
      stp <- len
      result <- c(result, substr(lst[i], strt, stp))
      
    }
  }
  write(result, file=tpath) 
  
  return(result)
}


# Funktion zum zerlegen eines Ensemble Peptid Headers
# Eingabe: Header-string
# Ausgabe: Named Vektor mit ENSP,ENST,ENSG IDs
enheaderToDetail <- function(header)
{
  # library(stringr)
  
  # schneide ENSP aus
  result <- c(str_extract(header, "ENSP[0-9]*"))  
  # schneide ENST aus
  result <- c(result, str_extract(header, "ENST[0-9]*"))
  # schneide ENSG aus
  result <- c(result, str_extract(header, "ENSG[0-9]*"))
  
  names(result) <- c("ENSP", "ENST", "ENSG")
  
  return(result)
}