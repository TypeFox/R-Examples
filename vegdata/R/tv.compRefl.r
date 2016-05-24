if(getRversion() >= "2.15.1")  utils::globalVariables(c("write.dbf"))

tv.compRefl <- function (refl1, refl2, tv_home, check.nr = FALSE, simplify = TRUE, verbose = FALSE, Sink = TRUE, filter.1, filter.2, new = FALSE, file="compRefl.txt", ...)  {
  if (missing(tv_home)) tv_home <- tv.home()
    refl.a <- if(is.character(refl1)) read.dbf(file.path(tv_home, "Species", refl1, "species.dbf")) else refl1
    refl.b <- if(is.character(refl2)) read.dbf(file.path(tv_home, "Species", refl2, "species.dbf")) else refl2
    names(refl.a) <- TCS.replace(names(refl.a))
    names(refl.b) <- TCS.replace(names(refl.b))
  
    refl1 <- deparse(substitute(refl1))
    refl2 <- deparse(substitute(refl2))
  	refl.a$TaxonNameOriginal <- refl.a[, "TaxonName"]; refl.b$TaxonNameOriginal <- refl.b[, "TaxonName"]
    refl.a[, "TaxonName"] <- taxname.abbr(refl.a[, "TaxonName"])
    refl.b[, "TaxonName"] <- taxname.abbr(refl.b[, "TaxonName"])
  if(simplify) {
    refl.a[, "TaxonName"] <- taxname.simplify(refl.a[, "TaxonName"])
    refl.b[, "TaxonName"] <- taxname.simplify(refl.b[, "TaxonName"])  	
  }
  diff.A <- sort(as.character(refl.b[!refl.b[, "TaxonName"] %in% refl.a[, "TaxonName"], "TaxonName"]))
  diff.B <- sort(as.character(refl.a[!refl.a[, "TaxonName"] %in% refl.b[, "TaxonName"], "TaxonName"]))
    if (check.nr) {
      merged.df <- merge(refl.a, refl.b, by = "TaxonName", all.x = FALSE)
      selectedcolumns <- if('EDITSTATUS' %in% names(merged.df)) c("TaxonName", "TaxonUsageID.x", "TaxonUsageID.y", 'BEGRUEND', 'EDITSTATUS') else c("TaxonName", "TaxonUsageID.x", "TaxonUsageID.y")
      nonmatchingNumbers <- merged.df[as.character(merged.df$TaxonUsageID.x) != as.character(merged.df$TaxonUsageID.y), selectedcolumns]
      nonmatchingNumbers <- nonmatchingNumbers[!is.na(nonmatchingNumbers[, 1]), ]
      #   nonmatchingNumbers <-  if('EDITSTATUS' %in% names(merged.df)) nonmatchingNumbers[order(nonmatchingNumbers$EDITSTATUS, nonmatchingNumbers[, 2]), ] else nonmatchingNumbers[order(nonmatchingNumbers[, 2]), ]
      print(nrow(nonmatchingNumbers), 'lines with non-matching numbers.')
      tab <- table(nonmatchingNumbers$TaxonName)
      print(paste('But', sum(tab == 2), 'duplicated names after simplification.'))
      nonmatchingNumbers$TaxonName[which(nonmatchingNumbers$TaxonName %in% names(tab[tab == 1]))]
    }
      merged.df <- merge(refl.a, refl.b, by = 'TaxonUsageID', all.x = FALSE)
      selectedcolumns <- if('EDITSTATUS' %in% names(merged.df)) c("TaxonUsageID", 'TaxonName.x', "TaxonName.y", 'BEGRUEND', 'EDITSTATUS') else c("TaxonUsageID", 'TaxonName.x', "TaxonName.y")
      nonmatchingNames <- merged.df[merged.df$TaxonName.x != merged.df$TaxonName.y, selectedcolumns]
      nonmatchingNames <- nonmatchingNames[!is.na(nonmatchingNames[, 1]), ]
      nonmatchingNames <- if('EDITSTATUS' %in% names(merged.df)) nonmatchingNames[order(nonmatchingNames[, 'EDITSTATUS'], nonmatchingNames[, 2]), ] else  nonmatchingNames <- nonmatchingNames[order(nonmatchingNames[, 2]), ]

  if (check.nr) {
    if (nrow(nonmatchingNumbers) == 0 & nrow(nonmatchingNames) == 0) 
            cat("\n Hurray! All TaxNr <-> TaxName combinations are identical. Species lists are identical or can be used as a combined list. \n")
        else cat("\n###############################################\n!!! Reference lists are not congruent !!!\n###############################################\n")
        if (nrow(nonmatchingNumbers) > 0) {
            cat("\n", nrow(nonmatchingNumbers), "identical taxon names with different numbers \n")
            if (verbose) 
                print(nonmatchingNumbers, row.names=FALSE)
        }
  }
        if (nrow(nonmatchingNames) > 0) {
            cat("\n", nrow(nonmatchingNames), "identical taxon numbers with different names \n")
            if (verbose) 
                print(nonmatchingNames, row.names=FALSE)
        }

        reflmerge <- merge(refl.a, refl.b, by = "TaxonName", all = TRUE)
#        refl <- reflmerge[is.na(reflmerge$TaxonUsageID.x) | is.na(reflmerge$TaxonUsageID.y),]
        combnames <- reflmerge$TaxonName # cat(reflmerge$TaxonName, ' ', reflmerge$Author)
        auct <- data.frame(Taxname = sort(grep("auct.", combnames, value = TRUE, fixed = TRUE, useBytes = TRUE)))
        auct$to_check_against <- sub(" auct.", "", auct$Taxname)
        if (nrow(auct) > 0 & verbose) {
            cat("\n", "Warning: Critical Pseudonyms in dataset, please check","\n")
            print(auct, row.names=FALSE)
        }
        sl <- data.frame(Taxname = sort(grep("s. l.", combnames, value = TRUE, fixed = TRUE, useBytes = TRUE)))
        sl$to_check_against <- sub(" s. l.", "", sl$Taxname)
        sstr <- data.frame(Taxname = sort(grep("s. str.", combnames, value = TRUE, fixed = TRUE, useBytes = TRUE)))
        sstr$to_check_against <- sub(" s. str.", "", sstr$Taxname)
        ext <- rbind(sl, sstr)
        if (nrow(ext) > 0 & verbose) {
            cat("\n", "Warning: Critical names/concepts in the lists, please check", "\n")
            print(ext, row.names=FALSE)
        }
  
    if(length(diff.B) == 0 & length(diff.A) == 0)
        cat("\n Species names are identical \n")
    else {
      if(!missing(filter.1)) diff.B <- diff.B[!diff.B %in% filter.1]
      if (length(diff.B) > 0) {
            cat("\n", length(diff.B), "TaxNames of reflist nr 1 =", refl1, "not occurring in reflist nr 2 =", refl2, "\n")
            if (verbose) 
                print(diff.B, quote = FALSE, row.names=FALSE)
      }
      if(!missing(filter.2)) {
        diff.A <- diff.A[!diff.A %in% filter.2]
      }
    
      if (length(diff.A) > 0) {
          cat("\n", length(diff.A), "TaxNames of reflist nr 2 =", refl2, "not occurring in reflist nr 1 =", refl1, ": \n")
          if (verbose) 
              print(diff.A, quote = FALSE, row.names=FALSE)
      }
    }
    if (Sink) {
        tmp.wid = getOption("width")
        options(width = 5000)
        sink(file)
        print(paste(".x =", refl1, ".y =", refl2))
        if (check.nr) {
           print(paste(nrow(nonmatchingNumbers), "taxon names with different numbers"), quote = FALSE)
           print(nonmatchingNumbers, row.names=FALSE)
#         write.csv2(cbind(nonmatchingNumbers, refl.a[match(nonmatchingNumbers[,1], refl.a$TaxonName), c("BEGRUEND","EDITSTATUS")]), file='differentNumbers.csv')
           print(paste(nrow(nonmatchingNames), "taxon numbers with different names"), quote = FALSE)
           print(nonmatchingNames, row.names=FALSE)
        }
        options(width = tmp.wid)
        cat('\n', length(diff.B), "TaxNames of", refl1, "not occurring in", refl2, ':\n')
        print(paste(diff.B, collapse = ', '))
        cat('\n', length(diff.A), "TaxNames of", refl2, "not occurring in", refl1, ":\n")
        print(paste(diff.A, collapse = ', '))
        sink()
        cat("\n Report is written to file \"", file, " \n")
        if (check.nr) write.csv2(nonmatchingNumbers, file='differentNumbers.csv')
        if(!missing(filter.1)) write.csv2(nonmatchingNames[!nonmatchingNames[,2] %in% filter.1,], file='differentNames.csv')
#         write.csv2(diff.B[!diff.B %in% nonmatchingNames], file='noMatches_inRefl_2.csv')
#         write.csv2(diff.A[!diff.A %in% nonmatchingNames], file='noMatches_inRefl_1.csv')
        write.csv2(refl.a[refl.a[, "TaxonName"] %in% diff.B, c("TaxonNameOriginal","TaxonName")], file='noMatches_inRefl_2.csv')
        write.csv2(refl.b[refl.b[, "TaxonName"] %in% diff.A, c("TaxonNameOriginal","TaxonName")], file='noMatches_inRefl_1.csv')
    }
    if (new) {
      names(refl.a) <- replace(names(refl.a), names(refl.a)=='TaxonName','ABBREVIAT')
      names(refl.a) <- replace(names(refl.a), names(refl.a)=='TaxonUsageID','SPECIES_NR')
      names(refl.a) <- replace(names(refl.a), names(refl.a)=='TaxonConcept','VALID_NAME')
      names(refl.a) <- replace(names(refl.a), names(refl.a)=='TaxonConceptID','VALID_NR')
      names(refl.a) <- replace(names(refl.a), names(refl.a)=='IsChildTaxonOf','AGG_NAME')
      names(refl.a) <- replace(names(refl.a), names(refl.a)=='VernacularName','NATIVENAME')

      names(refl.b) <- replace(names(refl.b), names(refl.b)=='TaxonName','ABBREVIAT')
      names(refl.b) <- replace(names(refl.b), names(refl.b)=='TaxonUsageID','SPECIES_NR')
      names(refl.b) <- replace(names(refl.b), names(refl.b)=='TaxonConcept','VALID_NAME')
      names(refl.b) <- replace(names(refl.b), names(refl.b)=='TaxonConceptID','VALID_NR')
      names(refl.b) <- replace(names(refl.b), names(refl.b)=='IsChildTaxonOf','AGG_NAME')
      names(refl.b) <- replace(names(refl.b), names(refl.b)=='VernacularName','NATIVENAME')
      
      inter <- intersect(names(refl.a),names(refl.b))
      comb <- rbind(refl.a[,inter], refl.b[refl.b$TaxonName %in% diff.A, inter])
      comb$Attention <- comb$ABBREVIAT %in% auct | comb$ABBREVIAT %in% ext
      cat("\n New names in refl2 added to refl1. Reference list \"combrefl\" saved in TURBOVEG species directory. Please check for critical species names before use. \n")
      dir.create(file.path(tv_home, "/Species/combrefl"), showWarnings = TRUE)
      write.dbf(comb, file.path(tv_home, "/Species/combrefl/species.dbf"))
    }
}
