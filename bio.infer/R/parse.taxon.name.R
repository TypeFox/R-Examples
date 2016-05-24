# 1.10.2009
# parse vector of taxon names into distinct strings
#

load.itis <- function(get.tax.env) {
  data(itis.ttable, envir = get.tax.env)
  tlevs <- names(get.tax.env[["itis.ttable"]])
  tlevs <- tlevs[-match("TAXON", tlevs)]
  return(tlevs)
}

parse.taxon.name <- function(tname.orig) {
    tname <- toupper(tname.orig)
    substr <- as.list(rep(NA, times = 1))
    i <- 1
    w1 <- regexpr("\\(", tname)
    w2 <- regexpr("\\)", tname)
    incvec <- (w1 != -1) & (w2 != -1)
    tname[incvec] <- paste(substring(tname[incvec], 1, w1[incvec] -
                                     1), substring(tname[incvec], w2[incvec] + 1, nchar(tname[incvec])))
    repeat {
      w <- regexpr("[A-Z]+", tname)
      if (sum(w != -1) == 0)
        break
      substr[[i]] <- substring(tname, w, w + attributes(w)$match.length -
                               1)
      w3 <- w + attributes(w)$match.length
      tname <- substring(tname, w3, nchar(tname))
      if (sum(tname != "") == 0)
        break
      i <- i + 1
    }
    substr[[i]] <- rep("", times = length(tname)) # add vector of blanks at the end
    df.parse <- matrix("", ncol = length(substr) + 1, nrow = length(tname))
    df.parse[, 1] <- tname.orig
    exlist <- c("DUPLICATE", "SETAE", "CODE", "GROUP", "TYPE",
        "GENUS", "PANEL", "SAND", "TURRET", "CASE", "LARVAE",
        "ADULT", "SENSU", "TZING", "RIBAUD", "RPEL", "STRUP", "NAWQA",
                "LLER","KANSSON", "UMICH", "ALBE" )
    for (i in 1:length(substr)) {
      incvec <- substr[[i]] %in% exlist | ((nchar(substr[[i]]) <= 3) &
                                           (nchar(substr[[i]]) > 0))
      incvec[is.na(incvec)] <- F
      while(sum(incvec) > 0) {
        if (i < (length(substr)-1)) {
          for (j in i:(length(substr)-1)) {
            substr[[j]][incvec] <- substr[[j+1]][incvec]
          }
        }
        else {
          substr[[i]][incvec] <- ""
        }
        incvec <- substr[[i]] %in% exlist | ((nchar(substr[[i]]) <= 3) &
                                             (nchar(substr[[i]]) > 0))
        incvec[is.na(incvec)] <- F
      }
    }

    for (i in 1:length(substr)) {
      df.parse[, i + 1] <- substr[[i]]
    }
    df.parse <- data.frame(df.parse, stringsAsFactors = FALSE)
    return(df.parse)
}

in.ITIS <- function(df.parse,get.tax.env, col.sel = NULL) {
  if (is.vector(df.parse)) {
    in.mat <- df.parse %in% get.tax.env[["itis.ttable"]]$TAXON
  }
  else {
    if (is.null(col.sel)) {
      in.mat <- matrix(NA, ncol = ncol(df.parse)-1, nrow = nrow(df.parse))
      for (i in 2:ncol(df.parse)) {
        in.mat[,i-1] <- df.parse[,i] %in% get.tax.env[["itis.ttable"]]$TAXON
      }
    }
    else {
      in.mat <- rep(NA, times = nrow(df.parse))
      in.mat <- df.parse[, col.sel] %in% get.tax.env[["itis.ttable"]]$TAXON
    }
  }
  return(in.mat)
}


resolve.mult <- function(parse.list, get.tax.env) {

  if (ncol(parse.list[[2]]) > 2) {
    in.mat <- in.ITIS(parse.list[[2]], get.tax.env)
    selvec <- apply(in.mat,1, sum) > 1
    if (sum(selvec) > 0) {
      df.parse <- parse.list[[2]][selvec,]
      tlevs <- names(get.tax.env[["itis.ttable"]])
      imatch <- match("TAXON", tlevs)
      tlevs <- tlevs[-imatch]

      imatch <- match("SUBCLASS", toupper(tlevs))
      tlevs.loc <- rev(tlevs[length(tlevs):imatch])

      imatch1 <- match(df.parse[,2], get.tax.env[["itis.ttable"]]$TAXON)
      imatch2 <- match(df.parse[,3], get.tax.env[["itis.ttable"]]$TAXON)

      comp1 <- get.tax.env[["itis.ttable"]][imatch1, tlevs.loc]
      comp2 <- get.tax.env[["itis.ttable"]][imatch2, tlevs.loc]
      tlev.sav <- rep(0, times = nrow(df.parse))
      tlev.o1 <- rep(0, times = nrow(df.parse))
      tlev.o2 <- rep(0, times = nrow(df.parse))
      for (j in 1:length(tlevs.loc))  {
        incvec <- comp1[j] == comp2[j]
        incvec[is.na(incvec)] <- FALSE
        tlev.sav[incvec] <- j
        tlev.o1[!is.na(comp1[j])] <- j
        tlev.o2[!is.na(comp2[j])] <- j
      }

      incvec1 <- tlev.sav != 0
      incvec2 <- abs(tlev.o1 - tlev.sav) < 5
      incvec.all <- incvec1 & incvec2

      x <- tlevs.loc[tlev.sav[incvec.all]]
      ind <- imatch1[incvec.all]
      str.save <- rep("", times = length(x))
      for (i in 1:length(x)) {
        str.save[i] <- get.tax.env[["itis.ttable"]][ind[i], x[i]]
      }
      df.parse[incvec.all,2] <- str.save
      parse.list[[1]] <- rbind(parse.list[[1]], df.parse)
      parse.list[[2]] <- parse.list[[2]][! selvec,]
    }
  }
  return(parse.list)
}

# merge taxalist with full itis.ttable
make.fulltab1 <- function(df.parse, get.tax.env) {
  df.parse <- data.frame(df.parse, stringsAsFactors = FALSE)
  names0 <- paste("t", 1:ncol(df.parse), sep = "")
  names0[1] <- "taxaname.orig"
  names0[2] <- "taxaname.merge"
  names(df.parse) <- names0
  df1 <- data.frame(TAXON = I(unique(df.parse$taxaname.merge)))
  fulltab <- merge(get.tax.env[["itis.ttable"]], df1,by = "TAXON",all.y = TRUE)
  return(fulltab)
}

# make species names
make.species <- function(df.parse, fulltab) {
    df.parse <- data.frame(df.parse, stringsAsFactors = FALSE)
    names0 <- names(df.parse)
    names0 <- paste("t", 1:length(names0), sep = "")
    names0[1] <- "taxaname.orig"
    names0[2] <- "TAXON"
    names(df.parse) <- names0
    df1 <- merge(df.parse, fulltab, by = "TAXON")
    incvec <- !is.na(df1$GENUS)
    df1$SPECIES <- rep(NA, times = nrow(df1))
    incvec2 <- df1$t3 != ""
    incvec.a <- incvec & incvec2
    df1$SPECIES[incvec.a] <- paste(df1$GENUS[incvec.a], df1$t3[incvec.a],
        sep = ".")
    df1$TAXON[incvec.a] <- df1$SPECIES[incvec.a]
    npos <- sum((nchar(names(df1)) == 2) & (substring(names(df1),
        1, 1) == "t"))
    tname.orig.cap <- toupper(df1$taxaname.orig)
    if (npos > 1) {
        for (i in 2:npos) {
            fname <- paste("t", i + 2, sep = "")
            incvec2 <- df1[, fname] != ""
            w <- rep(-1, times = length(tname.orig.cap))
            for (j in 1:length(tname.orig.cap)) {
              w[j] <- regexpr(df1[j, fname], tname.orig.cap[j])
            }
            # check to see if substring is capitalized in taxaname.orig
            # if so, do not add as a species name
            incvec3 <- regexpr("[A-Z]+", substring(df1$taxaname.orig,w,w)) == -1
            incvec.a <- incvec & incvec2 & incvec3
            df1$SPECIES[incvec.a] <- paste(df1$SPECIES[incvec.a],
                df1[incvec.a, fname], sep = "/")
            df1$TAXON[incvec.a] <- df1$SPECIES[incvec.a]
        }
    }
    return(df1)
}

# identify duplicates and return row numbers and summary string
locate.dupes <- function(fulltab) {
  t.sel <- unique(fulltab$TAXON[duplicated(fulltab$TAXON)])
  sumstr <- character(0)
  isav <- numeric(0)
  if (length(t.sel) > 0) {
    for (i in 1:length(t.sel)) {
      isav <- c(isav, (1:nrow(fulltab))[t.sel[i] == fulltab$TAXON])
    }
    for (i in 1:length(isav)) {
      stringtemp <- fulltab[isav[i], c("TAXON", "FAMILY", "ORDER",
                                                  "CLASS", "PHYLUM")]
      incvec <- ! is.na(stringtemp)
      stringtemp <- stringtemp[incvec]
      sumstr <- c(sumstr,paste(stringtemp, collapse = "-"))
    }
  }
  return(list(sumstr = sumstr, isav = isav))
}

get.dupe.sel <- function(sumstr) {
  isel <- numeric(0)
  if (length(sumstr > 0)) {
    w <- regexpr("-", sumstr)
    TAXON <- substring(sumstr, 1, w - 1)
    selvec <- TAXON == ""
    TAXON <- TAXON[!selvec]
    ntax <- length(unique(TAXON))
    repeat {
#      a <- tklist.modal("Select appropriate taxon", sumstr,
#                        selectmode = "multiple")
        a <- select.list(sumstr, multiple = TRUE,
                         title = "Select appropriate taxon")
      if (length(a) == ntax) {
        for (i in 1:length(a)) {
          isel <- c(isel, match(a[i], sumstr))
        }
        if (sum(duplicated(TAXON[isel])) > 0) {
          cat("Please select only one choice per taxon name\n")
          flush.console()
          isel <- numeric(0)
        }
        else break
      }
      else {
        cat("Please select one choice per taxon name\n")
        flush.console()
      }
    }
  }
  return(a)
}

remove.dupes <- function(fulltab, dupe.list, dupe.sel ) {

  isav <- dupe.list$isav
  if (length(dupe.sel) > 0) {
    for (i in 1:length(dupe.sel)) {
      isel <- match(dupe.sel[i], dupe.list$sumstr)
      isav <- isav[-isel]
    }
    fulltab <- fulltab[-isav,]
  }
  return(fulltab)
}


get.taxon.names <- function(bcnt) {
  names0 <- names(bcnt)
  f.tname <- names0[2]
  if (is.factor(bcnt[, 2])) {
    tname.a <- sort(unique(levels(bcnt[, f.tname])[bcnt[, f.tname]]))
  }
  else {
    if (is.character(bcnt[, f.tname])) {
      tname.a <- sort(unique(bcnt[, f.tname]))
    }
    else {
      tkmessageBox(message = "2nd field is neither factor nor character",
                   icon = "error", type = "ok")
    }
  }
  tname <- unique(tname.a)
  return(tname)
}


get.valid.names <- function(df.parse, get.tax.env)  {
  in.mat <- in.ITIS(df.parse, get.tax.env)
  tempmat <- in.mat
  if (ncol(tempmat) > 1) {
    for (i in 2:ncol(tempmat)) {
      tempmat[,i] <- ! tempmat[,i]
    }
    incvec <- apply(tempmat,1, all)
  }
  else {
    incvec <- tempmat[,1]
  }
  incvec <- !incvec # reverse trues and falses to get parselist ordered correctly
  parse.list <- split(df.parse, incvec)
#  df.parse.sav <- data.frame(df.parse[incvec,])
#  df.parse <- data.frame(df.parse[!incvec,])
  return(parse.list)
}

correct.taxanames <- function(tname.old, get.tax.env) {
  # prompt for taxaname corrections
  # This version checks each name for presence in ITIS
  # The java version will not
  tname.new <- tname.old
  dfedit <- data.frame(Orig.Taxon.Name = tname.old, Rev.Taxon.Name = tname.new,
                       stringsAsFactors = F)
  repeat {
#    tname.new <- modalDialog("Enter corrections", tname.old,
#                              tname.new, returnValOnCancel = tname.old)
      dfedit <- fix(dfedit)
      tname.new <- toupper(as.character(dfedit[,2]))
#    tname.new <- toupper(tname.new)
    in.mat.correct <- in.ITIS(tname.new, get.tax.env)
    if (sum(! in.mat.correct) > 0) {
      cat("The following taxa are not in ITIS:\n")
      cat(tname.new[!in.mat.correct], sep = "\n")
      flush.console()
    }
    done0 <- tkmessageBox(message = "Are you done editing?",
                          icon = "question", type = "yesno", default = "yes")
    if (as.character(done0) == "yes")  break
  }
  return(tname.new)
}


output.tax.table <- function(finaltab, tlevs) {
  tlevs.loc <- c(tlevs, "SPECIES", "TAXON", "taxaname.orig")

  df1 <- finaltab[, tlevs.loc]
  df1 <- df1[do.call(order, df1), ]
  write.table(df1, sep = "\t", file = "sum.tax.table.txt", row.names = FALSE)
}

incorp.correct <- function(tname.new, parse.list) {
  tname.new2 <- toupper(tname.new)
  tname.old <- sort(unique(parse.list[[2]][,2]))
  for (i in 1:length(tname.old)) {
    incvec <- tname.old[i] == parse.list[[2]][,2]
    parse.list[[2]][incvec,2] <- tname.new2[i]
  }
  parse.list[[1]] <- rbind(parse.list[[1]], parse.list[[2]])
  parse.list[[2]] <- parse.list[[2]][rep(F,times = nrow(parse.list[[2]])),]
  return(parse.list)
}


