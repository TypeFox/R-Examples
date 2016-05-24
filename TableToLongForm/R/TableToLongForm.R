##----------------------------------------------------------------------
## The code in this .R file is machine generated from the literate
##  program, TableToLongForm.Rnw
## Documentation can be found in the literate description for this
##  program, TableToLongForm.pdf
##----------------------------------------------------------------------
TableToLongForm =
  function(Table, IdentResult = NULL,
           IdentPrimary = "combound",
           IdentAuxiliary = "sequence",
           ParePreRow = NULL,
           ParePreCol = c("mismatch", "misalign", "multirow"),
           fulloutput = FALSE,
           diagnostics = FALSE, diagnostics.trim = TRUE){
    if(is.data.frame(Table)){
      warning("Table supplied is a data.frame.\n",
              "TableToLongForm is designed for a character matrix.\n",
              "The data.frame is being coerced to a matrix but this\n",
              "may lead to unexpected results.",
              immediate. = TRUE)
      Table = as.matrix(Table)
    }
    if(!is.matrix(Table))
      stop("Table argument must be a matrix or a data.frame")
    if(diagnostics != FALSE){
      if(!is.character(diagnostics))
        diagnostics = deparse(substitute(Table))
      assign("TCRunout", file(paste0(diagnostics, ".TCRunout"), "w"),
             envir = TTLFBaseEnv)
      assign("TCtrim", diagnostics.trim, envir = TTLFBaseEnv)
      on.exit({
        with(TTLFBaseEnv, {
          close(TCRunout)
          rm(TCRunout)
          rm(TCtrim)
        })
      })
    }
    fullout = ReconsMain(matFull = Table, IdentResult,
      IdentPrimary, IdentAuxiliary, ParePreRow, ParePreCol)
    if(fulloutput) fullout else fullout$datafr
  }
rbinddf =
  function(..., deparse.level = 0){
    bindlist = list(...)
    nameunion = NULL
    for(j in 1:length(bindlist))
      nameunion = union(nameunion, colnames(bindlist[[j]]))
    for(j in 1:length(bindlist)){
      curdf = bindlist[[j]]
      namediff = setdiff(nameunion, colnames(curdf))
      matdummy = matrix(NA, nrow = nrow(curdf), ncol = length(namediff),
        dimnames = list(NULL, namediff))
      bindlist[[j]] = cbind(curdf, matdummy)
    }
    outdf = do.call(rbind,
      c(bindlist, list(deparse.level = deparse.level)))
    for(j in 1:ncol(outdf))
      if(mode(outdf[,j]) == "character") outdf[,j] = factor(outdf[,j])
    outdf
  }
print.plist =
  function(x, ...){
    plistC = function(plist){
      pLoc = attr(plist, "Loc")
      if(is.list(plist)){
        namevec = names(plist)
        if(!is.null(pLoc))
          namevec = paste0(names(plist),
            " (", pLoc[,"rows"], ", ", pLoc[,"cols"], ")")
        namelist = as.list(namevec)
        for(i in 1:length(namelist))
          namelist[[i]] =
            c(paste("+", namelist[[i]]),
              paste("-", plistC(plist[[i]])))
        do.call(c, namelist)
      } else{
        if(!is.null(names(plist))){
          namevec = names(plist)
          if(!is.null(pLoc))
            namevec = paste0(names(plist),
              " (", plist, ", ", pLoc[,"cols"], ")")
          paste("+", namevec)
        } else paste(plist, collapse = " ")
      }
    }
  
    cat(plistC(x), sep = "\n")
  }
attrLoc =
  function(plist, rows = NULL, cols = NULL){
    attr(plist, "Loc") = cbind(rows, cols)
    class(plist) = "plist"
    plist
  }
TCRsink =
  function(ID, ...)
  if(exists("TCRunout", envir = TTLFBaseEnv)){
    varlist = list(...)
    names(varlist) = gsub(" ", "", as.character(match.call()[-(1:2)]))
    TCtrim = get("TCtrim", envir = TTLFBaseEnv)
    with(TTLFBaseEnv, sink(TCRunout))
    for(i in 1:length(varlist)){
      cat("###TCR", ID, names(varlist)[i], "\n")
      curvar = varlist[[i]]
      if(TCtrim == TRUE){
        curvar = head(curvar)
        if(is.matrix(curvar) || is.matrix(curvar))
          if(ncol(curvar) > 6)
            curvar = curvar[,1:6]
      }
      print(curvar)
    }
    sink()
  }
TTLFBaseEnv = new.env()
with(TTLFBaseEnv, {aliasmat = NULL})
TTLFaliasAdd =
  function(Type, Fname, Falias, Author = "", Description = "")
  assign("aliasmat",
         rbind(get("aliasmat", envir = TTLFBaseEnv),
               c(Type = Type, Name = Fname, Alias = Falias,
                 Author = Author, Description = Description)),
         envir = TTLFBaseEnv)

TTLFaliasGet =
  function(Type, Falias){
    aliasmat = get("aliasmat", envir = TTLFBaseEnv)
    matchRow = which(aliasmat[,"Type"] == Type &
      aliasmat[,"Alias"] == Falias)
    if(length(matchRow) == 1)
      aliasmat[matchRow,"Name"]
    else stop("Invalid algorithm specified for ", Type)
  }

TTLFaliasList =
  function(){
    aliasmat = get("aliasmat", envir = TTLFBaseEnv)
    Types = unique(aliasmat[,"Type"])
    for(Type in Types){
      cat("==Type: ", Type, "==\n", sep = "")
      Algos = aliasmat[aliasmat[,"Type"] == Type,,drop=FALSE]
      for(i in 1:nrow(Algos))
        cat("Name: ", Algos[i, "Name"], "\n",
            "Alias: ", Algos[i, "Alias"], "\n",
            "Author: ", Algos[i, "Author"], "\n",
            "Description: ", Algos[i, "Description"], "\n\n",
            sep = "")
    }
  }
  
IdentbyMostCommonBoundary =
  function(matFull){
    rowNonempty = (1:nrow(matFull))[IdentNonEmpty(matFull, 1)]
    colNonempty = (1:ncol(matFull))[IdentNonEmpty(matFull, 2)]
    rowData = IdentMostCommonBoundary(matFull, 2)
    colData = IdentMostCommonBoundary(matFull, 1)
    TCRsink("CIMCB", rowData, colData)
    rowslist = list(label = rowNonempty[rowNonempty < rowData[1]],
                    data = rowNonempty[(rowNonempty >= rowData[1]) &
                                       (rowNonempty <= rowData[2])])
    colslist = list(label = colNonempty[colNonempty < colData[1]],
                    data = colNonempty[(colNonempty >= colData[1]) &
                                       (colNonempty <= colData[2])])
    TCRsink("CRAC", rowslist, colslist)
    matRowLabel = matFull[rowslist$data, colslist$label,drop=FALSE]
    if(!all(is.na(matRowLabel)) && ncol(matRowLabel) > 1){
      RowLabelNonempty = IdentNonEmpty(matRowLabel, 2)
      if(max(RowLabelNonempty) < ncol(matRowLabel)){
        toshift = (max(RowLabelNonempty) + 1):ncol(matRowLabel)
        colslist$data = c(colslist$label[toshift], colslist$data)
        colslist$label = colslist$label[-toshift]
      }
    }
    list(rows = rowslist, cols = colslist)
  }
TTLFaliasAdd("IdentPrimary", "IdentbyMostCommonBoundary", "combound",
             "Base Algorithm", "Default IdentPrimary algorithm")
IdentbySequence =
  function(matFull, IdentResult)
  with(IdentResult, {
    matRowLabel = matFull[rows$data, cols$label]
    if(all(is.na(matRowLabel))){
      cols$label = cols$data[1]
      cols$data = cols$data[-1]
      IdentbySequence(matFull, list(rows = rows, cols = cols))
    }
    else{
      matRowLabel = suppressWarnings(as.numeric(matRowLabel))
      if(length(unique(matRowLabel)) > 1 &&
         length(unique(diff(matRowLabel))) == 1)
        list(rows = rows, cols = cols)
      else IdentResult
    }
  })
TTLFaliasAdd("IdentAuxiliary", "IdentbySequence", "sequence",
             "Base Algorithm", paste("Search for fully numeric row",
             "labels (e.g. Years) that were misidentified as data"))
IdentNonEmpty =
  function(mat, margin, emptyident = is.na){
    isnonempty = apply(mat, margin, function(x) !all(emptyident(x)))
    which(isnonempty)
  }
IdentPattern =
  function(vec){
    matchvec = match(vec, unique(vec))
    for(i in 1:length(unique(matchvec))){
      repind = unique(diff(which(matchvec == i)))
      if(length(repind) == 0)
        repind = length(vec)
      if(length(repind) == 1)
        break
    }
    curseg = paste0("^(", paste(vec[1:repind], collapse = ""), ")+$")
    if(length(grep(curseg, paste(vec, collapse = ""))) > 0)
      repind else length(vec)
  }
IdentMostCommonBoundary =
  function(matFull, margin){
    isnumber = suppressWarnings(apply(matFull, margin,
      function(x) which(!is.na(as.numeric(x)))))
    nstarts = table(sapply(isnumber,
      function(x) if(length(x) > 0) min(x) else NA))
    nends = table(sapply(isnumber,
      function(x) if(length(x) > 0) max(x) else NA))
    as.numeric(names(c(which.max(nstarts), which.max(rev(nends)))))
  }
## Empty
ParePreColMismatch =
  function(matData, matColLabel){
    colsData = IdentNonEmpty(matData, 2)
    colsLabels = IdentNonEmpty(matColLabel, 2)
    if(length(colsData) == length(colsLabels))
      if(ncol(matData) != length(colsData)){
        matColLabel = matColLabel[,colsLabels,drop=FALSE]
        matData = matData[,colsData,drop=FALSE]
      }
    list(matData = matData, matColLabel = matColLabel)
  }
TTLFaliasAdd("ParePreCol", "ParePreColMismatch", "mismatch",
             "Base Algorithm", paste("Correct for column labels",
             "not matched correctly over data (label in a",
             "different column to data)"))
ParePreColMisaligned =
  function(matData, matColLabel){
    TCRsink("MCPBefore", matColLabel)
    for(i in 1:nrow(matColLabel)){
      currow = matColLabel[i,]
      curPattern =
        if(all(is.na(currow))) NA
        else if(any(is.na(currow))) IdentPattern(is.na(currow))
        else IdentPattern(currow)
      if(!is.na(curPattern)){
        nParents = length(currow)/curPattern
        for(j in 1:nParents){
          curcols = 1:curPattern + curPattern * (j - 1)
          cursub = currow[curcols]
          currow[curcols] = c(cursub[!is.na(cursub)], cursub[is.na(cursub)])
          TCRsink("ACP", cursub, currow[curcols])
        }
        matColLabel[i,] = currow
      }
    }
    TCRsink("MCPAfter", matColLabel)
    list(matData = matData, matColLabel = matColLabel)
  }
TTLFaliasAdd("ParePreCol", "ParePreColMisaligned", "misalign",
             "Base Algorithm", paste("Correct for column labels",
             "not aligned correctly over data (parents not",
             "positioned on the far-left, relative to their",
             "children in the row below)"))
ParePreColMultirow =
  function(matData, matColLabel){
    fullrows = apply(matColLabel, 1, function(x) all(!is.na(x)))
    if(any(diff(fullrows) > 1))
      warning("full rows followed by not full rows!")
    if(any(fullrows)){
      pastestring = ""
      pasterows = which(fullrows)
      for(i in 1:length(pasterows))
        pastestring[i] = paste0("matColLabel[", pasterows[i],
                     ",,drop=FALSE]")
      collapsedlabels =
        eval(parse(text = paste0("paste(",
                     paste(pastestring, collapse = ", "), ")")))

      TCRsink("MCLBefore", matColLabel)
      matColLabel = rbind(matColLabel[!fullrows,,drop=FALSE],
        collapsedlabels, deparse.level = 0)
      TCRsink("MCLAfter", matColLabel)
    }
    list(matData = matData, matColLabel = matColLabel)
  }
TTLFaliasAdd("ParePreCol", "ParePreColMultirow", "multirow",
             "Base Algorithm", paste("Merge long column labels",
             "that were physically split over multiple rows",
             "back into a single label"))
PareFront =
  function(matLabel)
  PareMain(matSub = matLabel, plist =
           list(rows = 1:nrow(matLabel), cols = 1:ncol(matLabel)))
PareMain =
  function(matSub, plist){
    if(length(plist$cols) == 1){
      res = structure(plist$rows, .Names = matSub[plist$rows, plist$cols])
      res = attrLoc(res, cols = plist$col)
      TCRsink("IOOC", plist, res)
    }
    else if(all(is.na(matSub[plist$rows, plist$cols[1]]))){
      plist$cols = plist$cols[-1]
      res = PareMain(matSub, plist)
    }
    else if(length(plist$rows) == 1){
      res = structure(plist$rows,
        .Names = matSub[plist$rows, plist$cols[length(plist$cols)]])
      res = attrLoc(res, cols = plist$cols[length(plist$cols)])
      for(i in (length(plist$cols) - 1):1){
        res = list(res)
        names(res) = matSub[plist$rows, plist$cols[i]]
        res = attrLoc(res, rows = plist$rows, cols = plist$cols[i])
      }
      TCRsink("IOOR", plist, res)
    }
    else if(is.na(matSub[plist$rows[1], plist$cols[1]])){
      warning("cell[1, 1] is empty")
      print(plist)
      res = NA
    }
    else{
      res = PareByEmptyRight(matSub, plist)
      if(any(is.na(res)))
        res = PareByEmptyBelow(matSub, plist)
      for(i in 1:length(res))
        res[[i]] = PareMain(matSub, res[[i]])
      res
    }
    class(res) = "plist"
    res
  }
PareByEmptyRight =
  function(matSub, plist)
  with(plist,
       if(all(is.na(matSub[rows[1], cols[-1]]))){
         emptyrights = apply(matSub[rows, cols[-1],drop=FALSE], 1,
           function(x) all(is.na(x)))
         rowemptyright = rows[emptyrights]
         if(length(rowemptyright) == 1){
           res = list(list(rows = rows[-1], cols = cols))
           names(res) = matSub[rows[1], cols[1]]
           res = attrLoc(res, rows = rows[1], cols = cols[1])
           TCRsink("CSER", res)
         }
         else{
           rowdiff = diff(rowemptyright)
           if(any(rowdiff == 1))
             rowemptyright = rowemptyright[c(rowdiff == 1, FALSE)]
           
           rowstart = pmin(rowemptyright + 1, max(rows))
           rowend = c(pmax(rowemptyright[-1] - 1, min(rows)), max(rows))
           
           res = list()
           for(i in 1:length(rowstart))
             res[i] = list(list(rows = rowstart[i]:rowend[i], cols = cols))
           names(res) = matSub[rowemptyright, cols[1]]
           res = attrLoc(res, rows = rowemptyright, cols = cols[1])
           TCRsink("CMER", res)
         }
         res
       } else NA)
PareByEmptyBelow =
  function(matSub, plist)
  with(plist, {
    emptybelow = is.na(matSub[rows, cols[1]])
    rowstart = rows[!emptybelow]
    rowend = c(rowstart[-1] - 1, max(rows))
    res = list()
    for(i in 1:length(rowstart))
      res[i] = list(list(rows = rowstart[i]:rowend[i], cols = cols[-1]))
    names(res) = matSub[rowstart, cols[1]]
    res = attrLoc(res, rows = rowstart, cols = cols[1])
    TCRsink("PBEB", res)
    res
  })
ReconsMain =
  function(matFull, IdentResult,
           IdentPrimary, IdentAuxiliary,
           ParePreRow, ParePreCol){
    if(is.null(IdentResult)){
      IdentPrimary = TTLFaliasGet("IdentPrimary", IdentPrimary)
      IdentResult = do.call(IdentPrimary, list(matFull = matFull))
      if(!is.null(IdentAuxiliary))
        for(AuxAlgo in IdentAuxiliary){
          AuxAlgo = TTLFaliasGet("IdentAuxiliary", AuxAlgo)
          IdentResult = do.call(AuxAlgo,
            list(matFull = matFull, IdentResult = IdentResult))
        }
    }
    matData = with(IdentResult,
      matFull[rows$data, cols$data,drop=FALSE])
    matRowLabel = with(IdentResult,
      matFull[rows$data, cols$label,drop=FALSE])
    if(!is.null(ParePreRow))
      for(PreAlgo in ParePreRow){
        PreAlgo = TTLFaliasGet("ParePreRow", PreAlgo)
        PreOut = do.call(PreAlgo,
          list(matData = matData, matRowLabel = matRowLabel))
        matData = PreOut$matData
        matRowLabel = PreOut$matRowLabel
      }
    rowplist = PareFront(matRowLabel)
    rowvecs = ReconsRowLabels(rowplist)
    TCRsink("RRL", rowplist, rowvecs)
    matColLabel = with(IdentResult,
      matFull[rows$label, cols$data,drop=FALSE])
    if(!is.null(ParePreCol))
      for(PreAlgo in ParePreCol){
        PreAlgo = TTLFaliasGet("ParePreCol", PreAlgo)
        PreOut = do.call(PreAlgo,
          list(matData = matData, matColLabel = matColLabel))
        matData = PreOut$matData
        matColLabel = PreOut$matColLabel
      }
    colplist = PareFront(t(matColLabel))
    matDataReduced = matData[unlist(rowplist),,drop=FALSE]
    res = ReconsColLabels(colplist, matDataReduced, rowvecs)
    TCRsink("RCL", colplist, res)
    list(datafr = res, oriTable = matFull, IdentResult = IdentResult,
         rowplist = rowplist, colplist = colplist)
  }
ReconsRowLabels =
  function(plist)
  if(is.list(plist)){
    rowvecs = as.list(names(plist))
    for(i in 1:length(rowvecs))
      rowvecs[[i]] = cbind(rowvecs[[i]], ReconsRowLabels(plist[[i]]))
    do.call(rbind, rowvecs)
  } else as.matrix(names(plist))
ReconsColLabels =
  function(plist, matData, rowvecs){
    if(is.list(plist)){
      colvecs = as.list(names(plist))
      for(i in 1:length(colvecs)){
        colvecs[[i]] = cbind(colvecs[[i]],
                 ReconsColLabels(plist[[i]], matData, rowvecs))
        colnames(colvecs[[i]])[1] = "UNKNOWN"
      }
      datfr = do.call(rbinddf, colvecs)
    }
    else{
      datbit = matData[,plist,drop=FALSE]
      TCRsink("RCC", plist, matData, datbit)
      datlist = NULL
      for(j in 1:ncol(datbit)){
        asnumer = suppressWarnings(as.numeric(datbit[,j]))
        if(all(is.na(datbit[,j])) || !all(is.na(asnumer)))
          datlist[[j]] = asnumer
        else
          datlist[[j]] = datbit[,j]
      }
      datbit = do.call(cbind, datlist)
      ## Specify row.names to avoid annoying warnings
      datfr =
        cbind(as.data.frame(rowvecs, row.names = 1:nrow(rowvecs)), datbit)
      colnames(datfr) =
        c(rep("UNKNOWN", length = ncol(rowvecs)), names(plist))
    }
    datfr
  }
