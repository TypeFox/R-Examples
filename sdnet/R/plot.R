#########################################################################
# Plotting Methods  

setMethod("cnDot", "catNetwork", function(object, file="", format="ps", nodestyle=NULL, edgestyle=NULL) {
  ## style format is a list of node[shape, shape.color, edge.color] 
  if(nchar(object@meta)>0) 
    str <- sprintf("\"%s, \\nComplexity %d, \\nLogLikelihood %5.3f\"[shape=plaintext]\n", 
                   as.character(object@meta), object@complx, object@loglik) 
  else
    str <- ""
  pmat <- cnMatParents(object)
  if(!is.null(nodestyle)) {
    strout <- sapply(1:object@numnodes, function(n) {
      if(sum(pmat[n,]) + sum(pmat[,n]) <= 0)
        return("")
      if(!is.null(nodestyle) && length(nodestyle)>=n && length(nodestyle[[n]])>=3)
        paste("\"", object@nodes[[n]], "\"[shape=", nodestyle[[n]][1], ", style=filled, fillcolor=", nodestyle[[n]][2], ", fontcolor=", nodestyle[[n]][3], "];\n", collapse="", sep="")
    })
    str <- paste(str, paste(strout, collapse="", sep=""))
  }
  noedges <- TRUE
  strout <- sapply(seq(1, length(object@pars)), function(n) {
    if(is.null(object@pars[[n]]))
        return("")
    if(!is.null(object@pars[[n]])) {
      noedges <- FALSE
      paste(sapply(object@pars[[n]], function(j) { 
        if(!is.null(edgestyle) && length(edgestyle)>=n && length(edgestyle[[n]])>=3) {
          if(length(object@pars[[j]]) > 0 && length(which(object@pars[[j]] == n)) > 0 )
            paste("\"", object@nodes[[n]], "\"[shape=\"", edgestyle[[n]][1], "\", color=\"", edgestyle[[n]][2], "\"];\n",
                  "edge[color=\"", edgestyle[[n]][3], "\"];\n", 
                  "\"", object@nodes[j], "\" -> \"", object@nodes[n], "\" [style=dashed];\n", collapse="", sep="")
          else
            paste("\"", object@nodes[[n]], "\"[shape=\"", edgestyle[[n]][1], "\", color=\"", edgestyle[[n]][2], "\"];\n",
                  "edge[color=\"", edgestyle[[n]][3], "\"];\n",
                  "\"", object@nodes[j], "\" -> \"", object@nodes[n], "\";\n", collapse="", sep="")
        }
        else {
          if(length(object@pars[[j]]) > 0 && length(which(object@pars[[j]] == n)) > 0) 
            paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\" [style=dashed];\n", collapse="", sep="") 
          else 
            paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\";\n", collapse="", sep="")
        }
      }), collapse="", sep="") 
    } 
  }) 
  strout <- paste(str, paste(strout, collapse="", sep="")) 
  str <- paste("digraph G {\n", strout, "}\n", collapse="", sep="")

  if(is.null(format))
    format <- ""
  
  if(!missing(file) && !is.null(file)) { 
 
    ## get the full path to the file 
    file <- paste(getwd(), "/", file, sep="") 
 
    write(str, file=paste(file,".dot",sep="")) 
 
    dotviewer <- as.character(Sys.getenv("R_DOTVIEWER")) 
    if(dotviewer != "" && (format == "ps" || format == "pdf")) {
      if(format == "ps") 
        strdotcall<-paste(dotviewer, " -Tps \"", file, ".dot\"", " -o \"", file, ".ps\"", sep="")
      else
        strdotcall<-paste(dotviewer, " -Tpdf \"", file, ".dot\"", " -o \"", file, ".pdf\"", sep="") 
      try(system(strdotcall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
 
      pdfviewer <- as.character(Sys.getenv("R_PDFVIEWER")) 
      if(pdfviewer != "") {
        if(format == "ps") 
          strevincecall<-paste(pdfviewer, " \"", file, ".ps\"", sep="")
        else
          strevincecall<-paste(pdfviewer, " \"", file, ".pdf\"", sep="") 
        try(system(strevincecall, intern=FALSE, wait=FALSE, ignore.stderr=TRUE), silent = TRUE) 
      } 
    } 
  } 
  if(format == "dot")
    cat(str) 
  }) 
 
 
setMethod("cnDot", "list", function(object, file="", format="ps",  nodestyle=NULL, edgestyle=NULL) { 
  if(!is.list(object)) 
    return("") 
  objectlist <- object 
  liststr <- "" 
  i <- 1 
  for(object in objectlist) {

    if(is(object, "catNetwork")) {
      
      str <- sprintf("\"%s, \\nComplexity %d, \\nLogLikelihood %5.3f\"[shape=plaintext]", 
                     as.character(object@meta), object@complx, object@loglik) 
      strout <- sapply(seq(1, length(object@pars)), function(n) { 
        if(is.null(object@pars[[n]])) {
          warning("network without edges")
          return("")
        }
        else{ 
          paste(sapply(object@pars[[n]], function(j) { 
            if(length(object@pars[[j]]) > 0 && length(which(object@pars[[j]] == n)) > 0) { 
              paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\" [style=dashed];\n", collapse="", sep="") 
            } 
            else 
              paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\";\n", collapse="", sep="") 
          }), collapse="", sep="") 
        } 
      }) 
      strout <- paste(str, paste(strout, collapse="", sep="")) 
      str <- paste("digraph ", sprintf("G%d", i), "{\n", strout, "};\n", collapse="", sep="")
    } ## catNetwork

    if(is.matrix(object)) {
      
      medges <- as.matrix(object) 
      if(dim(medges)[1] != dim(medges)[2] || dim(medges)[1] < 2) {
        warning("Wrong matrix")
        next
      }
      rnames <- rownames(medges) 
      if(is.null(rnames)) 
        rnames <- 1:dim(medges)[1] 
      nnodes <- dim(medges)[1] 
      strout <- "" 
      for(row in 1:nnodes) { 
        for(col in 1:nnodes) { 
          if(medges[row,col] <= 0) 
            next 
          if(medges[col,row] > 0)  
            strout <- paste(strout, rnames[col], " -> ", rnames[row], " [style=dashed];\n", collapse="", sep="") 
          else 
            strout <- paste(strout, rnames[col], " -> ", rnames[row], ";\n", collapse="", sep="") 
        } 
      }
      str <- paste("digraph G {\n", strout, "}\n", collapse="", sep="")
    }
    
    liststr <- paste(liststr, str, "", sep="") 
    i <- i + 1 
  }

  if(is.null(format))
    format <- ""
  
  if(!missing(file) && !is.null(file)) { 
 
    ## get the full path to the file 
    file <- paste(getwd(), "/", file, sep="") 
 
    write(liststr, file=paste(file,".dot",sep="")) 
 
    dotviewer <- as.character(Sys.getenv("R_DOTVIEWER")) 
    if(dotviewer != "" && (format == "ps" || format == "pdf")) {
      if(format == "ps") 
        strdotcall<-paste(dotviewer, " -Tps \"", file, ".dot\"", " -o \"", file, ".ps\"", sep="")
      else
        strdotcall<-paste(dotviewer, " -Tpdf \"", file, ".dot\"", " -o \"", file, ".pdf\"", sep="") 
      try(system(strdotcall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
     
      pdfviewer <- as.character(Sys.getenv("R_PDFVIEWER")) 
      if(pdfviewer != "") {
        if(format == "ps") 
          strevincecall<-paste(pdfviewer, " \"", file, ".ps\"", sep="")
        else
          strevincecall<-paste(pdfviewer, " \"", file, ".pdf\"", sep="") 
        try(system(strevincecall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
      } 
    } 
  } 
  if(format == "dot")
    cat(liststr) 
}) 
 
setMethod("cnDot", "matrix", function(object, file="", format="ps",  nodestyle=NULL, edgestyle=NULL) { 
  if(!is(object, "matrix")) 
    stop("Specify a valid square matrix.") 
  medges <- as.matrix(object) 
  if(dim(medges)[1] != dim(medges)[2] || dim(medges)[1] < 2) 
    stop("Specify a valid square matrix.") 
  rnames <- rownames(medges) 
  if(is.null(rnames)) 
    rnames <- 1:dim(medges)[1] 
  nnodes <- dim(medges)[1] 
  strout <- "" 
  for(row in 1:nnodes) { 
    for(col in 1:nnodes) { 
      if(medges[row,col] <= 0) 
        next
      ##cat(rnames[col], " -> ", rnames[row], "\n") 
      if(medges[col,row] > 0)  
        ## double-edge in both directions 
        ##strout <- paste(strout, rnames[row], " -> ", rnames[col], " [dir=both, style=dashed];\n", collapse="", sep="")
        strout <- paste(strout, rnames[col], " -> ", rnames[row], " [style=dashed];\n", collapse="", sep="") 
      else 
        strout <- paste(strout, rnames[col], " -> ", rnames[row], ";\n", collapse="", sep="") 
    } 
  }
  str <- paste("digraph G {\n", strout, "}\n", collapse="", sep="")

  if(is.null(format))
    format <- ""
  
  if(!missing(file) && !is.null(file)) { 
 
    ## get the full path to the file 
    file <- paste(getwd(), "/", file, sep="") 
 
    write(str, file=paste(file,".dot",sep="")) 
 
    dotviewer <- as.character(Sys.getenv("R_DOTVIEWER"))
    if(dotviewer != "" && (format == "ps" || format == "pdf")) {
      if(format == "ps")
        strdotcall<-paste(dotviewer, " -Tps \"", file, ".dot\"", " -o \"", file, ".ps\"", sep="")
      else
        strdotcall<-paste(dotviewer, " -Tpdf \"", file, ".dot\"", " -o \"", file, ".pdf\"", sep="")
      try(system(strdotcall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
     
      pdfviewer <- as.character(Sys.getenv("R_PDFVIEWER")) 
      if(pdfviewer != "") {
        if(format == "ps") 
          strevincecall<-paste(pdfviewer, " \"", file, ".ps\"", sep="")
        else
          strevincecall<-paste(pdfviewer, " \"", file, ".pdf\"", sep="")
        try(system(strevincecall, intern=TRUE, ignore.stderr=TRUE)) 
      } 
    } 
  } 
  if(format == "dot")
    cat(str) 
  }) 
 
setMethod("cnPlot", "catNetwork", 
        function(object, file = NULL) { 
              if(is.null(file) || file == "") 
                return(cnDot(object, "unknown", "pdf")) 
              else 
                return(cnDot(object, file, "pdf")) 
          }) 

