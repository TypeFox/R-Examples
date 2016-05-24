#' Graphical User Interface for Crossover Designs
#' 
#' Starts a graphical user interface for accessing and creating crossover designs.
#' 
#' See the vignette of this package for further details, since describing a GUI
#' interface is better done with some nice pictures.
#' 
#' @return The function itself returns nothing of interest. But from the GUI designs
#' and objects can be created or edited that will be available in R under the
#' specified variable name after saving.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @keywords misc
#' @examples
#' 
#' 
#' \dontrun{
#' CrossoverGUI()
#' }
#' 
#' 
#' @export CrossoverGUI
CrossoverGUI <- function() {
	invisible(.jnew("org/mutoss/gui/CrossoverGUI"))	
}

getTable <- function(d, type="HTML", forceInteger=TRUE, digits=4, names=TRUE) {
  if (forceInteger) {
    d <- design2integer(d)
    rownames(d) <- paste("p", 1:dim(d)[1], sep="")
    colnames(d) <- paste("s", 1:dim(d)[2], sep="")
  } else {
    d <- round(d, digits)
    rownames(d) <- paste("t", 1:dim(d)[1], sep="")
    colnames(d) <- paste("t", 1:dim(d)[2], sep="")
  }
  if (type=="ASCII") {
    if (names) {
      return(paste("<pre>",paste(capture.output(print(d)), collapse="\n"),"</pre>"))
    } else {
      return(paste("<pre>",paste(capture.output(write.table(format(d, justify="right"),
                  row.names=F, col.names=F, quote=F)), collapse="\n"),"</pre>"))
    }
  } else if (type=="HTML") {
    return(paste(capture.output(print(xtable(d, digits=digits), comment=FALSE, type="html", include.rownames=names, include.colnames=names)),collapse="\n"))
  } else if (type=="R") {
    return(paste("<pre>",dputMatrix(d),"</pre>"))
  }
}

getDesignText <- function(d, model=1, type="HTML", carryover=TRUE, digits=4, var=TRUE, eff=TRUE, names=TRUE, model.param=list()) {  
  result <- ""
  if (var) {
    m <- general.carryover(d, model=model)$Var.trt.pair
    result <- paste(result, "<b>Var.trt.pair:</b><br>", getTable(m, type, forceInteger=FALSE, digits=digits, names=names))  
  }
  if (eff) {
    warn <- ""
    #if (model!=1) warn <- "(Warning: efficiency is calculated for model 1)"
    m <- design.efficiency(d, model=model, model.param=model.param)$eff.trt.pair.adj
    result <- paste(result, "<b>Eff.trt.pair",warn,":</b><br>", getTable(m, type, forceInteger=FALSE, digits=digits, names=names),sep="")  
  }
  if (carryover) {
    gco <- general.carryover(d, model=model)
    i <- 2
    while (i<length(gco)) {
      if (is.matrix(gco[[i]])) {
        result <- paste(result, "<b>",names(gco)[i],":</b><br>", getTable(gco[[i]], type, forceInteger=FALSE, digits=digits, names=names))  
      }
      i <- i + 1
    }
  }
  return(result)
}

dputMatrix <- function(m, name, indent=6, rowNames=FALSE) {
  s <- "rbind("
  if (!missing(name)) s <- paste(name,"<- rbind(") 
  for (i in 1:(dim(m)[1])) {
    nameLater <- FALSE
    if (any(make.names(row.names(m))!=row.names(m))) {
      rowNames <- FALSE
      nameLater <- TRUE
    }
    rName <- ifelse(rowNames, paste(row.names(m)[i],"=",sep=""), "")
    s <- paste(s, 
               ifelse(i==1,"",paste(rep(" ",indent),collapse="")),
               rName,
               dputS(unname(m[i,])),
               ifelse(i==dim(m)[1],")\n",",\n"),
               sep="")            
  }
  if (nameLater) {
    if (missing(name)) {
      warning("Can set row names if no name for matrix is given.")
      return(s)
    }
    s <- paste(s, 
               "row.names(",name,") <- ", dputS(row.names(m)), "\n", sep="")
  }
  return(s)
}

dputS <- function(x) {
  paste(capture.output(dput(x)), collapse="\n")
}

isRName <- function(x) {
  return(x==make.names(x))
}