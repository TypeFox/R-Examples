
## require(gtools)


##' Make a template that feeds into JASPAR databases
##'
##' NA
##' @title Make a template that feeds into JASPAR databases
##' @param x matrix, the pfm
##' @param PARAM a list, the PARAM(s)
##' @param TAG a list, the TAG(s)
##' @param sep a string, the delimiter
##' @param outFpre a string, a file path to save
##' @return A string of the template, and save it in output
##' format of `.template' and `.matrix' if `outFpre' specified.
##' @examples
##' 
##' x <-
##'   rbind(
##'         c(3, 0, 0, 0, 0, 0),
##'         c(8, 0, 23, 0, 0, 0),
##'         c(2, 23, 0, 23, 0, 24),
##'         c(11, 1, 1, 1, 24, 0)
##'         )
##' 
##' PARAM <-
##'   list(
##'        INT_ID=NULL,
##'        BASE_ID="MA0006",
##'        COLLECTION="CORE",
##'        VERSION=1,
##'        NAME="Arnt-Ahr",
##'        SPECIES="10090")
##' TAG <-
##'   list(
##'        class="bHLH",
##'        medline="7592839",
##'        tax_group="vertebrate",
##'        sysgroup="vertebrate",
##'        acc="P30561",
##'        acc="P53762",
##'        comment="dimer",
##'        type="SELEX",
##'        newest=1
##'        )
##' cat(make_template(x=x,PARAM=PARAM,TAG=TAG))
##' 
##' @author Xiaobei Zhao
make_template <-
  function(x,
           PARAM=NA,
           TAG=NA,
           sep="\t",
           outFpre=NULL
           )
{
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  if (!is.matrix(x)){
    stop("x must either be a matrix or a data.frame.")
  }
  if (nrow(x)!=4){
    stop("x must have 4 rows.")
  }
  
  rnames <- c("A","C","G","T")
  if (is.null(rownames(x))){
    rownames(x) <- rnames
    message(sprintf("`x' is a matrix without specific rownames, treating the rows as 
in order of %s.",paste(rnames,sep="",collapse=", ")))
  } else {
    rownames(x) <- toupper(rownames(x))
    if (all(rownames(x) %in% rnames)){
      x <- x[rnames,]
    } else {
      stop(sprintf("`x' must only specified by rownames in 
%s.",paste(rnames,sep="",collapse=", ")))
    }
  }

  if (!identical(NA,PARAM) | !identical(NA,TAG)){
    ## set default PARAM
    pnames <-
      ##c("INT_ID","BASE_ID","COLLECTION","VERSION","NAME","SPECIES")
      c("BASE_ID","COLLECTION","VERSION","NAME","SPECIES")
    pnames_na <- pnames[!pnames %in% names(PARAM)]

    if (length(names(PARAM))==0)
      PARAM <- list()
    
    ## if (! "INT_ID" %in% pnames_na) PARAM$ <- NULL
    if ("BASE_ID" %in% pnames_na) PARAM$BASE_ID <- 1
    if ("COLLECTION" %in% pnames_na) PARAM$COLLECTION <- "CORE"
    if ("VERSION" %in% pnames_na) PARAM$VERSION <- 1
    if ("NAME" %in% pnames_na) PARAM$NAME <- NA
    if ("SPECIES" %in% pnames_na) PARAM$SPECIES <- NA

    
    ## print(pnames); print(names(PARAM));
    for (a in pnames){
      if ( invalid(PARAM[[a]]) ) PARAM[[a]] <- "-"
    }

    ## set default TAG
    tnames <-
      c("class", "family", "acc", "pazar_tf_id", "medline", "tax_group", "type", "newest", "comment")

    if (length(names(TAG))!=0)
      tnames <- unique(c(tnames[tnames %in% names(TAG)],names(TAG)))
    else
      TAG <- list()

    tnames_na <- tnames[!tnames %in% names(TAG)]
    if ("class" %in% tnames_na) TAG$class <- NA
    if ("family" %in% tnames_na) TAG$family <- NA
    if ("acc" %in% tnames_na) TAG$acc <- NA
    if ("pazar_tf_id" %in% tnames_na) TAG$pazar_tf_id <- NA
    if ("medline" %in% tnames_na) TAG$medline <- NA
    if ("tax_group" %in% tnames_na) TAG$tax_group <- NA
    if ("type" %in% tnames_na) TAG$type <- NA
    if ("newest" %in% tnames_na) TAG$newest <- 1
    if ("comment" %in% tnames_na) TAG$comment <- NA
    ##if ("sysgroup" %in% tnames_na) TAG$sysgroup <- NA


    ##print(tnames); print(names(TAG));

    for (a in tnames){
      if ( invalid(TAG[[a]]) ) TAG[[a]] <- "-"
    }
  
    if (! all(pnames %in% names(PARAM))){
      stop(sprintf('PARAM must contain: %s',paste(pnames,sep="",collapse=",")))
    }
    for (i in 1:length(PARAM)) {
      if(identical(NA,PARAM[[i]])){
        PARAM[[i]] <- "#TBD"
      }
    }
    
    for (j in 1:length(TAG)){
      if(identical(NA,TAG[[i]])){
        TAG[[j]] <- "-"
      }
    }

  
    ##output
    names(TAG) <- paste("TAG",names(TAG),sep=sep)
    tnames <- paste("TAG",tnames,sep=sep)
    ##print(PARAM)

    ret0 <-
      paste(
            list2table.string(sortList(PARAM,pnames),sep=sep),
            list2table.string(sortList(TAG,tnames),sep=sep),
            sep="")
  } else {
    ret0 <- ""
  }

  ##
  ret0.matrix <- matrix2table.string(x)
  ret1 <-
    paste(ret0, "\n", ret0.matrix, "\n//\n", sep="")

  
  ret0.header <- sprintf("> %s.%s %s", PARAM$BASE_ID, PARAM$VERSION, PARAM$NAME)
  ret2 <-
    paste(ret0.header, "\n", ret0.matrix, "\n", sep="")

  
  if (!is.null(outFpre)){
    writeString(sprintf("%s.template",outFpre),ret1)
    writeString(sprintf("%s.matrix",outFpre),ret2)
  }
  
  invisible(ret1)
}



## template2matrix <- function(){}
## template2pfm <- function(){}
