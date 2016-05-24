write.mitmlSPSS <- function(x, filename, sep="\t", dec=".", na.value=-999, syntax=TRUE, locale=NULL){
# write text file to be read into SPSS

  if(!"mitml"%in%class(x) & !"mitml.list"%in%class(x)) stop("'x' must be of class 'mitml' or 'mitml.list'.")
  if(!dec%in%c(",",".")) stop("Only a dot '.' or a comma ',' may be specified as decimal separator.")

  if("mitml"%in%class(x)){
    il <- mitmlComplete(x,0:x$iter$m)
  }else{ 
    il <- x
  }

  for(ii in 1:length(il)){
    il[[ii]] <- cbind(ii-1,il[[ii]])
    colnames(il[[ii]])[1] <- "Imputation_"
  }
  out <- do.call(rbind,il)
  num <- sapply(out,is.numeric)
  chr <- sapply(out,is.character)
  fac <- sapply(out,is.factor)
  ord <- sapply(out,is.ordered)
  # convert factors
  conv <- as.list(which(fac))
  for(ff in which(fac)){
    out[,ff] <- as.factor(out[,ff])
    conv[[colnames(out)[ff]]] <- matrix(c(levels(out[,ff]),1:nlevels(out[,ff])),ncol=2)
    out[,ff] <- as.numeric(out[,ff])
  }

  ds <- paste(filename, ".dat", sep="")
  out[is.na(out)] <- na.value
  write.table(out, file=ds, sep=sep, dec=dec, col.names=T, row.names=F, quote=F)

  # gerate syntax
  if(syntax){

    sf <- paste(filename, ".sps", sep="")
    if(dec==".") d <- "DOT" else d <- "COMMA"
    cat(file=sf,"SET DECIMAL", d, ".\n")
    if(!is.null(locale)) cat(file=sf, "SET LOCALE",locale,".\n", append=T)
    cat(file=sf, "\n", append=T)
  
    cat(file=sf, append=T,
      "GET DATA\n",
      "/TYPE=TXT\n",
      paste("/FILE=\"",ds,"\"\n",sep=""),
      "/DELCASE=LINE\n",
      paste("/DELIMITERS=\"",sub("\t","\\\\t",sep),"\"\n",sep=""),
      "/ARRANGEMENT=DELIMITED\n",
      "/FIRSTCASE=2\n",
      "/IMPORTCASE=ALL\n",
      "/VARIABLES="
    )

    # class specific format
    width <- sapply(as.matrix(out)[1,], nchar, type="width")
    width[chr] <- sapply(out[,chr,drop=FALSE], function(z) max(nchar(z,type="width")))
    fmt <- data.frame(v=colnames(out),f=character(ncol(out)),stringsAsFactors=F)
    fmt[num|fac|ord,"f"] <- paste("F",width[num|fac|ord]+3,".2",sep="")
    fmt[chr,"f"] <- paste("A",width[chr],sep="")
    fmt[num,"l"] <- "SCALE"
    fmt[fac|chr,"l"] <- "NOMINAL"
    fmt[ord,"l"] <- "ORDINAL"
    fmt[1,"l"] <- "NOMINAL"

    cat(file=sf, "\n ", append=T)
    cat(file=sf, paste(fmt$v,fmt$f, collapse="\n "), ".\n\n", append=T)

    cat(file=sf, append=T, sep="",
      "CACHE .\n",
      "EXECUTE .\n",
      "DATASET NAME panImpute1 WINDOW=FRONT .\n\n"
    )

    # value labels
    cat(file=sf, "VALUE LABELS", append=T)
    for(cc in 1:length(conv)){
      cat(file=sf, "\n", paste("/",names(conv)[cc],sep=""), append=T)
      for(rr in 1:nrow(conv[[cc]])){
        cat(file=sf, "\n", conv[[cc]][rr,2], paste("\'",conv[[cc]][rr,1],"\'",sep=""), append=T)
      }
    }
    cat(file=sf, " .\n\n", append=T)

    # missing values
    cat(file=sf, append=T,
      "MISSING VALUES\n",
      paste(fmt$v[num|fac|ord], collapse=" "), paste("(",na.value,")",sep=""),"\n",
      paste(fmt$v[chr], collapse=" "), paste("(\"",na.value,"\")",sep=""),
      ".\n"
    )

  }

  invisible() 

}

