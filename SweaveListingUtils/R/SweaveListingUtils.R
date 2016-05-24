###############################################################################
# some R helpers for Sweave / listings
# will be included into some package some day...
#
###############################################################################

# .onLoad<-function(lib,pkg){require(methods)}

.onAttach <- function(library, pkg)
{
#  if (is.null(library)) 
#            library <- .libPaths()

     unlockBinding(".keywordsR", asNamespace("SweaveListingUtils"))
     unlockBinding(".alreadyDefinedPkgs", asNamespace("SweaveListingUtils"))
     unlockBinding(".numberofRstyleDefs", asNamespace("SweaveListingUtils"))
     unlockBinding(".tobeDefinedPkgs", asNamespace("SweaveListingUtils"))
     unlockBinding(".CacheFiles", asNamespace("SweaveListingUtils"))
     unlockBinding(".CacheLength", asNamespace("SweaveListingUtils"))
     unlockBinding(".SweaveListingOptions", asNamespace("SweaveListingUtils"))
     msga <- gettext(
    "NOTE: Support for this package will stop soon.\nPackage 'knitr' is providing the same functionality in a better way.\n"
                   )
     msgb <- gettext(
    "Some functions from package 'base' are intentionally masked ---see SweaveListingMASK().\n"
                   )
    msgc <- gettext(
    "Note that global options are controlled by SweaveListingoptions() ---c.f. ?\"SweaveListingoptions\"."
                   )

     buildStartupMessage(pkg = "SweaveListingUtils", msga, msgb, msgc,
                         library = library, packageHelp = TRUE,
                    VIGNETTE = gettext(
"There is a vignette to this package; try vignette(\"ExampleSweaveListingUtils\")."
                                      )
                         )

  invisible()
} 

SweaveListingMASK <- function(library = NULL) 
{
    infoShow(pkg = "SweaveListingUtils", filename="MASKING", library = library)
}



#------------------------------------------------------------------------------
# SweaveListingPreparations
#------------------------------------------------------------------------------

SweaveListingPreparations <- function(
      withOwnFileSection = FALSE,
      withVerbatim = FALSE,
      withSchunkDef = TRUE,
      gin = TRUE,
      ae = TRUE,
      LineLength = getOption("width"),
      Rset = getSweaveListingOption("Rset"),
      Rdset = getSweaveListingOption("Rdset"),
      Rin = getSweaveListingOption("Rin"),
      Rout = getSweaveListingOption("Rout"),
      Rcode = getSweaveListingOption("Rcode"),
      Rcolor = getSweaveListingOption("Rcolor"),
      RRecomdcolor = getSweaveListingOption("RRecomdcolor"),
      Rbcolor = getSweaveListingOption("Rbcolor"),
      Routcolor = getSweaveListingOption("Routcolor"),
      Rcommentcolor = getSweaveListingOption("Rcommentcolor"),
      pkg = getSweaveListingOption("pkg"),
      pkv = getSweaveListingOption("pkv"),
      fileCommand = getSweaveListingOption("fileCommand"),
      pkgCommand = getSweaveListingOption("pkgCommand"),
      lib.loc = NULL){

### begin of code to SweaveListingPreparations

   sws <- .SweaveListingOptions
   sws$inSweave <- TRUE
   if(getRversion()>"2.15.1")
      assignInMyNamespace(".SweaveListingOptions", sws) else{
   assignInNamespace(".SweaveListingOptions", sws, "SweaveListingUtils")
   }
   
   withVerbatim <- rep(withVerbatim, length.out=3)
   if(is.null(names(withVerbatim)))
      names(withVerbatim) <- c("Sinput", "Soutput", "Scode")


   line <- paste("%",paste(rep("-",LineLength-2),collapse=""),"%\n", sep="")



   cat(line,"%Preparations for Sweave and Listings\n",line,"%\n", sep = "")

   cat("\\RequirePackage{color}\n")
   cat("\\definecolor{Rcolor}{rgb}{",
       paste(Rcolor,collapse=", "),"}\n", sep = "")
   cat("\\definecolor{RRecomdcolor}{rgb}{",
       paste(RRecomdcolor,collapse=", "),"}\n", sep = "")
   cat("\\definecolor{Rbcolor}{rgb}{",
       paste(Rbcolor,collapse=", "),"}\n", sep = "")
   cat("\\definecolor{Routcolor}{rgb}{",
       paste(Routcolor,collapse=", "),"}\n", sep = "")
   cat("\\definecolor{Rcommentcolor}{rgb}{",
       paste(Rcommentcolor,collapse=", "),"}\n", sep = "")
   cat(line)
   writeLines(readLines(file.path(system.file(package = "SweaveListingUtils",
                                         lib.loc = lib.loc),
                             "TeX", "Rdlisting.sty")))
   cat(line)
   ### (first) definition of Rstyle
   lstsetR(Rset = Rset, LineLength = LineLength,
           startS ="\\lstdefinestyle{RstyleO1}{")
   cat("\\lstdefinestyle{Rstyle}{style=RstyleO1}\n")
   ###  definition of Rdstyle
   lstsetRd(Rdset=Rdset, LineLength=LineLength,
           startS ="\\lstdefinestyle{Rdstyle}{")
   cat(line)
   if(!withOwnFileSection)
       SweaveListingoptions("addRset" = FALSE, "addRdset" = FALSE)
   cat("\\global\\def\\Rlstset{\\lstset{style=Rstyle}}%\n")
   cat("\\global\\def\\Rdlstset{\\lstset{style=Rdstyle}}%\n")
   cat(line)
   cat("\\global\\def\\Rinlstset{\\lstset{style=Rinstyle}}%\n")
   cat("\\global\\def\\Routlstset{\\lstset{style=Routstyle}}%\n")
   cat("\\global\\def\\Rcodelstset{\\lstset{style=Rcodestyle}}%\n")
   cat(line)
   if(!withOwnFileSection)
      cat("\\Rlstset\n")
   cat(line,"%copying relevant parts of Sweave.sty\n",line,"%\n", sep = "")

   cat("\\RequirePackage{graphicx,fancyvrb}%\n")
   cat("\\IfFileExists{upquote.sty}{\\RequirePackage{upquote}}{}%\n\n")
   cat("\\RequirePackage{ifthen}%\n")
   ### you might still want to have the boolean TeX
   #   variables available in your code
   if(gin){
     cat("\\newboolean{Sweave@gin}%\n")
     cat("\\setboolean{Sweave@gin}{true}%\n")
     cat("\\setkeys{Gin}{width=0.8\\textwidth}%\n")
   }

   if(ae){
      cat("\\newboolean{Sweave@ae}\n")
      cat("\\setboolean{Sweave@ae}{true}%\n")
      cat("\\RequirePackage[T1]{fontenc}\n",
          "\\RequirePackage{ae}\n%\n", sep ="")
   }

   if(withSchunkDef)
      cat("\\newenvironment{Schunk}{}{}\n\n")

   cat("\\newcommand{\\Sconcordance}[1]{% \n",
     "\\ifx\\pdfoutput\\undefined% \n",
     "\\csname newcount\\endcsname\\pdfoutput\\fi% \n",
     "\\ifcase\\pdfoutput\\special{#1}% \n",
     "\\else\\immediate\\pdfobj{#1}\\fi} \n\n", sep ="")

   #% ---- end of parts of Sweave.sty
   cat(line,"% ---- end of parts of Sweave.sty\n",line,"%\n", sep = "")

   #% definition of Rstyle, Rdstyle, Rinstyle, Routstyle, Rcodestyle
   if(!withOwnFileSection){
      if(withVerbatim["Sinput"]){
         cat("% ---- input \n")
         cat("\\DefineVerbatimEnvironment{Sinput}{Verbatim}")
         cat("%\n  {formatcom=\\color{Rcolor}",
             "\\lstset{fancyvrb=true,escapechar='}}\n", sep = "")
         cat("%\n")
      }else{#### Thanks to Andrew Ellis !!
         cat("% ---- input \n")
         lstset(taglist(list=Rin), LineLength=LineLength,
                startS = "\\lstdefinestyle{RinstyleO}{")
         cat("\\lstdefinestyle{Rinstyle}{style=RinstyleO}\n")
         cat("\\lstnewenvironment{Sinput}{\\Rinlstset}{\\Rlstset}\n")
         cat("%\n")
      }
      if(withVerbatim["Soutput"]){
         cat("% ---- output \n")
         cat("\\DefineVerbatimEnvironment{Soutput}{Verbatim}")
         cat("%\n  {formatcom=\\color{Routcolor}\\small",
             "\\lstset{fancyvrb=false}}\n", sep = "")
         cat("%\n")
      }else{#### Thanks to Andrew Ellis !!
         cat("% ---- output \n")
         lstset(taglist(list=Rout), LineLength=LineLength,
                startS = "\\lstdefinestyle{RoutstyleO}{")
         cat("\\lstdefinestyle{Routstyle}{style=RoutstyleO}\n")
         cat("\\lstnewenvironment{Soutput}{\\Routlstset}{\\Rlstset}\n")
         cat("%\n")
      }
      if(withVerbatim["Scode"]){
         cat("% ---- code \n")
         cat("\\DefineVerbatimEnvironment{Scode}{Verbatim}")
         cat("%\n  {fontshape=sl,formatcom=\\color{Rcolor}",
             "\\lstset{fancyvrb=true}}\n", sep = "")
         cat("%\n")
      }else{#### Thanks to Andrew Ellis !!
         cat("% ---- code \n")
         lstset(taglist(list=Rcode), LineLength=LineLength,
                startS = "\\lstdefinestyle{RcodestyleO}{")
         cat("\\lstdefinestyle{Rcodestyle}{style=RcodestyleO}\n")
         cat("\\lstnewenvironment{Scode}{\\Rcodelstset}{\\Rlstset}\n")
         cat("%\n")
      }
   }
   cat(line)
   cat("\\let\\code\\lstinline\n")
   cat("\\def\\Code#1{{\\tt\\color{Rcolor} #1}}\n")
   cat(fileCommand,"\n")
   cat(pkgCommand,"\n")
   if(missing(pkv)){
        if(nzchar(rpkv <- readPkgVersion(package = pkg, lib.loc=lib.loc)))
           pkv <- rpkv
      }
   cat("\\newcommand{\\pkgversion}{{\\tt ",pkv,"}}\n", sep = "")
   cat(line)
   lstsetLanguage()
   cat(line,"%\n%\n",sep="")
   cat("\n")
   return(invisible())
}

#------------------------------------------------------------------------------
# readSourceFromRForge
#------------------------------------------------------------------------------
readSourceFromRForge <- function(PKG, TYPE, FILENAME, PROJECT,
                                 fromRForge = getSweaveListingOption("fromRForge"),
                                 base.url = getSweaveListingOption("base.url")){
  base.URL <- if(fromRForge) paste(base.url, PKG, "/", TYPE,"/", FILENAME,
                                   "?root=", PROJECT, sep ="") else base.url
  if(is.null(.CacheFiles[[base.URL]])){
    .CacheLength <<- .CacheLength + 1
    url.connection <- url(base.URL) 
    RL <- try(readLines(url.connection),silent=TRUE)
    if(is(RL,"try-error")) return(character(0))
    close(url.connection)
    .CacheFiles[[base.URL]] <<- RL
  }
  .CacheFiles[[base.URL]]
}

#------------------------------------------------------------------------------
# copySourceFromRForge
#------------------------------------------------------------------------------
copySourceFromRForge <- function(PKG, TYPE, FILENAME, PROJECT, from, to,
                                 offset.before = 0, offset.after = 0,
                                 fromRForge = getSweaveListingOption("fromRForge"),
                                 base.url = getSweaveListingOption("base.url"), ... ){
   RL <- readSourceFromRForge(PKG, TYPE, FILENAME, PROJECT, 
                              fromRForge = fromRForge, base.url = base.url)
   lR <- length(RL)
   if(lR>0){
      from <- if(missing(from)) 1 else {if(is.numeric(from))
                                            max(from-offset.before,1)
                                        else {if(length(gr0 <- grep(from,RL, ...)))
                                            max(gr0[1]-offset.before,1) else lR
                                          }
                                    }
      to <- if(missing(to)) lR else {if(is.numeric(to))
                                        min(to+offset.after,lR)
                                     else {if(length(gr1<-grep(to,RL[from:lR], ...)))
                                              min(from+gr1[1]-1+offset.after,lR)
                                           else 0
                                           }
                                    }
      if(to>=from) return(list(text=RL[from:to], lines=c(from,to)))
   }
   return(invisible(NULL))
}

#------------------------------------------------------------------------------
# lstinputSourceFromRForge
#------------------------------------------------------------------------------
lstinputSourceFromRForge <- function(PKG, TYPE, FILENAME, PROJECT, from, to,
                                 offset.before = 0, offset.after = 0,
                                 LineLength = getOption("width"),
                                 withLines = ifelse(TYPE=="R", TRUE, FALSE),
                                 fromRForge = getSweaveListingOption("fromRForge"),
                                 base.url = getSweaveListingOption("base.url"), 
                                 ...){
   line <- paste("%",paste(rep("-",LineLength-2),collapse=""),"%\n", sep="")
   dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
   argL0 <- list(PKG, TYPE, FILENAME, PROJECT)
   totl <- 1
   if(!missing(from)) totl <- max(totl, length(from))
   if(!missing(to)) totl <- max(totl, length(to))
   if(totl>1){
      argL <- vector("list",totl)
      fromL <- NULL; toL <- NULL
      offs.from <- rep(offset.before, length.out = totl)
      offs.to <- rep(offset.after, length.out = totl)
      fromRForge <- rep(fromRForge, length.out = totl)
      base.url <- rep(base.url, length.out = totl)
      if(!missing(from))
          fromL <- rep(as.list(from),length.out = totl)
      if(!missing(to))
          toL <- rep(as.list(to),length.out = totl)
      for (j in 1 : totl){
          argL[[j]] <- argL0
          if(!missing(from))
              argL[[j]] <- c(argL[[j]], from = fromL[[j]])
          if(!missing(to))
              argL[[j]] <- c(argL[[j]],to = toL[[j]])
          argL[[j]] <- c(argL[[j]], offset.before = offs.from[j],
                         offset.after = offs.to[j], 
                         fromRForge = fromRForge[j], base.url = base.url[j],dots)
      }
   }else {
      argL <- argL0
      if(!missing(from)) argL <- c(argL, from = from)
      if(!missing(to)) argL <- c(argL, to = to)
      argL <- list(c(argL, offset.before = offset.before,
                        offset.after = offset.after, 
                        fromRForge = fromRForge, base.url = base.url,dots))
   }
   erg <- lapply(argL, function(x)  do.call(copySourceFromRForge, args = c(x)))
   if(length(erg)==0) return(invisible(NULL))
   RL <- lapply(erg, function(x) x$text)
   lineNr <- lapply(erg, function(x) x$lines)
   lR <- lapply(RL, length)
   lE <- length(erg)
   if(withLines){
      for(k in 1:lE){
        if( !is.null(lineNr[[k]])){
          if( k == 1 ) {
               if( ( lineNr[[k]][1] < lineNr[[k]][2] ) || ( lE>1 ) )
                    cat("lines ")
               else cat("line ")
          }else{
               if( k < lE )
                    cat(", \n")
               else cat(", and\n")
               }
          if(lineNr[[k]][1] < lineNr[[k]][2])
             cat(lineNr[[k]][1], "--", lineNr[[k]][2], sep = "")
          else cat(lineNr[[k]][1])
        }
      }
      cat("\n")
   }
   for(k in 1:length(erg)){
     if(!is.null(lR[[k]])){ if(lR[[k]]){
        todo <- NULL
        if(TYPE=="man"){
          ex.from <- if(length(gr <- grep("\\\\examples\\{",RL[[k]])))
                     gr[1] else lR[[k]]
          ex.to <- if(length(gr <- grep("\\}",RL[[k]][ex.from:lR[[k]]])))
                   ex.from+gr[1]-1 else 1
          cat(line)
          cat("\\begin{lstlisting}[style = Rdstyle]\n")
          if(ex.from<=ex.to){
             writeLines(RL[[k]][1:(ex.from)])
             cat("\\end{lstlisting}\\vspace{-2ex}\n")
             cat(line)
             cat("\\begin{lstlisting}[style = Rstyle,")
             cat("basicstyle = \\scriptsize\\color{Rcolor},",
                 " xleftmargin = 2em]\n", sep = "")
             writeLines(RL[[k]][(ex.from+1):(ex.to-1)])
             if(ex.to <lR){
                cat("\\end{lstlisting}\\vspace{-3ex}\n")
                cat(line)
                cat("\\begin{lstlisting}[style = Rdstyle]\n")
                writeLines(RL[[k]][(ex.to):lR[[k]]])
                }
          }else writeLines(RL[[k]])
        }else{
          cat(line,"\\Rlstset\n",sep="")
          cat("\\begin{lstlisting}\n")
          writeLines(RL[[k]])
        }
        cat("\\end{lstlisting}\n",line,"%\n\n",sep="")
        }}
      }
   return(invisible())
}

readPkgVersion <- function(package, lib.loc = NULL){
       Dfile <- system.file("DESCRIPTION", package=package, lib.loc=lib.loc)
       return( if(nzchar(Dfile)) read.dcf(Dfile, fields="Version") else "")}
       
