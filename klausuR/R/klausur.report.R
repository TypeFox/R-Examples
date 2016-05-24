# Copyright 2009-2015 Meik Michalke <meik.michalke@hhu.de>
#
# This file is part of the R package klausuR.
#
# klausuR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# klausuR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with klausuR.  If not, see <http://www.gnu.org/licenses/>.


#' Generate individual reports on multipe choice test results
#'
#' \code{klausur.report} takes (at least) an object of class klausuR (or klausuR.mult) and a matriculation number to generate personal test results
#' in LaTeX and/or PDF format.
#'
#' The report contains, next to the individual results, a table with all given and correct  answers (using \code{\link[xtable]{xtable}}),
#' as well as nice histograms showing the distribution of the test results (points and/or marks are supportet). If the matriculation numer
#' is set to "all", reports for all subjects are produced. Setting it to "anon" will get you a printable version of the anonymized results.
#' 
#' By default output is sent to standard out. To save them to disk in LaTeX format a "save" parameter is provided. Alternatively, the reports
#' can be converted to PDF format as well. \code{klausur.report} is calling \code{\link[tools]{texi2dvi}} from the \code{tools} package for that.
#'
#' If the object is of class klausuR.mult, only the global results for tests with several test forms are evaluated. In case you'd rather like
#' reports on each test form, call \code{klausur.report} with the single slots from that object accordingly.
#'
#' @param klsr An object of class klausuR or klausuR.mult. To create reports from more than one object with the same configuration, you can
#'    also give them in one list here, which will cause the function to call itself recursively.
#' @param matn Matriculation number, "all" (produces individuall documents for all subjects), "anon" (produces anonymous feedback)
#'  or "glob" (produces a global results document).
#' @param save Logical: If TRUE, files are saved to disk (scheme: "\code{path}/\code{matn}.tex").
#' @param pdf Logical: If TRUE, LaTeX reports will be converted to PDF automatically, using \code{\link[tools]{texi2dvi}}.
#'  If \code{save} is FALSE, a temporary directory is used, that is only the PDF files will be saved.
#' @param path Path for \code{save} and \code{hist} files.
#' @param file.name File name scheme for the reports, either "matn" (matriculation number) or "name" (name and firstname).
#' @param hist A list with the logical elements \code{points} and \code{marks}: If TRUE, the reports will include histograms
#'  of the distribution of points and/or marks. The needed PDF files will be created by \code{\link[klausuR:plot]{plot}} and saved as well.
#'  (see \code{path}, \code{hist.points} and \code{hist.marks}).
#' @param hist.merge If you need/want to combine results from several \code{klausuR} class objects for the histograms, provide them all in a list here.
#' @param hist.points File name for the histogram of points.
#' @param hist.marks File name for the histogram of marks.
#' @param descr Details on the test: List with the elements \code{title} (title of the test), \code{name} (your name) and \code{date}.
#' @param marks.info A list with the logical elements \code{points} and \code{percent}: If TRUE, the reports will include a table showing
#'  how marks were assigned to points achieved and/or percent solved, respectively.
#' @param lang Set to "de" for reports in German, English is the default.
#' @param alt.candy If TRUE, a comma will be inserted for items with multiple alternatives ("235" becomes "2, 3, 5" in the printout)
#' @param anon.glob.file If \code{matn="anon"} or \code{matn="glob"}, you can specify a filename for this particular report.
#' @param decreasing Logical, whether sorting of output should be done increasing or decreasing (only relevant for \code{matn="anon"} or
#'  \code{matn="glob"}).
#' @param sort.by Character string naming a variable to sort the results by. Defaults to \code{"Marks"} (only relevant for \code{matn="anon"} or
#'  \code{matn="glob"}).
#' @param NRET.legend Logical, If ET/NRET data is reported, you can demand a legend in the table caption by setting this to true.
#' @param table.size Character string to shrink the tables, must be one of \code{"auto"}, \code{"normalsize"}, \code{"small"},
#'    \code{"footnotesize"}, \code{"scriptsize"} or \code{"tiny"}. The default \code{table.size="auto"} tries to decide between
#'    \code{"normalsize"} and \code{"footnotesize"} to avoid pages with only one or two rows. If that fails, try to manually set the size.
#' @param merge Logical, if \code{TRUE} no individual PDFs will be saved, but one large file with all reports. Uses the "pdfpages" package,
#'    and only useful if \code{pdf=TRUE} as well.
#' @param quiet Logical, if \code{TRUE} no feedback messages on the current status are given.
#' @aliases klausur.report
#' @keywords IO file
#' @return One or several LaTeX and/or PDF documents. If defined two histograms will be plotted.
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @seealso \code{\link[klausuR:klausur]{klausur}}, \code{\link[xtable]{xtable}}, \code{\link[tools]{texi2dvi}}
#' @import xtable graphics tools
#' @export
#' @examples
#' data(antworten)
#' 
#' # vector with correct answers:
#' richtig <- c(Item01=3, Item02=2, Item03=2, Item04=2, Item05=4,
#'  Item06=3, Item07=4, Item08=1, Item09=2, Item10=2, Item11=4,
#'  Item12=4, Item13=2, Item14=3, Item15=2, Item16=3, Item17=4,
#'  Item18=4, Item19=3, Item20=5, Item21=3, Item22=3, Item23=1,
#'  Item24=3, Item25=1, Item26=3, Item27=5, Item28=3, Item29=4,
#'  Item30=4, Item31=13, Item32=234)
#' 
#' # vector with assignement of marks:
#' notenschluessel <- c()
#' # scheme of assignments: marks[points_from:to] <- mark
#' notenschluessel[0:12]  <- 5.0
#' notenschluessel[13:15] <- 4.0
#' notenschluessel[16:18] <- 3.7
#' notenschluessel[19:20] <- 3.3
#' notenschluessel[21]    <- 3.0
#' notenschluessel[22]    <- 2.7
#' notenschluessel[23]    <- 2.3
#' notenschluessel[24]    <- 2.0
#' notenschluessel[25:26] <- 1.7
#' notenschluessel[27:29] <- 1.3
#' notenschluessel[30:32] <- 1.0
#' 
#' data.obj <- klausur.data(answ=antworten, corr=richtig, marks=notenschluessel)
#' klsr.obj <- klausur(data.obj)
#'
#' \dontrun{
#' klausur.report(klsr=klsr.obj, matn="all", descr=list(title="Klausur Tatort",
#'   name="Dr. T. Aeter", date="24.09.2010"))
#' }

klausur.report <- function(klsr, matn, save=FALSE, pdf=FALSE, path=NULL, file.name="matn",
          hist=list(points=FALSE, marks=FALSE), hist.merge=list(), hist.points="hist_points.pdf", hist.marks="hist_marks.pdf",
          descr=list(title=NULL, name=NULL, date=NULL), marks.info=list(points=FALSE, percent=FALSE),
          lang="en", alt.candy=TRUE, anon.glob.file="anon.tex", decreasing=TRUE, sort.by="Points", NRET.legend=FALSE, table.size="auto", merge=FALSE, quiet=FALSE){
  # to avoid NOTEs from R CMD check:
  marks.information <- NULL

  # before we start let's look at klsr
  # if klsr is a list, iterate through it recusively
  if(is.list(klsr)){
    for(this.klsr in klsr){
      klausur.report(klsr=this.klsr, matn=matn, save=save, pdf=pdf, path=path, file.name=file.name,
        hist=hist, hist.points=hist.points, hist.marks=hist.marks,
        descr=descr, marks.info=marks.info, lang=lang, alt.candy=alt.candy, anon.glob.file=anon.glob.file,
        decreasing=decreasing, sort.by=sort.by, NRET.legend=NRET.legend, table.size=table.size)
    }
    return("done")
  } else {}
  # check if output needs to be sorted
  if(identical(matn, "anon") | identical(matn, "glob")){
    klsr <- sort(klsr, decreasing=decreasing, sort.by=sort.by)
  } else {}
  # if it's of class "klausuR.mult", extract global results and drop the rest
  if(inherits(klsr, "klausuR.mult")){
    klsr <- klsr@results.glob
  } else{
    # check whether klsr is an object of class "klausuR" instead
    if(!inherits(klsr, "klausuR")){
      stop(simpleError("The given object is not of class \"klausuR\"!"))
    } else {}
  }

  if(!is.numeric(matn) && !identical(matn, "all") && !identical(matn, "anon") && !identical(matn, "glob")){
    stop(simpleError("Value assigned to matn must be numeric, \"all\", \"anon\" or \"glob\"!"))
  } else {}

  if((isTRUE(save) || isTRUE(pdf) || isTRUE(hist$points) || isTRUE(hist$marks)) && is.null(path)){
    stop(simpleError("Files have to be saved, but path is empty!"))
  } else {}

  ## path handling
  # check if path exists
  if(!isTRUE(file.info(path)$isdir)) {
    stop(simpleError(paste(path,"is not a valid path!")))
  } else {}
  # PDF creation will be done in an temporal directory if "save" is FALSE
  # we'll make "path" "path.orig" and override it with that tempdir internally
  if(!isTRUE(save) && isTRUE(pdf)){
    path.orig <- path
    path <- tempfile("klausuR")
    if(!dir.create(path, recursive=TRUE)){
      stop(simpleError("Couldn't create temporary directory! Try with save=TRUE"))
    } else {}
    # if the function is done, remove the tempdir
    on.exit(
      if(!identical(path, path.orig)){
        unlink(path, recursive=TRUE)
      } else {}
    )
  } else {}

  # check value for table.size
  if(!table.size %in% c("auto", "normalsize", "small", "footnotesize", "scriptsize", "tiny")){
    warning(paste("Invalid value for 'tabe.size':\n  ", table.size,"\n  Reverted to \"normal\"."), call.=FALSE)
    table.size <- "auto"
  } else {}

  # define the text of the LaTeX document...
  if(identical(lang, "de")){
    text <- list(Auswertung="Einzelauswertung",
        DozentIn="Dozent",
        MatrikelNr="Matrikel-Nr.",
        Datum="Datum der Klausur",
        Antwort="Antwort",
        Korrekt="L\u00f6sung",
        Punkte="Punkte",
        Ergebnisse="Ergebnisse",
        Erreicht="Erreichte Punktzahl",
        Prozent=" der erreichbaren Punkte",
        Note="Daraus ergibt sich die \\textbf{Note",
        VertErgs="Verteilung der Klausurergebnisse",
        Pseudonym="Pseudonym",
        AProzent="Prozent",
        ANote="Note",
        LfdNr="Nr.",
        GMatNr="MatrNr.",
        Name="Name",
        Vorname="Vorname",
        Notenschluessel="Notenschl\"ussel",
         NRET.expl=". \\emph{Erl\"auterung:} >>+<< -- richtig; >>-<< -- falsch; >>0<< -- keine Angabe; >>*<< -- fehlerhafte Angabe."
    )
    # ... and that of the plots
    hist.text <- list(P.xlab="Punkte",P.ylab="H\u00e4ufigkeit",P.main="Verteilung nach Punkten",
        N.xlab="Note",N.ylab="H\u00e4ufigkeit",N.main="Verteilung nach Noten"
    )
  } ## end of german l10n
  else {
    text <- list(Auswertung="Individual Report",
        DozentIn="Docent",
        MatrikelNr="Matriculation No.",
        Datum="Date of test",
        Antwort="Answer",
        Korrekt="Solution",
        Punkte="Points",
        Ergebnisse="Results",
        Erreicht="Achieved score",
        Prozent=" of achievable points",
        Note="This score gives \\textbf{mark",
        VertErgs="Distribution of test results",
        Pseudonym="Pseudonym",
        AProzent="Percent",
        ANote="Mark",
        LfdNr="No.",
        GMatNr="MatrNo.",
        Name="Name",
        Vorname="First name",
        Notenschluessel="Marks defined",
        NRET.expl=". \\emph{Explaination:} >>+<< -- right; >>-<< -- wrong; >>0<< -- not answered; >>*<< -- errenous answer."
    )
    # ... and that of the plots
    hist.text <- list(P.xlab="Points",P.ylab="Frequency",P.main="Distribution by points",
        N.xlab="Marks",N.ylab="Frequency",N.main="Distribution by marks"
    )
  }

  klsr.to.plot <- klsr
  if(is.list(hist.merge) & length(hist.merge) > 0){
    klsr.to.plot@results <- plot.merger(hist.merge)
  } else {}

  if(hist$points | hist$marks) {
    if(hist$points) {
      pdf(file=file.path(path, hist.points),
      width=10, height=10,
      pointsize=22, bg="white")
      plot(klsr.to.plot, xlab=hist.text$P.xlab, ylab=hist.text$P.ylab, main=hist.text$P.main)
      dev.off()
    } else {}

    if(hist$marks) {
      pdf(file=file.path(path, hist.marks),
      width=10, height=10,
      pointsize=22, bg="white")
      plot(klsr.to.plot, marks=TRUE, xlab=hist.text$N.xlab, ylab=hist.text$N.ylab, main=hist.text$N.main)
      dev.off()
    } else {}
  } else {}

  ## let's grab some info out of the klausuR-object for code readability
  results <- klsr@results
  res.points <- klsr@points
  truefalse <- klsr@trfls
  wght <- klsr@wght
  answers <- klsr@answ
  correct <- klsr@corr

  # for a nice printout, check numer of needed digits for points.
  # e.g, if you can get 1/2 points, you'd need one digit. but we won't allow more than two!
  if(identical(round(res.points[,-1], digits=0), res.points[,-1])){
    print.digits <- 0
  } else if(identical(round(res.points[,-1], digits=1), res.points[,-1])){
    print.digits <- 1
  } else {
    print.digits <- 2
  }

  ## function latex.head()
  latex.head <- function(one.file=FALSE, individual=TRUE, hist.and.marks=FALSE, marks.hist.stuff=NULL, person=list()){
    # use paste to create the LaTeX head, that is up to the table
    full.head <- paste("\\documentclass[a4paper,ngerman]{scrartcl}\n",
      ifelse(isTRUE(hist.and.marks), "      \\usepackage[a4paper,hmargin={2cm,2cm}]{geometry}\n", ""),
      "      \\usepackage{longtable}
      \\usepackage{mathptmx}
      \\usepackage{helvet}
      \\usepackage{courier}
      \\usepackage[T1]{fontenc}
      \\usepackage[latin9]{inputenc}
        \\setlength{\\parskip}{\\medskipamount}
        \\setlength{\\parindent}{0pt}
      \\usepackage{amsmath}
      \\usepackage{graphicx}\n",
      ifelse(isTRUE(one.file), "      \\usepackage{pdfpages}\n", ""),
      "      \\usepackage{amssymb}
      \\usepackage{thumbpdf}
      \\usepackage{color}
        \\definecolor{dunkelgrau}{gray}{.5}
        \\definecolor{hellgrau}{gray}{.7}
        \\definecolor{dunkelblau}{rgb}{.1,.1,.6}
      \\usepackage[
        a4paper,
        pdftex,
        bookmarksopen,
        bookmarksopenlevel=1,
        colorlinks,
        anchorcolor=dunkelgrau,
        linkcolor=dunkelgrau,
        urlcolor=hellgrau,
        citecolor=dunkelblau]{hyperref}
      \\usepackage{babel}
      \\title{",latex.umlaute(descr$title),"}\n",
      if(isTRUE(individual)){
        # this is for individual reports, so each subject's name is printed
        paste("      \\date{",latex.umlaute(text$DozentIn),": ",latex.umlaute(descr$name),"\\\\",text$Datum,": ",descr$date,"}
      \\author{",person$vorname," ",person$name,",\\\\",text$MatrikelNr," ",person$matn,"}\n", sep="")
      } else {
        # anonymous or global results
        paste("      \\date{",text$Datum,": ",descr$date,"}
      \\author{",latex.umlaute(text$DozentIn),": ",latex.umlaute(descr$name),"}\n", sep="")
      },
      if(isTRUE(individual)){
        # this is for individual reports
        paste("\\begin{document}\n\\maketitle\n\\section*{",text$Ergebnisse,"}\n",
        text$Erreicht,": \\textbf{",person$punkte,"} (\\textbf{",person$prozent,"\\%}",text$Prozent,").
        \\\\",text$Note,": ",format(person$note, nsmall=1),"}
        ",
        marks.hist.stuff,"
          \\newpage\n", sep="")
      } else if(isTRUE(one.file)){
        # combine all individual reports into one file
        paste("\\begin{document}\n\\maketitle\n", sep="")
      } else {
        # anonymous or global results
        paste("\\begin{document}\n\\maketitle\n\\section*{",text$Ergebnisse,"}\n", sep="")
      },
    sep="")
    return(full.head)
  } ## end function latex.head()

  ## function create.pdf()
  # this function will convert LaTeX to PDF
  # it is called in global.report() and tabellenbau()
  # the suppress option is for mergeing files to skip copying the individual reports
  create.pdf <- function(file, path, path.orig=path, suppress=FALSE){
    # save current working directory; unfortuneately, texi2dvi() doesn't seem to be able
    # to create PDF files anywhere but in the WD, so we'll have to cd there and back, afterwards
    current.wd <- getwd()
    # change to destined directory
    setwd(file.path(path))
    texi2dvi(file.path(file), pdf=TRUE, clean=TRUE)
    # in case save and merge were FALSE, move the PDFs to the actual destination
    if(!isTRUE(save) && !isTRUE(suppress) && isTRUE(file.info(path.orig)$isdir)){
      file.copy(gsub(".tex", ".pdf", file), path.orig, overwrite=TRUE)
    } else {}
    # get back to where we came from
    setwd(file.path(current.wd))
  } ## end function create.pdf()

  ## function tabellenbau()
  # this is the main function for individual reports
  tabellenbau <- function(matn){
    points.mtrx <- res.points[res.points$MatrNo==matn,]
    einzelergebnis <- results[results$MatrNo==matn,]
    geg.items <- grep("Item([[:digit:]]{1,3})",names(points.mtrx))
    geg.points <- as.numeric(points.mtrx[,geg.items])
    # the indices of all answers
    items <- grep("Item([[:digit:]]{1,3})",names(answers))
    if(isTRUE(alt.candy)){
      # some eye-candy: put commata between multiple answer alternatives, if applicable
      # calls the internal function answ.alternatives()
      geg.antw1 <- answ.alternatives(answers[answers$MatrNo==matn,items], latex=TRUE)
      loesungen <- answ.alternatives(correct, latex=TRUE)
    } else {
      geg.antw1 <- answers[answers$MatrNo==matn,items]
      loesungen <- correct
    }

    # define table size
    if(identical(table.size, "auto")){
      if(length(items) > 37 & length(items) < 43){
        # to avoid ugly tables with few lines on one page, shrink by heuristics
        table.size <- "\\footnotesize\n"
      } else {
        table.size <- "\\normalsize\n"
      }
    } else {
      table.size <- paste0("\\",table.size,"\n")
    }

    # name and first name from the answer matrix
    name <- latex.umlaute(einzelergebnis$Name)
    vorname <- latex.umlaute(einzelergebnis$FirstName)
    # points and percent of the results
    punkte <- einzelergebnis$Points
    prozent <- einzelergebnis$Percent
    note <- einzelergebnis$Mark

    if (!isTRUE(quiet)){
      # give some feedback on current status
      message(paste("Processing: ", einzelergebnis$FirstName, " ", einzelergebnis$Name, " (", matn, ")", sep=""))
    } else {}

    # check for file name scheme
    if(identical(file.name, "name")){
      name.scheme <- paste(file.umlaute(gsub("[[:space:]]", "_", paste(einzelergebnis$Name, einzelergebnis$FirstName))),".tex", sep="")
    } else {
      name.scheme <- paste(matn,".tex", sep="")
    }
    # create filename from name scheme
    if(isTRUE(save) || isTRUE(pdf)){
      dateiname <- file.path(path, name.scheme)
    } else {
      dateiname <- ""
    }

    # prepare histograms and/or marks info
    if(any(unlist(marks.info))){
      # if informatin on marks is wanted, only grab the intended stuff
      marks.information <- as.matrix(klsr@marks.sum[,unlist(marks.info)])
      colnames(marks.information) <- c(text$Punkte, text$AProzent)[unlist(marks.info)]
      # summary on the defined marks
      # use capture.output() to get rid of the printout; the table tags need to be removed,
      # we need only the tabular environment
      unused.garbage <- capture.output(
          marks.info.table <- print(xtable(marks.information),
          floating=FALSE,
          sanitize.text.function=function(x){latex.umlaute(x)})
        )
      marks.info.tabular <- paste(
        gsub("\\\\end\\{center\\}\\n\\\\end\\{table\\}\\n","",
        gsub("\\\\begin\\{table\\}\\[ht\\]\\n\\\\begin\\{center\\}\\n", "", marks.info.table))
      )
    } else {
      marks.info.table <- ""
      marks.info.tabular <- ""
    }
    # first cases including histograms
    if(any(unlist(hist))){
      if(any(unlist(marks.info))){
        if(sum(unlist(hist)) == 1){
          hist.and.marks <- FALSE
          # one graph and marks
          marks.hist.stuff <- paste("
          \\begin{table}[ht]
            \\begin{minipage}[b]{0.5\\linewidth}
            \\centering
            \\begin{tabular}{c}
              \\includegraphics[width=9cm]{",
              if(hist$points){
                paste(hist.points)
              } else {
                paste(hist.marks)
              },"}
            \\end{tabular}
            \\end{minipage}
            \\hspace{0.8cm}
            \\begin{minipage}[b]{0.5\\linewidth}
              \\small
              \\centering
              ",marks.info.tabular,"
            \\end{minipage}
            \\caption{", latex.umlaute(text$VertErgs), ", ", latex.umlaute(text$Notenschluessel),"}
          \\end{table}%", sep="")
        } else {
          hist.and.marks <- TRUE
          # two graphs and marks
          marks.hist.stuff <- paste("
          \\begin{table}[ht]
            \\begin{minipage}[b]{0.6\\linewidth}
            \\centering
            \\begin{tabular}{c@{\\hskip 0cm}c}
              \\includegraphics[width=0.6\\linewidth]{",hist.points,"}&\\includegraphics[width=0.6\\linewidth]{",hist.marks,"}
            \\end{tabular}
            \\end{minipage}
            \\hspace{0.8cm}
            \\begin{minipage}[b]{0.4\\linewidth}
              \\footnotesize
              \\centering
              ",marks.info.tabular,"
            \\end{minipage}
            \\caption{", latex.umlaute(text$VertErgs), ", ", latex.umlaute(text$Notenschluessel),"}
          \\end{table}%", sep="")
        }
      } else {
        # just one or two graphs
        hist.and.marks <- FALSE
        marks.hist.stuff <- paste("
        \\begin{figure}[h]
        \\centerline{",
        if(hist$points){
          paste("\\mbox{\\includegraphics[width=9cm]{",hist.points,"}}", sep="")
        } else {},
        if(hist$marks){
          paste("\\mbox{\\includegraphics[width=9cm]{",hist.marks,"}}", sep="")
        } else {},
        "}
          \\caption{",text$VertErgs,"}
        \\end{figure}", sep="")
      }
    } else if(any(unlist(marks.info))){
      hist.and.marks <- FALSE
      # cases without histograms but marks
      marks.hist.stuff <- marks.info.table
    } else {
      hist.and.marks <- FALSE
      # cases with neither histograms nor marks
      marks.hist.stuff <- ""
    }

    # here comes the foot, that is after the table
    latex.foot <- paste("
      \\end{document}\n",
    sep="")

    # combine parts to a document
    write(paste(
      latex.head(one.file=FALSE, individual=TRUE, hist.and.marks=hist.and.marks,
      person=list(vorname=vorname, name=name, matn=matn, punkte=punkte, prozent=prozent, note=note),
      marks.hist.stuff=marks.hist.stuff),
      table.size), file=dateiname)
    # create table
    pre.erg.tabelle <- rbind(geg.antw1,loesungen,geg.points)
    rownames(pre.erg.tabelle) <- c(text$Antwort,text$Korrekt,text$Punkte)
    colnames(pre.erg.tabelle) <- names(answers)[items]
    if(isTRUE(NRET.legend)){
      cap.extra <- text$NRET.expl
    } else {
      cap.extra <- ""
    }
    print(xtable(t(pre.erg.tabelle), digits=c(0,0,0,print.digits),
    caption=paste(text$Auswertung," ",vorname," ",name," (",text$MatrikelNr," ",matn,")", cap.extra, sep="")),
    file=dateiname, append=TRUE, sanitize.text.function=function(x){latex.umlaute(x)}, tabular.environment="longtable", floating=FALSE)
    write(latex.foot, file=dateiname, append=TRUE)

    # check if PDF creation is demanded
    if(isTRUE(pdf) && is.character(name.scheme)){
      if(isTRUE(merge)){
        suppress=TRUE
      } else {
        suppress=FALSE
      }
      create.pdf(file=name.scheme, path=path, path.orig=path.orig, suppress=suppress)
    } else {}
  } ## end function tabellenbau()

  ## function global.report()
  global.report <- function(form){
    # set the file name
    if((isTRUE(save) || isTRUE(pdf)) && is.character(anon.glob.file)){
      dateiname <- file.path(path, anon.glob.file)
    } else {
      dateiname <- ""
    }

      # prepare the table
      if(identical(form, "anon")){
        anon.glob.table <- klsr@anon
        colnames(anon.glob.table) <- c(text$Pseudonym,text$Punkte,text$AProzent,text$ANote)
        anon.glob.digits <- c(0,0,print.digits,1,1)
        if (!isTRUE(quiet)){
          # give some feedback on current status
          message(paste("Processing: Anonymous feedback...", sep=""))
        } else {}
      } else {
      anon.glob.table <- klsr@results
      colnames(anon.glob.table) <- c((if(!is.null(anon.glob.table$No)) text$LfdNr),text$Name,text$Vorname,text$GMatNr,text$Punkte,text$AProzent,text$ANote,(if(!is.null(anon.glob.table$Pseudonym)) text$Pseudonym))
      anon.glob.digits <- c(0,(if(!is.null(anon.glob.table$No)) 0),0,0,0,print.digits,1,1,(if(!is.null(anon.glob.table$Pseudonym)) 0))
      if (!isTRUE(quiet)){
        # give some feedback on current status
        message(paste("Processing: Global results...", sep=""))
      } else {}
      }

    # define table size
    if(identical(table.size, "auto")){
      if((nrow(anon.glob.table) > 25 & nrow(anon.glob.table) < 31)
        | (nrow(anon.glob.table) > 68 & nrow(anon.glob.table) < 78)){
        # to avoid ugly tables with few lines on one page, shrink by heuristics
        table.size <- "\\footnotesize\n"
      } else {
        table.size <- "\\normalsize\n"
      }
    } else {
      table.size <- paste0("\\",table.size,"\n")
    }

    latex.foot <- paste(
    if(hist$points | hist$marks){
      paste("
      \\begin{figure}[h]
      \\centerline{",
      if(hist$points){
        paste("\\mbox{\\includegraphics[width=9cm]{",hist.points,"}}", sep="")
      } else {},
      if(hist$marks){
        paste("\\mbox{\\includegraphics[width=9cm]{",hist.marks,"}}", sep="")
      } else {},
      "}
      \\caption{",text$VertErgs,"}
      \\end{figure}", sep="")
    } else {},"
    \\end{document}\n",
    sep="")
    # combine parts to a document
    write(paste(latex.head(one.file=FALSE, individual=FALSE), table.size), file=dateiname)
    # create table with anonymous feedback
    print(xtable(anon.glob.table, digits=anon.glob.digits,
    caption=paste(text$Ergebnisse,": ",latex.umlaute(descr$title)," (",latex.umlaute(descr$name),", ",descr$date,")", sep="")),
    file=dateiname, append=TRUE, sanitize.text.function=function(x){latex.umlaute(x)}, tabular.environment="longtable", floating=FALSE)
    if(sum(unlist(marks.info)) > 0){
      # summary on the defined marks
      print(xtable(marks.information, caption=latex.umlaute(text$Notenschluessel)),
      file=dateiname, append=TRUE, sanitize.text.function=function(x){latex.umlaute(x)}, floating=TRUE)
    } else {}
    write(latex.foot, file=dateiname, append=TRUE)

    # check if PDF creation is demanded
    if(isTRUE(pdf) && is.character(anon.glob.file)){
      create.pdf(file=anon.glob.file, path=path, path.orig=path.orig)
    } else {}
  } ## end function global.report()

  ## function merge.reports()
  # creates one PDF file from the individual reports
  merge.reports <- function(){
    merge.file <-  file.path(path, "individual_reports.tex")
    if(identical(file.name, "name")){
      name.scheme <- sapply(results$MatrNo, function(matn){
          einzelergebnis <- results[results$MatrNo==matn,]
          paste(file.umlaute(gsub("[[:space:]]", "_", paste(einzelergebnis$Name, einzelergebnis$FirstName))),".pdf", sep="")
        })
    } else {
      name.scheme <- paste(results$MatrNo,".pdf", sep="")
    }
    # create filename from name scheme
    all.pdf.files <- file.path(path, name.scheme)

    # here comes the foot
    latex.foot <- paste("
      \\end{document}\n",
    sep="")

    # combine parts to a document
    write(paste(
      latex.head(one.file=TRUE, individual=FALSE),
      paste("      \\includepdf[pages=-]{", all.pdf.files, "}", sep="", collapse="\n"),
      "\n",
      latex.foot, sep=""),
      file=merge.file)

    if (!isTRUE(quiet)){
      # give some feedback on current status
      message(paste("Merging individual reports into one file...", sep=""))
    } else {}

    # check if PDF creation is demanded
    if(isTRUE(pdf)){
      create.pdf(file="individual_reports.tex", path=path, path.orig=path.orig)
    } else {}
  } ## end function merge.reports()

  if(identical(matn, "all")){
    for(i in results$MatrNo) tabellenbau(matn=i)
    if(isTRUE(merge) & isTRUE(pdf)){
      merge.reports()
    } else {}
  } else if(identical(matn, "anon")){
    global.report(form="anon")
  } else if(identical(matn, "glob")){
    global.report(form="global")
  } else {
    tabellenbau(matn)
  }
}
