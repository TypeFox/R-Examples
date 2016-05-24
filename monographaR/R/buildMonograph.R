buildMonograph <-
function(headings, tableToDescription.data, examinedSpecimens.data=NULL, collectorList.data=NULL, output="Word", title="Taxonomic treatment") {
   
  if (class(headings) != "data.frame") {
    stop("headings must be a data.frame")
  }
  if (ncol(headings) != 3) {
    stop("headings must have 3 columns, see help(\"buildMonograph\")")
  }
  if (class(tableToDescription.data) != "data.frame") {
    stop("tableToDescription.data must be a data.frame")
  }
  if (is.null(examinedSpecimens.data)) {
    include.specimens = F
  } else {
    if (class(examinedSpecimens.data) != "data.frame") {
      stop("examinedSpecimens.data must be a data.frame")
    }
    if (ncol(examinedSpecimens.data) != 8) {
      stop("examinedSpecimens.data must have 8 columns, see help(\"examinedSpecimens\")")
    }
    include.specimens = T
    sub(" ", "_", examinedSpecimens.data[,1]) -> examinedSpecimens.data[,1]
  }
  if (is.null(collectorList.data)) {
    include.collector = F
  } else {
    if (class(collectorList.data) != "data.frame") {
      stop("collectorList.data must be a data.frame")
    }
    if (ncol(collectorList.data) < 3 | ncol(collectorList.data) == 4 | ncol(collectorList.data) > 5) {
      stop("collectorList.data must have 3 or 5 columns, see help(\"collectorList\")")
    }
    include.collector = T
    as.factor(sub(" ", "_", collectorList.data[,1])) -> collectorList.data[,1]
  }
  
  sub(" ", "_", headings[,1]) -> headings[,1]
  headings[,1] -> spp
  sort(spp) -> spp
  sub("\\.", "_", spp) -> spp
  sub("\\.", "_", colnames(tableToDescription.data)) -> colnames(tableToDescription.data)
  
  cat(  "---\n", "title: ", "\"", title, "\"\n", file = "temp.Rmd", sep="")
  if (output == "html") {
    cat("output: rmarkdown::html_vignette\n", file = "temp.Rmd", append=T)
  }  
  if (output == "Word") {
    cat("output: rmarkdown::word_document\n", file = "temp.Rmd", append=T)
  }
  cat(
    "vignette: >\n",
    " %\\VignetteIndexEntry{Monograph skeleton}\n",
    " %\\VignetteEngine{knitr::rmarkdown}\n",
    " %\\VignetteEncoding{UTF-8}\n",
    "---\n",
    file = "temp.Rmd", sep="", append = T)
  cat("&nbsp;\n", fill=T, file = "temp.Rmd", sep="", append=T)
  cat("### Taxonomic treatment\n", "&nbsp;\n", fill=T, file = "temp.Rmd", sep="", append=T)
  if (include.specimens) {
    cat("```{r, echo=FALSE, results='asis'}\n",
        "for (i in 1:length(spp)) {\n", 
        "spp[i] -> sp0\n", 
        "as.character(headings[match(sp0, headings[,1]),2]) -> head0\n",
        "head1 <- strsplit(head0, \" \")[[1]][-c(1:2)]\n", 
        "paste(head1, collapse=\" \") -> head1\n",
        "head0 <- strsplit(head0, \" \")[[1]][c(1:2)]\n", 
        "head0 <- paste(head0, collapse=\" \")",
        "cat(\"*\", head0, \"* \", head1, \"\\n\", sep=\"\", file=\"\")\n", 
        "cat(\" \", file=\"\", fill=T)\n",
        "kc <- match(sp0,colnames(tableToDescription.data))\n",
        "desc.data <- tableToDescription.data[,c(1:3,kc)] \n", 
        "tableToDescription2(desc.data)\n", 
        "cat(\"&nbsp;\", file=\"\", fill=T)\n",
        "cat(\" \", file=\"\", fill=T)\n",
        "cat(\"**Comments:** \")\n", 
        "cat(as.character(headings[match(sp0, headings[,1]),3]), \"\\n\", sep=\"\", file=\"\")\n", 
        "cat(\"&nbsp;\\n\", file=\"\", fill=T)\n",
        "exam.data <- examinedSpecimens.data[which(examinedSpecimens.data[,1] == sp0),]\n",
        "cat(\"**Examined specimens:** \", file=\"\", fill=F)\n",
        "cat(examinedSpecimens2(exam.data), file=\"\", fill=T)\n",
        "cat(\"&nbsp;\\n\", file=\"\", fill=T)\n",
        "}\n", "```\n",
        file = "temp.Rmd", sep="", append = T, fill=T)
  } else {
    cat("```{r, echo=FALSE, results='asis'}\n",
        "for (i in 1:length(spp)) {\n", 
        "spp[i] -> sp0\n", 
        "cat(as.character(headings[match(sp0, headings[,1]),2]), \"\\n\", sep=\"\", file=\"\")\n", 
        "cat(\" \", file=\"\", fill=T)\n",
        "kc <- match(sp0,colnames(tableToDescription.data))\n",
        "desc.data <- tableToDescription.data[,c(1:3,kc)] \n", 
        "tableToDescription2(desc.data)\n", 
        "cat(\"&nbsp;\", file=\"\", fill=T)\n",
        "cat(\" \", file=\"\", fill=T)\n",
        "cat(\"**Comments:** \")\n", 
        "cat(as.character(headings[match(sp0, headings[,1]),3]), \"\\n\", sep=\"\", file=\"\")\n", 
        "cat(\"&nbsp;\\n\", file=\"\", fill=T)\n",
        "cat(\"&nbsp;\\n\", file=\"\", fill=T)\n",
        "}\n", "```\n",
        file = "temp.Rmd", sep="", append = T, fill=T)
  }
  if (include.collector) {
    collectorList.data <- collectorList.data[which(is.na(match(collectorList.data[,1], spp)) == F),]
    cat("### List of all specimens by collector\n", fill=T, file = "temp.Rmd", sep="", append=T)
    cat("```{r, echo=FALSE, results='asis'}\n",
        "collectorList2(data=collectorList.data)\n", "```\n",
        fill=T, file = "temp.Rmd", sep="", append=T)
  }
  if (output == "Word") {
    suppressWarnings(rmarkdown::render("temp.Rmd", rmarkdown::word_document(), output_file = "Monograph_skeleton.docx", quiet=T))
    try(system("open Monograph_skeleton.docx"), silent=T)
  }
  if (output == "html") {
    suppressWarnings(rmarkdown::render("temp.Rmd", rmarkdown::html_vignette(), output_file = "Monograph_skeleton.html", quiet=T))
    try(system("open Monograph_skeleton.html"), silent=T)
  }
  unlink("temp.Rmd")
  cat("The monograph skeleton was saved in:")
  cat("\n", getwd())
}
