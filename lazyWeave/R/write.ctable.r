#' @rdname WritePrintCtable
#' @export write.ctable
#' 

write.ctable <- function(x, round = 2, percent = TRUE,
                         quartile=TRUE, cwidth=NULL,
                         caption=NULL, footnote = NULL, 
                         byVarN=FALSE, size="\\normalsize", 
                         descripCombine = TRUE,
                         oddsCombine = TRUE, markSignificant = FALSE, 
                         statHeader="Statistics",
                         name = FALSE, var.label = TRUE, level = TRUE,
                         total = TRUE, descriptive = TRUE, missing=FALSE, 
                         missing.perc=FALSE, testStat = TRUE,
                         odds = FALSE, pval = TRUE, oneLine = FALSE, 
                         pvalFormat="default", pvalArgs=list(), 
                         cat=getOption("lazyWeave_cat"), ...){
  
  reportFormat <- getOption("lazyReportFormat")
  pm <- if (reportFormat %in% "latex") "$\\pm$" else if (reportFormat == "html") "&plusmn" else "$\\pm$"
  ln.break <- if (reportFormat %in% "latex") "\\\\" else if (reportFormat == "html") "<br>" else "  \n"
  prc <- if (reportFormat %in% "latex") "\\%" else "%"
  na.char <- if (reportFormat %in% c("latex", "markdown")) "" else " &nbsp; "
  
  if (!is.factor(attributes(x)$byVar))
    attributes(x)$byVar <- factor(attributes(x)$byVar)
  nlev <- nlevels(attributes(x)$byVar)
  lev <- levels(attributes(x)$byVar)
  
  x <- as.data.frame(x)
  
  if (oneLine){  #*** Set to print only one line for binary variables
    binary.var <- diff(c(which(!is.na(x$name)), nrow(x) + 1)) == 3
    binary.row <- which(x$name %in% x$name[!is.na(x$name)][binary.var])
    delete.row <- sort(c(binary.row + 1, binary.row+2))
    
    upOneRow <- binary.row + 1
    
    oneRowStay <- c("name", "label", "total", "test", "test.mark", "test.stat", "pvalue", "significant", "type")
    oneRowMove <- colnames(x)[!colnames(x) %in% oneRowStay]
    
    x[binary.row, oneRowMove] <- x[upOneRow, oneRowMove]
    x[binary.row, "level"] <- paste(" -", x[binary.row, "level"])
    x <- x[-delete.row, ]
  }
  
  p <- if (percent) 100 else 1
  
  if (is.null(caption)) 
    caption <- paste("Listing of ", latexTranslate(Hmisc::label(attributes(x)$byVar)),
                     "vs. Selected Continuous Measures.")
  
  if (name && var.label){
    name <- FALSE
    warning("Both 'name' and 'label' were TRUE--'name' has been set to FALSE")
  } 
  
  summary.names <- function(r){
    k <- x$type[r]
    if (k == "Bootstrap Mean")       n <- c("boot",  "lowerb", "upperb")
    else if (k == "Parametric Mean") n <- c("mean",   "sd",     "prop")
    else if (k == "Median")
      n <- if (quartile) c("median", "p25",    "p75")
    else          c("median", "min",    "max")
    else                             n <- c("prop",   "mean",   "sd")
    
    n <- as.vector(t(sapply(n, grep, names(x))))
    summ <- x[r, n]
    names(summ) <- paste("stat", rep(LETTERS[1:3], nlev), sep="")
    names(summ) <- paste(names(summ), ".", rep(lev, each=3), sep="")
    summ
  }
  
  combine.stats <- function(r, d=descrip, s=statnames, pm){
    if (descrip$type[r] == "Bootstrap Mean")
      paste(d[r, s[,1]], " (", d[r, s[,2]], ", ", d[r, s[,3]], ")", sep="")
    else if (descrip$type[r] == "Parametric Mean")
      paste(d[r, s[,1]], pm, d[r, s[,2]])
    else if (descrip$type[r] == "Median")
      paste(d[r, s[,1]], " [", d[r, s[,2]], ", ", d[r, s[,3]], "]", sep="")
    else d[r, s[,1]]
  }
  
  #x = as.data.frame(x)
  x[, grep("prop.", names(x))] <-
    x[, grep("prop.", names(x))] * p
  
  descrip <- do.call("rbind", lapply(1:nrow(x), summary.names))
  descrip <- round(descrip, round)
  statnames <- matrix(names(descrip), ncol=3, byrow=TRUE)
  descrip$type <- x$type
  
  combine <- do.call("rbind", lapply(1:nrow(descrip), combine.stats, pm=pm))
  combine <- as.data.frame(combine, stringsAsFactors=FALSE)
  names(combine) <- gsub("statA", "combine", names(combine))
  rownames(combine) <- rownames(descrip) 
  
  odd.var <- c("odds", "odds.lower", "odds.upper", "odds.scale")
  x[, odd.var] <- lapply(x[, odd.var], round, round)
  
  x$pvalue <- ifelse(is.na(x$pvalue), na.char, do.call("pvalString", 
                                                       c(list(p=x$pvalue, 
                                                              format=pvalFormat), 
                                                         pvalArgs)))
  if (reportFormat %in% "latex") x$pvalue <- latexTranslate(x$pvalue)
  x$test.stat <- round(x$test.stat, round)
  x$missing.perc <- ifelse(!is.na(x$missing.perc), format(x$missing.perc, digits=1), x$missing.perc)
  x[is.na(x)] <- na.char
  descrip[is.na(descrip)] <- na.char
  combine[is.na(combine)] <- na.char
  
  descrip$type <- NULL
  
  if (oddsCombine)
    odds.sec <- cbind(ifelse(x$odds != "",
                             paste(x$odds, " (",
                                   x$odds.lower, ", ",
                                   x$odds.upper, ")", sep=""),
                             ""),
                      paste(x$odds.scale, x$odds.unit))
  else odds.sec <- x[, odd.var]
  odds.sec[, 1] <- gsub("1 [(], [)]", "REF", odds.sec[,1])
  
  if (descripCombine) output <- combine else output <- descrip
  tmp <- cbind(output, x[, paste("count.", lev, sep=""), drop=FALSE]) 
  output <- tmp[, c(rbind(matrix(paste("count.", lev, sep=""), ncol=nlev),
                          matrix(names(output), ncol=nlev)))]
  if (total) output <- cbind(x$total, output)
  if (var.label) output <- 
    cbind(ifelse(x$label!="",
                 if (reportFormat %in% "latex") latexTranslate(paste(x$label, x$level, sep="")) 
                 else paste(x$label, x$level, sep=""),
                 if (reportFormat %in% "latex") latexTranslate(paste(x$label, x$level, sep="\\hspace{.2in}"))
                 else if (reportFormat == "html") paste(x$label, x$level, sep=" &nbsp&nbsp&nbsp&nbsp ")
                 else paste(x$label, x$level, sep = " - ")), 
          output)
  if (name) output <- 
    cbind(ifelse(x$name!="",
                 if (reportFormat %in% "latex") latexTranslate(paste(x$name, x$level, sep="")) 
                 else paste(x$name, x$level, sep=""),
                 if (reportFormat %in% "latex") latexTranslate(paste(x$name, x$level, sep="\\hspace{.2in}"))
                 else if (reportFormat == "html") paste(x$name, x$level, sep=" &nbsp&nbsp&nbsp&nbsp ")
                 else paste(x$name, x$level, sep = " - ")), 
          output)
  if (missing) output <- cbind(output, x$missing)
  if (missing.perc) output <- cbind(output, x$missing.perc)
  if (odds) output <- cbind(output, odds.sec)
  if (testStat) output <- cbind(output, x$test.stat)
  if (pval) output <- 
    cbind(output,
          if (reportFormat %in% "latex") 
            paste(x$pvalue, ifelse(!x$test.mark %in% "", paste("$^{", x$test.mark, "}$", sep=""), x$test.mark), sep="")
          else if (reportFormat == "html") 
            paste(x$pvalue, ifelse(!x$test.mark %in% "", paste("<sup>", x$test.mark, "</sup>", sep=""), x$test.mark), sep="")
          else paste(x$pvalue, ifelse(!x$test.mark %in% "", paste("^", x$test.mark, "^", sep=""), x$test.mark), sep=""))
          
  
  denote <- ifelse(rownames(output) %in% attributes(x)$vars,
                   as.character(x$type), NA)
  
  denote <- ifelse(denote %in% c("Chi-Square", "Fisher", "Logistic", "CMH"),
                   "Proportion", denote)
  
  s.type <- levels(x$type)
  s.type <- ifelse(s.type %in% c("Chi-Square", "Fisher", "Logistic", "CMH"),
                   "Proportion", s.type)
  
  type.mark <- if (reportFormat %in% "latex") paste("$^", letters[1:length(s.type)], "$", sep="")
               else if (reportFormat == "html") paste("<sup>", letters[1:length(s.type)], "</sup>", sep="")
               else paste("^", letters[1:length(s.type)], "^", sep="")
  names(type.mark) <- s.type
  
  
  if (descripCombine){
    type.note <- ifelse(s.type == "Bootstrap Mean", paste("Mean (95", prc, " Bootstrap CI)", sep=""),
                        ifelse(s.type == "Parametric Mean", paste("Mean", pm, "SD"),
                               ifelse(s.type == "Median" & quartile, "Median [P25, P75]",
                                      ifelse(s.type == "Median" & !quartile, "Median [Min, Max]",
                                             ifelse(percent, "Percentage", "Proportion")))))
    type.note <- paste(paste(type.mark, type.note, sep="  "), collapse="; ")
  }
  else{
    type.note <- ifelse(s.type == "Bootstrap Mean", paste("Mean, 95", prc, " Boostrap CI", sep=""),
                        ifelse(s.type == "Parametric Mean", "Mean, SD",
                               ifelse(s.type == "Median" & quartile, "Median, P25, P75",
                                      ifelse(s.type == "Median" & !quartile, "Median, Min, Max",
                                             ifelse(percent, "Percentage", "Proportion")))))
    type.note <- paste(paste(type.mark, type.note, sep="  "), collapse="; ")
  }
  
  tmark <- paste(x$test.mark, ": ", x$test, sep="")
  tmark <- sort(tmark[!duplicated(x$test.mark) & x$test.mark != na.char])
  tmark <- paste(tmark, collapse=ln.break)
  
  fnote <- paste(type.note, tmark, sep=ln.break)
  
  
  
  output[,1] <- ifelse(!is.na(type.mark[denote]),
                       paste(output[,1], type.mark[denote]),
                       as.character(output[,1]))
  output <- as.data.frame(lapply(output, as.character), stringsAsFactors=FALSE)
  
  if(markSignificant){
    x$significant <- as.logical(x$significant)
    x$significant[is.na(x$significant)] <- FALSE
    output[x$significant, ] <- 
      lapply(output[x$significant, ], 
             function(x){ if (reportFormat == "latex") paste("\\textbf{", x, "}", sep="") 
                          else if (reportFormat == "html") paste("<b>", x, "</b>", sep="")
                          else if (reportFormat == "markdown") paste("**", x, "**", sep="")})
  }
  
  #******************************************************************************
  #* For LaTeX output (PDF)
  #* a. Part 1
  #* b. Part 2
  #* c. Part 3
  #******************************************************************************
  
  #   if (reportFormat %in% "latex"){
  cspan <- c( rep(1, name), rep(1, var.label), rep(1, total),
              {if (descripCombine) rep(2, nlev) else rep(4, nlev)},
              rep(1, missing), rep(1, missing.perc), rep(1, testStat),
              if (odds){ if (oddsCombine) 2 else if(oddsCombine) 4},
              rep(1, pval))
  if (reportFormat == "markdown") cspan <- rep(1, sum(cspan))
  cwidth1 <- rep("", length(cspan))
  
  align.lev <- if (reportFormat == "markdown") {if (descripCombine) nlev * 2 else nlev * 4}
               else nlev
  odds.lev <- if (reportFormat == "markdown") {if (oddsCombine) 2*odds else 4*odds} else odds

  align <- c( rep("l", name + var.label),                
              rep("c", total + align.lev + missing + missing.perc + testStat + odds.lev ),
              rep("r", pval))
  head <- c( if(name) na.char,
             if(var.label) na.char,
             if(total) na.char,
             if (reportFormat == "markdown") unlist(lapply(lev, function(x) c(x, rep("", as.numeric(!descripCombine*2) + 1)))) else lev,
             if(missing) na.char,
             if(missing.perc) na.char,
             if(testStat) na.char,
             if (odds) {if (reportFormat == "markdown") {if (oddsCombine) 2 else 4} else na.char},
             if(pval) na.char)
  

  if (options("lazyReportFormat") == "latex") head <- latexTranslate(head)
  
  #    return(name + var.label + total + 1 + nlev)
  part1 <- lazy.table(head, align=align, justify="left", 
                      cspan=cspan, cwidth=cwidth1,
                      rborder=if(!byVarN) c(0, 0, 1) else c(0, 0), 
                      rbspan=if (options("lazyReportFormat") == "latex"){ 
                        c(name+var.label+total+1,
                          name+var.label+total+2 * nlev + missing + missing.perc + (1 + (nlev - 1) * 3)* !descripCombine)
                      }
                      else c((name + var.label + total + 1) : (name + var.label + total + nlev)),
                      caption=caption, size=size, 
                      close=FALSE, translate=FALSE, cborder=NULL, 
                      cat = FALSE, ...)

  if(byVarN){ 
    Nline <-   head <- c( if(name) "",
                          if(var.label) "",
                          if(total) "",
                          paste("(N = ", table(attributes(x)$byVar), ")", sep=""),
                          if(missing) "",
                          if(missing.perc) "",
                          if(testStat) "",
                          if(odds) "",
                          if(pval) "")
    Nline <- lazy.table(Nline, align=align, justify="left",
                        cspan=cspan, cwidth=cwidth1,
                        rborder=1, 
                        rbspan=if (options("lazyReportFormat") == "latex"){ 
                          c(name+var.label+total+1,
                            name+var.label+total+2 * nlev + missing + missing.perc + (1 + (nlev - 1) * 3)* !descripCombine)
                        }
                        else c(name + var.label + total + 1, name + var.label + total + 1 + nlev),
                        open=FALSE, close=FALSE, translate=FALSE, cborder=NULL,
                        cat = FALSE)
    part1 <- paste(part1, Nline)
  }
  
  #******************************************************************************
  #*** Part 2
  
  cspan <- c( rep(1, name), rep(1, var.label), rep(1, total),
{if (descripCombine) rep(c(1, 1), nlev) else rep(c(1, 3), nlev)},
              rep(1, missing), rep(1, missing.perc),
              if (odds){ if (oddsCombine) 2 else if(oddsCombine) 4},
              rep(1, testStat),
              rep(1, pval))
  
  align <- c( rep("l", name + var.label),
              rep("c", total + 2 * nlev + missing + missing.perc + testStat + odds ),
              rep("r", pval))
  
  if (is.null(cwidth)){
    cwidth <- rep("", length(cspan))
    if (var.label) cwidth[name + var.label] <- 1.5
  }
  if (length(cspan) != length(cwidth)) 
    stop(paste("With specified arguments, cwidth must have length",
               length(cspan)))
  head <- c( if(name) "Factor",
             if(var.label) "Factor",
             if(total) "Total",
             rep(c("N", statHeader), nlev),
             if(missing) "Missing",
             if(missing.perc) latexTranslate("% Missing"),
             if(odds) "Odds Ratio",
             if(testStat) "Test Statistic",
             if(pval) "p-value")
  
  part2 <- lazy.table(head, align=align, cspan=cspan, 
                      cwidth=cwidth, open=FALSE, close=FALSE, translate=FALSE, 
                      rborder=1, cat = FALSE, ...)
  
  
  #******************************************************************************
  #*** Part 3
  
  align <- c(rep("l", name + var.label),
             rep("c", total + (2 * nlev + (nlev * 2 * !descripCombine)) + missing + missing.perc +
                   testStat + (odds + 1*odds + (2 * !oddsCombine))),
             rep("r", pval))
  
  part3 <- lazy.table(output, align=align, open=FALSE, justify="left",
                      cwidth=c(cwidth, rep("", odds)),
                      rborder=c(0, nrow(output)),
                      footnote=paste(fnote, footnote, sep=ln.break),
                      translate=FALSE, 
                      cat = FALSE, ...)
  
  if (cat) cat(paste(part1, part2, part3, sep="\n"))
  else return(paste(part1, part2, part3, sep="\n"))
}
