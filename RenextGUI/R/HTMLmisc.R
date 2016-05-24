##*******************************************************************
## Copy of R2HTML generic function.
##
##********************************************************************

HTML <- function (x, ...) { UseMethod("HTML") }


##*******************************************************************
## Function to HTML-ize summary results for objects with (old) class
## "Renouv" 
##
##********************************************************************

HTML.summary.Renouv <-
  function(x,
           file = get(".HTML.file"),
           append,
           coef = TRUE,
           pred = TRUE,
           probT = FALSE,
           digits = max(3, getOption("digits") - 3),
           symbolic.cor = x$symbolic.cor,
           signif.stars = getOption("show.signif.stars"),
           ...) {

    
    cat("<div class=\"summary\">\n", file = file, append = append)
    cat("<ul>\n", file = file, append = append)

    ## sample info
    cat("<li>Main sample 'Over Threshold'\n", file = file, append = append)

    cat("   <ul>\n", file = file, append = append)
    cat(sprintf(paste("      <li>Threshold        <span class=\"res\">%8.2f</span></li>\n",
                      "      <li>Effect. duration <span class=\"res\">%8.2f</span></li>\n",
                      "      <li>Nb. of exceed.   <span class=\"res\">%5d</span></li>\n",
                      collapse = ""),
                x$threshold, x$effDuration, length(x$y.OT)),
        file = file, append = append)
    
    cat("   </ul>\n", file = file, append = append)
    cat("</li>\n\n", file = file, append = append)
    
    ## rate info
    cat(sprintf(paste("<li>Estimated rate 'lambda' for Poisson process (events):",
                      "<span class=\"res\">%5.2f</span> evt/year.</li>\n\n"),
                x$estimate["lambda"]),
        file = file, append = append)
    ## dist info
    cat("<li>\n", file = file, append = append)
    cat(sprintf(paste("Distribution for exceedances y: <span class=\"res\">\"%s\"</span>,",
                      "with <span class=\"res\">%d</span> par. "),
                x$distname.y, x$p.y),
        file = file, append = append)

    ## print on same same line or not
    if (length(x$parnames.y) > 2) {
      cat("<br/>\n", file = file, append = append)
    }
    cat(paste(sprintf("<span class = \"res\">\"%s\"</span>", x$parnames.y), collapse = ", "),
          file = file, append = append)
    
    cat("\n</li>\n", file = file, append = append)
    
    if (coef) {
      cat("<li>Coefficients\n", file = file, append = append)
      
      if (!probT) {
        HTML(x$coefficients[ , 1:3],
             align="left",
             file = file, append = append)
        
      } else {
        HTML(x$coefficients,
             digits = digits,
             signif.stars = signif.stars,
             na.print = "NA", ...,
             file = file, append = append)
      }
      cat("   <p>\n", file = file, append = append)
      cat(sprintf(paste("Degrees of freedom: <span class=\"res\">%d</span> (param.)",
                        "and <span class=\"res\">%d</span> (obs)\n"),
                  x$df["par"], x$df["obs"]),
          file = file, append = append)
      cat("   </p>\n", file = file, append = append)
      cat("</li>\n", file = file, append = append)
      
    }
    if (any(x$fixed)) {
      cat("The following coef. were fixed\n")
      print(names(x$fixed)[x$fixed])
      cat("\n")
    }
 
    if (pred) {
      cat("<li>Return levels\n", file = file, append = append)
      M <- as.matrix(format(roundPred(x$pred)))
      rownames(M) <- NULL
      
      HTML(x = M,
           file = file,
           row.names = FALSE,
           align = "left",
           nsmall = c(3, 5, 2, rep(2, 4)),
           append = TRUE)
      
      cat("</li>", file = file, append = append)
    }
    if (x$history.MAX$flag) {
      cat(sprintf(paste("<li> 'MAX' historical info: <span class=\"res\">%d</span> blocks,",
                        "<span class=\"res\">%d</span> obs.,",
                        "total duration = <span class=\"res\">%5.2f</span> years\n"),
                  nlevels(x$history.MAX$block),
                  length(unlist(x$history.MAX$data)),
                  sum(x$history.MAX$effDuration)),
          file = file, append = append)

      cat("<br></br>\n<ul>\n", file = file, append = append)
      
      for ( i in 1:nlevels(x$history.MAX$block) ) {

        Zi <- x$history.MAX$data[[i]]
        ri <- length(Zi)

        if ( ri <= 12L ) {
          obs.str <- paste(paste("<span class=\"res\">", format(Zi), "</span>", sep = ""),  collapse = ", ")
        } else {
          obs.str <- paste(paste(paste("<span class= \"res\">", format(Zi[1L:3L]), "</span>", sep = ""),
                                 collapse = ", "),
                           "...",
                           paste(paste("<span class =\"res\">", format(Zi[(ri-2L):ri]), "</span>", sep = ""),
                                 collapse = ", "),
                           sep = ", ")
        }
        cat(sprintf(paste("<li> block <span class=\"res\">%d</span>:",
                          "<span class=\"res\">%5.2f</span> years,",
                          "<span class=\"res\">%d</span> obs.\n<br></br>\n",
                          "%s  \n</li>"),
                    i, x$history.MAX$effDuration[[i]], ri, obs.str),
            file = file, append = append)
        
      }
      
      cat("</ul>\n</li>\n", file = file, append = append)
      
    } else {
      cat("<li> no 'MAX' historical data</li>\n",
          file = file, append = append)
    }
    
    if (x$history.OTS$flag) {
      
      cat(sprintf(paste("<li> 'OTS' historical info: <span class=\"res\">%d</span> blocks,",
                        "<span class=\"res\">%d</span> obs.,",
                        "total duration = <span class=\"res\">%5.2f</span> years\n"),
                  nlevels(x$history.OTS$block),
                  length(unlist(x$history.OTS$data)),
                  sum(x$history.OTS$effDuration)),
          file = file, append = append)


      cat("<br></br>\n<ul>\n", file = file, append = append)
      
      for ( i in 1:nlevels(x$history.OTS$block) ) {
        
        Zi <- x$history.OTS$data[[i]]
        ri <- length(Zi)

        if ( (ri > 1L) && (ri <= 12L) ) {
          obs.str <- paste(paste("<span class=\"res\">", format(Zi), "</span>", sep = ""),  collapse = ", ")
        } else {
          obs.str <- paste(paste(paste("<span class= \"res\">", format(Zi[1L:3L]), "</span>", sep = ""),
                                 collapse = ", "),
                           "...",
                           paste(paste("<span class =\"res\">", format(Zi[(ri-2L):ri]), "</span>", sep = ""),
                                 collapse = ", "),
                           sep = ", ")
        }
        cat(sprintf(paste("<li> block <span class=\"res\">%d</span>:",
                          "<span class=\"res\">%5.2f</span> years,",
                          "<span class=\"res\">%d</span> obs.\n<br></br>\n",
                          "%s  \n</li>"),
                    i, x$history.OTS$effDuration[[i]], ri, obs.str),
            file = file, append = append)
      }
      
      cat("</ul>\n</li>\n", file = file, append = append)

      
    } else  {
      cat("<li> no 'OTS' historical data</li>\n",
          file = file, append = append)
    }


    ##==================================================================
    ## KOLMOGOROV TEST
    ##
    ## NB. Not dood, because of unsuitable tags <h2/>, ...
    ## HTML(x$KS.test, file = file, append = append)
    ##===================================================================
    
    cat("<li> Kolmogorov-Smirnov test\n", file = file, append = append)
    cat("   <ul>\n", file = file, append = append)
    cat(sprintf("      <li>method <span class=\"res\">%s</span></li>\n", x$KS.test$method),
        file = file, append = append)
    cat(sprintf("      <li>data <span class=\"res\">%s</span></li>\n", x$KS.test$data.name),
        file = file, append = append)
    cat(sprintf("      <li>alternative hyp. <span class=\"res\">%s</span></li>\n",
                x$KS.test$alternative),
        file = file, append = append)
    cat(sprintf(paste("      <li>D = <span class=\"res\">%6.4f</span>, ",
                      "p-value = <span class=\"res\">%6.4f</span></li>\n"),
                x$KS.test$stat, x$KS.test$p.value),
        file = file, append = append)  

    cat("   </ul>\n", file = file, append = append)
    cat("</li>\n\n", file = file, append = append)

    ##==================================================================
    ## BARTLETT's TEST
    ##
    ## NB. Not dood, because of unsuitable tags <h2/>, ...
    ## HTML(x$KS.test, file = file, append = append)
    ##===================================================================

    if (x$distname == "exponential") {
      cat("<li>Test for exponentiality\n", file = file, append = append)
      cat("   <ul>\n", file = file, append = append)
      cat(sprintf("      <li>method <span class=\"res\">%s</span></li>\n", x$expon.test$method),
          file = file, append = append)
      cat(sprintf(paste("      <li>stat = <span class=\"res\">%6.2f</span>, ",
                        "df = <span class=\"res\">%d</span>, ",
                        "p-value = <span class=\"res\">%6.4f</span></li>\n"),
                  x$expon.test$statistic, x$expon.test$df, x$expon.test$p.value),
          file = file, append = append)  
      cat("   </ul>\n", file = file, append = append)
      cat("</li>\n\n", file = file, append = append)
    }

    
    ## copied from 'print.summary.lm' in the 'stats' package
    correl <- x$correlation
    if (!is.null(correl)) {
      p <- NCOL(correl)
      if (p > 1L) {
        cat("<li> Correlation of Coefficients:\n",
            file = file, append = append)
        if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
          HTML(symnum(correl, abbr.colnames = NULL),
               file = file, append = append)
          cat("\n", file = file, append = append)
        } else {
          correl <- format(round(correl, 2), nsmall = 2, digits = digits)
          correl[!lower.tri(correl)] <- ""
          HTML(correl[-1, -p, drop = FALSE], quote = FALSE,
               file = file, append = append)
          cat("\n", file = file, append = append)
        }
        cat("</li>\n", file = file, append = append)
      }
    }
    cat("</ul>\n", file = file, append = append)
    cat("</div>\n", file = file, append = append)
  }



