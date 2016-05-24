###     -*- Coding: utf-8 -*-          ###
### Analyste Charles-Édouard Giguère   ###

### Creation of an object to manage a two-dimensional table.
cross <- function(x, ...){
    ## if x is not a table we first create the table
    ## with t and ... options .
    if( !( "table" %in% class(x) )){
        if( ("formula" %in% class(x)) ){
            T <- xtab(x, ...)
        }
        else{
            T <- table(x, ...)
        }
    }
    else{
        T <- x
    }
    cr.T <- list()
    cr.T$T <- T
    cr.T$PCT.ROW.T <- addmargins(prop.table(T,1), 2, Total)
    cr.T$PCT.COL.T <- addmargins(prop.table(T,2), 1, Total)
    cr.T$NAMES <- names(dimnames(cr.T$T))
    names(dimnames(cr.T$T)) <- NULL
    names(dimnames(cr.T$PCT.ROW.T)) <- NULL
    names(dimnames(cr.T$PCT.COL.T)) <- NULL
    class(cr.T) <- "cross"
    cr.T
}

### Methods to print a cross object.
### A list of tests can be passed as parameters.
### A chi-square test is displayed by default.
### Results can be exported to pdf (latex is needed to do that).
###  or to an Excel spreadsheet (xlsx).
print.cross <- function(x, ..., test = "chisq.test", export = NULL){
    p <- paste(x$NAMES[1], x$NAMES[2], sep = "*")
    cat("Contingency table for ", p, "\n\n", sep = "")
    print.table(round(addmargins(x$T, FUN = Total, quiet = TRUE), 1),
                zero.print = ".")

    cat("\n\n", "Percentages by row for ", p, "\n\n", sep = "")
    print.table(round(x$PCT.ROW.T * 100, 2), zero.print = ".")

    cat("\n\n", "Percentages by column for ", p, "\n\n", sep = "")
    print.table(round(x$PCT.COL.T*100, 2), zero.print = ".")

    cat("\n\n", "Statistics for ", p, "\n\n", sep = "")
    for(i in test){
        FUN <- match.fun(i)
        print(FUN(x$T))
    }
    if("pdf" %in% export){
        wd <- setwd(Sys.getenv("temp"))
        f1 <- file("cross_output.tex",
                   open="w")
        writeLines(con = f1,
                   text = c("\\documentclass{article}",
                            "\\begin{document}",
                            "\\begin{table}[htp!]",
                            "\\centering",
                            "\\pagenumbering{gobble}"))
        
        xt1 <- xtable(round(addmargins(x$T, FUN = Total, quiet = TRUE), 1),
                      align = paste("l|",paste(rep("r",dim(x$T)[2]),collapse=""),
                                    "|r",sep=""),
                      digits= 0)
        print.xtable(xt1,file=f1,append=TRUE,
                     hline.after = c(0,0,dim(x$T)[1]),
                     floating=FALSE)
        cat("\\caption{Contingency table for ",
            sub("[*]","$\\\\times$",p),"}",file=f1,"")
        cat("\\end{table}", file = f1)
        writeLines(con = f1,
                   text = c("\\begin{table}[htp!]",
                            "\\centering"))
        xt2 <- xtable(round(x$PCT.ROW.T*100, 2),
                      align = paste("l|",
                                    paste(rep("r",dim(x$PCT.ROW.T)[2]-1),
                                          collapse=""),"|r",
                                    sep=""),
                      digits= 2)
        print.xtable(xt2,file=f1,append=TRUE,
                     hline.after = c(0,0,dim(x$PCT.ROW.T)[1]),
                     floating=FALSE)
        cat("\\caption{Percentages by row for ",
            sub("[*]","$\\\\times$",p),"}",file=f1,"")
        cat("\\end{table}", file = f1)
                writeLines(con = f1,
                   text = c("\\begin{table}[htp!]",
                            "\\centering"))
        xt2 <- xtable(round(x$PCT.COL.T*100, 2),
                      align = paste("l|",
                                    paste(rep("r",dim(x$PCT.COL.T)[2]),
                                          collapse=""),"|",
                                    sep=""),
                      digits= 2)
        print.xtable(xt2,file=f1,append=TRUE,
                     hline.after = c(0,0,dim(x$PCT.COL.T)[1]-1),
                     floating=FALSE)
        cat("\\caption{Percentages by column for ",
            sub("[*]","$\\\\times$",p),"}",file=f1,"",fill=TRUE)
        cat("\\end{table}", file = f1,fill=TRUE)
        if("chisq.test" %in% test){
            test <- setdiff(test,"chisq.test")
            tst <- chisq.test(x$T)
            writeLines(c("\\begin{table}",
                         "\\centering",
                         "\\begin{tabular}{rrr}",
                         "\\hline",
                         "$\\chi^2$ & df & p.value\\\\",
                         "\\hline\\hline",
                         sprintf("%5.2f & %d & %6.4f\\\\",tst$statistic,
                                 tst$parameter,tst$p.value),
                         "\\hline",
                         "\\end{tabular}",
                         paste("\\caption{$\\chi^2$ test for ",
                               sub("[*]","$\\\\times$",p),"}"),
                         "\\end{table}"),
                         con=f1)
        }
        writeLines(c("\\begin{table}",
                         "\\centering",
                         "\\begin{verbatim}"),
                       con = f1)
        for(i in test){

            FUN <- match.fun(i)
            writeLines(capture.output(FUN(x$T)), con = f1)


        }
        writeLines(c("\\end{verbatim}",
                     paste("\\caption{Other tests for ",
                           sub("[*]","$\\\\times$",p),"}"),
                     "\\end{table}"),
                   con = f1)
        cat("\\end{document}", file = f1)
        close(f1)
        out <- shell("pdflatex cross_output.tex",intern = TRUE)
        shell("start cross_output.pdf")
        setwd(wd)
    }
    if("xlsx" %in% export){
        wd <- setwd(Sys.getenv("temp"))
        
        wb <- createWorkbook("cross_output.xlsx")
        addWorksheet(wb, "Frequencies")
        addWorksheet(wb, "Row percentages")
        addWorksheet(wb, "Column percentages")
        writeData(wb, 1, paste("Contingency table for",p,
                               sep = " "),2,2)
        writeData(wb, 1, addmargins(x$T, FUN = Total, quiet = TRUE), 2,3)
        writeData(wb, 2, paste("Row percentages for",p,
                               sep = " "),2,2)
        writeData(wb, 2, round(x$PCT.ROW.T * 100,2),2,3)
        writeData(wb, 3, paste("Col percentages for",p,
                               sep = " "),2,2)
        writeData(wb, 3, round(x$PCT.COL.T * 100,2),2,3)
        addWorksheet(wb, "Statistics")
        stat.output <- character()
        for(i in test){
            FUN <- match.fun(i)
            stat.output <- c(stat.output,
                             capture.output(FUN(x$T)))

        }
        writeData(wb, 4, stat.output,2,3)
        saveWorkbook(wb, file = "cross_output.xlsx", overwrite = TRUE)
        shell("start cross_output.xlsx")
        setwd(wd)
    }

}


### Fonction xtab overloads xtabs with more parameters to handle
### missing variables in categorical variables.
xtab <- function(formula, data = parent.frame(), useNA = FALSE,
                 exclude = c(NA,NaN), miss.char = "-", na.action = na.exclude,
                 subset = NULL, sparse = FALSE, drop.unused.levels = FALSE){
    dtaNA <- model.frame(formula,data,na.action = na.pass)

    if( useNA ){
        for(i in names(dtaNA)){
            if( is.factor(dtaNA[[i]]) ){
                dtaNA[[i]] <- addNA(dtaNA[[i]],ifany=TRUE)
            }
        }
        exclude <- setdiff(exclude,c(NA,NaN))
        if(length(exclude)==0){
            xt <- do.call("xtabs",list(formula=formula,data = dtaNA,
                                       exclude = NULL, na.action = na.pass,
                                       subset=subset, sparse = sparse,
                                       drop.unused.levels = drop.unused.levels))
            for(i in seq(along = dimnames(xt)))
                dimnames(xt)[[i]][is.na(dimnames(xt)[[i]])] <- miss.char
        }
        else{
            xt <- do.call("xtabs",list(formula =formula,data=dtaNA,
                                       exclude = exclude, na.action = na.pass,
                                       subset=subset, sparse = sparse,
                                       drop.unused.levels = drop.unused.levels))
            for(i in seq(along = dimnames(xt)))
                dimnames(xt)[[i]][is.na(dimnames(xt)[[i]])] <- miss.char

        }
    }
    else{
        xt <- do.call("xtabs",list(formula, data = data, subset = subset,
                                   sparse = sparse, na.action = na.action,
                                   exclude = exclude,
                                   drop.unused.levels = drop.unused.levels))
        for(i in seq(along = dimnames(xt)))
            dimnames(xt)[[i]][is.na(dimnames(xt)[[i]])] <- miss.char

    }
    xt
}


### Function that overloads sum with na.rm=TRUE and replaces NA with 0.
Total = function(x){
    sx <- sum(x, na.rm = TRUE)
    ifelse(is.na(sx), 0, sx)
}
