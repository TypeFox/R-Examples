cat("   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   +++                     Begin Functions               +++
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

redToBlue <- function(n,seqColor=matrix(c(1,0,0, 1,0,1, 0,0,1),3)){
    c1 <- seqColor[,1]
    c2 <- seqColor[,2]
    c3 <- seqColor[,3]
    cat("C1=",c1," C2=",c2," C3=",c3,"\n")
    s <- seq(0,1,,n)
    return(c(rgb(c1[1]*(1-s)+c2[1]*s,c1[2]*(1-s)+c2[2]*s,c1[3]*(1-s)+c2[3]*s),
             rgb(c2[1]*(1-s)+c3[1]*s,c2[2]*(1-s)+c3[2]*s,c2[3]*(1-s)+c3[3]*s)))
}

redToBlue <- function(n){
    sequ <- seq(0,1,,n)
    return(rgb(pmin(1,2-sequ*2),0,pmin(1,sequ*2)))
}


######################################################################
###                           High level                           ###
######################################################################

setClass(Class="continuous",contains=c("numeric"))
setClass(Class="discrete",contains=c("numeric"))


changeClass <- function(x,limDiscreteX){
    if(class(x)[1]=="character"){x <- as.factor(x)}else{}
    if(class(x)[1]=="logical"){x <- as.factor(x)}else{}
    if(length(table(x))==2){class(x) <- c("logical",class(x))}
    if(class(x)[1] %in% c("numeric","integer")){
        if(length(table(x))<=limDiscreteX) {
            class(x) <- c("discrete","numeric")
        }else{
            class(x) <- c("continuous","numeric")
        }
    }else{}
    return(x)
}



###############################################################
###      Utility functions to write out pieces of text      ###
###                    Both univ and biv                    ###
###############################################################

### This is the ratio used when converting between centimeters and pixels.
### 80 seems to be a reasonable choice. HTML <IMG> tag expects pixels,
### LaTeX \includegraphics is in centimeters.
r2lPixelsPerCm <- function() return(80)


r2lEol <- function(str, LF=TRUE) {
    if (LF) {
        return(paste(str, "\n", sep=""))
    } else {
        return(str)
    }
}


r2lComment <- function(x, out="latex", LF=TRUE) {
    str <- switch(out,
        "html" = paste("<!--", as.character(x), "-->"),
	"latex" = paste("%", as.character(x))
    )
    return(r2lEol(str, LF=LF))
}


r2lStartTable <- function(out="latex", hline=TRUE, LF=TRUE, attrs=NULL) {
	str <- switch(out,
		"html" = paste("<TABLE", attrs, ">", sep=""),
		"latex" = paste("\\begin{tabular}{", attrs, "}", "\n",if(hline){"\\hline"}else{""}, sep="")
	)
	return(r2lEol(str, LF=LF))
}


r2lEndTable <- function(out="latex", LF=TRUE) {
	str <- switch(out,
		"html" = "</TABLE>",
		"latex" = "\\end{tabular}"
	)
	return(r2lEol(str, LF=LF))
}


r2lBuildRow <- function(x, span=NULL, hline=TRUE, out="latex", LF=TRUE,border="|") {
    x <- as.character(x)
    len <- length(x)
    start <- switch(out,
                    "html" = "<TR ALIGN='center'>",
                    "latex" = ""
                    )
    rowTxt <- character()
    if (is.null(span)) {
        span = rep(1,len)
    }

    for (i in 1:length(span)) {
        if (span[i]>1) {
            rowTxt <- switch(out,
                             "html" = paste(rowTxt,"<TD COLSPAN=",span[i], ">", x[i], "</TD>", sep=""),
                             "latex" = paste(rowTxt,"\\multicolumn{",span[i],"}{",border,"c",border,"}{", x[i], "}",if(i<len){"&"}else{""}, sep="")
                             )
        }else{
            rowTxt <- switch(out,
                             "html" = paste(rowTxt, "<TD>", x[i], "</TD>", sep=""),
                             "latex" = paste(rowTxt, x[i], if(i<len){"&"}else{""})
                             )
        }
    }

    end <- switch(out,
                  "html" = "</TR>",
                  "latex" = if(hline){"\\\\\n\\hline"}else{"\\\\\n"}
                  )
    str <- paste(start, rowTxt, end, sep="")
    return(r2lEol(str, LF=LF))
}

r2lBuildColumnTitle <- function(x, span=NULL, hline=TRUE, out="latex", LF=TRUE,border="|") {
	x <- as.character(x)
	len <- length(x)
	for (i in 1:len) {
            x[i] <- r2lBold(x[i],out=out)
	}
	return(r2lBuildRow(x, span=span, hline=hline, out=out, LF=LF,border=border))
}

r2lPercentRow <- function(x, out="latex") {
	x <- round(x,digits=4)*100
	percent <- switch(out,
		"html" = paste("<small>&rarr;",x,"%</small>"),
		"latex" = paste("\\small$\\rightarrow$",x,"\\%\\normalsize")
	)
	return(percent)
}


### On entry, width is in pixels
###
###
r2lIncludeGraphics <- function(fileName, width, out="latex", LF=TRUE) {
	if (width == 0) {width <- 300}
	ppc <- r2lPixelsPerCm()
	str <- switch(out,
		"html" = paste("<IMG SRC='",fileName,"' WIDTH=", width,">", sep=""),
		"latex" = paste("\\parbox{",width/ppc,"cm}{\\includegraphics[width=",width/ppc,"cm]{",fileName,"}}", sep="")
	)
	return(r2lEol(str, LF=LF))
}


r2lBold <- function(str, out="latex") {
	txt <- switch(out,
		"html" = paste("<B>", str, "</B>", sep=""),
		"latex" = paste("\\textbf{", str, "}", sep=""))
	return(txt)
}

r2lUnderline <- function(str, out="latex") {
	txt <- switch(out,
		"html" = paste("<U>", str, "</U>", sep=""),
		"latex" = paste("\\underline{", str, "}", sep=""))
	return(txt)
}


r2lItal <- function(str, out="latex") {
	txt <- switch(out,
		"html" = paste("<I>", str, "</I>", sep=""),
		"latex" = paste("\\textit{", str, "}", sep=""))
	return(txt)
}


r2lStartBlock <- function(title, out="latex", hline=FALSE, LF=TRUE) {
	attrs <- switch(out,
		"html" = " align='center'",
		"latex" = "@{}c@{}"
	)
	str <- r2lStartTable(out=out, hline=hline, LF=LF, attrs=attrs)
	str <- paste(str, r2lBuildRow(c(r2lBold(title, out=out)), out=out, hline=hline, LF=LF))
	return(r2lEol(str, LF=LF))
}


r2lEndBlock <- function(out="latex", LF=TRUE) {
	return(r2lEndTable(out=out, LF=LF))
}


################# r2lEndStruct #################
### close r2lUnivBeginStruct and r2lBivBeginStruct

r2lEndStruct <- function(out="latex") {
	txt <- r2lEndTable(out=out)
	if (out == "latex") {
		txt <- paste(txt, "\\end{center}", sep="")
	}
	return(txt)
}




###############################################################
###                    Functions for Univ                   ###
###############################################################

################# r2lUnivBeginStruct #################

r2lUnivBeginStruct <- function(x,tabTitle,nbColumn=1,tabSpec="c", out="latex") {
	txt <- character()
	if (out == "latex") {
		txt <- paste(txt, "\t\\begin{center}\n\t\\addtolength{\\leftskip}{-4cm}\\addtolength{\\rightskip}{-4cm}\n\n")
	}

	attrs <- switch(out,
		"html"= " BORDER=1",
		"latex" = tabSpec#paste("|", paste(rep("c", nbColumn), sep="", collapse="|"), "|", sep="")
	)
	txt <- paste(txt, r2lStartTable(attrs=attrs,hline=TRUE,out=out), sep="")

        titleLine <- r2lBold(r2lUnderline(tabTitle,out=out),out=out)

	headerLine <- paste(r2lBivBannerStr(x),
            switch(out, "html"= " - ", "latex" = " \\hfill "),
            "N=", length(x), " ; NA=", sum(is.na(x)), " (", round(sum(is.na(x))/length(x)*100, digits=2),
	    switch(out, "html"= "%)", "latex" = "\\%)"),
            sep="")

	txt <- paste(txt,
                     r2lBuildRow(titleLine, span=nbColumn, hline=FALSE, out=out),
                     r2lBuildRow(headerLine, span=nbColumn, hline=TRUE, out=out),
                     sep="")
	return(txt)
}


#############################
### Print each Frequency, its size and frequency
### Used in: r2lNominal, r2lOrdinal, r2lDiscrete

r2lUnivFrequency <- function(x, out="latex") {
    attrs <- switch(out,
                    "html"= "",
                    "latex" = "@{}l@{ : }cl@{}"
                    )

    txt <- r2lStartTable(out=out, attrs=attrs,hline=FALSE)

    tableVar <- table(x)
    sumVar <- sum(tableVar)
    for (i in 1:nrow(tableVar)) {
        rowTxt <- r2lBuildRow( c(labels(tableVar[i]),
                              tableVar[i][[1]],
                              paste("(", round(tableVar[i][[1]]/sumVar*100, digits=2), "\\%)",sep="")
                              ), span=c(1,1,1),hline=FALSE,out=out)
        txt <- paste(txt, rowTxt)
    }
    txt <- paste(txt, r2lEndTable(out=out))
    return(txt)
}


r2lUnivSummary <- function(x, out="latex") {
    attrs <- switch(out,
        "html"= "",
        "latex" = "@{}l@{ : }c@{}"
    )

    txt <- r2lStartTable(out=out, attrs=attrs,hline=FALSE)

    if (class(x)[1]=="ordered") {
        quartile <- summary(as.integer(x))
        niveaux <- levels(x)
        txt <- paste(txt, r2lBuildRow( c("        Min.   ", niveaux[quartile[1]]),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Q1     ", niveaux[round(quartile[2])]),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Median ", niveaux[round(quartile[3])]),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Q3     ", niveaux[round(quartile[5])]),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Max.   ", niveaux[quartile[6]]),hline=FALSE, out=out))
    } else {
        quartile <- summary(x)
        txt <- paste(txt, r2lBuildRow( c("        Mean   ", round(mean(na.omit(x)), digits=3)),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Var.   ", round(var(na.omit(x)), digits=3)),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        SD     ", round(sd(na.omit(x)), digits=3)), hline=TRUE,out=out))
        txt <- paste(txt, r2lBuildRow( c("        Min.   ", round(quartile[1], digits=3)),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Q1     ", round(quartile[2], digits=3)),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Median ", round(quartile[3], digits=3)),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Q3     ", round(quartile[5], digits=3)),hline=FALSE, out=out))
        txt <- paste(txt, r2lBuildRow( c("        Max.   ", round(quartile[6], digits=3)),hline=FALSE, out=out))
    }

    txt <- paste(txt, r2lEndTable(out=out))
    return(txt)
}



###############################################################
###                     Functions for Biv                   ###
###############################################################

################# r2lBivBeginStruct #################

r2lBivBeginStruct <- function(y,x,tabTitle,nbColumn=1,tabSpec="c",out="latex") {
	txt <- character()
	if (out == "latex") {
		txt <- paste(txt, "\t\\begin{center}\n\t\\addtolength{\\leftskip}{-4cm}\\addtolength{\\rightskip}{-4cm}\n\n")
	}

	numNA <- sum(is.na(x) | is.na(y))
	total <- length(y)

	attrs <- switch(out,
		"html" = " BORDER=1",
		"latex" = tabSpec
	)
	txt <- paste(txt, r2lStartTable(attrs=attrs, out=out), sep="")

        titleLine <- r2lBold(r2lUnderline(tabTitle,out=out),out=out)

	formula <- paste(r2lBivBannerStr(y),
		switch(out, "html" = " ~ ", "latex" = " $\\sim$ "),
		r2lBivBannerStr(x))
	headerLine <- paste(formula,
		switch(out, "html" = " - ", "latex" = " \\hfill "),
		"N=", length(y), " ; NA=", numNA, " (", round(numNA/total*100, digits=2),
		switch(out, "html" = "%)", "latex" = "\\%)"),
		sep="")
	txt <- paste(txt,
                     r2lBuildRow(titleLine, span=nbColumn, hline=FALSE, out=out),
                     r2lBuildRow(headerLine, span=nbColumn, hline=TRUE, out=out),
                     sep="")
	return(txt)
}


r2lBivBannerStr <- function(y) {
	if (class(y)[1]=="discrete") {
		ylen <- length(levels(as.factor(y)))
		ytype <- "Discrete"
	} else if (class(y)[1]=="ordered") {
		ylen <- length(levels(y))
		ytype <- "Ordered"
	} else if (class(y)[1]=="factor" || class(y)[1]=="character" || class(y)[1]=="logical") {
		ylen <- length(levels(as.factor(y)))
		ytype <- "Nominal"
	} else if (class(y)[1]=="continuous" || class(y)[1]=="numeric") {
		ylen <- 0
		ytype <- "Continuous"
	} else {
		ylen <- 0
		ytype <- ""
	}

	str <- ytype
	if (ylen>0) {
		str <- paste(str,"(",ylen,")",sep="")
	}
	return(str)
}


################# r2lBivSummaryArray #################
### y is continuous, x is a factor or is continuous.
### Return an array with 10 rows corresponding to the statistics "N", "NA",
### "Mean", "Var", "SD", "Min", "Q1", "Median", "Q3", "Max"
### and as many columns as there are modalities in x.

r2lBivSummaryArray <- function(y, x) {
	lvl <- levels(x)
	ln <- length(lvl)
	arr <- array(dim=c(10,ln))
	for (i in 1:ln) {
		data <- y[x==lvl[i]]
		arr[1,i] <- length(data)
		arr[2,i] <- sum(is.na(data))
		arr[3,i] <- round(mean(na.omit(data)), digits=2)
		arr[4,i] <- round(var(na.omit(data)), digits=2)
		arr[5,i] <- round(sd(na.omit(data)), digits=2)
		smry <- summary(data)
		arr[6,i] <- round(smry[1][[1]], digits=2)
		arr[7,i] <- round(smry[2][[1]], digits=2)
		arr[8,i] <- round(smry[3][[1]], digits=2)
		arr[9,i] <- round(smry[5][[1]], digits=2)
		arr[10,i] <- round(smry[6][[1]], digits=2)
	}

	return(arr)
}


r2lBivSummary <- function(y,x,out="latex") {
    names <- c("N", "NA", "Mean", "Var", "SD", "Min", "Q1", "Median", "Q3", "Max")
    hline <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)

    if (!is.factor(x) || class(x)[1] == "numeric") {
        len <- length(y)
        y <- c(y,x)
        x <- factor(rep(c(1:2),c(len,len)), labels=c("y","x"))
    }
    arr <- r2lBivSummaryArray(y,x)
    nbCol <- dim(arr)[2]

    attrs <- switch(out,
        "html" = " BORDER=1 ALIGN='center'",
        "latex" = paste("@{}l@{ : }", paste(rep("c", nbCol), sep="", collapse="|"), "@{}", sep="")
    )

    tbl <- r2lStartTable(attrs=attrs,out=out,hline=FALSE)
    tbl <- paste(tbl, r2lBuildRow( r2lBold(c("Modal.", levels(x)), out=out), hline=TRUE, out=out))
    for (i in 1:length(names)) {
        tbl <- paste(tbl, r2lBuildRow( c(names[i], arr[i,]), hline=hline[i], out=out))
    }
    tbl <- paste(tbl, r2lEndTable(out=out))
    return(tbl)
}
#cat(r2lBivSummary(y,x1))



################# r2lBivQuartilesArray #################
### y is ordered, x is a factor.
### Return an array with 5 rows corresponding to the quartiles "Q0", "Q1", "Q2", "Q3", "Q4".
### and as many columns as there are modalities in x.

r2lBivQuartilesArray <- function(y, x) {
	x <- as.factor(x)
	# y is ordered, x is a factor
	lvlx <- levels(x)
	ln <- length(lvlx)
	arr <- array(dim=c(5,ln))

	lvly <- levels(y)
	for (i in 1:ln) {
		data <- y[x==lvlx[i]]
		smry <- summary(as.numeric(data))
		arr[1,i] <- lvly[smry[1]]
		arr[2,i] <- lvly[smry[2]]
		arr[3,i] <- lvly[smry[3]]
		arr[4,i] <- lvly[smry[5]]
		arr[5,i] <- lvly[smry[6]]
	}
	return(arr)
}


r2lBivQuartilesTable <- function(y,x,out="latex") {
    x <- as.factor(x)
    arr <- r2lBivQuartilesArray(y,x)
    nbLi <- dim(arr)[1]
    nbCol <- dim(arr)[2]

    attrs <- switch(out,
                    "html" = " BORDER=1 ALIGN='center'",
                    "latex" = paste("@{}l@{ : }", paste(rep("c", nbCol), sep="", collapse="|"), "@{}", sep=""))

    # Build the table
    tbl <- r2lStartTable(attrs=attrs,out=out,hline=FALSE)
    tbl <- paste(tbl, r2lBuildRow( r2lBold(c("Modal. ", levels(x)), out=out), hline=TRUE, out=out))
    for (i in 1:nbLi) {
        name <- paste("Q",i,sep="")
        tbl <- paste(tbl, r2lBuildRow( c(name, arr[i,]), hline=FALSE, out=out))
    }
    tbl <- paste(tbl, r2lEndTable(out=out))

	return(tbl)
}


################# r2lBivContingencyTable #################
### Display the contingency table between x (nominal(2)) and y (nominal,
### ordered, discrete).

r2lBivContingencyTable <- function(y,x,out="latex") {
    tt <- table(as.factor(x),as.factor(y))
    attrs <- switch(out,
        "html" = " BORDER=1 ALIGN='center'",
	"latex" = paste("l|",paste(rep("c|",ncol(tt)),collapse=""),"|c",sep="")
    )

    rowNames <- rownames(tt)
    colNames <- colnames(tt)
    att <- addmargins(tt)
    ptt <- prop.table(tt, margin=1)
    dim <- dim(tt)

	# Build the table
    tbl <- r2lStartTable(attrs=attrs,out=out, hline=FALSE)
    tbl <- paste(tbl, r2lBuildRow( c("", colNames, "Total"), out=out))
    for (i in 1:dim[1]) {
        tbl <- paste(tbl, r2lBuildRow( c(rowNames[i], att[i,]), out=out, hline=FALSE))
        tbl <- paste(tbl, r2lBuildRow( c("", r2lPercentRow(c(ptt[i,],1), out=out)), out=out))
    }

    tbl <- switch(out,
		"html" = tbl,
		"latex" = paste(tbl,"\\hline")
	)

    tbl <- paste(tbl, r2lBuildRow( c("Total", att[dim[1]+1,]), hline=FALSE,out=out))
    tbl <- paste(tbl, r2lBuildRow( c("", r2lPercentRow(att[dim[1]+1,]/sum(tt), out=out)), hline=FALSE,out=out))

    tbl <- paste(tbl, r2lEndTable(out=out))

	return(tbl)
}



###############################################################
###                           tests                         ###
###############################################################

r2lSignificance <- function(pvalue, out="latex") {
    str <- character()

    if(is.na(pvalue)){
        str="NA"
    }else{
        if(pvalue < 0.001){
            str <- "***"
        }else{
            if(pvalue < 0.01){
                str <- "**"
            }else{
                if(pvalue < 0.05){
                    str <- "*"
                }else{
                    if(pvalue < 0.1){
                        str <- "."
                    }else{
                        str <- " "
                    }
                }
            }
        }

        if(out=="latex"){str <- paste("$",str,"$",sep="")}else{}
    }
    return(str)
}


r2lPValueStr <- function(pvalue, out="latex") {
	str <- character()
	if (pvalue < 1e-6) {
		str <- paste(" < 1e",trunc(log(pvalue,base=10)),sep="")
	} else if (pvalue < 0.001) {
		str <- format(pvalue,scientific=TRUE,digits=4)
	} else {
		str <- as.character(round(pvalue,digits=4))
	}
	return(str)
}


r2lBivTestCorPearson <- function(y,x,out="latex") {
    x <- as.numeric(x)
    y <- as.numeric(y)

    str <- switch(out,
        "html" = "&rho;<sub>P</sub>",
	"latex" = "$\\rho_{P}$"
    )

    corP <- cor.test(y,x, na.action="na.omit", method="pearson")
    txt <- r2lBuildRow( c("Cor Pearson",
            paste(str, "=", round(corP$estimate,digits=4)),
            paste("p =", round(corP$p.value,digits=4)),
            r2lSignificance(corP$p.value,out=out)
        ),hline=FALSE,out=out)
    return(txt)
}

r2lBivTestCorSpearman <- function(y,x,out="latex") {
    x <- as.numeric(x)
    y <- as.numeric(y)

    str <- switch(out,
        "html" = "&rho;<sub>S</sub>",
	"latex" = "$\\rho_{S}$"
    )

    corS <- cor.test(y,x, na.action="na.omit", method="spearman",exact=FALSE)
    txt <- r2lBuildRow( c("Cor Spearman",
            paste(str, "=", round(corS$estimate,digits=4)),
            paste("p =", round(corS$p.value,digits=4))
        ), hline=FALSE,out=out
    )
    return(txt)
}



################# r2lBivTestStudent #################

r2lBivTestStudent <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.numeric(y)

    reg <- t.test(y~x, na.action="na.omit")
    txt <- r2lBuildRow( c("T test",
                          paste("T =", round(reg$statistic,digits=4)),
                          paste("p =", round(reg$p.value,digits=4)),
                          r2lSignificance(reg$p.value, out=out)
                          ),hline=FALSE, out=out
                       )
    return(txt)
}


################# r2lBivTestWilcoxon #################

r2lBivTestWilcoxon <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.numeric(y)

    reg <- wilcox.test(y~x, na.action="na.omit")
    txt <- r2lBuildRow( c("Wilcoxon",
                          paste("W =", round(reg$statistic,digits=4)),
                          paste("p =", round(reg$p.value,digits=4)),
                          r2lSignificance(reg$p.value, out=out)
                          ), hline=FALSE,out=out
                       )
    return(txt)
}


################# r2lBivTestAnova #################
### Anova (with Fisher)

r2lBivTestAnova <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.numeric(y)

    anv <- anova(lm(y~x, na.action="na.omit"))
    F <- anv$"F value"[1]
    pvalue <- anv$"Pr(>F)"[1]
    txt <- r2lBuildRow( c("ANOVA",
            paste("F =", round(F,digits=4)),
            paste("p =", round(pvalue,digits=4)),
            r2lSignificance(pvalue, out=out)
        ),hline=FALSE,out=out
    )
    return(txt)
}


################# r2lBivTestKruskalWallis #################

r2lBivTestKruskalWallis <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.numeric(y)

    kwt <- kruskal.test(y~x, na.action="na.omit")
    K <- kwt$statistic
    pvalue <- kwt$p.value
    txt <- r2lBuildRow( c("Kruskal-Wallis (y~x)",
                          paste("K =", round(K,digits=4)),
                          paste("p =", round(kwt$p.value,digits=4)),
                          r2lSignificance(kwt$p.value, out=out)), hline=FALSE,out=out)
    return(txt)
}


r2lBivTestKruskalWallisInv <- function(y,x,out="latex") {
    y <- as.factor(y)
    x <- as.numeric(x)

    kwt <- kruskal.test(x~y, na.action="na.omit")
    K <- kwt$statistic
    pvalue <- kwt$p.value
    txt <- r2lBuildRow( c("Kruskal-Wallis (x~y)",
                          paste("K =", round(K,digits=4)),
                          paste("p =", round(kwt$p.value,digits=4)),
                          r2lSignificance(kwt$p.value, out=out)), hline=FALSE,out=out)
    return(txt)
}


################# r2lBivTestKhi2 #################
### Chi-Square and Fisher's Exact Test

r2lBivTestKhi2 <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.factor(y)

    smry <- summary(table(x,y))
    str <- switch(out,
                  "html" = "&chi;<sup>2</sup>",
                  "latex" = "$\\chi^{2}$"
                  )

    pvalue <- smry$p.value
    txt <- r2lBuildRow( c(paste(str," test"),
                paste(str, "=", round(smry$statistic,digits=4)),
#		paste("dl =", smry$parameter),
		paste("p =", round(pvalue,digits=4)),
		r2lSignificance(pvalue, out=out)), hline=FALSE,out=out)
    return(txt)
}

r2lBivTestFisherExact <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.factor(y)

    ft <- fisher.test(y,x,simulate.p.value=TRUE)
    pvalue <- ft$p.value
    txt <- r2lBuildRow(c("Fisher's Exact Test","",
                         paste("p =", round(pvalue,digits=4)),
                         r2lSignificance(pvalue, out=out)), hline=FALSE,out=out)
    return(txt)
}


################# r2lBivTestOddRR #################
### Computation of Odds ratio and Relative Risk (RR).
#         \begin{tabular}{lllll}
#             Odds Ratio & ODD = 0.34 & & $p < 0.001$ & ** \\
#             Relative Risk & RR = 0.35 & & $p < 0.002$ & ** \\
#         \end{tabular}
# 		   | x1 |  x0
# 		_____________
# 		y1 | a  |  b |
# 		_____________
# 		y0 | c  |  d |
# 		_____________
# The odds ratio is defined as
#    OR = ad/bc
#    Standard error for the log of the odd ratio:
#    LSE = sqrt(1/a + 1/b + 1/c + 1/d)
# The relative risk is defined as
#         a/(a+b)
#    RR = _______
#         c/(c+d)
#    Standard error for the log of the relative risk:
#    LSE = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))

r2lBivTestOddsRatio <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.factor(y)

    tt <- table(x,y)
	# Odds Ratio
    odd <- (tt[1,1]*tt[2,2]) / (tt[1,2]*tt[2,1])
	# Compute the standard error LSE for the log of the odd ratio
    lse <- sqrt(sum(1/tt))
	# log(ODD)/LSE has a normal distribution N(0,1).
	# The p-value is 2*P( Z < -|log(ODD)|/LSE ).
    pvalue <- 2*pnorm(-abs(log(odd))/lse)
    txt <- r2lBuildRow( c("Odds Ratio",
                          paste("ODD =", round(odd,digits=4)),
                          paste("p =",r2lPValueStr(pvalue)),
                          r2lSignificance(pvalue, out=out)), hline=FALSE, out=out)
    return(txt)
}

r2lBivTestRelativeRisk <- function(y,x,out="latex") {
    x <- as.factor(x)
    y <- as.factor(y)

    tt <- table(x,y)
	# Relative Risk
    rr <- (tt[1,1]*(tt[2,1]+tt[2,2])) / (tt[2,1]*(tt[1,1]+tt[1,2]))
	# Compute the standard error LSE for the log of the relative risk
    lse <- sqrt(1/tt[1,1] - 1/(tt[1,1]+tt[1,2]) + 1/tt[2,1] - 1/(tt[2,1]+tt[2,2]))
	# log(RR)/LSE has a normal distribution N(0,1).
	# The p-value is 2*P( Z < -|log(RR)|/LSE ).
    pvalue <- 2*pnorm(-abs(log(rr))/lse)
    txt <- r2lBuildRow( c("Relative Risk",
                          paste("RR =", round(rr,digits=4)),
                          paste("p =",r2lPValueStr(pvalue)),
                          r2lSignificance(pvalue, out=out)), hline=FALSE,out=out)
    return(txt)
}


r2lBivTest <- function(y,x,test,line,out="latex"){
    attrs <- switch(out,
        "html" = " align='center' cellpadding=5",
        "latex" = "llll"
	)

    #    txt <- r2lStartBlock("Tests", out=out)
    horizLine <- switch(out,
                        "html" = paste("<TR ALIGN='center'><TD COLSPAN=4><HR></TD></TR>", sep=""),
                        "latex" = "\\hline"
                        )

	# Build the table
    tbl <- r2lStartTable(attrs=attrs,out=out,hline=FALSE)
    nbLine <- 1
    for (i in test){
        tbl <- switch(out,
                      "html" = paste(tbl,"", sep=""),
                      "latex" = paste(tbl,"\n", sep="")
        )

        tbl <- paste(tbl, horizLine[line[nbLine]],do.call(paste("r2lBivTest",i,sep=""),list(y,x,out=out)))
        nbLine <- nbLine + 1
    }

    tbl <- paste(tbl,horizLine[line[nbLine]],r2lEndTable(out=out))
	# Put the table in the block
	#txt <- paste(txt, r2lBuildRow(tbl, out=out))
	#txt <- paste(txt, r2lEndBlock(out=out))

    return(tbl)
}



###############################################################
###                          graphs                         ###
###############################################################

### r2lDensities is a basic function (like boxplot, barplot,...)
r2lDensities <- function(y,x,...) {
    if (is.factor(x) || class(x)[1] == "discrete") {
        myColor <- redToBlue(length(table(x)))
    # This is the case continuous ~ (nominal|ordered|discrete)
    # Clean up the NAs
        df <- na.omit(data.frame(y,x))
        y <- df$y
        x <- df$x

        if (!is.factor(x)){
            x <- as.factor(x)
        }

        lv <- levels(x)
        len <- length(lv)
#       colors <- rev(rainbow(len,end=0.17))

        dd <- list()
        # First pass to calcul the densities
        for (i in 1:len) {
            sub <- y[x==lv[i]]
            if(length(sub)>1){
                dd <- c(dd,list(density(sub)))
            }else{}
        }

        # Second pass to calcul the top of ylim
        top <- 0
        for (i in 1:length(dd)) {
            if (top < max(dd[[i]]$y)) {
                top <- max(dd[[i]]$y)
            }else{}
        }

        # Third pass to draw the densities
        plot(dd[[1]],ylim=c(0,top), col=myColor[1],type="l",lwd=6,...)
        for (i in 2:length(dd)) {
            lines(dd[[i]], col=myColor[i],lwd=6)
        }
    } else {
        # This is the case continuous ~ continuous
        # Clean up the NAs
        x <- na.omit(x)
        y <- na.omit(y)

        # Calculate the densities
        ddx <- density(x)
        ddy <- density(y)
        top <- max(ddx$y,ddy$y)

        # Transmit the arguments to plot
        plot(ddx, ylim=c(0,top), col="red",...)
        lines(ddy, ylim=c(0,top), col="orange",...)
    }
}


################# r2lMakePlot #################
### Compute a plot, save it in a file on disk and return the proper string
### to include the graphic in the LaTeX or Html document. The argument
### 'kind' can be: "boxplot", "barplot", "hist", "mosaicplot", "qqplot",
### "plot". If x is NULL, y can be either a vector or a formula. If x is a
### vector, then y must be a vector too. Used by the r2lGraphXXX functions.

r2lMakePlot <- function(kind, arguments, graphDir, graphName, type, out="latex", wd = 0) {
	# Close the device
	on.exit(dev.off())

	if (type == "postscript") {
		ext <- "eps"
	} else {
		ext <- tolower(type)
	}
	# Build the full path
	plotName <- paste(graphName, "-", kind, sep="")
	if (graphDir=="") {
		plotDir <- getwd()
	} else {
		plotDir <- paste(getwd(), graphDir, sep="/")
		if (!file.exists(plotDir)) {
			dir.create(plotDir)
		} else {}
	}
	plotName <- paste(plotName, ext, sep=".")
	plotPath <- paste(plotDir, plotName, sep="/")

	# Open the device
	devArg <- list(plotPath)
	if (type == "postscript") {
		devArg <- c(devArg, horizontal=FALSE, onefile=FALSE)
	}
	if (wd > 0) {
		if (type == "postscript") {
			# The postscript device expects width and height in inches.
			# The default device size is 7 inches square.
			ppc <- r2lPixelsPerCm()
			plotwd <- (wd*7)/(ppc*2.54)
			# Hmm, setting the width does not seem to work...
			# (temporarily back to the default)
			devArg <- c(devArg, height=7, width=7)
		} else {
			plotwd <- wd
			devArg <- c(devArg, width=plotwd)
		}
	}
	do.call(type, devArg)

	# Execute the plotting command
	do.call(kind, arguments)

	# Build the string to include the graphics
	nch <- nchar(paste(getwd(),"/",sep=""))
	txt <- r2lIncludeGraphics( substring(plotPath, nch+1), width=wd, out=out)
	return(txt)
}


### Among the r2lGraphXXX functions below, r2lGraphQQPlot and r2lGraphHist
### expect one var argument.

r2lGraphQQPlot <- function(x,graphDir,graphName,type,out="latex") {
    ## To exchange axes -> datax=T
    arguments <- list(y=x,main="", xlab="", ylab="")
    ppc <- r2lPixelsPerCm()
    txt <- r2lMakePlot(kind="qqnorm", wd=2*ppc, arguments, graphDir=graphDir, graphName=graphName, type=type, out=out)
    return(txt)
}


r2lGraphHist <- function(x,graphDir,graphName,type,out="latex") {
	arguments <- list(x=x, main="", xlab="", ylab="", col=redToBlue(length(hist(x,plot=FALSE)$counts)))
	txt <- r2lMakePlot(kind="hist", arguments, graphDir=graphDir, graphName=graphName, type=type, out=out)
	return(txt)
}


### The others have two var arguments y and x. If
### x is missing, the function was called from the univariate case.
### They are all implemented via r2lMakePlot.

r2lGraphBoxplot <- function(y,x,graphDir,graphName,type,out="latex") {
    arguments <- list(main="", xlab="", ylab="")
    if(missing(x)){
        num <- 2
        arguments <- list(y,main="", xlab="", ylab="",col="orange")
    }else{
        num <- length(levels(as.factor(x)))
        arguments <- list(y~x,main="", xlab="", ylab="",col=redToBlue(num))
    }
    ppc <- r2lPixelsPerCm()
    txt <- r2lMakePlot(kind="boxplot", wd=min(2+num/2,9)*ppc, arguments=arguments, graphDir=graphDir, graphName=graphName, type=type, out=out)
    return(txt)
}


r2lGraphMosaicPlot <- function(y,x,graphDir,graphName,type,out="latex") {
	arguments <- list(x~y,main="", xlab="", ylab="",dir=c("h","v"), col=redToBlue(length(levels(as.factor(y)))))
	txt <- r2lMakePlot(kind="mosaicplot", arguments=arguments, graphDir=graphDir, graphName=graphName, type=type, out=out)
	return(txt)
}


r2lGraphScatterPlot <- function(y,x,graphDir,graphName,type,out="latex") {
	arguments <- list(y,x,main="", xlab="", ylab="")
	# jitter(y) ?
	txt <- r2lMakePlot(kind="plot", arguments=arguments, graphDir=graphDir, graphName=graphName, type=type, out=out)
	return(txt)
}


r2lGraphBarplot <- function(y,x,graphDir,graphName,type,out="latex") {
    if(missing(x)){
        if (class(y)[1]=="discrete") {
            tabVar <- table(c(y,range(na.omit(y))[1]:range(na.omit(y))[2]))-1
        } else {
            tabVar <- table(y)
        }
        numCol <- num <- dim(tabVar)

    }else{
        tabVar <- table(y, x)
        num <- prod(dim(tabVar))
        numCol <- dim(tabVar)[1]
    }
    ppc <- r2lPixelsPerCm()
    arguments <- list(tabVar, main="", xlab="", ylab="", beside=TRUE, col=redToBlue(numCol))
    txt <- r2lMakePlot(kind="barplot", wd=min(2+num/2,9)*ppc, arguments=arguments, graphDir=graphDir, graphName=graphName, type=type, out=out)
    return(txt)
}


r2lGraphDensity <- function(y,x,graphDir,graphName,type,out="latex") {
    arguments <- list(y=y,x=x,main="", xlab="", ylab="")
    ppc <- r2lPixelsPerCm()
    txt <- r2lMakePlot(kind="r2lDensities", wd=3*ppc,arguments=arguments, graphDir=graphDir, graphName=graphName, type=type, out=out)
    return(txt)
}



cat("   ---------------------------------------------------------
   ---                       Fin Functions               ---
   ---------------------------------------------------------\n")
