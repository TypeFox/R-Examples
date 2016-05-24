
#Load required packages. (If R commander is started with an soption of
#R_DEFAULT_PACKAGES="Rcmdr", other packages will not be not loaded.

#library(methods, quietly=TRUE)
#library(datasets, quietly=TRUE)
#library(multcomp, quietly=TRUE)
#library(mvtnorm, quietly=TRUE)
#library(abind, quietly=TRUE)
#library(aplpack, quietly=TRUE)
#library(foreign, quietly=TRUE)
#library(survival, quietly=TRUE)
#library(cmprsk, quietly=TRUE)
#library(aod, quietly=TRUE)

require("datasets")
#requireNamespace("car") 
require("methods")
window.type <- "width=7, height=7"
par.option <- 'lwd=1, las=1, family="sans", cex=1, mgp=c(3.0,1,0)'
#The first parameter of mpg defines the distance between axis and axis labels.
#This number is changed to 2.5 only for survival plots to avoid overlapping 
#with "Number at risk".

par.lwd <- "lwd=1"
par.cex <- "1"

#assign("window.type", "width=7, height=7", envir=.GlobalEnv)
#assign("par.option", 'lwd=1, las=1, family="sans", cex=1', envir=.GlobalEnv)

currentFields <- NULL	#A variable to send diaglog memory to Formula
#currentFields$subset <- ""
#Rcmdr=list(dialog.memory=TRUE)

#cat("\n")
#cat(gettext(domain="R-RcmdrPlugin.EZR","EZR comes with ABSOLUTELY NO WARRANTY just like R itself.", "\n"))
#cat(gettext(domain="R-RcmdrPlugin.EZR","Conditions for redistribution are also the same with R and R commander.", "\n"))
#cat(gettext(domain="R-RcmdrPlugin.EZR","Changes made from the original R commander include", "\n"))
#cat(gettext(domain="R-RcmdrPlugin.EZR","1. Replacing Rcmdr-menus.txt in //Rcmdr//etc with a file of the same name for EZR (menu file of R commander).", "\n"))
#cat(gettext(domain="R-RcmdrPlugin.EZR","2. Adding EZR.R, the main script of EZR written by Y.Kanda to //Rcmdr//etc.", "\n"))
#cat(gettext(domain="R-RcmdrPlugin.EZR","3. Replacing R-Rcmdr.mo in //Rcmdr//po//ja//LC_MESSAGES with a file of the same name for EZR (for translation in EZR).", "\n"))
#cat(gettext(domain="R-RcmdrPlugin.EZR","4. Replacing R-Rcmdr.po in //Rcmdr//po//ja//LC_MESSAGES with a file of the same name for EZR (for translation in EZR).", "\n"))
#cat(gettext(domain="R-RcmdrPlugin.EZR","5. Minimally modifying Commander.R in Rcmdr package.", "\n"))
cat("\n")
cat("-----------------------------------\n")
cat(gettext(domain="R-RcmdrPlugin.EZR","Starting EZR...", "\n"))
cat("   Version 1.32", "\n")
cat(gettext(domain="R-RcmdrPlugin.EZR","Use the R commander window.", "\n"))
cat("-----------------------------------\n")
cat("\n")


# for assignments to the global environment, from Rcmdr_1.9-3

#gassign <- function(x, value){
#    if (!(is.valid.name(x))) stop("argument x not a valid R name")
#    G <- .GlobalEnv
#    assign(x, value, envir=G)
#}


ifelse2 <- function (test, yes, no) 	#Treat the condition of NA as FALSE.
{
    storage.mode(test) <- "logical"
#    if(is.factor(yes)) yes <- as.character(yes)
#    if(is.factor(no)) yes <- as.character(no)
    	nas <- is.na(test)

	test[nas] <- FALSE
    	ans <- test
      ans[test] <- rep(yes, length.out = length(ans))[test]
      ans[!test] <- rep(no, length.out = length(ans))[!test]
    	ans
}

###modified from hist(), add one group below the lowest group, change default from "Sturges" to "Scott"
hist2 <- function(x, breaks="scott", plot=TRUE, ...){
	res <- hist(x, plot=F, breaks=breaks)
	if(res$breaks[1]==min(x, na.rm=TRUE)){
		breaks <- c(res$breaks[1]*2-res$breaks[2], res$breaks) #add a group below the lowest group
	}
	hist(x, breaks=breaks, plot=plot, ...)
}

###modified from original Hist() to use hist2() instead of hist()
HistEZR <- function (x, scale = c("frequency", "percent", "density"), xlab = deparse(substitute(x)), 
    ylab = scale, main = "", ...) 
{
    xlab
    x <- na.omit(x)
    scale <- match.arg(scale)
    if (scale == "frequency") 
        hist2(x, xlab = xlab, ylab = ylab, main = main, ...)
    else if (scale == "density") 
        hist2(x, freq = FALSE, xlab = xlab, ylab = ylab, main = main, 
            ...)
    else {
        n <- length(x)
        hist2(x, axes = FALSE, xlab = xlab, ylab = ylab, main = main, 
            ...)
        axis(1)
        max <- ceiling(10 * par("usr")[4]/n)
        at <- if (max <= 3) 
            (0:(2 * max))/20
        else (0:max)/10
        axis(2, at = at * n, labels = at * 100)
    }
    box()
    abline(h = 0)
    invisible(NULL)
}


nchar.ZenToHan <- function(x) {
	if(length(x)==1){
		return(length(charToRaw(x)))
	} else {
		x2 <- NULL	
		for(i in 1:length(x)){
			x2[i] <- length(charToRaw(x[i]))
		}		
		return(x2)
	}
}


###Print dataframe with ruled lines.
dataframe_print <- function(x, printrow=1) {

	row.number <- length(x[,1])
	col.number <- length(colnames(x))

	group.name.max.nchar <- max(nchar.ZenToHan(colnames(x)[1:col.number]))
	group.data.max.nchar <- 0
	for (i in 1:(col.number)){
		if(max(nchar.ZenToHan(as.character(x[,i]))) > group.data.max.nchar){
			group.data.max.nchar <- max(nchar.ZenToHan(as.character(x[,i])))
		}
	}
	group.nchar <- max(group.name.max.nchar, group.data.max.nchar)
	for (i in 1:(col.number)){
		margin <- group.nchar - nchar.ZenToHan(colnames(x)[i])
		colnames(x)[i] <- paste(paste(rep(" ", floor(margin/2)), collapse=""), colnames(x)[i], paste(rep(" ", ceiling(margin/2)), collapse=""), sep="")
	}

	if(printrow==1){
		rownames.nchar <- max(nchar.ZenToHan(row.names(x)))
		line.nchar <- rownames.nchar + group.nchar * col.number
		line.nchar <- line.nchar + 3 * (col.number - 1)
		table.line <- NULL
		table.line[1] <- paste(rep("-", line.nchar), collapse="")
		table.line.1 <- paste(rep(" ", rownames.nchar), collapse="")
		table.line.2 <- paste(colnames(x), collapse=" | ")
		table.line[2] <- paste(table.line.1, table.line.2, sep=" | ")
		table.line[3] <- table.line[1]
		substring(table.line[3], rownames.nchar + 2) <- "+"
		for(i in 1:col.number - 1){
			substring(table.line[3], rownames.nchar + 3 + (group.nchar +3 ) * (i - 1) + group.nchar + 2) <- "+"	
		}
		for(i in 1:row.number){
			table.line[3+i] <- paste(rep(" ", rownames.nchar - nchar.ZenToHan(row.names(x)[i])), collapse="")
			table.line[3+i] <- paste(row.names(x)[i], table.line[3+i], sep="")
			for(j in 1:col.number){
				cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[i,j]))), collapse="")
				cell <- paste(cell, x[i,j], sep="")			
				table.line[3+i] <- paste(table.line[3+i], " | ", cell, sep="")
			}
		}
		table.line[4+row.number] <- table.line[1]
	} else {
		line.nchar <- group.nchar * col.number
		line.nchar <- line.nchar + 3 * (col.number - 1)
		table.line <- NULL
		table.line[1] <- paste(rep("-", line.nchar), collapse="")
		table.line[2] <- paste(colnames(x), collapse=" | ")
		table.line[3] <- table.line[1]
		for(i in 1:(col.number-1)){
			substring(table.line[3], (group.nchar +3 ) * (i - 1) + group.nchar + 2) <- "+"	
		}
		for(i in 1:row.number){
 			cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[i,1]))), collapse="")
			table.line[3+i] <- paste(cell, x[i,1], sep="")			
			for(j in 2:col.number){
				cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[i,j]))), collapse="")
				cell <- paste(cell, x[i,j], sep="")			
				table.line[3+i] <- paste(table.line[3+i], " | ", cell, sep="")
			}
		}
		table.line[4+row.number] <- table.line[1]
	}
	cat(table.line, sep="\n")
}


###Print twoway dataframe with ruled lines.
twoway_dataframe_print <- function(x) {

	group.name <- x[1,3]
	row.number <- length(row.names(x))
	col.number <- length(colnames(x))
	x <- x[3:row.number,]
	row.number <- row.number-2

	group.name.max.nchar <- max(nchar.ZenToHan(colnames(x)[1:col.number]))
	group.name.max.nchar <- max(group.name.max.nchar, nchar.ZenToHan(group.name))
	group.data.max.nchar <- 0
	for (i in 1:(col.number)){
		if(max(nchar.ZenToHan(as.character(x[,i]))) > group.data.max.nchar){
			group.data.max.nchar <- max(nchar.ZenToHan(as.character(x[,i])))
		}
	}
	group.nchar <- max(group.name.max.nchar, group.data.max.nchar)
	for (i in 1:(col.number)){
		margin <- group.nchar - nchar.ZenToHan(colnames(x)[i])
		colnames(x)[i] <- paste(paste(rep(" ", floor(margin/2)), collapse=""), colnames(x)[i], paste(rep(" ", ceiling(margin/2)), collapse=""), sep="")
	}
	margin <- group.nchar - nchar.ZenToHan(group.name)
	group.name <- paste(paste(rep(" ", floor(margin/2)), collapse=""), group.name, paste(rep(" ", ceiling(margin/2)), collapse=""), sep="")

	rownames.nchar <- max(nchar.ZenToHan(row.names(x)))

	line.nchar <- group.nchar * col.number
	line.nchar <- line.nchar + 3 * (col.number - 1)
	table.line <- NULL
	table.line[1] <- paste(rep("-", line.nchar), collapse="")
	dummy.colname <- paste(rep(" ", group.nchar), collapse="")
	table.line[2] <- paste(paste(rep(dummy.colname, 2), collapse=" | "), paste(group.name, paste(rep(dummy.colname, col.number-4), collapse="   "), sep="   "), dummy.colname, sep=" | ")
	table.line[3] <- paste(colnames(x)[1], colnames(x)[2], paste(colnames(x)[3:(col.number-1)], collapse="   "), colnames(x)[col.number], sep=" | ")
	table.line[4] <- table.line[1]
	for(i in c(1, 2, col.number-1)){
		substring(table.line[4], (group.nchar +3 ) * (i - 1) + group.nchar + 2) <- "+"	
	}
	for(i in 1:row.number){
		cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[i,1]))), collapse="")
		table.line[4+i] <- paste(cell, x[i,1], sep="")			
		for(j in 2:col.number){
			cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[i,j]))), collapse="")
			cell <- paste(cell, x[i,j], sep="")			
			table.line[4+i] <- paste(table.line[4+i], " | ", cell, sep="")
		}
	}
	table.line[5+row.number] <- table.line[1]

	cat(table.line, sep="\n")
}


finaltable_dataframe_print <- function(x) {

	grouping=1
	if(x[1,1]!=""){	#No grouping
		grouping=0
		flag <- x[1,2]
		group.name <- "Overall"
	} else {		#Summary with grouping
		flag <- x[2,2]
		if(flag=="Group"){	#At least one categorical variable
			group.name <- x[1,3]
		} else {			#Only continuoue variables
			group.name <- x[1,2]
		}
	}

	row.number <- length(row.names(x))
	col.number <- length(colnames(x))
	if (grouping==1){
		x <- x[3:row.number,]
		row.number <- row.number-2
	}

	group.name.max.nchar <- max(nchar.ZenToHan(colnames(x)[1:col.number]))
	group.name.max.nchar <- max(group.name.max.nchar, nchar.ZenToHan(group.name))
	group.data.max.nchar <- 0
	for (i in 1:(col.number)){
		if(max(nchar.ZenToHan(as.character(x[,i]))) > group.data.max.nchar){
			group.data.max.nchar <- max(nchar.ZenToHan(as.character(x[,i])))
		}
	}
	group.nchar <- max(group.name.max.nchar, group.data.max.nchar)
	for (i in 1:(col.number)){
		margin <- group.nchar - nchar.ZenToHan(colnames(x)[i])
		colnames(x)[i] <- paste(paste(rep(" ", floor(margin/2)), collapse=""), colnames(x)[i], paste(rep(" ", ceiling(margin/2)), collapse=""), sep="")
	}
	margin <- group.nchar - nchar.ZenToHan(group.name)
	group.name <- paste(paste(rep(" ", floor(margin/2)), collapse=""), group.name, paste(rep(" ", ceiling(margin/2)), collapse=""), sep="")

	rownames.nchar <- max(nchar.ZenToHan(row.names(x)))

	line.nchar <- group.nchar * col.number
	line.nchar <- line.nchar + 3 * (col.number - 1)
	table.line <- NULL
	table.line[1] <- paste(rep("-", line.nchar), collapse="")
	dummy.colname <- paste(rep(" ", group.nchar), collapse="")
	if (grouping==0){
		table.line[2] <- ""
		table.line[3] <- ""
		table.line[4] <- ""
	} else {
		if(flag=="Group"){	#At least one categorical variable
			table.line[2] <- paste(paste(rep(dummy.colname, 2), collapse=" | "), paste(group.name, paste(rep(dummy.colname, col.number-4), collapse="   "), sep="   "), dummy.colname, sep=" | ")
			table.line[3] <- paste(colnames(x)[1], colnames(x)[2], paste(colnames(x)[3:(col.number-1)], collapse="   "), colnames(x)[col.number], sep=" | ")
		} else {			#Only continuoue variables
			table.line[2] <- paste(dummy.colname, paste(group.name, paste(rep(dummy.colname, col.number-3), collapse="   "), sep="   "), dummy.colname, sep=" | ")
			table.line[3] <- paste(colnames(x)[1], paste(colnames(x)[2:(col.number-1)], collapse="   "), colnames(x)[col.number], sep=" | ")
		}
		table.line[4] <- table.line[1]
		for(i in c(1, 2, col.number-1)){
			substring(table.line[4], (group.nchar +3 ) * (i - 1) + group.nchar + 2) <- "+"	
		}
	}
	cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[1,1]))), collapse="")
	table.line[5] <- paste(cell, x[1,1], sep="")			
	for(j in 2:col.number){
			margin <- group.nchar - nchar.ZenToHan(as.character(x[1,j]))
			cell1 <- paste(rep(" ", floor(margin/2)), collapse="")
			cell2 <- paste(rep(" ", ceiling(margin/2)), collapse="")
			cell <- paste(cell1, x[1,j], cell2, sep="")			
			table.line[5] <- paste(table.line[5], " | ", cell, sep="")
	}
	cell <- paste(rep("-", group.nchar), collapse="")
	table.line[6] <- paste(rep(cell, col.number), collapse="-+-")
#	table.line[6] <- paste(table.line[6], paste(rep("-", 3 * (col.number)), collapse=""), sep="")
	for(i in 2:row.number){
		cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[i,1]))), collapse="")
		table.line[5+i] <- paste(cell, x[i,1], sep="")			
		for(j in 2:col.number){
			cell <- paste(rep(" ", group.nchar - nchar.ZenToHan(as.character(x[i,j]))), collapse="")
			cell <- paste(cell, x[i,j], sep="")			
			table.line[5+i] <- paste(table.line[5+i], " | ", cell, sep="")
		}
	}
	table.line[6+row.number] <- table.line[1]
	if(grouping==0) table.line <- c(table.line[1], table.line[5:7], table.line[6], table.line[8:length(table.line)])
	cat(table.line, sep="\n")
}


###Output the results of multivariate analysis to clipboard and files.
w.multi <- function (table = cox.table, filename = "clipboard", CI = 0, signif = 0, en = 1) {
	#Jan 2016 modified to work correctly when the number of independent covariate is only one

    table[, 4] <- as.numeric(table[, 4])
    if (signif > 0) {
        table <- signif(table, digits = signif)
    }
    if (signif == 0) {
        table[, 1] <- floor(table[, 1] * 100 + 0.5)/100
        table[, 2] <- floor(table[, 2] * 100 + 0.5)/100
        table[, 3] <- floor(table[, 3] * 100 + 0.5)/100
        table[, 1] <- formatC(table[, 1], format = "f", digits = 2)
        table[, 2] <- formatC(table[, 2], format = "f", digits = 2)
        table[, 3] <- formatC(table[, 3], format = "f", digits = 2)
        table[, 4] <- signif(as.numeric(table[, 4]), digits = 2)
        table[, 4] <- formatC(as.numeric(table[, 4]), format = "fg")
    }

    if(length(rownames(table))==1){
	  table <- rbind(table, c(" ", " ", " ", " "))
    }

    table2 <- table
    if (CI == 0) {
        table2[, 1] <- paste(table[, 1], " (", table[, 2], "-", 
            table[, 3], ")", sep = "")
        table2[, 2] <- table[, 4]
        table2 <- table2[, 1:2]
    }

    if(table2[2,1]=="  ( - )") table2[2,1] <- " "
    table2 <- cbind(rownames(table), table2)
    colnames(table2)[1] <- ifelse(en == 1, "Factor", gettext(domain="R-RcmdrPlugin.EZR", 
        "Factor"))

    rownames(table2) <- NULL
    if (en == 1 & colnames(table2)[2] == gettext(domain="R-RcmdrPlugin.EZR", 
        "Hazard ratio")) 
        colnames(table2)[2] <- "Hazard ratio"
    if (en == 1 & colnames(table2)[2] == gettext(domain="R-RcmdrPlugin.EZR", 
        "odds ratio")) 
        colnames(table2)[2] <- "Odds ratio"
    if (en == 1 & CI == 1) 
        colnames(table2)[3:4] <- c("Lower 95%CI", "Upper 95%CI")
    if (CI == 0) 
        colnames(table2)[3] <- ifelse(en == 1, "p.value", gettext(domain="R-RcmdrPlugin.EZR", 
            "p.value"))
    if (CI == 1) 
        colnames(table2)[5] <- ifelse(en == 1, "p.value", gettext(domain="R-RcmdrPlugin.EZR", 
            "p.value"))
     #print(data.frame(table2), quote=FALSE, row.names=FALSE, col.names=TRUE)
	dataframe_print(table2, printrow=0)

#    print(table2)
#    print(paste("Write to ", filename, sep = ""))
    if (filename == "clipboard") {
        write.table(data.frame(table2), "clipboard", sep = "\t", 
            row.names = FALSE)
    }
    else {
        write.csv(data.frame(table2), file = as.character(filename), row.names = FALSE)
    }
}


w.multireg <- function (table = multireg.table, filename = "clipboard", CI = 0, signif = 0, en = 1) {
    if (signif > 0) {
        table <- signif(table, digits = signif)
    }
    if (signif == 0) {
        table[, 1] <- floor(table[, 1] * 100 + 0.5)/100
        table[, 2] <- floor(table[, 2] * 100 + 0.5)/100
        table[, 3] <- floor(table[, 3] * 100 + 0.5)/100
        table[, 1] <- formatC(table[, 1], format = "f", digits = 2)
        table[, 2] <- formatC(table[, 2], format = "f", digits = 2)
        table[, 3] <- formatC(table[, 3], format = "f", digits = 2)
        table[, 4] <- signif(as.numeric(table[, 4]), digits = 2)
        table[, 4] <- formatC(as.numeric(table[, 4]), format = "fg")
    }
    table <- cbind(rownames(table), table)
	if(en==1){
		colnames(table) <- c("Factor", "Estimate", "Std. Error", "t value", "p.value")
	} else {
		colnames(table) <- gettext(domain="R-RcmdrPlugin.EZR", c("Factor", "Estimate", "Std. Error", "t value", "p.value"))
	}
    rownames(table) <- NULL
#	print(data.frame(table), quote=FALSE, row.names=FALSE, col.names=TRUE)
	dataframe_print(table, printrow=0)
#    print(paste("Write to ", filename, sep = ""))
    if (filename == "clipboard") {
        write.table(data.frame(table), "clipboard", sep = "\t", 
            row.names = FALSE)
    }
    else {
        write.csv(data.frame(table), file = as.character(filename), row.names = FALSE)
    }
}


###Output two-way table to clipboard and files.
w.twoway <- function (table = Fisher.summary.table, filename = "clipboard", en = 1) {
    table <- as.matrix(table)

    rows <- length(table[, 1])
    columns <- length(table)
    Factor <- substring(row.names(table), 1, regexpr("=", row.names(table)) - 
        1)
    Group <- substring(row.names(table), regexpr("=", row.names(table)) + 
        1)
    for (i in 1:(rows - 1)) {
        j <- 1
        while (Factor[i] == Factor[i + j]) {
            Factor[i + j] <- ""
            j <- j + 1
            if ((i + j) > rows) 
                break
        }
    }
    StratifyFactor <- substring(colnames(table), 1, regexpr("=", colnames(table)) - 
        1)
    StratifyGroup <- substring(colnames(table), regexpr("=", colnames(table)) + 
        1)
    colnames(table) <- StratifyGroup
    table <- cbind(Factor, Group, table)
#    rownames(table) <- NULL
    colnames(table)[length(colnames(table))] <- ifelse(en == 
        1, "p.value", gettextRcmdr("p.value"))
    if (en == 0) colnames(table) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(table))
#    print(data.frame(table), quote=F, row.names = FALSE, col.names = FALSE)
    ncol <- length(colnames(table))
    row1 <- colnames(table)
    row1 <- matrix(row1, ncol=ncol)
    table <- rbind(row1, table)
    row0 <- rep(" ", ncol)
    row0[3] <- StratifyFactor[1]
    table <- rbind(row0, table)
	twoway_dataframe_print(table)
	#    print(paste("Write to ", filename, sep = ""))
    if (filename == "clipboard") {
        write.table(data.frame(table), "clipboard", sep = "\t", 
            row.names = FALSE, col.names = FALSE)
    }
    else {
        write.table(data.frame(table), file = as.character(filename), sep = ",", row.names = FALSE, col.names = FALSE)
    }
}


###Output the results of t-test to clipboard and files.
w.ttest <- function (table = summary.ttest, filename = "clipboard", en = 1) {
    rows <- length(table[, 1])
    columns <- length(table)
    Factor <- substring(row.names(table), 1, regexpr("=", row.names(table)) - 
        1)
    Group <- substring(row.names(table), regexpr("=", row.names(table)) + 
        1)
    for (i in 1:(rows - 1)) {
        j <- 1
        while (Factor[i] == Factor[i + j]) {
            Factor[i + j] <- ""
            j <- j + 1
            if ((i + j) > rows) 
                break
        }
    }
    table[, 3] <- as.numeric(as.character(data.frame(table)[, 
        3]))
    table <- signif(data.frame(table), digits = 3)
    table[, 3] <- ifelse(is.na(table[, 3]), "", table[, 3])
    table[, 1] <- paste(table[, 1], " +- ", table[, 2], sep = "")
    table <- table[, c(1, 3)]
    colnames(table)[1] <- "mean +- SD"
    table <- cbind(Factor, Group, table)
    rownames(table) <- NULL
    colnames(table)[4] <- ifelse(en == 1, "p.value", gettext(domain="R-RcmdrPlugin.EZR", 
        "p.value"))
    if (en == 0) 
        colnames(table)[1:3] <- gettext(domain="R-RcmdrPlugin.EZR", 
            c("Factor", "Group", "mean +- SD"))
#    print(data.frame(table), quote=FALSE, row.names=FALSE)
	dataframe_print(table, printrow=0)
#    print(table)
#    print(paste("Write to ", filename, sep = ""))
    if (filename == "clipboard") {
        write.table(data.frame(table), "clipboard", sep = "\t", 
            row.names = FALSE)
    }
    else {
        write.csv(data.frame(table), file = as.character(filename), row.names = FALSE)
    }
}


w.survival <- function (table = km.summary.table, filename = "clipboard", en = 1) {
    rows <- length(table[, 1])
    columns <- length(table)
    Factor <- substring(row.names(table), 1, regexpr("=", row.names(table)) - 
        1)
    Group <- substring(row.names(table), regexpr("=", row.names(table)) + 
        1)
    for (i in 1:(rows - 1)) {
        j <- 1
        while (Factor[i] == Factor[i + j]) {
            Factor[i + j] <- ""
            j <- j + 1
            if ((i + j) > rows) 
                break
        }
    }
    if (colnames(table)[2] == gettext(domain="R-RcmdrPlugin.EZR", 
        "median survival")) {
        table[, 2] <- paste(table[, 2], " (", table[, 3], ")", 
            sep = "")
        table <- table[, c(1, 2, 4)]
        if (en == 1) {
            colnames(table)[1:3] <- c("n", "median survival", 
                "p.value")
        }
        else {
            colnames(table)[1:3] <- gettext(domain="R-RcmdrPlugin.EZR", 
                c("n", "median survival", "p.value"))
        }
    }
    if (colnames(table)[2] == gettext(domain="R-RcmdrPlugin.EZR", 
        "survival rate")) {
        table[, 2] <- paste(table[, 2], " ", table[, 3], sep = "")
        table[, 4] <- paste(table[, 4], " (", table[, 5], ")", 
            sep = "")
        table <- table[, c(1, 2, 4, 6)]
        if (en == 1) {
            colnames(table)[1:4] <- c("n", "survival rate", "median survival", 
                "p.value")
        }
        else {
            colnames(table)[1:4] <- gettext(domain="R-RcmdrPlugin.EZR", 
                c("n", "survival rate", "median survival", "p.value"))
        }
    }
    table <- cbind(Factor, Group, table)
    rownames(table) <- NULL
    if (en == 0) 
        colnames(table)[1:2] <- gettext(domain="R-RcmdrPlugin.EZR", 
            c("Factor", "Group"))
#    print(table, quote=FALSE, row.names=FALSE, col.names=TRUE)
	dataframe_print(table, printrow=0)
#    print(table)
#    print(paste("Write to ", filename, sep = ""))
    if (filename == "clipboard") {
        write.table(data.frame(table), "clipboard", sep = "\t", 
            row.names = FALSE)
    }
    else {
        write.csv(data.frame(table), file = as.character(filename), row.names = FALSE)
    }
}


w.ci <- function (table = ci.summary.table, filename = "clipboard", en = 1) {
    rows <- length(table[, 1])
    columns <- length(table)
    Group_Factor <- substring(row.names(table), regexpr(",", 
        row.names(table)) + 2)
    Factor <- substring(Group_Factor, 1, regexpr("=", Group_Factor) - 
        1)
    Group <- substring(Group_Factor, regexpr("=", Group_Factor) + 
        1)
    Event <- substring(row.names(table), 1, regexpr(",", row.names(table)) - 
        1)
    for (i in 1:(rows - 1)) {
        j <- 1
        while (Factor[i] == Factor[i + j]) {
            Factor[i + j] <- ""
            j <- j + 1
            if ((i + j) > rows) 
                break
        }
    }
    if (colnames(table)[2] == gettext(domain="R-RcmdrPlugin.EZR", 
        "incidence")) {
        table[, 2] <- paste(table[, 2], " ", table[, 3], sep = "")
        table <- table[, c(1, 2, 4, 5)]
        if (en == 1) {
            colnames(table)[1:4] <- c("n", "incidence", "median time", 
                "p.value")
        }
        else {
            colnames(table)[1:4] <- gettext(domain="R-RcmdrPlugin.EZR", 
                c("n", "incidence", "median time", "p.value"))
        }
    }
    else {
        if (en == 1) {
            colnames(table)[1:3] <- c("n", "median time", "p.value")
        }
        else {
            colnames(table)[1:3] <- gettext(domain="R-RcmdrPlugin.EZR", 
                c("n", "median time", "p.value"))
        }
    }
    table <- cbind(Factor, Group, Event, table)
    rownames(table) <- NULL
    if (en == 0) 
        colnames(table)[1:3] <- gettext(domain="R-RcmdrPlugin.EZR", 
            c("Factor", "Group", "Event"))
#    print(table, quote=FALSE, row.names=FALSE, col.names=TRUE)
	dataframe_print(table, printrow=0)
#    print(table)
#    print(paste("Write to ", filename, sep = ""))
    if (filename == "clipboard") {
        write.table(data.frame(table), "clipboard", sep = "\t", 
            row.names = FALSE)
    }
    else {
        write.csv(data.frame(table), file = as.character(filename), row.names = FALSE)
    }
}


ChrToFactor <- function(dataset){
	for (i in 1:length(dataset)){
		if (is.character(dataset[,i])==TRUE){
			dataset[,i] <- factor(dataset[,i])
			cat(paste(colnames(dataset[i]), " ", gettext(domain="R-RcmdrPlugin.EZR","was converted to a factor."), "\n", sep=""))				
		}
	}
	return(dataset)
}


.funincrisk <- function(cdat, conf.level) {
### from epiR package, required for epi.tests()
        N. <- 1 - ((1 - conf.level)/2)
        a <- cdat[, 1]
        n <- cdat[, 2]
        b <- n - a
        p <- a/n
        a. <- ifelse(a == 0, a + 1, a)
        b. <- ifelse(b == 0, b + 1, b)
        low <- a./(a. + (b. + 1) * (1/qf(1 - N., 2 * a., 2 * 
            b. + 2)))
        up <- (a. + 1)/(a. + 1 + b./(1/qf(1 - N., 2 * b., 2 * 
            a. + 2)))
        low <- ifelse(a == 0, 0, low)
        up <- ifelse(a == n, 1, up)
        rval <- data.frame(est = p, lower = low, upper = up)
        rval
}


epi.tests <- function (dat, conf.level = 0.95, verbose = FALSE) 
{
### from epiR package 0.9-45
    N. <- 1 - ((1 - conf.level)/2)
    z <- qnorm(N., mean = 0, sd = 1)
    .funincrisk <- function(cdat, conf.level) {
        N. <- 1 - ((1 - conf.level)/2)
        a <- cdat[, 1]
        n <- cdat[, 2]
        b <- n - a
        p <- a/n
        a. <- ifelse(a == 0, a + 1, a)
        b. <- ifelse(b == 0, b + 1, b)
        low <- a./(a. + (b. + 1) * (1/qf(1 - N., 2 * a., 2 * 
            b. + 2)))
        up <- (a. + 1)/(a. + 1 + b./(1/qf(1 - N., 2 * b., 2 * 
            a. + 2)))
        low <- ifelse(a == 0, 0, low)
        up <- ifelse(a == n, 1, up)
        rval <- data.frame(est = p, lower = low, upper = up)
        rval
    }
    a <- dat[1]
    b <- dat[3]
    c <- dat[2]
    d <- dat[4]
    M1 <- a + c
    M0 <- b + d
    N1 <- a + b
    N0 <- c + d
    total <- a + b + c + d
    tdat <- as.matrix(cbind(M1, total))
    trval <- .funincrisk(tdat, conf.level)
    tp <- trval$est
    tp.low <- trval$lower
    tp.up <- trval$upper
    tprev <- data.frame(est = tp, lower = tp.low, upper = tp.up)
    tdat <- as.matrix(cbind(N1, total))
    trval <- .funincrisk(tdat, conf.level)
    ap <- trval$est
    ap.low <- trval$lower
    ap.up <- trval$upper
    aprev <- data.frame(est = ap, lower = ap.low, upper = ap.up)
    tdat <- as.matrix(cbind(a, M1))
    trval <- .funincrisk(tdat, conf.level)
    se <- trval$est
    se.low <- trval$lower
    se.up <- trval$upper
    sensitivity <- data.frame(est = se, lower = se.low, upper = se.up)
    tdat <- as.matrix(cbind(d, M0))
    trval <- .funincrisk(tdat, conf.level)
    sp <- trval$est
    sp.low <- trval$lower
    sp.up <- trval$upper
    specificity <- data.frame(est = sp, lower = sp.low, upper = sp.up)
    tdat <- as.matrix(cbind(a, N1))
    trval <- .funincrisk(tdat, conf.level)
    ppv <- trval$est
    ppv.low <- trval$lower
    ppv.up <- trval$upper
    pv.positive <- data.frame(est = ppv, lower = ppv.low, upper = ppv.up)
    tdat <- as.matrix(cbind(d, N0))
    trval <- .funincrisk(tdat, conf.level)
    npv <- trval$est
    npv.low <- trval$lower
    npv.up <- trval$upper
    pv.negative <- data.frame(est = npv, lower = npv.low, upper = npv.up)
    lrpos <- (a/M1)/(1 - (d/M0))
    lrpos.low <- exp(log(lrpos) - z * sqrt((1 - se)/(M1 * se) + 
        (sp)/(M0 * (1 - sp))))
    lrpos.up <- exp(log(lrpos) + z * sqrt((1 - se)/(M1 * se) + 
        (sp)/(M0 * (1 - sp))))
    lr.positive <- data.frame(est = lrpos, lower = lrpos.low, 
        upper = lrpos.up)
    lrneg <- (1 - (a/M1))/(d/M0)
    lrneg.low <- exp(log(lrneg) - z * sqrt((se)/(M1 * (1 - se)) + 
        (1 - sp)/(M0 * (sp))))
    lrneg.up <- exp(log(lrneg) + z * sqrt((se)/(M1 * (1 - se)) + 
        (1 - sp)/(M0 * (sp))))
    lr.negative <- data.frame(est = lrneg, lower = lrneg.low, 
        upper = lrneg.up)
    tdat <- as.matrix(cbind((a + d), total))
    trval <- .funincrisk(tdat, conf.level)
    da <- trval$est
    da.low <- trval$lower
    da.up <- trval$upper
    diag.acc <- data.frame(est = da, lower = da.low, upper = da.up)
    dOR.p <- (a * d)/(b * c)
    lndOR <- log(dOR.p)
    lndOR.var <- 1/a + 1/b + 1/c + 1/d
    lndOR.se <- sqrt(1/a + 1/b + 1/c + 1/d)
    lndOR.l <- lndOR - (z * lndOR.se)
    lndOR.u <- lndOR + (z * lndOR.se)
    dOR.se <- exp(lndOR.se)
    dOR.low <- exp(lndOR.l)
    dOR.up <- exp(lndOR.u)
    diag.or <- data.frame(est = dOR.p, lower = dOR.low, upper = dOR.up)
    ndx <- 1/(se - (1 - sp))
    ndx.1 <- 1/(se.low - (1 - sp.low))
    ndx.2 <- 1/(se.up - (1 - sp.up))
    ndx.low <- min(ndx.1, ndx.2)
    ndx.up <- max(ndx.1, ndx.2)
    nnd <- data.frame(est = ndx, lower = ndx.low, upper = ndx.up)
    c.p <- se - (1 - sp)
    c.1 <- se.low - (1 - sp.low)
    c.2 <- se.up - (1 - sp.up)
    c.low <- min(c.1, c.2)
    c.up <- max(c.1, c.2)
    youden <- data.frame(est = c.p, lower = c.low, upper = c.up)
    if (verbose == TRUE) {
        rval <- list(aprev = aprev, tprev = tprev, se = sensitivity, 
            sp = specificity, diag.acc = diag.acc, diag.or = diag.or, 
            nnd = nnd, youden = youden, ppv = pv.positive, npv = pv.negative, 
            plr = lr.positive, nlr = lr.negative)
        return(rval)
    }
    if (verbose == FALSE) {
        r1 <- c(a, b, N1)
        r2 <- c(c, d, N0)
        r3 <- c(M1, M0, M0 + M1)
        tab <- as.data.frame(rbind(r1, r2, r3))
        colnames(tab) <- gettext(domain="R-RcmdrPlugin.EZR",c("Disease positive", "Disease negative", "Total"))
        rownames(tab) <- gettext(domain="R-RcmdrPlugin.EZR",c("Test positive", "Test negative", "Total"))
        tab <- format.data.frame(tab, digits = 3, justify = "right")
        print(tab)
        cat("\n", gettext(domain="R-RcmdrPlugin.EZR","Point estimates and"), conf.level * 100, "%", gettext(domain="R-RcmdrPlugin.EZR","CIs:"))
        cat("\n---------------------------------------------------------\n")
		res.table <- c(aprev$est, aprev$lower, aprev$upper)
		res.table <- rbind(res.table, c(tprev$est, tprev$lower, tprev$upper))
		res.table <- rbind(res.table, c(sensitivity$est, sensitivity$lower, sensitivity$upper))
		res.table <- rbind(res.table, c(specificity$est, specificity$lower, specificity$upper))
		res.table <- rbind(res.table, c(pv.positive$est, pv.positive$lower, pv.positive$upper))
		res.table <- rbind(res.table, c(pv.negative$est, pv.negative$lower, pv.negative$upper))
		res.table <- rbind(res.table, c(diag.acc$est, diag.acc$lower, diag.acc$upper))
		res.table <- rbind(res.table, c(lr.positive$est, lr.positive$lower, lr.positive$upper))
		res.table <- rbind(res.table, c(lr.negative$est, lr.negative$lower, lr.negative$upper))
		res.table <- round(res.table, digits=3)
		colnames(res.table) <- gettext(domain="R-RcmdrPlugin.EZR", c("Estimation", "Lower CI", "Upper CI"))
		rownames(res.table) <- gettext(domain="R-RcmdrPlugin.EZR", c("Apparent prevalence", "True prevalence", "Sensitivity", "Specificity",
		"Positive predictive value", "Negative predictive value", "Diagnstic accuracy", "Likelihood ratio of a positive test",
		"Likelihood ratio of a negative test"))
		print(res.table)
        cat("---------------------------------------------------------")
        cat("\n")
    }
}


epi.kappa <- function (dat, conf.level = 0.95) 
{
### from epiR package 0.9-27. In this version, mcNemar test is pweformed.
    a <- dat[1]
    b <- dat[3]
    c <- dat[2]
    d <- dat[4]
    N. <- 1 - ((1 - conf.level)/2)
    z <- qnorm(N., mean = 0, sd = 1)
    lower <- "lower"
    upper <- "upper"
    n <- a + b + c + d
    pO <- (a + d)/n
    pE.pos <- ((a + b) * (a + c))/n^2
    pE.neg <- ((c + d) * (b + d))/n^2
    pE <- pE.pos + pE.neg
    kappa <- (pO - pE)/(1 - pE)
    se.kappa <- sqrt((pO * (1 - pO))/(n * (1 - pE)^2))
    kappa.low <- kappa - (z * se.kappa)
    kappa.up <- kappa + (z * se.kappa)
    mcnemar <- (b - c)^2/(b + c)
    p.chi2 <- 1 - pchisq(mcnemar, df = 1)
    kappa <- as.data.frame(cbind(kappa, kappa.low, kappa.up))
    names(kappa) <- c("est", lower, upper)
    mcnemar <- as.data.frame(cbind(test.statistic = mcnemar, 
        df = 1, p.value = p.chi2))
    rval <- list(kappa = kappa, mcnemar = mcnemar)
    return(rval)
}


dot.plot <- function(x, y, accu=0, stp=0, log.flag=FALSE, simple=FALSE, symmetrical=TRUE, ...) {                                           
       #modified from http://aoki2.si.gunma-u.ac.jp/R/dot_plot.html
        OK <- complete.cases(x, y)                                      
        x <- x[OK]
        x <- as.factor(x)
        y <- y[OK]
		x.name <- unique(x)
		if (is.factor(x)) {                                             
                x <- as.integer(x)                                      
        }
        if (log.flag == TRUE) {                                         
                y0 <- y                                                 
                y <- log10(y)                                           
        }
        if (accu == 0) {                                                
                accu <- diff(range(y))/100                              
        }
        if(stp == 0) {                                                   
                stp <- (diff(range(x))+1)/100                           
		}
        y <- round(y/accu)*accu                                         
        x1 <- unique(x)                                                 
        for (i in seq(along=x1)) {                                      
                freq <- table(y[x==x1[i]])                              
                for (j in seq(along=freq)) {                            
                        if (freq[j] >= 2) {                             
                                offset <- ifelse(symmetrical, (freq[j]-1)/2*stp, 0)
                                for (k in seq(along=y)) {
                                        if (abs(y[k]-as.numeric(names(freq)[j])) < 1e-10 && abs(x[k]-x1[i]) < 1e-10) {
                                                freq[j] <- freq[j]-1
                                                x[k] <- x[k]-offset+freq[j]*stp
                                        }
                                }
                        }
                }
        }
        if (log.flag) {                                       
                plot(x, y, type="n", xaxt="n", yaxt="n", xlim=c(min(x)-0.5, max(x)+0.5), ...)
                options(warn=-1)
                points(x, y, ...)
                options(warn=0)
                y0 <- floor(log10(y0))
                log.min <- min(y0)
                y2 <- 1:10*10^log.min
                n <- max(y0)-log.min
                y1 <- rep(y2, n+1)*10^rep(0:n, each=10)
                if (simple) {
                        y2 <- y1[abs(log10(y1)-round(log10(y1))) < 1e-6]
                        axis(2, at=log10(y1), labels=FALSE)
                        axis(2, at=log10(y2), labels=y2)
                }
                else {
                        axis(2, at=log10(y1), labels=y1)
                }
        }
        else {
                plot(x, y, xaxt="n", xlim=c(min(x)-0.5, max(x)+0.5), ...)
        }
        if (length(x.name)>1) {
			axis(1, at=x1, labels=as.character(x.name))
		}	
}


OrderedPlot <- function(y, group=NULL, type="line", xlab="", ylab="Value", ylog=FALSE, lowlim=NULL, uplim=NULL, decreasing=FALSE){
	#For waterfall plot, ordered chart
	if (is.null(group)){		
		cc <- complete.cases(y)
	} else {
		cc <- complete.cases(y, group)
	}
	y <- y[cc]

	if (is.null(lowlim)) lowlim <- min(y)
	if (type=="box" & ylog==FALSE & lowlim>0) lowlim <- 0
	if (is.null(uplim)) uplim=max(y)
	ylim=c(lowlim, uplim)

	ylog <- ifelse(ylog==TRUE, "y", "")

	if (is.null(group) | type=="box"){
		if (type=="line"){
			Order <- order(y, decreasing=decreasing)
			plot(x=seq(from=0, to=1, length.out=length(y)), y=y[Order], xaxp=c(0,1,10), type="l", ylim=ylim, log=ylog, xlab=xlab, ylab=ylab)
		} else {
			Order <- order(y, decreasing=decreasing)
#			names.arg <- c("0", rep("", length(y)-2), "1")
			names.arg=NULL
			barplot(y[Order], names.arg=names.arg, axis.lty=1, ylim=ylim, log=ylog, axisnames=TRUE)
		}
	}

	if (!is.null(group) & type=="line"){
		group <- group[cc]
		group <- factor(group)
		levels <- levels(group)
		j <- 1
		for (i in levels){
			Order <- order(y[group==i], decreasing=decreasing)
			axt <- "s"
			if (j>1) {par(new=T);axt <- "n"; xlab <- ""; ylab <- ""}
			plot(x=seq(from=0, to=1, length.out=length(y[group==i])), y=y[group==i][Order], xaxp=c(0,1,10), type="l", ylim=ylim, log=ylog, xlab=xlab, ylab=ylab, xaxt=axt, yaxt=axt, lty=j)
			j <- j+1
		}
		if (decreasing==FALSE){
			legend("topleft", levels, col=1, lty=1:32, lwd=1,  box.lty=0)
		} else {
			legend("topright", levels, col=1, lty=1:32, lwd=1,  box.lty=0)			
		}
	}
}


BarplotFor3Factors <- function(First, Second, Third, prop=0, col=0, data){

	dataset <- eval(parse(text=data))	
	legend <- eval(parse(text=paste("levels(factor(", data, "$", First, "))", sep="")))
	groups <- eval(parse(text=paste("levels(factor(", data, "$", Second, "))", sep="")))
	levels <- eval(parse(text=paste("levels(factor(", data, "$", Third, "))", sep="")))
	num <- length(levels)
	
	colors <- gray(2:(length(legend)+1) / (length(legend)+2))
	if (col==1) colors <- 2:(1+length(legend))

	res <- eval(parse(text=paste("xtabs(~", First, "+", Second, ", data=dataset, subset=", Third, "=='", levels[1], "')", sep="")))

	if (prop==0){
		barplot.table <- res
	}else{
		barplot.table <- prop.table(res,2)
	}

	dummy <- rep(0, length(barplot.table[,1]))

	for (i in 2:num){
		res <- eval(parse(text=paste("xtabs(~", First, "+", Second, ", data=dataset, subset=", Third, "=='", levels[i], "')", sep="")))
		if (prop==0){
			barplot.table <- cbind(barplot.table, " "=dummy, res)
		}else{
			barplot.table <- cbind(barplot.table, " "=dummy, prop.table(res,2))
		}
	}

	mar <- par("mar")
	mar[1] <- mar[1] + 2.5
	mar[3] <- mar[3] + 1.5
	par(mar=mar)
	opar <- par(mar = mar)
	on.exit(par(opar))

	if(prop==1){
		legend.y <- 1.2
	} else {
		max.height <- 0
		for(i in 1:length(barplot.table[1,])) {
			if (sum(barplot.table[,i]) > max.height) {max.height <- sum(barplot.table[,i])}
		}
		legend.y <- max.height * 1.2
	}

	(bplot <- barplot(barplot.table, beside=FALSE, xlab=NULL, ylab="Frequency", col=colors,
  	legend=legend, args.legend=list(y=legend.y, horiz=TRUE, title=First, box.lty=0), axis.lty=1))

	at <- NULL
	for (i in 1:num){
		at <- c(at, (bplot[(length(groups)+1)*(i-1)+1]+bplot[(length(groups)+1)*(i-1)+length(groups)])/2)
	}
	center <- (bplot[1] + bplot[length(bplot)])/2

	axis(1, at = center, labels = Second, line = 2, tick = FALSE, las=0)
#	axis(1, at = at, labels = rep(Third, length(levels)), line = 4, tick = FALSE, las=0)
	axis(1, at = center, labels = Third, line = 4, tick = FALSE, las=0)
	axis(1, at = at, labels = levels, line = 5, tick = FALSE, las=0)

}


nrisk <- function (x, times = pretty(x$time)) {
		#Function to count number at risk in KM plot.  Modified from survplot package to be applied for single group.
#    stopifnot(class(x) == "survfit")
    if (!is.null(x$strata)){
	ns <- length(x$strata)
    	idx <- rep.int(1:ns, x$strata)
    	str.n.risk <- split(x$n.risk, idx)
    	str.times <- split(x$time, idx)
    	m <- sapply(times, function(y) {
        sapply(1:ns, function(i) {
            w <- which(str.times[[i]] >= y)[1]
            ifelse(is.na(w), 0, str.n.risk[[i]][w])
        })
    	})
	rownames(m) <- names(x$strata)
    	colnames(m) <- times
 	} else {    
	str.n.risk <- x$n.risk
	str.times <- x$time
    	m <- sapply(times, function(y) {
            w <- which(str.times >= y)[1]
            ifelse(is.na(w), 0, str.n.risk[w])
    	})
#	rownames(m) <- names(x$strata)
#    	colnames(m) <- times
	}
    	m
}


prop.conf <- function(  r, n, conf){
       #modified from http://aoki2.si.gunma-u.ac.jp/R/p-conf.html
	p <- r/n 
    alpha <- 1-conf/100               
    if (p == 0) {                                   
		pl <- 0
        pu <- 1-alpha^(1/n)
    } else if (p == 1) {                            
        pl <- alpha^(1/n)
        pu <- 1
    } else {                                        
        nu1 <- 2*(n-r+1)
        nu2 <- 2*r
        Fv <- qf(alpha/2, nu1, nu2, lower.tail=FALSE)
        pl <- nu2/(nu1*Fv+nu2)
        nu1 <- 2*(r+1)
        nu2 <- 2*(n-r)
        Fv <- qf(alpha/2, nu1, nu2, lower.tail=FALSE)
        pu <- nu1*Fv/(nu1*Fv+nu2)
    }
        print(paste(gettext(domain="R-RcmdrPlugin.EZR","Probability :"), " ", round(p,3), sep=""), quote=F)
        print(paste(conf, gettext(domain="R-RcmdrPlugin.EZR","% confidence interval :"), " ", round(pl,3), " - ", round(pu,3), sep=""), quote=F)
}


prop.diff.conf <- function(r1, n1, r2, n2, conf) {
	alpha <- 1-conf/100
	p1 <- r1/n1
	p2 <- r2/n2
	D <- p1-p2
	SE <- sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
	pl <- D-qnorm(1-alpha/2)*SE
	pu <- D+qnorm(1-alpha/2)*SE
	print(paste(gettext(domain="R-RcmdrPlugin.EZR","Difference :"), " ", round(D,3), sep=""), quote=F)
	print(paste(conf, gettext(domain="R-RcmdrPlugin.EZR","% confidence interval :"), " ", round(pl,3), " - ", round(pu,3), sep=""), quote=F)
}

		
prop.ratio.conf <- function(r1, n1, r2, n2, conf) {
	alpha <- 1-conf/100
	p1 <- r1/n1
	p2 <- r2/n2
	RR<- p1/p2
	SE <- sqrt((n1-r1)/r1/n1+(n2-r2)/r2/n2)
	pl <- exp(log(RR)-qnorm(1-alpha/2)*SE)
	pu <- exp(log(RR)+qnorm(1-alpha/2)*SE)
	print(paste(gettext(domain="R-RcmdrPlugin.EZR","Ratio : "), round(RR,3), sep=""), quote=F)
	print(paste(conf, gettext(domain="R-RcmdrPlugin.EZR","% confidence interval : "), round(pl,3), " - ", round(pu,3), sep=""), quote=F)
}

		
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop(gettext(domain="R-RcmdrPlugin.EZR","vectors must be same length"))
		arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


StatMedplotMeans <- function(response, factor1, factor2, error.bars = c("se", "sd", "conf.int", "none"),
	level=0.95, xlab=deparse(substitute(factor1)), ylab=paste("mean of", deparse(substitute(response))),
	legend.lab=deparse(substitute(factor2)), main="Plot of Means",
	pch=1:n.levs.2, lty=1:n.levs.2, lwd=1:n.levs.2, col=palette(), yrange=NULL){
	if (!is.numeric(response)) stop(gettext(domain="R-RcmdrPlugin.EZR","Argument response must be numeric."))
	xlab # force evaluation
	ylab
	legend.lab
	error.bars <- match.arg(error.bars)
	if (missing(factor2)){
		if (!is.factor(factor1)) stop(gettext(domain="R-RcmdrPlugin.EZR","Argument factor1 must be a factor."))
		valid <- complete.cases(factor1, response)
		factor1 <- factor1[valid]
		response <- response[valid]
		means <- tapply(response, factor1, mean)
		sds <- tapply(response, factor1, sd)
		ns <- tapply(response, factor1, length)
		if (error.bars == "se") sds <- sds/sqrt(ns)
		if (error.bars == "conf.int") sds <- qt((1 - level)/2, df=ns - 1, lower.tail=FALSE) * sds/sqrt(ns)
		sds[is.na(sds)] <- 0
		if (is.null(yrange)){
			yrange <-  if (error.bars != "none") c( min(means - sds, na.rm=TRUE), max(means + sds, na.rm=TRUE)) else range(means, na.rm=TRUE)
		}
		levs <- levels(factor1)
		n.levs <- length(levs)
		plot(c(1, n.levs), yrange, type="n", xlab=xlab, ylab=ylab, axes=FALSE, main=main)
		points(1:n.levs, means, type="b", pch=16, cex=2)
		box()
		axis(2)
		axis(1, at=1:n.levs, labels=levs)
		if (error.bars != "none") arrows(1:n.levs, means - sds, 1:n.levs, means + sds,
				angle=90, lty=2, code=3, length=0.125)
	}
	else {
		if (!(is.factor(factor1) | is.factor(factor2))) stop(gettext(domain="R-RcmdrPlugin.EZR","Arguments factor1 and factor2 must be factors."))
		valid <- complete.cases(factor1, factor2, response)
		factor1 <- factor1[valid]
		factor2 <- factor2[valid]
		response <- response[valid]
		means <- tapply(response, list(factor1, factor2), mean)
		sds <- tapply(response, list(factor1, factor2), sd)
		ns <- tapply(response, list(factor1, factor2), length)
		if (error.bars == "se") sds <- sds/sqrt(ns)
		if (error.bars == "conf.int") sds <- qt((1 - level)/2, df=ns - 1, lower.tail=FALSE) * sds/sqrt(ns)
		sds[is.na(sds)] <- 0
		if (is.null(yrange)){
			yrange <-  if (error.bars != "none") c( min(means - sds, na.rm=TRUE), max(means + sds, na.rm=TRUE)) else range(means, na.rm=TRUE)
		}
		levs.1 <- levels(factor1)
		levs.2 <- levels(factor2)
		n.levs.1 <- length(levs.1)
		n.levs.2 <- length(levs.2)
		if (length(pch) == 1) pch <- rep(pch, n.levs.2)
		if (length(col) == 1) col <- rep(col, n.levs.2)
		if (length(lty) == 1) lty <- rep(lty, n.levs.2)
		if (length(lwd) == 1) lwd <- rep(lwd, n.levs.2)
		if (n.levs.2 > length(col)) stop(sprintf(gettext(domain="R-RcmdrPlugin.EZR","Number of groups for factor2, %d, exceeds number of distinct colours, %d."), n.levs.2, length(col)))		
		plot(c(1, n.levs.1 * 1.2), yrange, type="n", xlab=xlab, ylab=ylab, axes=FALSE, main=main)
		box()
		axis(2)
		axis(1, at=1:n.levs.1, labels=levs.1)
		for (i in 1:n.levs.2){
			points(1:n.levs.1, means[, i], type="b", pch=pch[i], cex=2, col=col[i], lty=lty[i], lwd=lwd[i])
			if (error.bars != "none") arrows(1:n.levs.1, means[, i] - sds[, i],
					1:n.levs.1, means[, i] + sds[, i], angle=90, code=3, col=col[i], lty=lty[i], lwd=lwd[i], length=0.125)
		}
#		x.posn <- n.levs.1 * 1.4
		x.posn <- n.levs.1
		y.posn <- sum(c(0.1, 0.9) * par("usr")[c(3,4)])
#		text(x.posn, y.posn, legend.lab, adj=c(0, -.5))
#		legend(x.posn, y.posn, levs.2, pch=pch, col=col, lty=lty)
		legend("topright", levs.2, pch=pch, col=col, lty=lty, lwd=lwd, title=legend.lab, box.lty=0)
	}
	invisible(NULL)
}


skewness.kurtosis <- function(x){
	x <- x[!is.na(x)]
#	skewness <- signif(mean((x-mean(x))^3)/(sd(x)^3), digits=3) #sample skewness
#	kurtosis <- signif(mean((x-mean(x))^4)/(sd(x)^4)-3, digits=3) #sample kurtosis
	n <- length(x)
	m <- mean(x)
	sd <- sd(x)
	skewness <- signif({n/(n-1)/(n-2)} * sum((x-m)^3) / sd^3, digits=3) #population skewness, same as kurt(), skew() in excel
	kurtosis <- signif({n*(n+1)/(n-1)/(n-2)/(n-3)} * sum((x-m)^4) / sd^4 - 3*(n-1)^2/(n-2)/(n-3), digits=3) #population kurtosis
	res <- data.frame(c(gettext(domain="R-RcmdrPlugin.EZR","Skewness (0 for normal distribution)"), gettext(domain="R-RcmdrPlugin.EZR","Kurtosis (0 for normal distribution)")), c(skewness, kurtosis))
	rownames(res) <- c("", " ")
	colnames(res) <- c("", " ")
	return(res)
}

	
Cochran.Q.test <- function(x) {
#http://aoki2.si.gunma-u.ac.jp/R/Cochran-Q-test.html
        data.name <- deparse(substitute(x))
        method <- "Cochran's Q test"
        x <- subset(x, complete.cases(x))                       
        k <- ncol(x)                                            
        g <- colSums(x)                                         
        l <- rowSums(x)                                         
        Q <- ((k-1)*(k*sum(g^2)-sum(g)^2))/(k*sum(l)-sum(l^2))  
        df <- k-1                                               
        p <- pchisq(Q, df, lower.tail=FALSE)                    
        names(Q) <- "X-squared"
        names(df) <- "df"
        return(structure(list(statistic=Q, parameter=df, p.value=p,
                method=method, data.name=data.name), class="htest"))
}


pairwise.prop2.test <- function (x, n, p.adjust.method = p.adjust.methods, test.function=prop.test, ...){
#modified from http://aoki2.si.gunma-u.ac.jp/R/p_multi_comp2.html
#to extend for fisher.test() and to show the results with group names.
    p.adjust.method <- match.arg(p.adjust.method)
    METHOD <- deparse(substitute(test.function))
    DNAME <- deparse(substitute(x))
    if (is.matrix(x)) {
        if (ncol(x) < 2) 
            stop("'x' must have at least 2 columns")
    } else if (is.vector(x) && is.vector(n))
        x <- cbind(x, n-x)
    else
        stop("'x' must be a matrix, or 'x', and 'n' must be a vector")
    if (nrow(x) < 2) 
        stop("too few groups")
	group.names <- rownames(x)
    compare.levels <- function(i, j) {
        test.function(x[c(i, j),], ...)$p.value                                 
    }
    level.names <- names(x)
    if (is.null(level.names)) 
        level.names <- group.names[seq_along(1:nrow(x))]
    PVAL <- pairwise.table(compare.levels, level.names, p.adjust.method)        
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
        p.adjust.method = p.adjust.method)
    class(ans) <- "pairwise.htest"
    ans
}


pairwise.pairedt.test <- function (response, group=NULL, data.name, p.adjust.method = p.adjust.methods){
#modified from http://aoki2.si.gunma-u.ac.jp/R/p_multi_comp2.html
	if(!is.null(group)){
		group <- factor(group)
		contrasts(group) <- "contr.Sum"
	}
	
    	p.adjust.method <- match.arg(p.adjust.method)
	method <- "Paired t-test"
	time.names <- colnames(response)
	n <- length(time.names)	
	p <- NULL
	n.comp <- 0				
	for (i in 2:n){
		for (j in 1:(n-1)){
			if (j < i){
				pairwise.response <- response[, time.names==time.names[i] | time.names==time.names[j]]
				if(!is.null(group)){
					AnovaModel <- lm(pairwise.response ~ group, na.action=na.omit)
					time <- colnames(pairwise.response)
					time <- data.frame(Time = time)
					res <- Anova(AnovaModel, idata=time, idesign=~Time, type="III")
					res <- capture.output(summary(res, multivariate=FALSE))
				}else{
					AnovaModel <- lm(pairwise.response ~ 1, na.action=na.omit)
					time <- colnames(pairwise.response)
					time <- data.frame(Time = time)
					res <- Anova(AnovaModel, idata=time, idesign=~Time, type="III")
					res <- capture.output(summary(res, multivariate=FALSE))
				}
				###The results of Anova() cannot be obtained by summry(), and thus should be extracted from the output text.
				l <- 0	###Looking for a row that stat with "Time"		
				for(k in 1:length(res)){
					if(substr(res[k],1,4)=="Time"){
						res <- res[k]
						break
					}
				}
				res <- strsplit(res, split=" ")
				l <- 0	###Extract p value
				for(k in 1:length(res[[1]])){
					if(res[[1]][k]!="")l <- l+1
					if(l==7){
						p[j] <- res[[1]][k] 
						break
					}
				}
				n.comp <- n.comp+1
			} else {
				p[j] <- NA
			}
		}
		if (i==2){
			pairwise.table <- p
		} else {
			pairwise.table <- rbind(pairwise.table, p)
		}
	}
	pairwise.table <- matrix (p.adjust(pairwise.table, method=p.adjust.method, n.comp), n-1)
	rownames(pairwise.table) <- time.names[2:n]
	colnames(pairwise.table) <- time.names[1:n-1]
    	ans <- list(method=method, data.name=data.name, p.value = pairwise.table,  p.adjust.method = p.adjust.method)
   	class(ans) <- "pairwise.htest"
    	ans
}


pairwise.kruskal.test <- function (response, group, data.name, p.adjust.method = p.adjust.methods){
#modified from http://aoki2.si.gunma-u.ac.jp/R/p_multi_comp2.html
    	p.adjust.method <- match.arg(p.adjust.method)
	method <- "Mann-Whitney U test"
	group.names <- levels(factor(group))
	n <- length(group.names)	
	p <- NULL
	n.comp <- 0			
	for (i in 2:n){
		for (j in 1:(n-1)){
			if (j < i){
				pairwise.response <- response[group==group.names[i] | group==group.names[j]]
				pairwise.group <- group[group==group.names[i] | group==group.names[j]]
				res <- wilcox.test(pairwise.response ~ factor(pairwise.group))
				p[j] <- signif(res$p.value,digits=3)
				n.comp <- n.comp+1
			} else {
				p[j] <- NA
			}
		}
		if (i==2){
			pairwise.table <- p
		} else {
			pairwise.table <- rbind(pairwise.table, p)
		}
	}
	pairwise.table <- matrix (p.adjust(pairwise.table, method=p.adjust.method, n.comp), n-1)
	rownames(pairwise.table) <- group.names[2:n]
	colnames(pairwise.table) <- group.names[1:n-1]
    	ans <- list(method=method, data.name=data.name, p.value = pairwise.table,  p.adjust.method = p.adjust.method)
   	class(ans) <- "pairwise.htest"
    	ans
}


pairwise.friedman.test <- function (response, data.name, p.adjust.method = p.adjust.methods){
#modified from http://aoki2.si.gunma-u.ac.jp/R/p_multi_comp2.html
    	p.adjust.method <- match.arg(p.adjust.method)
	method <- "Wilcoxon signed rank test"
	time.names <- colnames(response)
	n <- length(time.names)	
	p <- NULL
	n.comp <- 0		
	for (i in 2:n){
		for (j in 1:(n-1)){
			if (j < i){
				pairwise.response1 <- response[, time.names==time.names[i]]
				pairwise.response2 <- response[, time.names==time.names[j]]
				res <- wilcox.test(pairwise.response1, pairwise.response2, alternative='two.sided', paired=TRUE)
				p[j] <- signif(res$p.value, digits=3)
				n.comp <- n.comp+1
			} else {
				p[j] <- NA
			}
		}
		if (i==2){
			pairwise.table <- p
		} else {
			pairwise.table <- rbind(pairwise.table, p)
		}
	}
	pairwise.table <- matrix (p.adjust(pairwise.table, method=p.adjust.method, n.comp), n-1)
	rownames(pairwise.table) <- time.names[2:n]
	colnames(pairwise.table) <- time.names[1:n-1]
    	ans <- list(method=method, data.name=data.name, p.value = pairwise.table,  p.adjust.method = p.adjust.method)
   	class(ans) <- "pairwise.htest"
    	ans
}


pairwise.logrank.test <- function (time, event, group, strata=NULL, data.name, p.adjust.method = p.adjust.methods, rho=0){
#modified from http://aoki2.si.gunma-u.ac.jp/R/p_multi_comp2.html
    	p.adjust.method <- match.arg(p.adjust.method)
	if (rho==0) method <- "logrank test" else method <- "Generalized Wilcoxon test"
	group.names <- levels(factor(group))
	n <- length(group.names)	
	p <- NULL
	n.comp <- 0				
	for (i in 2:n){
		for (j in 1:(n-1)){
			if (j < i){
				pairwise.time <- time[group==group.names[i] | group==group.names[j]]
				pairwise.event <- event[group==group.names[i] | group==group.names[j]]
				pairwise.group <- group[group==group.names[i] | group==group.names[j]]
				if(is.null(strata)){
					res <- survdiff(Surv(pairwise.time, pairwise.event==1)~pairwise.group, rho=rho)
				} else {
					pairwise.strata <- strata[group==group.names[i] | group==group.names[j]]
					res <- survdiff(Surv(pairwise.time, pairwise.event==1)~pairwise.group+strata(pairwise.strata), rho=rho)
				}
				p[j] <- signif(pchisq(c(res$chisq), df=1, lower.tail=FALSE),digits=3)
				n.comp <- n.comp+1
			} else {
				p[j] <- NA
			}
		}
		if (i==2){
			pairwise.table <- p
		} else {
			pairwise.table <- rbind(pairwise.table, p)
		}
	}
	pairwise.table <- matrix (p.adjust(pairwise.table, method=p.adjust.method, n.comp), n-1)
	rownames(pairwise.table) <- group.names[2:n]
	colnames(pairwise.table) <- group.names[1:n-1]
    	ans <- list(method=method, data.name=data.name, p.value = pairwise.table,  p.adjust.method = p.adjust.method)
   	class(ans) <- "pairwise.htest"
    	ans
}

		
pairwise.gray.test <- function (time, event, group, data.name, p.adjust.method = p.adjust.methods, endpoint=1){
#modified from http://aoki2.si.gunma-u.ac.jp/R/p_multi_comp2.html
    	p.adjust.method <- match.arg(p.adjust.method)
	method <- "Gray test"
	group.names <- levels(factor(group))
	n <- length(group.names)	
	p <- NULL
	n.comp <- 0		
	for (i in 2:n){
		for (j in 1:(n-1)){
			if (j < i){
				pairwise.time <- time[group==group.names[i] | group==group.names[j]]
				pairwise.event <- event[group==group.names[i] | group==group.names[j]]
				pairwise.group <- group[group==group.names[i] | group==group.names[j]]
				res <- cuminc(pairwise.time, pairwise.event, pairwise.group, cencode=0, na.action=na.omit)
				p[j] <- signif(res$Tests[endpoint, 2],digits=3)
				n.comp <- n.comp+1
			} else {
				p[j] <- NA
			}
		}
		if (i==2){
			pairwise.table <- p
		} else {
			pairwise.table <- rbind(pairwise.table, p)
		}
	}
	pairwise.table <- matrix (p.adjust(pairwise.table, method=p.adjust.method, n.comp), n-1)
	rownames(pairwise.table) <- group.names[2:n]
	colnames(pairwise.table) <- group.names[1:n-1]
    	ans <- list(method=method, data.name=data.name, p.value = pairwise.table,  p.adjust.method = p.adjust.method)
   	class(ans) <- "pairwise.htest"
    	ans
}


Steel.Dwass <- function(data, group){
#modified from http://aoki2.si.gunma-u.ac.jp/R/Steel-Dwass.html
        OK <- complete.cases(data, group)                               
        data <- data[OK]
        group <- factor(group[OK])									
        n.i <- table(group)                                             
        ng <- length(n.i)                                               
        t <- combn(ng, 2, function(ij) {
                i <- ij[1]
                j <- ij[2]
                r <- rank(c(data[group == levels(factor(group))[i]], data[group == levels(factor(group))[j]]))   
                R <- sum(r[1:n.i[i]])                                   
                N <- n.i[i]+n.i[j]                                      
                E <- n.i[i]*(N+1)/2                                     
                V <- n.i[i]*n.i[j]/(N*(N-1))*(sum(r^2)-N*(N+1)^2/4)     
                return(abs(R-E)/sqrt(V))                                
        })
        p <- ptukey(t*sqrt(2), ng, Inf, lower.tail=FALSE)               
        result <- cbind(t, p)                                           
        rownames(result) <- combn(levels(factor(group)), 2, paste, collapse=":")
        return(result)
}


Steel <- function(data, group) {
#modified from http://aoki2.si.gunma-u.ac.jp/R/Steel.html
        get.rho <- function(ni)                                 
        {
                k <- length(ni)
                rho <- outer(ni, ni, function(x, y) { sqrt(x/(x+ni[1])*y/(y+ni[1])) })
                diag(rho) <- 0
                sum(rho[-1, -1])/(k-2)/(k-1)
        }

        OK <- complete.cases(data, group)                       
        data <- data[OK]
        group <- factor(group[OK])                             
        ni <- table(group)                                     
        a <- length(ni)                                        
        control <- data[group == levels(factor(group))[1]]     
        n1 <- length(control)                                  
        t <- numeric(a)
        rho <- ifelse(sum(n1 == ni) == a, 0.5, get.rho(ni))    
	p.value <- numeric(a) 
       for (i in 2:a) {
                r <- rank(c(control, data[group == levels(factor(group))[i]]))     
                R <- sum(r[1:n1])                               
                N <- n1+ni[i]                                   
                E <- n1*(N+1)/2                                 
                V <- n1*ni[i]/N/(N-1)*(sum(r^2)-N*(N+1)^2/4)    
                t[i] <- abs(R-E)/sqrt(V)                        
			p.value[i] <-  pdunnett(t[i], a, df=0, r=rho)
       }
        result <- cbind(t, rho, p.value)[-1,]                   
        rownames(result) <- paste(levels(factor(group))[1], levels(factor(group))[2:a], sep=":")
        return(result)
}


pdunnett <- function(x, a, df, r) {        
# Used in Steel(). Originated from Dunnet()
                corr <- diag(a-1)
                corr[lower.tri(corr)] <- r
                1-pmvt(lower=-x, upper=x, delta=numeric(a-1), df=df, corr=corr, abseps=0.0001)
}

		
RemoveOutlier <- function(x, return){
	i <- 0
	repeat{
		x1 <- x[!is.na(x)]
	      n <- length(x1)
      	if(max(x1)-mean(x1)>=mean(x1)-min(x1)){
			p.value <- n*pt(sqrt((n-2)/((n-1)^2/((max(x1)-mean(x1))/sd(x1))^2/n-1)), n-2, lower.tail=FALSE)
			if(p.value < 0.05) {
				cat(gettext(domain="R-RcmdrPlugin.EZR","Identify data"), " ", max(x1), " ", gettext(domain="R-RcmdrPlugin.EZR","as an outlier. (Smirnov-Grubbs p-value="), p.value, ")\n", sep="")
				x[x==max(x1)] <- NA
				i <- i + 1
				}
		} else {
			p.value <- n*pt(sqrt((n-2)/((n-1)^2/((mean(x1)-min(x1))/sd(x1))^2/n-1)), n-2, lower.tail=FALSE)
			if(p.value < 0.05) {
				cat(gettext(domain="R-RcmdrPlugin.EZR","Identify data"), min(x1), gettext(domain="R-RcmdrPlugin.EZR","as an outlier. (Smirnov-Grubbs p-value="), p.value, ")\n", sep="")
				x[x==min(x1)] <- NA
				i <- i + 1
			}
		}
		if(p.value >= 0.05) break		
	}
	if (i==0) cat(gettext(domain="R-RcmdrPlugin.EZR","No outliers were identified."), "\n")
	if (return==1) return(x)
}

	
summary.table.twoway <- function(object, ..., table, res){
	p.value <- signif(res$p.value, digits=3)
	summary.table <- data.frame(cbind(table, p.value))
	groups1 <- length(levels(factor(data.frame(table)[,1])))
	groups2 <- length(levels(factor(data.frame(table)[,2])))
	for (i in 1:groups1){
		rownames(summary.table)[i] <- paste(names(data.frame(table))[1], "=", levels(factor(data.frame(table)[,1]))[i], sep="") 		
		if (i >=2) summary.table$p.value[i] <- ""
	}
	for (i in 1:groups2){	
		colnames(summary.table)[i] <- paste(names(data.frame(table))[2], "=", levels(factor(data.frame(table)[,2]))[i], sep="") 
	}
	if(res$method=="Fisher's Exact Test for Count Data"){
		colnames(summary.table)[length(summary.table)] <- "Fisher.p.value"
	} else {
		colnames(summary.table)[length(summary.table)] <- "Chisq.p.value"	
	}	
	return(summary.table)
}	


summary.table.MH <- function(object, ..., table, res){
	MH.p.value <- signif(res$p.value, digits=3)
	summary.table <- data.frame(cbind(table, MH.p.value))
	groups1 <- length(levels(factor(data.frame(table)[,1])))
	groups2 <- length(levels(factor(data.frame(table)[,2])))
	for (i in 1:groups1){
		rownames(summary.table)[i] <- paste(names(data.frame(table))[1], "=", levels(factor(data.frame(table)[,1]))[i], sep="") 		
		if (i >=2) summary.table$MH.p.value[i] <- ""
	}
	for (i in 1:groups2){	
		colnames(summary.table)[i] <- paste(names(data.frame(table))[2], "=", levels(factor(data.frame(table)[,2]))[i], sep="") 
	}
	return(summary.table)
}	

		
summary.km <- function (object, ..., survfit, survdiff=NULL, time=0){
	km <- survfit
	km.table <- summary(survfit)
	if (is.null(survdiff)){
		p.value <- NULL
	}else{
		p.value <- signif(pchisq(c(survdiff$chisq), df=length(survdiff$n)-1, lower.tail=FALSE),digits=3)
	}
	if (is.null(survdiff)){
		groups <- 1
		samples <- km.table$table[1]
		medians <- km.table$table[5]
		med.ci <- paste(km.table$table[6], "-", km.table$table[7], sep="")
		km$strata[1] <- samples
	}else{
		group.names <- row.names(km.table$table)
		groups <- length(group.names)
		samples <- km.table$table[,1]
		medians <- km.table$table[,5]
		med.ci <- paste(km.table$table[,6], "-", km.table$table[,7], sep="")
	}
	surv <- NULL
	surv.ci <- NULL	
	if (time > 0){	# show survival rate at time
		start <- 1
		for(i in 1:groups){
			numbers <- km$strata[i]
			stop <- start + numbers - 1
			timetoevent <- km$time[start:stop]
			if (max(timetoevent, na.rm=TRUE) >= time){
				point <- max((1:length(timetoevent))[timetoevent<=time], na.rm=TRUE)
				surv[i] <- formatC(km$surv[start+point-1], format="f", digits=3)
				surv.ci[i] <- paste("(", formatC(km$lower[start+point-1], format="f", digits=3), "-", formatC(km$upper[start+point-1], format="f", digits=3), ")", sep="")
			}else{
				surv[i] <- NA
				surv.ci[i] <- NA
			}
			start <- stop + 1
		}
	}
	if(groups==1){
		surv.table <- data.frame(t(c(samples, probability=surv, CI=surv.ci, median=medians, medianCI=med.ci)))
#		colnames(surv.table)[1:3] <- c("n", "median survival", "95% CI")
#		if(length(surv.table)==5){
#			colnames(surv.table)[1:5] <- c("n", "survival rate", "95% CI", "median survival", "95% CI")
#		}
		}else{
		for(i in 2:groups){
			p.value[i] <- ""
		}
		if(!is.null(surv.ci)){	
			surv.table <- data.frame(n=samples, probability=surv, CI=surv.ci, median=medians, medianCI=med.ci, p.value)
#		colnames(surv.table)[1:3] <- c("n", "median survival", "95% CI")
#		if(length(surv.table)==5){
#			colnames(surv.table)[1:5] <- c("n", "survival rate", "95% CI", "median survival", "95% CI")
#		}
		}else{
			surv.table <- data.frame(n=samples, median=medians, medianCI=med.ci, p.value)
#			colnames(surv.table)[1:3] <- c("n", "median survival", "95% CI")
		}
	}
	colnames(surv.table)[1:3] <- c("n", "median survival", "95% CI")
	if(length(surv.table)>=5){
			colnames(surv.table)[1:5] <- c("n", "survival rate", "95% CI", "median survival", "95% CI")
	}
	colnames(surv.table) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(surv.table))
	return(surv.table)
}


summary.ci <- function (object, ..., ci, res, event=1, time=0){
	ci.table <- summary(ci)
	if(is.null(ci$strata)){
		ngroups <- 1
		p.value <- NULL
	} else {
		groups <- levels(ci.table$strata)
		ngroups <- length(groups)
		p.value <- signif(res$Tests[event, 2],digits=3)	
	}

	nevents <- length(ci$surv[1,])
	samples <- ci.table$table[as.numeric(substring(row.names(ci.table$table),1,1))==event,1]
	medians <- ci.table$table[as.numeric(substring(row.names(ci.table$table),1,1))==event,5]
	surv <- NULL
	surv.ci <- NULL	
	if (time > 0){	# show survival rate at time
		for(i in 1:ngroups){
			survival <- timepoints(res, time)$est[ngroups*(event-1)+i]
			hazard <- log(survival)
			se <- sqrt(timepoints(res, time)$var[ngroups*(event-1)+i])			
			lower <- survival^exp(-qnorm(0.975)*se/(survival*hazard))	#log-log
#			lower <- survival*exp(-qnorm(0.975)*se/(survival))			#log
			if(is.nan(lower)) lower<-0
			if(is.na(lower)==FALSE & lower>1) lower<-1
			upper <- survival^exp(qnorm(0.975)*se/(survival*hazard))	#log-lgo
#			upper <- survival*exp(qnorm(0.975)*se/(survival))			#log
			if(is.nan(upper)) upper<-0
			if(is.na(upper)==FALSE & upper>1) upper<-1
			surv[i] <- formatC(survival, format="f", digits=3)
			surv.ci[i] <- paste("(", formatC(lower, format="f", digits=3), "-", formatC(upper, format="f", digits=3), ")", sep="")
		}
	}
	if(ngroups==1){
		if(!is.null(surv)){	
			surv.table <- data.frame(n=samples, incidence=surv, CI=surv.ci, median=medians)
			colnames(surv.table)[3] <- "95% CI"
		}else{
			surv.table <- data.frame(n=samples, median=medians)
			colnames(surv.table)[2] <- "median time"
		}
	}else{
		p.value[2:ngroups] <- ""
		if(!is.null(surv)){	
			surv.table <- data.frame(n=samples, incidence=surv, CI=surv.ci, median=medians, p.value)
			colnames(surv.table)[3] <- "95% CI"
			colnames(surv.table)[4] <- "median time"
		}else{
			surv.table <- data.frame(n=samples, median=medians, p.value)
			colnames(surv.table)[2] <- "median time"
		}
	}
	colnames(surv.table) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(surv.table))
	return(surv.table)
}


print.ci.summary <- function (x, ..., ci, res) {
    ngroups <- length(ci$n)
    group.names <- names(ci$strata)
    if (is.null(ci$surv)) 
        ci$surv <- 1 - ci$prev
    nevents <- length(ci$surv[1, ])	
#event column with no event left here in the new survival package
#event column with no event should be deleted
    zerocolumn <- NA
    for (i in 1:nevents) {
        zerocolumn[i] <- ifelse(sum(1-ci$surv[,i])==0, 0, 1)
    }
    ci$surv <- ci$surv[,zerocolumn==1]
    nevents <- sum(zerocolumn)

    start <- 1
    for (i in 1:ngroups) {
        if (ngroups == 1) {
            stop <- start + length(ci$time) - 1
        } else {
            stop <- start + ci$strata[i] - 1
        }
        ci.summary.table <- data.frame(time = ci$time[start:stop], 
            n.risk = ci$n.risk[start:stop], n.event = ci$n.event[start:stop])
        for (j in 1:nevents) {
		if (nevents >1){
	            ci.summary.table <- cbind(ci.summary.table, incidence = round(1 - ci$surv[start:stop, j], 3))
		} else {
	            ci.summary.table <- cbind(ci.summary.table, incidence = round(1 - ci$surv[start:stop], 3))
		}
            time <- res[[ngroups * (j - 1) + i]]$time
            var <- res[[ngroups * (j - 1) + i]]$var
            ci95 <- NULL
            for (k in start:stop) {
                point <- max((1:length(time))[time <= ci$time[k]])
   		    if (nevents >1){
                	est <- 1 - ci$surv[k, j]
		    } else {
                	est <- 1 - ci$surv[k]
		    }
                variance <- var[point]
                se <- sqrt(variance)
                hazard <- log(est)
                lower <- est^exp(-qnorm(0.975) * se/(est * hazard))
                if (is.nan(lower)) 
                  lower <- 0
                if (is.na(lower) == FALSE & lower > 1) 
                  lower <- 1
                upper <- est^exp(qnorm(0.975) * se/(est * hazard))
                if (is.nan(upper)) 
                  upper <- 0
                if (is.na(upper) == FALSE & upper > 1) 
                  upper <- 1
                ci95[k - start + 1] <- paste("(", formatC(lower, 
                  format = "f", digits = 3), "-", formatC(upper, 
                  format = "f", digits = 3), ")", sep = "")
            }
            ci.summary.table <- cbind(ci.summary.table, ci95)
            colnames(ci.summary.table)[2 + j * 2] <- paste("incidence-", 
                j, sep = "")
            colnames(ci.summary.table)[3 + j * 2] <- paste("95% CI-", 
                j, sep = "")
        }
        cat("\t\t", names(ci$strata[i]), "\n")
        print(ci.summary.table)
        cat("\n")
        start <- stop + 1
    }
}


StatMedTableOne  <- function(){
    Library("tableone")
	defaults <- list(group=NULL, cat=NULL, cont=NULL, contnonnormal=NULL, exact="auto", range="TRUE", explain="FALSE", output="clipboard", language="1", subset = "")
	dialog.values <- getDialog("StatMedTableOne", defaults)
	currentFields$subset <- dialog.values$subset	
	currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Summary table of sample characteristics"))
    groupBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), listHeight=10, initialSelection=varPosn(dialog.values$group, "all"))
    variableFrame <- tkframe(top)
    categoryBox <- variableListBox(variableFrame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Categorical variables"), listHeight=10, initialSelection=varPosn(dialog.values$cat, "all"))
    contBox <- variableListBox(variableFrame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Continuous variables (normal distribution)"), listHeight=10, initialSelection=varPosn(dialog.values$cont, "all"))
    contnonnormalBox <- variableListBox(variableFrame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Continuous variables (non-normal distribution)"), listHeight=10, initialSelection=varPosn(dialog.values$contnonnormal, "all"))
    optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="exact", buttons=c("chisq", "fisher", "auto"), values=c("chisq", "exact", "auto"),
        initialValue=dialog.values$exact, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Chi-square test with continuity correction", "Fisher's exact test",  "Automatic selection")), title=gettext(domain="R-RcmdrPlugin.EZR","Test for categorical variables"))
    radioButtons(optionsFrame, name="range", buttons=c("MinMax", "IQR"), values=c("TRUE", "FALSE"),
        initialValue=dialog.values$range, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Minimum and maximum values", "Interquartile ranges")), title=gettext(domain="R-RcmdrPlugin.EZR","Range for non-normal categorical variables"))			
    radioButtons(optionsFrame, name="explain", buttons=c("No", "Yes"), values=c("FALSE", "TRUE"),
        initialValue=dialog.values$explain, labels=gettext(domain="R-RcmdrPlugin.EZR",c("No", "Yes")), title=gettext(domain="R-RcmdrPlugin.EZR","Show explantation for continuous variables"))			
    options2Frame <- tkframe(top)
    radioButtons(options2Frame, name="output", buttons=c("Clipboard", "CSVfile"), values=c("clipboard", "CSVfile"),
        initialValue=dialog.values$output, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Clipboard", "CSV file")), title=gettext(domain="R-RcmdrPlugin.EZR","Output destination"))			
    radioButtons(options2Frame, name="language", buttons=c("Eng", "Local"), values=c("1", "0"),
        initialValue=dialog.values$language, labels=gettext(domain="R-RcmdrPlugin.EZR",c("English", "Local")), title=gettext(domain="R-RcmdrPlugin.EZR","Language"))			
	StatMedSubsetBox(model=TRUE)

	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Summary table of sample characteristics"), "#####", sep=""))
    group <- getSelection(groupBox)
    cat <- getSelection(categoryBox)
    cont <- getSelection(contBox)
    contnonnormal <- getSelection(contnonnormalBox)
	exact <- tclvalue(exactVariable)	
	range <- tclvalue(rangeVariable)	
	explain <- tclvalue(explainVariable)	
    output <- tclvalue(outputVariable)
    language <- tclvalue(languageVariable)
	dataSet <- activeDataSet()
    subset <- tclvalue(subsetVariable)

	putDialog("StatMedTableOne", list(group=group, cat=cat, cont=cont, contnonnormal=contnonnormal, exact=exact, range=range, explain=explain, output=output, language=language, subset = tclvalue(subsetVariable)))
		if(output=="Screen") output <- ""		
		if(output=="CSVfile") {
			output <- tclvalue(tkgetSaveFile(filetypes=
				gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"Text Files" {".txt" ".TXT" ".csv" ".CSV"}}'),
				defaultextension="csv", initialfile=paste("tableone.csv", sep=".")))
			if (output == "") return()
		}	
		if (.Platform$OS.type != 'windows' & output=="clipboard"){
            errorCondition(recall=StatMedTableOne, message=gettext(domain="R-RcmdrPlugin.EZR","Clipboard can be selected only in Windows."))
            return()			
		}

    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
		subdataSet <- dataSet
		}
    else{
		subdataSet <- paste("subset(", dataSet, ", ", subset, ")", sep="")
	}
	if (length(group==1)){
		levels <- eval(parse(text=paste("length(levels(factor(", subdataSet, "$", group, ")))", sep="")))	
	}
	if (exact=="auto" & length(group)==0) exact <- "exact"
	if (exact=="auto" & length(group)==1){
		if (levels>=3){
			exact <- "chisq"
		} else{
			exact <- "exact"
		}
	}	
    if (length(cat)+length(cont)+length(contnonnormal)==0){
            errorCondition(recall=StatMedTableOne, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable"))
            return()
	}

    closeDialog()
	
#	doItAndPrint("matCatTable <- NULL")
#	doItAndPrint("matContTable <- NULL")
#	doItAndPrint("matContnonnormalTable <- NULL")
	if(length(cat)>0){
		catVariables <- paste('c("', cat[1], '"', sep="")
		if(length(cat)>1){
			for (i in 2:length(cat)) {
				catVariables <- paste(catVariables, ', "', cat[i], '"', sep="")
			}
		}
		catVariables <- paste(catVariables, ")", sep="")
		if (length(group) == 0){
			doItAndPrint(paste("CatTable <- CreateCatTable(vars = ", catVariables, ', data=', subdataSet, ")", sep=""))			
		}else if(length(group)==1 & levels<2){
			doItAndPrint(paste("CatTable <- CreateCatTable(vars = ", catVariables, ', data=', subdataSet, ")", sep=""))			
		}else{
			doItAndPrint(paste("CatTable <- CreateCatTable(vars = ", catVariables, ', strata="', group, '", data=', subdataSet, ")", sep=""))
		}
		if (exact=="chisq"){
			doItAndPrint("matCatTable <- print(CatTable, printToggle = FALSE, showAllLevels = TRUE)")
		} else if (exact=="exact"){
			doItAndPrint(paste("matCatTable <- print(CatTable, printToggle = FALSE, showAllLevels = TRUE, exact=", catVariables, ")", sep=""))		
		}	
#		doItAndPrint("matCatTable <- data.frame(matCatTable)")
		doItAndPrint('if(colnames(matCatTable)[length(colnames(matCatTable))]=="test"){\nmatCatTable<-matCatTable[,1:length(colnames(matCatTable))-1]\n}')
		doItAndPrint("matCatTable <- cbind(Factor=row.names(matCatTable), matCatTable)")
	}
	if(length(cont)>0){
		contVariables <- paste('c("', cont[1], '"', sep="")
		if(length(cont)>1){
			for (i in 2:length(cont)) {
				contVariables <- paste(contVariables, ', "', cont[i], '"', sep="")
			}
		}
		contVariables <- paste(contVariables, ")", sep="")
		if (length(group) == 0){
			doItAndPrint(paste("ContTable <- CreateContTable(vars = ", contVariables, ', data=', subdataSet, ")", sep=""))			
		}else if(length(group)==1 & levels<2){
			doItAndPrint(paste("ContTable <- CreateContTable(vars = ", contVariables, ', data=', subdataSet, ")", sep=""))			
		}else{
			doItAndPrint(paste("ContTable <- CreateContTable(vars = ", contVariables, ', strata="', group, '", data=', subdataSet, ")", sep=""))
		}
		doItAndPrint(paste("matContTable <- print(ContTable, printToggle = FALSE, explain = ", explain, ")", sep=""))
#		doItAndPrint("matContTable <- data.frame(matContTable)")
		doItAndPrint('if(colnames(matContTable)[length(colnames(matContTable))]=="test"){\nmatContTable<-matContTable[,1:length(colnames(matContTable))-1]\n}')
		#Add a dummy column to ContTable, because CatTable has a grouping column
		if(length(cat)>0) doItAndPrint('matContTable <- cbind(level="", matContTable)')
		if(language==0 & explain=="TRUE") {
			doItAndPrint('row.names(matContTable)[2:length(row.names(matContTable))] <- paste(substring(row.names(matContTable)[2:length(row.names(matContTable))], 1, nchar(row.names(matContTable)[2:length(row.names(matContTable))])-11), gettext(domain="R-RcmdrPlugin.EZR", "(mean (sd))"), sep="")')
		}
		doItAndPrint("matContTable <- cbind(Factor=row.names(matContTable), matContTable)")
#		if(length(cat)>0) doItAndPrint("matContTable <- matContTable[2:length(rownames(matContTable)),]")
	}
	if(length(contnonnormal)>0){
		contnonnormalVariables <- paste('c("', contnonnormal[1], '"', sep="")
		if(length(contnonnormal)>1){
			for (i in 2:length(contnonnormal)) {
				contnonnormalVariables <- paste(contnonnormalVariables, ', "', contnonnormal[i], '"', sep="")
			}
		}
		contnonnormalVariables <- paste(contnonnormalVariables, ")", sep="")
		if (length(group) == 0) {
			doItAndPrint(paste("ContnonnormalTable <- CreateContTable(vars = ", contnonnormalVariables, ', data=', subdataSet, ")", sep=""))			
		}else if(length(group)==1 & levels<2){
			doItAndPrint(paste("ContnonnormalTable <- CreateContTable(vars = ", contnonnormalVariables, ', data=', subdataSet, ")", sep=""))						
		}else{
			doItAndPrint(paste("ContnonnormalTable <- CreateContTable(vars = ", contnonnormalVariables, ', strata="', group, '", data=', subdataSet, ")", sep=""))
		}
		doItAndPrint(paste("matContnonnormalTable <- print(ContnonnormalTable, printToggle = FALSE, nonnormal = TRUE, explain = ", explain, ", minMax=", range, ")", sep=""))
#		doItAndPrint("matContnonnormalTable <- data.frame(matContnonnormalTable)")
		doItAndPrint('if(colnames(matContnonnormalTable)[length(colnames(matContnonnormalTable))]=="test"){\nmatContnonnormalTable<-matContnonnormalTable[,1:length(colnames(matContnonnormalTable))-1]\n}')
		#Add a dummy column to ContTable, because CatTable has a grouping column
		if(length(cat)>0) doItAndPrint('matContnonnormalTable <- cbind(level="", matContnonnormalTable)')
		if(language==0 & explain=="TRUE") {
			if(range=="TRUE"){
				doItAndPrint('row.names(matContnonnormalTable)[2:length(row.names(matContnonnormalTable))] <- paste(substring(row.names(matContnonnormalTable)[2:length(row.names(matContnonnormalTable))], 1, nchar(row.names(matContnonnormalTable)[2:length(row.names(matContnonnormalTable))])-16), gettext(domain="R-RcmdrPlugin.EZR", "(median [range])"), sep="")')
			} else {
				doItAndPrint('row.names(matContnonnormalTable)[2:length(row.names(matContnonnormalTable))] <- paste(substring(row.names(matContnonnormalTable)[2:length(row.names(matContnonnormalTable))], 1, nchar(row.names(matContnonnormalTable)[2:length(row.names(matContnonnormalTable))])-14), gettext(domain="R-RcmdrPlugin.EZR", "(median [IQR])"), sep="")')			
			}
		}
		doItAndPrint("matContnonnormalTable <- cbind(Factor=row.names(matContnonnormalTable), matContnonnormalTable)")
#		if(length(cat)>0 | length(cont)>0) doItAndPrint("matContnonnormalTable <- matContnonnormalTable[2:length(rownames(matContnonnormalTable)),]")
	}
	if(length(cat)>0){
		doItAndPrint("FinalTable <- as.matrix(matCatTable)")
		ncol <- eval(parse(text=paste("length(colnames(FinalTable))")))	
		doItAndPrint("tempStrata <- attributes(FinalTable)[[2]][2]")
			if(length(cont>0)){
#				doItAndPrint(paste("FinalTable <- rbind(FinalTable, matrix(matContTable, ncol=", ncol, "))", sep=""))
#				doItAndPrint("FinalTable <- rbind(FinalTable, matContTable[2:length(rownames(matContTable)),])")
				doItAndPrint("FinalTable <- rbind(FinalTable, matContTable)")
			}
			if(length(contnonnormal>0)){
#				doItAndPrint(paste("FinalTable <- rbind(FinalTable, matrix(matContnonnormalTable, ncol=", ncol, "))", sep=""))				
#				doItAndPrint("FinalTable <- rbind(FinalTable, matContnonnormalTable[2:length(rownames(matContnonnormalTable)),])")
				doItAndPrint("FinalTable <- rbind(FinalTable, matContnonnormalTable)")
			}
	}
	if(length(cat)==0 & length(cont)>0){
		doItAndPrint("FinalTable <- as.matrix(matContTable)")
		ncol <- eval(parse(text=paste("length(colnames(FinalTable))")))	
		doItAndPrint("tempStrata <- attributes(FinalTable)[[2]][2]")
			if(length(contnonnormal>0)){
#				doItAndPrint(paste("FinalTable <- rbind(FinalTable, matrix(matContnonnormalTable, ncol=", ncol, "))", sep=""))				
#				doItAndPrint("FinalTable <- rbind(FinalTable, matContnonnormalTable[2:length(rownames(matContnonnormalTable)),])")
				doItAndPrint("FinalTable <- rbind(FinalTable, matContnonnormalTable)")
			}
	}
	if(length(cat)==0 & length(cont)==0 & length(contnonnormal>0)){
		doItAndPrint("FinalTable <- as.matrix(matContnonnormalTable)")
		doItAndPrint("tempStrata <- attributes(FinalTable)[[2]][2]")
	}	
	doItAndPrint("attributes(FinalTable) <- c(list(dim=attributes(FinalTable)[[1]]), list(dimnames=c(attributes(FinalTable)[[2]][1], tempStrata)))")
	if(length(cat)>0) doItAndPrint('colnames(FinalTable)[2] <- "Group"')
	if(length(group)==1) {if (levels>1) doItAndPrint('colnames(FinalTable)[length(colnames(FinalTable))] <- "p.value"')}
#	doItAndPrint("print(as.matrix(FinalTable), quote=FALSE)")
#	doItAndPrint("FinalTable <- cbind(Factor=row.names(FinalTable), FinalTable)")
	if(language==0) {
		doItAndPrint('colnames(FinalTable) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(FinalTable))')
#		doItAndPrint('colnames(FinalTable)[1] <- gettext(domain="R-RcmdrPlugin.EZR", "Factor")')
#		doItAndPrint('if(colnames(FinalTable)[2] == "Group") colnames(FinalTable)[2] <- gettext(domain="R-RcmdrPlugin.EZR", "Group")')		
#		if(length(group)==1) {if (levels>1) doItAndPrint('colnames(FinalTable)[length(colnames(FinalTable))] <- gettext(domain="R-RcmdrPlugin.EZR","p.value")')}
	}
	doItAndPrint("row0 <- colnames(FinalTable)")
	doItAndPrint("row1 <- FinalTable[1,]")
	doItAndPrint("row1 <- matrix(row1, nrow=1)")
	doItAndPrint("colnames(row1) <- row0")
	doItAndPrint('FinalTable <- FinalTable[which(rownames(FinalTable)!="n"),]')
	doItAndPrint("FinalTable <- rbind(n=row1, FinalTable)")
#	doItAndPrint("row.names(FinalTable) <- NULL")

	if(length(cat)>0){
#		doItAndPrint("print(FinalTable[,2:length(FinalTable[1,])], quote=FALSE)")
	} else if (length(group)==1){
		if(levels>1){
			if(length(cont)==1 && length(contnonnormal)==0){
				doItAndPrint(paste('rownames(FinalTable) <- c("n", "', cont[1], '")',sep="")) 
			}
			if(length(cont)==0 && length(contnonnormal)==1){
				doItAndPrint(paste('rownames(FinalTable) <- c("n", "', contnonnormal[1], '")',sep="")) 
			}
#			doItAndPrint("print(FinalTable[,2:length(FinalTable[1,])], quote=FALSE)")
		} else {		
			doItAndPrint('rownames(FinalTable) <- rep("", length(rownames(FinalTable)))')
#			doItAndPrint("print(FinalTable, quote=F)")
		}
	} else {
			doItAndPrint('rownames(FinalTable) <- rep("", length(rownames(FinalTable)))')
#			doItAndPrint("print(FinalTable, quote=F)")
	}	
	#	doItAndPrint("FinalTable <- cbind(Factor=row.names(FinalTable), FinalTable)")
	doItAndPrint("FinalTable <- rbind(row0, FinalTable)")
	if(length(group)==1) {
		if (levels>1) {
			doItAndPrint('row0 <- rep("", length(colnames(FinalTable)))')
			if(length(cat)==0){
				doItAndPrint(paste('row0[2] <- "', group, '"', sep=""))
			}else{
				doItAndPrint(paste('row0[3] <- "', group, '"', sep=""))			
			}
			doItAndPrint("FinalTable <- rbind(row0, FinalTable)")
		}
	}
	doItAndPrint("finaltable_dataframe_print(FinalTable)")
	if (output=="clipboard"){
		doItAndPrint('write.table(FinalTable, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)')
	} else {
		doItAndPrint(paste('write.table(FinalTable, file="', output, '", sep=",", row.names=FALSE, col.names=FALSE)', sep=""))	
	}
	tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="tableone", apply="StatMedTableOne", reset="StatMedTableOne")
    tkgrid(getFrame(groupBox), sticky="nw")
	tkgrid(labelRcmdr(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
    tkgrid(getFrame(categoryBox), labelRcmdr(variableFrame, text="    "), getFrame(contBox), labelRcmdr(variableFrame, text="    "), getFrame(contnonnormalBox), sticky="nw")
    tkgrid(variableFrame, sticky="nw")
	tkgrid(exactFrame, labelRcmdr(optionsFrame, text="   "), rangeFrame, labelRcmdr(optionsFrame, text="   "), explainFrame, sticky="nw")	
	tkgrid(optionsFrame, sticky="w")
	tkgrid(outputFrame, labelRcmdr(options2Frame, text="   "), languageFrame, sticky="nw")	
	tkgrid(options2Frame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Clipboard can be selected only in Windows."), fg="blue"), sticky="w")
	tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


objectCheck <- function(name, obj){
	#obj <- objects() should be performed before executing this function.
	#Used in StatMedSummaryResults
	present <- 0
	for(i in 1:length(obj)){
		if (name==obj[i]) present <- 1
	}
#	if (present==0) print(paste("Object ", name, " was not found.", sep=""))
	if (present==0) print(gettext(domain="R-RcmdrPlugin.EZR","You must perform analysis before outputting."))
	return(present)
}


StatMedSummaryResults <- function() {
	defaults <- list(analysis="twoway", output="clipboard", language="1")
	dialog.values <- getDialog("StatMedSummaryResults", defaults)
	currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Summary table of results"))
    optionsFrame <- tkframe(top)
	radioButtons(optionsFrame, name="analysis", buttons=c("twoway", "ttest", "survival", "ci", "logistic", "multireg", "cox", "finegray"), values=c("twoway", "ttest", "survival", "ci", "logistic", "multireg", "cox", "finegray"), initialValue=dialog.values$analysis, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-way table", "T-test",  "Survival test", "Cumulative incidence", "Multivariate logistic regression", "Multivariate linear regression", "Proportional hazard regression", "Fine-Gray regression")), title=gettext(domain="R-RcmdrPlugin.EZR","Test for outputting result"))
    radioButtons(optionsFrame, name="output", buttons=c("Clipboard", "CSVfile"), values=c("clipboard", "CSVfile"),
        initialValue=dialog.values$output, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Clipboard", "CSV file")), title=gettext(domain="R-RcmdrPlugin.EZR","Output destination"))			
    radioButtons(optionsFrame, name="language", buttons=c("Eng", "Local"), values=c("1", "0"),
        initialValue=dialog.values$language, labels=gettext(domain="R-RcmdrPlugin.EZR",c("English", "Local")), title=gettext(domain="R-RcmdrPlugin.EZR","Language"))			
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Summary table of results"), "#####", sep=""))
        analysis <- tclvalue(analysisVariable)
        output <- tclvalue(outputVariable)
        language <- tclvalue(languageVariable)
		putDialog("StatMedSummaryResults", list(analysis=analysis, output=output, language=language))
		if(analysis=="twoway") table <- "Fisher.summary.table"
		if(analysis=="ttest") table <- "summary.ttest"
		if(analysis=="survival") table <- "km.summary.table"
		if(analysis=="ci") table <- "ci.summary.table"
		if(analysis=="logistic") table <- "odds"
		if(analysis=="multireg") table <- "multireg.table"
		if(analysis=="cox") table <- "cox.table"
		if(analysis=="finegray") table <- "crr.table"		
		if(output=="Screen") output <- ""
		if(output=="CSVfile") {
			output <- tclvalue(tkgetSaveFile(filetypes=
				gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"Text Files" {".txt" ".TXT" ".csv" ".CSV"}}'),
				defaultextension="csv", initialfile=paste(table, "csv", sep=".")))
			if (output == "") return()
		}	
		if (.Platform$OS.type != 'windows' & output=="clipboard"){
            errorCondition(recall=StatMedSummaryResults, message=gettext(domain="R-RcmdrPlugin.EZR","Clipboard can be selected only in Windows."))
            return()			
		}
#		findobject <- eval(parse(text=paste('objectCheck("', table, '", objects())', sep="")))
#		doItAndPrint(paste('findobject <- objectCheck("', table, '", objects())', sep=""))
#		if(findobject==0){
#            errorCondition(recall=StatMedSummaryResults, message=gettext(domain="R-RcmdrPlugin.EZR","You must perform analysis before outputting."))
#            return()		
#		}
		if(analysis=="twoway") doItAndPrint(paste('if(objectCheck("Fisher.summary.table", objects())) w.twoway(Fisher.summary.table, filename="', output, '", en=', language, ")", sep=""))
		if(analysis=="ttest") doItAndPrint(paste('if(objectCheck("summary.ttest", objects())) w.ttest(summary.ttest, filename="', output, '", en=', language, ")", sep=""))
		if(analysis=="survival") doItAndPrint(paste('if(objectCheck("km.summary.table", objects())) w.survival(km.summary.table, filename = "', output, '", en=', language, ")", sep=""))
		if(analysis=="ci") doItAndPrint(paste('if(objectCheck("ci.summary.table", objects())) w.ci(ci.summary.table, filename = "', output, '", en=', language, ")", sep=""))
		if(analysis=="logistic") doItAndPrint(paste('if(objectCheck("odds", objects())) w.multi(odds, filename = "', output, '", en=', language, ")", sep=""))
		if(analysis=="multireg") doItAndPrint(paste('if(objectCheck("multireg.table", objects())) w.multireg(multireg.table, filename = "', output, '", en=', language, ")", sep=""))
		if(analysis=="cox") doItAndPrint(paste('if(objectCheck("cox.table", objects())) w.multi(cox.table, filename = "', output, '", en=', language, ")", sep=""))
		if(analysis=="finegray") doItAndPrint(paste('if(objectCheck("crr.table", objects())) w.multi(crr.table, filename = "', output, '", en=', language, ")", sep=""))			
        closeDialog()
        }
    OKCancelHelp(helpSubject="w.multi")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Clipboard can be selected only in Windows."), fg="blue"), sticky="w")
	tkgrid(analysisFrame, labelRcmdr(optionsFrame, text="   "), outputFrame, labelRcmdr(optionsFrame, text="   "), languageFrame, sticky="nw")	
	tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=7, columns=2)
}


stsplit <- function (dataframe, timetoevent, event, timeon, covariate, timeoff){
	Temp1 <- dataframe
	PatientsNumber <- length(Temp1[,1])
	Temp1$start_td <- 0
	Temp1$stop_td <- timetoevent
	Temp1$endpoint_td <- event
	Temp1$covariate_td <- covariate
	timeon <- ifelse(timeon<0, 0, timeon)
	timeoff <- ifelse(timeoff<0, 0, timeoff)
	for (i in 1:PatientsNumber){
		Temp1$patientsnumber_td[i] <- i + 0.1
		if ( is.na(timetoevent[i]) == TRUE || is.na(timeon[i])==TRUE || is.na(timeoff[i]) == TRUE || is.na(covariate[i]) == TRUE){
			Temp1$covariate_td[i] <- NA
		}
		else {
		if (covariate[i] == 1 && timetoevent[i] > timeon[i]) {
			Temp1$stop_td[i] <- timeon[i]
			Temp2 <- Temp1[i,]
			Temp2$start_td[1] <- Temp1$stop_td[i]
			Temp2$stop_td[1] <- timetoevent[i]
			Temp2$patientsnumber_td[1] <- i + 0.2
			if ( timeoff[i] > timeon[i] && timetoevent[i] > timeoff[i]){
				Temp3 <- Temp2
				Temp2$stop_td[1] <- timeoff[i]
				Temp3$start_td[1] <- Temp2$stop_td[1]
				Temp3$stop_td[1] <- timetoevent[i]
				Temp3$covariate_td[1] <- 0
				Temp3$patientsnumber_td[1] <- i + 0.3
				Temp2$endpoint_td[1] <- 0
				Temp1<- rbind(Temp1, Temp3)
			}
			Temp1<- rbind(Temp1, Temp2)
			Temp1$endpoint_td[i] <- 0
			}
			Temp1$covariate_td[i] <- 0
		}
	}
	sortlist <- order(Temp1$patientsnumber_td)
	return (Temp1[sortlist,])
}


Mantel.Byar <- function(Group=TempTD$covariate_td, Event=TempTD$endpoint_td, StartTime=TempTD$start_td,	StopTime=TempTD$stop_td, method=c("SAS", "Tominaga"), plot=0, landmark=0) {
	#modified from logrank test in http://aoki2.si.gunma-u.ac.jp/R/logrank.html
	#Reuire TempTD dataset created by Cox with TD variable in EZR
	method <- match.arg(method)
	data.name <- sprintf("StartTime: %s, StopTime: %s, Event: %s, Group: %s",
                deparse(substitute(StartTime)),
	deparse(substitute(StopTime)), deparse(substitute(Event)),
	deparse(substitute(Group)))
	OK <- complete.cases(Group, Event, StartTime, StopTime)
	Group <- Group[OK]
	Event <- Event[OK]
	StartTime <- StartTime[OK]
	StopTime <- StopTime[OK]

	Start <- pmin(StartTime, StopTime)					#for samples with StartTime>StopTime
	Stop <- pmax(StartTime, StopTime)
	StartTime <- Start
	StopTime <- Stop

	len <- length(Group)
	stopifnot(length(Event) == len, length(StopTime) == len)

	tg <- table(c(StopTime, rep(NA, 4)),                        
                    c(Group, 1, 1, 2, 2)*10+c(Event, 1, 0, 1, 0))
	k <- nrow(tg)
	nia <- table(Group)[1]
	nib <- len-nia
	na <- c(nia, (rep(nia, k)-cumsum(tg[,1]+tg[,2]))[-k])	
	nb <- c(nib, (rep(nib, k)-cumsum(tg[,3]+tg[,4]))[-k])	
											#following part is different from log-rank test
	minus <- NULL							
	for (i in 1:length(tg[,1])){
		if(as.integer(rownames(tg))[i]==0){
			minus[i] <- sum((as.integer(rownames(tg))[i] < StartTime))
		} else {
			minus[i] <- sum((as.integer(rownames(tg))[i] <= StartTime))
		}
	}
	nb <- nb - minus						
											#Following part is same wtih log-ranktest 
	da <- tg[,2]							
	db <- tg[,4]							
	dt <- da+db								
	nt <- na+nb								
	d <- dt/nt								
	O <- c(sum(da), sum(db))				
	ea <- na*d								
	eb <- nb*d								
	E <- c(sum(ea), sum(eb))				
	result <- data.frame(da, db, dt, na, nb, nt, d, ea, eb)
	if (method == "Tominaga") {                   
                method <- "Mantel Byar(Tominaga)"
                chi <- sum((O-E)^2/E)
	} else {                                      
                method <- "Mantel Byar test"
                v <- sum(dt*(nt-dt)/(nt-1)*na/nt*(1-na/nt), na.rm=TRUE)
                chi <- (sum(da)-sum(na*d))^2/v
#			print (paste("(O-E) = ", sum(da)-sum(na*d), ", V=", v, sep="") )
#			HR <- 1 / exp((sum(da)-sum(na*d))/v)
#			print (paste("HR = ", HR, sep="") )			
	}
	P <- pchisq(chi, 1, lower.tail=FALSE)
	if(plot>=1){							#If plot>=1, draw Simon Makuch plot with a landmark as specified.
		StartTime2 <- StartTime[StopTime>=landmark]
		StopTime2 <- StopTime[StopTime>=landmark]
		Event2 <- Event[StopTime>=landmark]
		Group2 <- Group[StopTime>=landmark]
		km <- survfit(Surv(StartTime2,StopTime2,Event2)~Group2, na.action = na.omit, conf.type="log-log")
		diff <- survdiff(Surv(StopTime2,Event2)~Group2)
		summary(km)
		km$n.risk[1] <- diff$n[1]			#To correct number at risk at zero point in no event group
		len <- nchar("Group2")
		legend <- substring(names(km$strata), len+2)
#		windows(width=7, height=7); par(lwd=1, las=1, family="sans", cex=1)
#		dev.new()
		if (.Platform$OS.type == 'windows'){
			justDoIt(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
		} else if (MacOSXP()==TRUE) {
			justDoIt(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
		} else {
			justDoIt(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
		}
		mar <- par("mar")
		mar[1] <- mar[1] + length(km$strata) + 0.5
		mar[2] <- mar[2] + 2
		par(mar=mar)
		opar <- par(mar = mar)
		on.exit(par(opar))
#		plot(km, ylab="Probability", bty="l", col=1:32, lty=1, lwd=1, conf.int=FALSE, mark.time=TRUE)
		if(plot==1) plot(km, ylab="Probability", bty="l", col=1:32, lty=1, lwd=1, conf.int=FALSE, mark.time=TRUE)
		if(plot==2) plot(km, ylab="Probability", bty="l", col=1, lty=1:32, lwd=1, conf.int=FALSE, mark.time=TRUE)
		if(plot>=3) plot(km, ylab="Probability", bty="l", col=1, lty=1, lwd=1:32, conf.int=FALSE, mark.time=TRUE)
		xticks <- axTicks(1)
		n.atrisk <- nrisk(km, xticks)
		for (i in 1:length(km$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}
		for (i in 1:length(km$strata)){mtext(legend[i], at=-(xticks[2]-xticks[1])/2, side=1, line=4+i, cex=1)}
		title(xlab = "Number at risk", line = 3.5, adj = 0)
#		legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title="Time-dependent covariate")
		if(plot==1) legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title="Time-dependent covariate")
		if(plot==2) legend ("topright", legend, col=1, lty=1:32, lwd=1,  box.lty=0, title="Time-dependent covariate")
		if(plot>=3) legend ("topright", legend, col=1, lty=1, lwd=1:32,  box.lty=0, title="Time-dependent covariate")
	}
	return(structure(list(statistic=c("X-squared"=chi), parameter=c(df=1), p.value=P, 
        method=method, data.name=data.name, result=result), class="htest"))
}


step.p.lm <- function (lm, dataframe.name, waldtest=0, subset=NULL){
	formula1 <- lm$terms[[2]]
	res <- summary(lm)
	reslist <- rownames(res$coefficients)[2:length(rownames(res$coefficients))]	
	var <- colnames(lm$model)[2:length(colnames(lm$model))]
	nvar <- length(var)

	dum <- NA
	fac <- NA
	for (i in 1:length(reslist)){
		dum[i] <- NA
		fac[i] <- NA
		if (regexpr(".Dummy.", reslist[i])>0) {		#Check dummy variables
			dum[i] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])+6)	
			next
		}
		for (j in 1:nchar(reslist[i])){				#Check factors
			if (substring(reslist[i], j, j)=="["){
			fac[i] <- substring(reslist[i], 1, j)
			next
			}				
		}
	}
	
	dum.list <- levels(factor(dum))		
	fac.list <- levels(factor(fac))		
	reslist <- rownames(res$coefficients)[2:length(rownames(res$coefficients))]	
	p.value <- res$coefficients[,4][2:length(rownames(res$coefficients))]
	subset <- ifelse(is.null(subset), "", paste(", subset=", subset, sep=""))	
	while(max(p.value) >= 0.05) {
		colnames(res$coefficients) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(res$coefficients))
		print(res$coefficients)
		if(length(dum.list)!=0){		#set the p values of dummy variables at minimum value
			for(i in 1:length(dum.list)){
				if (length(p.value[substring(reslist, 1, nchar(dum.list[i]))==dum.list[i]])>=2){
					wald <- wald.test(vcov(lm), lm$coef, which(substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]))
					p.value[substring(reslist, 1, nchar(dum.list[i]))==dum.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(length(fac.list)!=0){		#set the p values of factors at minimum value
			for(i in 1:length(fac.list)){
				if (length(p.value[substring(reslist, 1, nchar(fac.list[i]))==fac.list[i]])>=2){
					wald <- wald.test(vcov(lm), lm$coef, which(substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]))
					p.value[substring(reslist, 1, nchar(fac.list[i]))==fac.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(max(p.value) < 0.05) break	
		del <- reslist[p.value==max(p.value)]
		if(length(del)>1)del <- del[1]

		delete.flag <- 0
		if(length(dum.list)!=0){		
			for(i in 1:length(dum.list)){
				if (substring(del, 1, nchar(dum.list[i]))==dum.list[i]){
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", substring(del, 1, regexpr(".Dummy.", del)+6), " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test"), "\n\n", sep=""))
					var <- subset(var, substring(var, 1, nchar(dum.list[i]))!=dum.list[i])
					delete.flag <- 1
				}
			}
		}
		if(length(fac.list)!=0){		
			for(i in 1:length(fac.list)){
				if (substring(del, 1, nchar(fac.list[i]))==fac.list[i]){
					del <- substring(fac.list[i], 1, nchar(fac.list[i])-1)
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					var <- subset(var, var!=del)
					delete.flag <- 1
				}
			}
		}
		if(delete.flag==0){
			cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), ")\n\n", sep=""))
			var <- subset(var, var!=del)
		}

		nvar <- length(var)
		if (nvar==0) {
			cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----All variables were removed from the model."), "\n\n", sep="")
			nvar <- 0
			break
		}
		formula <- paste(formula1, " ~ ", var[1], sep="")
		if (nvar > 1){
			for(i in 2:nvar){
				formula <- paste(formula, "+", var[i])
			}
		}
		command <- paste("lm <- lm(", formula, ", data=", dataframe.name, subset, ")", sep="")
#		cat(command, "\n\n")
		eval(parse(text=command))
		res <- summary(lm)
		reslist <- rownames(res$coefficients)[2:length(rownames(res$coefficients))]	
		p.value <- res$coefficients[,4][2:length(rownames(res$coefficients))]
	}
	if(nvar>=1){
		cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Final model"), "\n\n", sep="")
		print(res$coefficients)
		if(waldtest==1) {waldtest(lm)}
	}
}


step.p.glm <- function (glm, dataframe.name, waldtest=0, subset=NULL){
	formula1 <- glm$terms[[2]]
	res <- summary(glm)
	reslist <- rownames(res$coefficients)[2:length(rownames(res$coefficients))]		
	var <- colnames(glm$model)[2:length(colnames(glm$model))]
	nvar <- length(var)

	dum <- NA
	fac <- NA
	for (i in 1:length(reslist)){
		dum[i] <- NA
		fac[i] <- NA
		if (regexpr(".Dummy.", reslist[i])>0) {		
			dum[i] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])+6)	
			next
		}
		for (j in 1:nchar(reslist[i])){			
			if (substring(reslist[i], j, j)=="["){
			fac[i] <- substring(reslist[i], 1, j)
			next
			}				
		}
	}
	
	dum.list <- levels(factor(dum))		
	fac.list <- levels(factor(fac))		
	reslist <- rownames(res$coefficients)[2:length(rownames(res$coefficients))]		
	p.value <- res$coefficients[,4][2:length(rownames(res$coefficients))]
	subset <- ifelse(is.null(subset), "", paste(", subset=", subset, sep=""))	
	while(max(p.value) >= 0.05) {
		odds <- data.frame(exp(res$coef[,1:2] %*% rbind(c(1,1,1), 1.96*c(0,-1,1))))
		odds <- cbind(odds, res$coefficients[,4])
		odds <- signif(odds, digits=3)
		names(odds) <- gettext(domain="R-RcmdrPlugin.EZR",c("odds ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		print(odds)

		if(length(dum.list)!=0){	
			for(i in 1:length(dum.list)){
				if (length(p.value[substring(reslist, 1, nchar(dum.list[i]))==dum.list[i]])!=0){
					wald <- wald.test(vcov(glm), glm$coef, which(substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]))
					p.value[substring(reslist, 1, nchar(dum.list[i]))==dum.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(length(fac.list)!=0){	
			for(i in 1:length(fac.list)){
				if (length(p.value[substring(reslist, 1, nchar(fac.list[i]))==fac.list[i]])!=0){
					wald <- wald.test(vcov(glm), glm$coef, which(substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]))
					p.value[substring(reslist, 1, nchar(fac.list[i]))==fac.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(max(p.value) < 0.05) break	
		del <- reslist[p.value==max(p.value)]
		if(length(del)>1)del <- del[1]
		delete.flag <- 0
		if(length(dum.list)!=0){	
			for(i in 1:length(dum.list)){
				if (substring(del, 1, nchar(dum.list[i]))==dum.list[i]){
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", substring(del, 1, regexpr(".Dummy.", del)+6), " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					var <- subset(var, substring(var, 1, nchar(dum.list[i]))!=dum.list[i])
					delete.flag <- 1
				}
			}
		}
		if(length(fac.list)!=0){	
			for(i in 1:length(fac.list)){
				if (substring(del, 1, nchar(fac.list[i]))==fac.list[i]){
					del <- substring(fac.list[i], 1, nchar(fac.list[i])-1)
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					var <- subset(var, var!=del)
					delete.flag <- 1
				}
			}
		}
		if(delete.flag==0){
			cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), ")\n\n", sep=""))
			var <- subset(var, var!=del)
		}

		nvar <- length(var)
		if (nvar==0) {
			cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----All variables were removed from the model."), "\n\n", sep="")
			nvar <- 0
			break
		}
		formula <- paste(formula1, " ~ ", var[1], sep="")
		if (nvar > 1){
			for(i in 2:nvar){
				formula <- paste(formula, "+", var[i])
			}
		}
		command <- paste("glm <- glm(", formula, ", data=", dataframe.name, subset, ", family=binomial(logit))", sep="")
#		cat(command, "\n\n")
		eval(parse(text=command))
		res <- summary(glm)
		reslist <- rownames(res$coefficients)[2:length(rownames(res$coefficients))]		
		p.value <- res$coefficients[,4][2:length(rownames(res$coefficients))]
	}
	if(nvar>=1){
		odds <- data.frame(exp(res$coef[,1:2] %*% rbind(c(1,1,1), 1.96*c(0,-1,1))))
		odds <- cbind(odds, res$coefficients[,4])
		odds <- signif(odds, digits=3)
		names(odds) <- gettext(domain="R-RcmdrPlugin.EZR",c("odds ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Final model"), "\n\n", sep="")
		print(odds)
		if(waldtest==1) {waldtest(glm)}
	}
}


step.p.cox <- function (cox, dataframe.name, waldtest=0, subset=NULL){
	formula1 <- cox$terms[[2]]
#	formula1 <- paste("Surv(", formula1[[2]], ", ", formula1[[3]], "==1)", sep="")
	formula1 <- paste("Surv(", formula1[[2]], ", ", as.character(formula1[[3]][2]), "==1)", sep="") #Change from EZR 1.20 according to the update of survival package

	res <- summary(cox)
	reslist <- rownames(res$coefficients)	
	
	dum <- NA
	fac <- NA
	k <- 1
	var <- NA

	for (i in 1:length(reslist)){
		dum[i] <- NA
		fac[i] <- NA

		if (regexpr(".Dummy.", reslist[i])>0) {		
			dum[i] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])+6)	
			var[k] <- reslist[i]
			k <- k+1
			next
		}
		for (j in 1:nchar(reslist[i])){			
			if (substring(reslist[i], j, j)=="["){
			fac[i] <- substring(reslist[i], 1, j)
			reslist[i] <- substring(fac[i], 1, nchar(fac[i])-1)
			next
			}				
		}
		if (k==1) {
			var[k] <- reslist[i]
			k <- k+1
		} else if (reslist[i]!=var[k-1]) {
			var[k] <- reslist[i]
			k <- k+1
		}
	}
	
	dum.list <- levels(factor(dum))		
	fac.list <- levels(factor(fac))		
	nvar <- length(var)

	res <- summary(cox)
	p.value <- res$coefficients[,5]
	subset <- ifelse(is.null(subset), "", paste(", subset=", subset, sep=""))	
	print(res$call)
	cat("\n")
	while(max(p.value) >= 0.05) {
#		if(nvar==1){
		if(length(res$coefficients[,5])==1){
			cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)
			rownames(cox.table) <- rownames(res$coefficients)
			colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		} else {
			cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)
			cox.table <- data.frame(cox.table)
			names(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		}
		print(cox.table)
		if(length(dum.list)!=0){		
			for(i in 1:length(dum.list)){
				if (length(p.value[substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]])>=2){
					wald <- wald.test(cox$var, cox$coef, which(substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]))
					p.value[substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(length(fac.list)!=0){		
			for(i in 1:length(fac.list)){
				if (length(p.value[substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]])>=2){
					wald <- wald.test(cox$var, cox$coef, which(substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]))
					p.value[substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(max(p.value) < 0.05) break	
		del <- rownames(res$coefficients)[p.value==max(p.value)]
		if(length(del)>1) del <- del[1]
		delete.flag <- 0
		if(length(dum.list)!=0){		
			for(i in 1:length(dum.list)){
				if (substring(del, 1, nchar(dum.list[i]))==dum.list[i]){
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", dum.list[i], " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					var <- subset(var, substring(var, 1, nchar(dum.list[i]))!=dum.list[i])
					delete.flag <- 1
				}
			}
		}
		if(length(fac.list)!=0){		
			for(i in 1:length(fac.list)){
				if (substring(del, 1, nchar(fac.list[i]))==fac.list[i]){
					del <- substring(fac.list[i], 1, nchar(fac.list[i])-1)
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					var <- subset(var, var!=del)
					delete.flag <- 1
				}
			}
		}
		if(delete.flag==0){
			cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), ")\n\n", sep=""))
			var <- subset(var, var!=del)
		}
		nvar <- length(var)
		if (nvar==0) {
			cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----All variables were removed from the model."), "\n\n", sep="")
			break
		}
		formula <- paste(formula1, " ~ ", var[1], sep="")
		if (nvar > 1){
			for(i in 2:nvar){
				formula <- paste(formula, "+", var[i])
			}
		}
		command <- paste("cox <- coxph(", formula, ", data=", dataframe.name, subset, ', method="breslow")', sep="")
		cat(command, "\n\n")
		eval(parse(text=command))
		res <- summary(cox)
		p.value <- res$coefficients[,5]
	}
#	if(nvar==1){
	if(length(res$coefficients[,5])==1){
		cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)
		rownames(cox.table) <- rownames(res$coefficients)
		colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
	} 
#	if (nvar>=2){
	if (length(res$coefficients[,5])>=2){
		cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)
		cox.table <- data.frame(cox.table)
		names(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
	}
	if(nvar>=1){
		cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Final model"), "\n\n", sep="")
		print(cox.table)
		if(waldtest==1) {waldtest(cox)}
	}
}


step.p.coxtd <- function (cox, dataframe.name, waldtest=0, subset=NULL){

	formula1 <- "Surv(start_td, stop_td, endpoint_td==1)"    #Only the different point from step.p.cox

	res <- summary(cox)
	reslist <- rownames(res$coefficients)	
	
	dum <- NA
	fac <- NA
	k <- 1
	var <- NA

	for (i in 1:length(reslist)){
		dum[i] <- NA
		fac[i] <- NA

		if (regexpr(".Dummy.", reslist[i])>0) {		
			dum[i] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])+6)	
			var[k] <- reslist[i]
			k <- k+1
			next
		}
		for (j in 1:nchar(reslist[i])){				
			if (substring(reslist[i], j, j)=="["){
			fac[i] <- substring(reslist[i], 1, j)
			reslist[i] <- substring(fac[i], 1, nchar(fac[i])-1)
			next
			}				
		}
		if (k==1) {
			var[k] <- reslist[i]
			k <- k+1
		} else if (reslist[i]!=var[k-1]) {
			var[k] <- reslist[i]
			k <- k+1
		}
	}
	
	dum.list <- levels(factor(dum))	
	fac.list <- levels(factor(fac))	
	nvar <- length(var)

	res <- summary(cox)
	p.value <- res$coefficients[,5]
	subset <- ifelse(is.null(subset), "", paste(", subset=", subset, sep=""))	
	print(res$call)
	cat("\n")
	while(max(p.value) >= 0.05) {
#		if(nvar==1){
		if(length(res$coefficients[,5])==1){
			cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)
			rownames(cox.table) <- rownames(res$coefficients)
			colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		} else {
			cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)
			cox.table <- data.frame(cox.table)
			names(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		}
		print(cox.table)
		if(length(dum.list)!=0){	
			for(i in 1:length(dum.list)){
				if (length(p.value[substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]])>=2){
					wald <- wald.test(cox$var, cox$coef, which(substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]))
					p.value[substring(rownames(res$coefficients), 1, nchar(dum.list[i]))==dum.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(length(fac.list)!=0){	
			for(i in 1:length(fac.list)){
				if (length(p.value[substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]])>=2){
					wald <- wald.test(cox$var, cox$coef, which(substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]))
					p.value[substring(rownames(res$coefficients), 1, nchar(fac.list[i]))==fac.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(max(p.value) < 0.05) break	
		del <- rownames(res$coefficients)[p.value==max(p.value)]
		if(length(del)>1) del <- del[1]
		delete.flag <- 0
		if(length(dum.list)!=0){	
			for(i in 1:length(dum.list)){
				if (substring(del, 1, nchar(dum.list[i]))==dum.list[i]){
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", dum.list[i], " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					var <- subset(var, substring(var, 1, nchar(dum.list[i]))!=dum.list[i])
					delete.flag <- 1
				}
			}
		}
		if(length(fac.list)!=0){	
			for(i in 1:length(fac.list)){
				if (substring(del, 1, nchar(fac.list[i]))==fac.list[i]){
					del <- substring(fac.list[i], 1, nchar(fac.list[i])-1)
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					var <- subset(var, var!=del)
					delete.flag <- 1
				}
			}
		}
		if(delete.flag==0){
			cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), ")\n\n", sep=""))
			var <- subset(var, var!=del)
		}
		nvar <- length(var)
		if (nvar==0) {
			cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----All variables were removed from the model."), "\n\n", sep="")
			break
		}
		formula <- paste(formula1, " ~ ", var[1], sep="")
		if (nvar > 1){
			for(i in 2:nvar){
				formula <- paste(formula, "+", var[i])
			}
		}
		command <- paste("cox <- coxph(", formula, ", data=", dataframe.name, subset, ', method="breslow")', sep="")
		cat(command, "\n\n")
		eval(parse(text=command))
		res <- summary(cox)
		p.value <- res$coefficients[,5]
	}
#	if(nvar==1){
	if(length(res$coefficients[,5])==1){
		cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)
		rownames(cox.table) <- rownames(res$coefficients)
		colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
	} 
#	if (nvar>=2){
	if (length(res$coefficients[,5])>=2){
		cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)
		cox.table <- data.frame(cox.table)
		names(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
	}
	if(nvar>=1){
		cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Final model"), "\n\n", sep="")
		print(cox.table)
		if(waldtest==1) {waldtest(cox)}
	}
}


step.p.crr <- function (crr, cov, dataframe.name, waldtest=0, subset=NULL){
	
	dataframe.name <- ifelse(is.null(subset), dataframe.name, paste("subset(", dataframe.name, ", ", subset, ")",sep=""))		
	command <- paste("cbind(", dataframe.name, "$", cov[1], sep="")
	if(length(cov)>1){
		for(i in 2:length(cov)){
			command <- paste(command, ", ", dataframe.name, "$", cov[i], sep="")
		}
	}
	command <- paste(command, ")", sep="")
	cov.matrix <- eval(parse(text=command))
	ncov <- length(cov)

	dum <- NA
	for (i in 1:ncov){
		dum[i] <- NA
		if (regexpr(".Dummy.", cov[i])>0) {	
			dum[i] <- substring(cov[i], 1, regexpr(".Dummy.", cov[i])+6)	
		}
	}
	dum.list <- levels(factor(dum))	

	call <- as.character(crr$call)
	command <- paste("with(", dataframe.name, ", crr(",  call[2], ", ", call[3], ", cov.matrix, failcode=", call[5], ", cencode=", call[6], ", na.action=na.omit))", sep="")
	res <- summary(crr)
	p.value <- res$coef[,5]
#	print(command)
	cat("\n")
	while(max(p.value) >= 0.05) {
		if(ncov==1){
			crr.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), res$coef[,5]), digits=4)
		} else {
			crr.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coef[,5]), digits=4)
		}	
		rownames(crr.table) <- cov
		colnames(crr.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		print(crr.table)
		if(length(dum.list)!=0){	
			for(i in 1:length(dum.list)){
				if (length(p.value[substring(cov, 1, nchar(dum.list[i]))==dum.list[i]])!=0){
					wald <- wald.test(crr$var, crr$coef, which(substring(cov, 1, nchar(dum.list[i]))==dum.list[i]))
					p.value[substring(cov, 1, nchar(dum.list[i]))==dum.list[i]] <- wald[[6]][[1]][3]	
				}
			}
		}
		if(max(p.value) < 0.05) break	
		del <- cov[p.value==max(p.value)]
		if(length(del)>1)del <- del[1]
		delete.flag <- 0
		if(length(dum.list)!=0){	
			for(i in 1:length(dum.list)){
				if (substring(del, 1, nchar(dum.list[i]))==dum.list[i]){
					cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", dum.list[i], " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), " ", gettext(domain="R-RcmdrPlugin.EZR","by Wald test)"), "\n\n", sep=""))
					cov.matrix <- cov.matrix[,substring(cov, 1, nchar(dum.list[i]))!=dum.list[i]]
					cov <- subset(cov, substring(cov, 1, nchar(dum.list[i]))!=dum.list[i])
					delete.flag <- 1
				}
			}
		}
		if(delete.flag==0){
			cat(paste("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Remove"), " ", del, " ", gettext(domain="R-RcmdrPlugin.EZR","from the model. (p="), signif(max(p.value),4), ")\n\n", sep=""))
			if(is.matrix(cov.matrix)){
				cov.matrix <- cov.matrix[,cov!=del]
			} else {
				cov.matrix <- cov.matrix[cov!=del]
			}
			cov <- cov[cov!=del]
		}
		ncov <- length(cov)
		if (ncov==0) {
			cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----All variables were removed from the model."), "\n\n", sep="")
			break
		}
		command <- paste("crr <- ", command, sep="")
		eval(parse(text=command))
		res <- summary(crr)
		p.value <- res$coef[,5]
	}
	if(ncov>0){
		if(ncov==1){
			crr.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), res$coef[,5]), digits=4)
		} else {
			crr.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coef[,5]), digits=4)
		}	
		rownames(crr.table) <- cov
		colnames(crr.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
		cat("\n", gettext(domain="R-RcmdrPlugin.EZR","-----Final model"), "\n\n", sep="")
		print(crr.table)
		if (waldtest==1) waldtest.crr(crr, rownames(crr.table))
	}
}


step.AIC.crr <- function (crr, cov, dataframe.name, BIC = 0, subset = NULL, waldtest=0) {
    method <- ifelse(BIC==0, "AIC", "BIC")
    dataframe.name <- ifelse(is.null(subset), dataframe.name, 
        paste("subset(", dataframe.name, ", ", subset, ")", sep = ""))
    command <- paste(paste(dataframe.name, "$", cov, sep=""), collapse=", ")
    command <- paste("cbind(", command, ")", sep = "")
    cov.matrix <- eval(parse(text = command))
    ncov <- length(cov)
    dum <- NA					#NA for non-dummy variables
    for (i in 1:ncov) {
        dum[i] <- NA
        if (regexpr(".Dummy.", cov[i]) > 0) {
            dum[i] <- substring(cov[i], 1, regexpr(".Dummy.", 
                cov[i]) + 6)
		}
	}
    dum.list <- levels(factor(dum))		#list of dummmy variables ("....Dummy.")
	dum.list.num <- NA
	if (length(dum.list)>=1){
		for (i in 1:length(dum.list)){
			dum.list.num[i] <- length(cov[substring(cov, 1, nchar(dum.list[i]))==dum.list[i]])
		}
	}
    if (length(dum.list) == 0){
	    var.list <- cov
    } else {
	    var.list <- c(cov[is.na(dum)], substring(dum.list, 1, nchar(dum.list)-7))
    }
    cov <- c(cov[is.na(dum)], cov[!is.na(dum)])			#rearrange cov according to the var.list
    dum <- c(dum[is.na(dum)], dum[!is.na(dum)])			#rearrange dum according to the var.list
	j <- length(cov[is.na(dum)])
    cov.to.var.list <- 1:j
	if (length(dum.list)>=1){
		for (i in 1:length(dum.list)){
			j <- j + 1
			cov.to.var.list <- c(cov.to.var.list, rep(j, dum.list.num[i]))
		}
	}
	
    in.model <- rep(1, length(var.list))				#1 if in model
    var.list.dum <- rep(1, length(var.list))			#1 if dummy
    if(length(cov[is.na(dum)])>=1)  var.list.dum[1:length(cov[is.na(dum)])] <- 0
    call <- as.character(crr$call)
    command <- paste("with(", dataframe.name, ", crr(", call[2], 
        ", ", call[3], ", cov.matrix, failcode=", call[5], ", cencode=", 
        call[6], ", na.action=na.omit))", sep = "")
    currentAIC <- crrAIC(crr, BIC)

    cat("\n\n", gettext(domain="R-RcmdrPlugin.EZR", "Current model:"), " ", paste(var.list[in.model==1], collapse=" + "), "\n", sep="")
    cat(method, " = ", currentAIC, "\n\n", sep="")

    cat("\n")
    flag <- 0

    while (flag==0) {		####while routine for forward/backward selection
	newAIC <- NA
	action <- NA
	target <- NA
	for (i in 1:length(var.list)){
		if (in.model[i]==0){
			action[i] <- "+"
			target[i] <- var.list[i]
			in.model[i] <- 1
		    	command <- "cbind("
			first.var <- 1
			for (j in 1:length(var.list)){
				if (in.model[j]==1){
				    if (first.var==0){
					    command <- paste(command, ",", sep="")
				    }	
				    command2 <- paste(paste(dataframe.name, "$", cov[cov.to.var.list==j], sep=""), collapse=",")
				    command <- paste(command, command2, sep="")
				    first.var <- 0
				}
			}
			command <- paste(command, ")", sep = "")
    			cov.matrix <- eval(parse(text = command))
			in.model[i] <- 0

    			command <- paste("with(", dataframe.name, ", crr(", call[2], 
        		", ", call[3], ", cov.matrix, failcode=", call[5], ", cencode=", 
        		call[6], ", na.action=na.omit))", sep = "")
    			crr2 <- eval(parse(text = command))
			newAIC[i] <- crrAIC(crr2, BIC)
		} else {
			action[i] <- "-"
			target[i] <- var.list[i]
			in.model[i] <- 0
		    if (sum(in.model)>0){
				command <- "cbind("
				first.var <- 1
				for (j in 1:length(var.list)){
					if (in.model[j]==1){
						if (first.var==0){
							command <- paste(command, ",", sep="")
						}	
						command2 <- paste(paste(dataframe.name, "$", cov[cov.to.var.list==j], sep=""), collapse=",")
						command <- paste(command, command2, sep="")
						first.var <- 0
					}
				}
			
			command <- paste(command, ")", sep = "")
    			cov.matrix <- eval(parse(text = command))

    			command <- paste("with(", dataframe.name, ", crr(", call[2], 
        		", ", call[3], ", cov.matrix, failcode=", call[5], ", cencode=", 
        		call[6], ", na.action=na.omit))", sep = "")
    			crr2 <- eval(parse(text = command))
			newAIC[i] <- crrAIC(crr2, BIC)
			} else {
				newAIC[i] <- ifelse(BIC==0, -2 * crr$loglik.null, -2 * crr$loglik.null)
			}
			in.model[i] <- 1
		}
	}

	action[length(var.list)+1] <- "<none>"
	target[length(var.list)+1] <- ""
	newAIC[length(var.list)+1] <- currentAIC
	res <- cbind(action, target, signif(newAIC, digits=7))
	res <- data.frame(res[order(newAIC),])
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR", c("action", "variable", method))
	print(res)

	min <- min(newAIC)
	if(currentAIC <= min){
		flag <- 1
	} else {
		change.var <- which(newAIC==min)
		currentAIC <- min
		if (in.model[change.var]==1){
			in.model[change.var] <- 0
	        cat("\n", gettext(domain="R-RcmdrPlugin.EZR", "-----Variable"), " ", var.list[change.var], " ", gettext(domain="R-RcmdrPlugin.EZR", 
            "removed from the model."), "(", method, "=", newAIC[change.var], ")\n\n", sep = "")
		} else {
			in.model[change.var] <- 1
	        cat("\n", gettext(domain="R-RcmdrPlugin.EZR", "-----Variable"), " ", var.list[change.var], " ", gettext(domain="R-RcmdrPlugin.EZR", 
            "removed from the model."), "(", method, "=", newAIC[change.var], ")\n\n", sep = "")
		}
		cat(gettext(domain="R-RcmdrPlugin.EZR", "Next model:"), " ", paste(var.list[in.model==1], collapse=" + "), "\n", sep="")
		cat(method, " = ", currentAIC, "\n\n", sep="")
	}
	}

#final model

     if (sum(in.model) == 0) {
            cat("\n", gettext(domain="R-RcmdrPlugin.EZR", "-----All variables were removed from the model."), 
                "\n\n", sep = "")
     } else {
		command <- "cbind("
		first.var <- 1
		final.cov <- NULL
		for (j in 1:length(var.list)){
			if (in.model[j]==1){
				if (first.var==0){
					command <- paste(command, ",", sep="")
				}	
				final.cov <- c(final.cov, cov[cov.to.var.list==j])
				command2 <- paste(paste(dataframe.name, "$", cov[cov.to.var.list==j], sep=""), collapse=",")
				command <- paste(command, command2, sep="")
				first.var <- 0
			}
		}
		command <- paste(command, ")", sep = "")
			cov.matrix <- eval(parse(text = command))
			command <- paste("with(", dataframe.name, ", crr(", call[2], 
      		", ", call[3], ", cov.matrix, failcode=", call[5], ", cencode=", 
        		call[6], ", na.action=na.omit))", sep = "")
			crr <- eval(parse(text = command))
		res <- summary(crr)

		ncov <- length(cov.matrix[1,])
		if (ncov > 0) {
			if (ncov == 1) {
				crr.table <- signif(cbind(t(res$conf.int[, c(1, 3, 
                4)]), res$coef[, 5]), digits = 4)
			} else {
				crr.table <- signif(cbind(res$conf.int[, c(1, 3, 
					4)], res$coef[, 5]), digits = 4)
			}
#		rownames(crr.table) <- cov[in.model[cov.to.var.list]==1]
		rownames(crr.table) <- final.cov
		colnames(crr.table) <- gettext(domain="R-RcmdrPlugin.EZR", 
            c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))
			cat("\n", gettext(domain="R-RcmdrPlugin.EZR", "-----Final model"), 
            "\n\n", sep = "")
			print(crr.table)
			if (waldtest == 1) 
				waldtest.crr(crr, rownames(crr.table))
		}
	}
}


crrAIC <- function(crr, BIC=0){
	AIC <- ifelse(BIC==0, -2 * crr$loglik + 2 * length(crr$coef),
		-2 * crr$loglik + log(crr$n) * length(crr$coef))
	return(AIC)	
}


waldtest <- function (cox){
#	This function can be used not only for cox but also for lm and  glm.

	res <- summary(cox)
	reslist <- rownames(res$coefficients)	
	
	dum <- NA	
	fac <- NA	
	k <- 1
	for (i in 1:length(reslist)){
		if (regexpr(".Dummy.", reslist[i])>0) {		
			if(k==1){
				dum[k] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])-1)
				k <- k + 1
			} else if(substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])-1)!=dum[k-1]){
				dum[k] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])-1)
				k <- k + 1
			}
		}
	}
	
	k <- 1
	for (i in 1:length(reslist)){
		for (j in 1:nchar(reslist[i])){			
			if (substring(reslist[i], j, j)=="["){
				if(k==1){
					fac[k] <- substring(reslist[i], 1, j-1)
					k <- k + 1
					next
				} else if (substring(reslist[i], 1, j-1)!=fac[k-1]){
					fac[k] <- substring(reslist[i], 1, j-1)
					k <- k + 1
					next
				} else {
					next
				}
			}
		}
	}

	dum <- levels(factor(dum))	
	fac <- levels(factor(fac))	

	if(length(dum)>0){
		for(i in 1:length(dum)){
			terms <- which(substring(rownames(res$coefficients), 1, nchar(dum[i]))==dum[i])
			if (length(terms)>=2){
				wald <- wald.test(vcov(cox), cox$coef, terms)
				cat(gettext(domain="R-RcmdrPlugin.EZR","\nOverall p value for"), dum[i], ": ", wald[[6]][[1]][3], "\n")			
			}
		}
	}

	if(length(fac)>0){
		for(i in 1:length(fac)){
			terms <- which(substring(rownames(res$coefficients), 1, nchar(fac[i]))==fac[i])
			if (length(terms)>=2){
				wald <- wald.test(vcov(cox), cox$coef, terms)
				cat(gettext(domain="R-RcmdrPlugin.EZR","\nOverall p value for"), fac[i], ": ", wald[[6]][[1]][3], "\n")			
			}
		}
	}

}


waldtest.crr <- function (crr, cov){

	reslist <- cov	
	
	dum <- NA	
	fac <- NA	
	k <- 1
	for (i in 1:length(reslist)){
		if (regexpr(".Dummy.", reslist[i])>0) {		
			if(k==1){
				dum[k] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])-1)
				k <- k + 1
			} else if(substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])-1)!=dum[k-1]){
				dum[k] <- substring(reslist[i], 1, regexpr(".Dummy.", reslist[i])-1)
				k <- k + 1
			}
		}
	}
		
	dum <- levels(factor(dum))	

	if(length(dum)>0){
		for(i in 1:length(dum)){
			terms <- which(substring(cov, 1, nchar(dum[i]))==dum[i])
			if (length(terms)>=2){
				wald <- wald.test(crr$var, crr$coef, terms)
				cat(gettext(domain="R-RcmdrPlugin.EZR","\nOverall p value for"), dum[i], ": ", wald[[6]][[1]][3], "\n")			
			}
		}
	}
}


logrank.trend <- function(survdiff.res, W = 1:length(survdiff.res[[1]])){
	#Calculation method from http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3311/handout4.pdf, Newcastle University
	#consistent with the results by MedCalc software
	# W = score for each group
	group.names <- survdiff.res[[1]]
	O <- survdiff.res$obs
	E <- survdiff.res$exp
	V <- survdiff.res$var
	Wupperbar <- sum (W * E) / sum(E)

	WOE <- W * (O - E)		
	UT <- sum(WOE)
	VT <- sum((W - Wupperbar)^2*E) 
	WT <- UT^2 / VT
	P <- pchisq(WT, df=1, lower.tail=FALSE)
	res <- data.frame(c(formatC(WT, format="g", digits=3), formatC(1, format="d"), formatC(P, format="g", digits=2)))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Logrank trend test")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("Chi square", "DF", "p-value"))
	return(res)
}


stackcuminc <- function(timetoevent, event, xlim=NULL, ylim=c(0,1), xlab=NULL, ylab=NULL, atrisk=1, main=""){
	num <- length(levels(factor(event)))
	max <- max(timetoevent, na.rm=TRUE)

	if(min(event, na.rm=TRUE)==0){	#for censored events
		censor <- 1
		num <- num-1				#Type of event will be num-1
	} else {
		censor <- 0
	}

	if(num==0)stop("No event")

	if (atrisk==1){
		doItAndPrint('mar <- par("mar")')
		doItAndPrint("mar[1] <- mar[1] + 1 + 0.5")
		doItAndPrint("par(mar=mar)")
		doItAndPrint("opar <- par(mar = mar)")
		doItAndPrint("on.exit(par(opar))")
	}

	if(num <= 1){	#Error occurs when there is only one event type observed
		ci <- survfit(Surv(timetoevent, event>0)~1, na.action=na.omit)
	}else{
#		ci <- survfit(Surv(timetoevent, event>0)~1, na.action=na.omit, etype=event)	
#If there are no censoring, an event with a smallest event number will be
#treated as censoring in the new survival package. To avoid this, make the smallest 
#level as "0".
	if(censor==0){
			event <- factor(event, levels=c("0", levels(as.factor(event))))
		}else{
			event <- as.factor(event)
		}
		ci <- survfit(Surv(timetoevent, event, type="mstate")~1, na.action=na.omit)	
	}

	time <- rep(ci$time[1], 2)
	for (i in 2:length(ci$time)){
		time <- c(time, rep(ci$time[i], 2))
	}
	time <- c(time, rep(max, 2))	
	if (is.null(ci$surv)) ci$surv <- 1-ci$prev	#added from EZR ver 1.11 
	ci$surv <- 1-ci$surv
	y <- rep(0, num)
	for (i in 1:length(ci$time)){	
		next.y <- NULL
		for(j in 1:num){
			if (num==1){
				next.y[j] <- sum(ci$surv[i])		
			} else {
				next.y[j] <- sum(ci$surv[i, j:num])
			}
		}
		y <- rbind(y, next.y)
		y <- rbind(y, next.y)
	}
	y <- rbind(y, rep(0, num))

	for(i in 1:num){
		if (i==1) {
			plot(ci, fun="event", col=0, bty="l", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
			if(atrisk==1){
				xticks <- axTicks(1)
				n.atrisk <- nrisk(ci, xticks)			
				axis(1, at = xticks, labels = n.atrisk, line = 3, tick = FALSE)
				title(xlab = "Number at risk", line = 3, adj = 0)			
			}
		}
		if (num==1){
			polygon(c(0, time, max), c(0, y, 0), col=gray(1-0.1*i))
		}else{
			polygon(c(0, time, max), c(0, y[ ,i], 0), col=gray(1-0.1*i))
		}
	}
	legend("topleft", legend=levels(factor(event))[(censor+1):(censor+num)],  col=gray(seq(0.9, 1-0.1*num, by=-0.1)), bty="n", lty=1, lwd=10)
}


#roc.best <- function (..., roc){

#	n <- length(roc[[4]])
#	min.distance <- 1
#	sensitivity <- roc[[4]][1]
#	specificity <- roc[[5]][1]
#	threshold <- roc[[6]][1]	
#	for (i in 1:n){
#		distance <- (1-roc[[4]][i])^2 + (1-roc[[5]][i])^2
#		if (distance < min.distance){
#			min.distance <- distance
#			sensitivity <- roc[[4]][i]			
#			specificity <- roc[[5]][i]
#			threshold <- roc[[6]][i]	
#		}		
#	}
#	res <- c(threshold=threshold, specificity=specificity, sensitivity=sensitivity)
#	return(res)
#}


SampleProportionSingleArm <- function (p1, p2, alpha, power, method, continuity) {
	#method = 2 for two sided, 1 for one sided
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	#Continuity correction method is from S-PLUS manual; binomial.sample.size()
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha/method
	ZA <- qnorm(1-alpha2)
	ZB <- qnorm(power)
	N <- ceiling((ZA*sqrt(p1*(1-p1))+ZB*sqrt(p2*(1-p2)))^2 / (p2-p1)^2 )
	if(continuity==1){
		N <- ceiling(N + (2 / abs(p2-p1))) #from S-PLUS manual
	} 
	res <- data.frame(c(p1, p2, alpha, side, power, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P in the population", "Alternative P", "Alpha", " ", "Power", "  ", "   ", "Required sample size"))
	y <- seq(0.2, 1, 0.05) 
	if(continuity==1){
		plot((ZA*sqrt(p1*(1-p1))+qnorm(y)*sqrt(p2*(1-p2)))^2 / (p2-p1)^2 + (2 / abs(p2-p1)), y, ylim=c(0,1), type="l", ylab="Power", xlab="N")
	} else {
		plot((ZA*sqrt(p1*(1-p1))+qnorm(y)*sqrt(p2*(1-p2)))^2 / (p2-p1)^2, y, ylim=c(0,1), type="l", ylab="Power", xlab="N")
	}
	abline(h=power, lty=2)
	return(res)
}
		

PowerProportionSingleArm <- function (p1, p2, alpha, n, method, continuity) {
	#method = 2 for two sided, 1 for one sided
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	#Continuity correction method is from S-PLUS manual; binomial.sample.size()
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha/method
	N <- n
	if(continuity==1){
		N <- N - (2 / abs(p2-p1)) #from S-PLUS manual
	} 
	ZA <- qnorm(1-alpha2)
	ZB <- (sqrt(N) * (p2-p1) - qnorm(1-alpha2)*sqrt(p1*(1-p1))) / sqrt(p2*(1-p2))
	power <- signif(pnorm(ZB), digits=3)
	res <- data.frame(c(p1, p2, alpha, side, n, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), power))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P in the population", "Alternative P", "Alpha", " ", "Sample size", "  ", "   ", "Power"))
	y <- seq(0.2, 1, 0.05) 
	
	if(continuity==1){
		plot((ZA*sqrt(p1*(1-p1))+qnorm(y)*sqrt(p2*(1-p2)))^2 / (p2-p1)^2 + (2 / abs(p2-p1)), y, ylim=c(0,1), type="l", ylab="Power", xlab="N")
	} else {
		plot((ZA*sqrt(p1*(1-p1))+qnorm(y)*sqrt(p2*(1-p2)))^2 / (p2-p1)^2, y, ylim=c(0,1), type="l", ylab="Power", xlab="N")
	}
#	plot((ZA*sqrt(p1*(1-p1))+qnorm(y)*sqrt(p2*(1-p2)))^2 / (p2-p1)^2, y, ylim=c(0,1), type="l", ylab="Power", xlab="N")
	abline(v=n, lty=2)
	return(res)
}
		

SampleProportionCI <- function (p, delta, ci) {
	#From Igakuteki Kenkyuno Design
	alpha <- (100 - ci) / 100
	ZA <- qnorm(1-alpha/2)
	N <- ceiling((4*ZA^2*p*(1-p)) / (delta^2))
	res <- data.frame(c(p, delta, ci/100, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P", "Confidence interval", "Confidence level", " ", "  ", "Required sample size"))
	y <- seq(delta/2, delta*2, length=20) 
	plot((4*ZA^2*p*(1-p)) / (y^2), y, ylim=c(0,1), type="l", ylab="Confidence interval", xlab="N")
	abline(h=delta, lty=2)
	return(res)		
}
		
	
SampleMeanCI <- function (sd, delta, ci) {
	#From Igakuteki Kenkyuno Design
	alpha <- (100 - ci) / 100
	ZA <- qnorm(1-alpha/2)
	N <- ceiling((4*ZA^2*sd^2) / (delta^2) )
	res <- data.frame(c(sd, delta, ci/100, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("Standard deviation", "Confidence interval", "Confidence level", " ", "  ", "Required sample size"))
	y <- seq(delta/2, delta*2, length=20) 
	plot((4*ZA^2*sd^2) / (y^2), y, ylim=c(0,delta*2.2), type="l", ylab="Delta", xlab="N")
	abline(h=delta, lty=2)
	return(res)
}


SamplePhaseII <- function (p1, p2, alpha, power) {
	ZA <- qnorm(1-alpha)
	ZB <- qnorm(power)
	N <- ceiling(((ZA*sqrt(p1*(1-p1))+ZB*sqrt(p2*(1-p2)))^2)/((p2-p1)^2))
	res <- data.frame(c(p1, p2, alpha, power, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P0 (not worth studying further)", "P1 (worth studying further", "Alpha", "Power", "  ", "   ", "Required sample size"))
	return(res)
}		


SampleMean <- function (difference, sd, alpha, power, method, r) {
	#method = 2 for two sided, 1 for one sided, r for group2/group1 ratio
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
	ZB <- qnorm(power)
	N1 <- ceiling((1+1/r)*((ZA+ZB)^2)*((sd/difference)^2))
	N2 <- N1 * r
	res <- data.frame(c(difference, sd, alpha, side, power, r, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N1, N2))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("Difference in means", "Standard deviation", "Alpha", " ", "Power", "N2/N1", "  ", "Required sample size", "N1", "N2"))
	y <- seq(0.2, 1, 0.05) 
	plot((1+1/r)*((ZA+qnorm(y))^2)*((sd/difference)^2), y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(h=power, lty=2)
	return(res)
}


SampleMeanPaired <- function (difference, sd, alpha, power, method) {
	#method = 2 for two sided, 1 for one sided
	side <- ifelse(method == 1, "one.sided", "two.sided")
	n <- power.t.test(power=power, delta=difference, sd=sd, sig.level=alpha, alternative=side, type="paired")
	res <- data.frame(c(difference, sd, alpha, side, power, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), ceiling(n$n)))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("Difference in means", "Standard deviation", "Alpha", " ", "Power", "  ", "Required sample size", "N"))
	x <- NULL
	y <- NULL
	for (i in 1:16){
		y[i] <- 0.15 + i * 0.05
		x[i] <- (power.t.test(power=y[i], delta=difference, sd=sd, sig.level=alpha, alternative=side, type="paired"))$n
	}
	plot(x, y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(h=power, lty=2)
	return(res)
}


PowerMean <- function (difference, sd, alpha, n, method, r) {
	#method = 2 for two sided, 1 for one sided, r for group2/group1 ratio
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
	N <- n
	ZB <- difference/sd*(1/sqrt((1+1/r)/N))-ZA
	power <- signif(pnorm(ZB), digits=3)
	res <- data.frame(c(difference, sd, alpha, side, " ", N, round(N*r, 0), " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), power))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("Difference in means", "Standard deviation", "Alpha", " ", "Sample size", "N1", "N2", "  ", "   ", "Power"))
	y <- seq(0.2, 1, 0.05) 
	plot((1+1/r)*((ZA+qnorm(y))^2)*((sd/difference)^2), y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(v=N, lty=2)
	return(res)
}


PowerMeanPaired <- function (difference, sd, alpha, n, method) {
	#method = 2 for two sided, 1 for one sided
	side <- ifelse(method == 1, "one.sided", "two.sided")
	power <- power.t.test(n=n, delta=difference, sd=sd, sig.level=alpha, alternative=side, type="paired")
	res <- data.frame(c(difference, sd, alpha, side, n, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), signif(power$power, digits=3)))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("Difference in means", "Standard deviation", "Alpha", " ", "Sample size", "  ", "    ", "Power"))
	x <- NULL
	y <- NULL
	for (i in 1:16){
		y[i] <- 0.15 + i * 0.05
		x[i] <- (power.t.test(power=y[i], delta=difference, sd=sd, sig.level=alpha, alternative=side, type="paired"))$n
	}
	plot(x, y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(v=n, lty=2)
	return(res)
}


SampleProportion <- function (group1, group2, alpha, power, method, r, continuity) {
	#method = 2 for two sided, 1 for one sided, r for group2/group1 ratio
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
	ZB <- qnorm(power)
	WeightedMean <- (group1 + group2 * r) / (1 + r)
	Delta <- abs(group1-group2)
	Ndash <- (1/Delta^2)*(ZA*sqrt((1+r)*WeightedMean*(1-WeightedMean))+ZB*sqrt(r*group1*(1-group1)+group2*(1-group2)))^2
	if(continuity==1){
		N1 <- ceiling(Ndash/r + (1+r)/(r*Delta)) #from S-PLUS manual
		N2 <- N1 * r
	} else {
		N1 <- ceiling(Ndash/r)
		N2 <- N1 * r
	}
	res <- data.frame(c(group1, group2, alpha, side, power, r, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N1, N2))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P1", "P2", "Alpha", " ", "Power", "N2/N1", "  ", "Required sample size", "N1", "N2"))
	y <- seq(0.2, 1, 0.05) 
	if(continuity==1){	
		plot((1/Delta^2)*(ZA*sqrt((1+r)*WeightedMean*(1-WeightedMean))+qnorm(y)*sqrt(r*group1*(1-group1)+group2*(1-group2)))^2 / r + (1+r)/(r*Delta), y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	} else {
		plot((1/Delta^2)*(ZA*sqrt((1+r)*WeightedMean*(1-WeightedMean))+qnorm(y)*sqrt(r*group1*(1-group1)+group2*(1-group2)))^2 / r, y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	}
	abline(h=power, lty=2)
	return(res)
}		


PowerProportion <- function (group1, group2, alpha, n, method, r, continuity) {
	#method = 2 for two sided, 1 for one sided, r for group2/group1 ratio
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
	N <- n
	WeightedMean <- (group1 + group2 * r) / (1 + r)
	Delta <- abs(group1-group2)
	if(continuity==1){
		Ndash <- (N-(1+r)/(r*Delta)) * r
	} else {
		Ndash <- N * r
	}
	ZB <- (sqrt(Ndash/(1/Delta^2))-ZA*sqrt((1+r)*WeightedMean*(1-WeightedMean))) / (sqrt(r*group1*(1-group1)+group2*(1-group2)))  #from S-PLUS manual
	power <- signif(pnorm(ZB), digits=3)
	res <- data.frame(c(group1, group2, alpha, side, " ", N, round(N*r, 0), " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), power))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P1", "P2", "Alpha", " ", "Sample size", "N1", "N2", "  ", "   ", "Power"))
	y <- seq(0.2, 1, 0.05) 
	if(continuity==1){	
		plot((1/Delta^2)*(ZA*sqrt((1+r)*WeightedMean*(1-WeightedMean))+qnorm(y)*sqrt(r*group1*(1-group1)+group2*(1-group2)))^2 / r + (1+r)/(r*Delta), y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	} else{
		plot((1/Delta^2)*(ZA*sqrt((1+r)*WeightedMean*(1-WeightedMean))+qnorm(y)*sqrt(r*group1*(1-group1)+group2*(1-group2)))^2 / r, y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")	
	}
	abline(v=N, lty=2)
	return(res)
}


SampleProportionNonInf <- function (group1, group2, delta, alpha, power, method) {
	#From Musakui Waritsuke Hikaku Rinsho Shiken Page 66
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
	ZB <- qnorm(power)
	Mean <- (group1 + group2) / 2
	N <- ceiling(((ZA*sqrt(2*Mean*(1-Mean))+ZB*sqrt(group1*(1-group1)+group2*(1-group2)))^2) / ((group1-group2-delta)^2))
	res <- data.frame(c(group1, group2, delta, alpha, side, power, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N, N))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P1", "P2", "Delta", "Alpha", " ", "Power", "  ", "Required sample size", "N1", "N2"))
	y <- seq(0.2, 1, 0.05) 
	plot(((ZA*sqrt(2*Mean*(1-Mean))+qnorm(y)*sqrt(group1*(1-group1)+group2*(1-group2)))^2) / ((group1-group2-delta)^2), y, ylim=c(0,1), type="l", ylab="Power", xlab="N")
	abline(h=power, lty=2)
	return(res)		
}		


SampleHazard <- function (enrol, observe, followup, group1, group2, alpha, power, method, ratio) {
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
	ZB <- qnorm(power)
	L1 <- -log(group1) / followup
	L2 <- -log(group2) / followup
	LBER <- (L1 + L2) / 2
	Q1 <- 1 / (1 + ratio)
	Q2 <- ratio / (1 + ratio)
	P00 <- LBER ^2
	P01 <- L1 ^ 2
	P02 <- L2 ^ 2
	P10 <- LBER ^ 2 / (1 - exp(-LBER * observe))
	P11 <- L1 ^ 2 / ( 1 - exp(-L1 * observe))
	P12 <- L2 ^ 2 / (1 - exp(-L2 * observe))
	P20 <- LBER ^ 3 * followup / (LBER * observe - 1 + exp(-LBER * observe))
	P21 <- L1 ^ 3 * followup / (L1 * observe - 1 + exp(-L1 * observe))
	P22 <- L2 ^ 3 * followup / (L2 * observe - 1 + exp(-L2 * observe))
	if (enrol>0){
		P30 <- LBER ^ 2 * (1 - (exp(-LBER * (observe - enrol)) - exp(-LBER * observe)) / (LBER * enrol)) ^ (-1)
		P31 <- L1 ^ 2 * (1 - (exp(-L1 * (observe - enrol)) - exp(-L1 * observe)) / (L1 * enrol)) ^ (-1)
		P32 <- L2 ^ 2 * (1 - (exp(-L2 * (observe - enrol)) - exp(-L2 * observe)) / (L2 * enrol)) ^ (-1)
	}
	if (enrol == 0) {
		N <- ((ZA * sqrt(P10 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P11/Q1 + P12/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
	} else if (enrol > 0 && observe > enrol) {
		N <- ((ZA * sqrt(P30 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P31/Q1 + P32/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
	} else if (enrol == observe) {
#		N <- ((ZA * sqrt(P20 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P21/Q1 + P22/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
		N <- ((ZA * sqrt(P30 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P31/Q1 + P32/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
	}
	N1 <- ceiling(N)
	N2 <- ceiling(N * ratio)
	res <- data.frame(c(group1, group2, followup, enrol, observe, alpha, side, power, ratio, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N1, N2))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P1", "P2", "(Follow-up duration for P1, P2)", "Enrollment duration", "Total study duration", "Alpha", " ", "Power", "N2/N1", "  ", "Required sample size", "N1", "N2"))

	x <- NULL	
	y <- seq(0.2, 1, 0.05) 
	for (i in 1: length(y)){
		ZB <- qnorm(y[i])
		if (enrol == 0) {
			x[i] <- ((ZA * sqrt(P10 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P11/Q1 + P12/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
		} else if (enrol > 0 && observe > enrol) {
			x[i] <- ((ZA * sqrt(P30 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P31/Q1 + P32/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
		} else if (enrol == observe) {
#			x[i] <- ((ZA * sqrt(P20 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P21/Q1 + P22/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
			x[i] <- ((ZA * sqrt(P30 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P31/Q1 + P32/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
		}
	}
	plot(x, y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(h=power, lty=2)

	return(res)
}		


PowerHazard <- function (enrol, observe, followup, group1, group2, alpha, sample, method, ratio) {
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
	L1 <- -log(group1) / followup
	L2 <- -log(group2) / followup
	LBER <- (L1 + L2) / 2
	Q1 <- 1 / (1 + ratio)
	Q2 <- ratio / (1 + ratio)
	P00 <- LBER ^2
	P01 <- L1 ^ 2
	P02 <- L2 ^ 2
	P10 <- LBER ^ 2 / (1 - exp(-LBER * observe))
	P11 <- L1 ^ 2 / ( 1 - exp(-L1 * observe))
	P12 <- L2 ^ 2 / (1 - exp(-L2 * observe))
	P20 <- LBER ^ 3 * followup / (LBER * observe - 1 + exp(-LBER * observe))
	P21 <- L1 ^ 3 * followup / (L1 * observe - 1 + exp(-L1 * observe))
	P22 <- L2 ^ 3 * followup / (L2 * observe - 1 + exp(-L2 * observe))
	if (enrol>0){
		P30 <- LBER ^ 2 * (1 - (exp(-LBER * (observe - enrol)) - exp(-LBER * observe)) / (LBER * enrol)) ^ (-1)
		P31 <- L1 ^ 2 * (1 - (exp(-L1 * (observe - enrol)) - exp(-L1 * observe)) / (L1 * enrol)) ^ (-1)
		P32 <- L2 ^ 2 * (1 - (exp(-L2 * (observe - enrol)) - exp(-L2 * observe)) / (L2 * enrol)) ^ (-1)
	}
	if (enrol == 0) {
		ZB <- (sqrt(sample * (1 + ratio)) * abs(L1 - L2) - ZA * sqrt(P10 * (1 / Q2 + 1 / Q1))) / sqrt(P11/Q1 + P12/Q2)
#			N <- ((ZA * sqrt(P10 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P11/Q1 + P12/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
	} else if (enrol > 0 && observe > enrol) {
		ZB <- (sqrt(sample * (1 + ratio)) * abs(L1 - L2) - ZA * sqrt(P30 * (1 / Q2 + 1 / Q1))) / sqrt(P31/Q1 + P32/Q2)
#			N <- ((ZA * sqrt(P30 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P31/Q1 + P32/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
	} else if (enrol == observe) {
#		ZB <- (sqrt(sample * (1 + ratio)) * abs(L1 - L2) - ZA * sqrt(P20 * (1 / Q2 + 1 / Q1))) / sqrt(P21/Q1 + P22/Q2)		
		ZB <- (sqrt(sample * (1 + ratio)) * abs(L1 - L2) - ZA * sqrt(P30 * (1 / Q2 + 1 / Q1))) / sqrt(P31/Q1 + P32/Q2)
#			N <- ((ZA * sqrt(P20 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P21/Q1 + P22/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
	}
	N1 <- sample
	N2 <- round(N1 * ratio, 0)
	power <- signif(pnorm(ZB), digits=3)
	res <- data.frame(c(group1, group2, followup, enrol, observe, alpha, side, " ", N1, N2, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), power))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P1", "P2", "(Follow-up duration for P1, P2)", "Enrollment duration", "Total study duration", "Alpha", " ", "Sample size", "N1", "N2", "  ", "   ", "Power"))

	x <- NULL	
	y <- seq(0.2, 1, 0.05) 
	for (i in 1: length(y)){
		ZB <- qnorm(y[i])
		if (enrol == 0) {
			x[i] <- ((ZA * sqrt(P10 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P11/Q1 + P12/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
		} else if (enrol > 0 && observe > enrol) {
			x[i] <- ((ZA * sqrt(P30 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P31/Q1 + P32/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
		} else if (enrol == observe) {
#			x[i] <- ((ZA * sqrt(P20 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P21/Q1 + P22/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
			x[i] <- ((ZA * sqrt(P30 * (1 / Q2 + 1 / Q1)) + ZB * sqrt(P31/Q1 + P32/Q2)) / abs(L1 - L2)) ^ 2 / (1 + ratio)
		}
	}
	plot(x, y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(v=sample, lty=2)

	return(res)
}		


SampleMeanNonInf <- function (difference, delta, sd, alpha, power, method) {
	#method = 2 for two sided, 1 for one sided, r for group2/group1 ratio
	#from Jitsuyo SAS Seibutsu Tokei Handbook
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
#	alpha <- alpha / method
#	ZA <- qnorm(1-alpha)
	ZB <- qnorm(power)
	N1 <- ceiling(2*(((ZA+ZB)/((difference+delta)/sd))^2))
	N2 <- N1

	res <- data.frame(c(difference, delta, sd, alpha, side, power, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N1, N2))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("Difference in means", "Delta", "Standard deviation", "Alpha", " ", "Power", "  ", "Required sample size", "N1", "N2"))
	y <- seq(0.2, 1, 0.05) 
	plot(2*(((ZA+qnorm(y))/((difference+delta)/sd))^2), y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(h=power, lty=2)
	return(res)
}


SampleHazardNonInf <- function (enrol, observe, followup, group1, group2, lowerlimit, alpha, power, method, ratio) {
	#From SWOG https://stattools.crab.org/
	side <- ifelse(method == 1, "one-sided", "two-sided")
	side <- gettext(domain="R-RcmdrPlugin.EZR",side)
	alpha2 <- alpha / method
	ZA <- qnorm(1-alpha2)
#	alpha <- alpha / method
#	ZA <- qnorm(1-alpha)
	ZB <- qnorm(power)
	L1 <- -log(group1) / followup		#L1 hazard rate 1 (control)
	L2 <- -log(group2) / followup		#L2 hazard rate 2 (experimental), thus L2/L1 is the hazard ratio
	L3 <- -log(lowerlimit) / followup	#L3 hazard of lower limit

	Q1 <- 1 / (1 + ratio)
	Q2 <- ratio / (1 + ratio)
	E11 <- 1 - exp(-L1 * observe)									#Event rate in 1 if (enrol == 0)
	E12 <- 1 - exp(-L2 * observe)									#Event rate in 1 if (enrol == 0)
	E21 <- (L1 * observe - 1 + exp(-L1 * observe)) / (L1*followup)			#Event rate in 1 if (enrol == observe)
	E22 <- (L2 * observe - 1 + exp(-L2 * observe)) / (L2*followup)			#Event rate in 1 if (enrol == observe)
	E31 <- 1 - (exp(-L1 * (observe - enrol)) - exp(-L1 * observe)) / (L1 * enrol)	#Event rate in 1 if (enrol > 0 && observe > enrol)
	E32 <- 1 - (exp(-L2 * (observe - enrol)) - exp(-L2 * observe)) / (L2 * enrol)	#Event rate in 2 if (enrol > 0 && observe > enrol)
	Delta <- log((L2/L1) / (L3/L1))
	if (enrol == 0) {
		E1 <- E11; E2 <- E12
	} else if (enrol > 0 && observe > enrol) {
		E1 <- E31; E2 <- E32
	} else if (enrol == observe) {
#		E1 <- E21; E2 <- E22
		E1 <- E31; E2 <- E32
	}

	N <- (ZA + ZB)^2 * ((1/(Q1*E1)) + (1/(Q2*E2))) / (Delta^2)
	N1 <- ceiling(N * Q1)
	N2 <- ceiling(N * Q2)

	res <- data.frame(c(group1, group2, lowerlimit, followup, enrol, observe, alpha, side, power, ratio, " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), N1, N2))
	colnames(res) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")
	rownames(res) <- gettext(domain="R-RcmdrPlugin.EZR",c("P1", "P2", "Non-inferiority lower limit", "(Follow-up duration for P1, P2)", "Enrollment duration", "Total study duration", "Alpha", " ", "Power", "N2/N1", "  ", "Required sample size", "N1", "N2"))

	x <- NULL	
	y <- seq(0.2, 1, 0.05) 
	for (i in 1: length(y)){
		ZB <- qnorm(y[i])
		x[i] <- ((ZA + ZB)^2 * ((1/(Q1*E1)) + (1/(Q2*E2))) / (Delta^2)) * Q1
	}
	plot(x, y, ylim=c(0,1), type="l", ylab="Power", xlab="N1")
	abline(h=power, lty=2)

	return(res)	
}


StatMedGroupsBox <- defmacro(recall=NULL, label=gettext(domain="R-RcmdrPlugin.EZR","Plot by:"), initialLabel=gettext(domain="R-RcmdrPlugin.EZR","Plot by groups"),
	plotLinesByGroup=FALSE, positionLegend=FALSE, plotLinesByGroupsText=gettext(domain="R-RcmdrPlugin.EZR","Plot lines by group"),
	expr={
		env <- environment()
		.groups <- FALSE
		.linesByGroup <- FALSE
		.groupsLabel <- tclVar(paste(initialLabel, "...", sep=""))
		.factors <- Variables() 
		onGroups <- function(){
			if (length(.factors) == 0){
				errorCondition(recall=recall, message=gettext(domain="R-RcmdrPlugin.EZR","There are no factors in the active data set."))
				return()
			}
			initializeDialog(subdialog, title=gettext(domain="R-RcmdrPlugin.EZR","Groups"))
			groupsBox <- variableListBox(subdialog, .factors, title=gettext(domain="R-RcmdrPlugin.EZR","Groups variable (pick one)"))
			if (plotLinesByGroup){
				linesByGroupFrame <- tkframe(subdialog)
				linesByGroup <- tclVar("1")
				linesCheckBox <- tkcheckbutton(linesByGroupFrame, variable=linesByGroup)
				tkgrid(labelRcmdr(linesByGroupFrame, text=plotLinesByGroupsText), linesCheckBox, sticky="w")
			}
			onOKsub <- function() {
				groups <- getSelection(groupsBox)
				if (length(groups) == 0){
					assign(".groups", FALSE, envir=env)
					tclvalue(.groupsLabel) <- paste(initialLabel, "...", sep="")
					tkconfigure(groupsButton, foreground="black")
					if (GrabFocus()) tkgrab.release(subdialog)
					tkdestroy(subdialog)
					tkwm.deiconify(top)
					if (GrabFocus()) tkgrab.set(top)
					tkfocus(top)
					tkwait.window(top)
					return()
				}
				assign(".groups", groups, envir=env)
				tclvalue(.groupsLabel) <- paste(label, groups)
				tkconfigure(groupsButton, foreground="blue")
				if (plotLinesByGroup) {
					lines <- as.character("1" == tclvalue(linesByGroup))
					assign(".linesByGroup", lines, envir=env)
				}
				if (GrabFocus()) tkgrab.release(subdialog)
				tkdestroy(subdialog)
				tkwm.deiconify(top)
				if (GrabFocus()) tkgrab.set(top)
				tkfocus(top)
				tkwait.window(top)
			}
			subOKCancelHelp()
			tkgrid(getFrame(groupsBox), sticky="nw")
			if (plotLinesByGroup) tkgrid(linesByGroupFrame, sticky="w")
			tkgrid(subButtonsFrame, sticky="w")
			if (positionLegend) tkgrid(labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Position legend with mouse click"), fg="blue"))
			dialogSuffix(subdialog, onOK=onOKsub, rows=3+plotLinesByGroup+positionLegend, columns=2, focus=subdialog, force.wait=TRUE)
		}
		groupsFrame <- tkframe(top)
		groupsButton <- tkbutton(groupsFrame, textvariable=.groupsLabel, command=onGroups, borderwidth=3)
		tkgrid(labelRcmdr(groupsFrame, text="    "), groupsButton, sticky="w")
	})

	
StatMedSubsetBox <- defmacro(window=top, model=FALSE,
	expr={
		subsetVariable <- if (model){
				if (currentModel && currentFields$subset != "")
					tclVar(currentFields$subset) else tclVar(gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>"))
			}
			else tclVar(gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>"))

#			subsetVariable <- ifelse (!is.null(currentFields$subset) & currentFields$subset != "", 
#					tclVar(currentFields$subset),
#					tclVar(gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")))
			
			subsetFrame <- tkframe(window)
		subsetEntry <- ttkentry(subsetFrame, width="60", textvariable=subsetVariable)
		subsetScroll <- ttkscrollbar(subsetFrame, orient="horizontal",
			command=function(...) tkxview(subsetEntry, ...))
		tkconfigure(subsetEntry, xscrollcommand=function(...) tkset(subsetScroll, ...))
		tkgrid(labelRcmdr(subsetFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Condition to limit samples for analysis. Ex1. age>50 & Sex==0  Ex2. age<50 | Sex==1"), foreground="blue"), sticky="w")
		tkgrid(subsetEntry, sticky="w")
		tkgrid(subsetScroll, sticky="ew")
	})


StatMedModelFormula <- defmacro(frame=top, hasLhs=TRUE, expr={
		checkAddOperator <- function(rhs){
			rhs.chars <- rev(strsplit(rhs, "")[[1]])
			if (length(rhs.chars) < 1) return(FALSE)
			check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
					rhs.chars[1] else rhs.chars[2]
			!is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
		}
		.variables <- Variables()
		word <- paste("\\[", gettext(domain="R-RcmdrPlugin.EZR","factor"), "\\]", sep="")
		variables <- paste(.variables,
			ifelse(is.element(.variables, Factors()), paste("[", gettext(domain="R-RcmdrPlugin.EZR","factor"), "]", sep=""), ""))
		xBox <- variableListBox(frame, variables, selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables (double-click to formula)"), listHeight=10)
		onDoubleClick <- if (!hasLhs){
				function(){
					var <- getSelection(xBox)
					tkselection.clear(xBox$listbox, "0", "end")					
					if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
					tkfocus(rhsEntry)
					rhs <- tclvalue(rhsVariable)
					rhs.chars <- rev(strsplit(rhs, "")[[1]])
					check.char <- if (length(rhs.chars) > 0){
							if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
								rhs.chars[1] else rhs.chars[2]
						}
						else ""
					tclvalue(rhsVariable) <- if (rhs == "" ||
							is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
							paste(rhs, var, sep="")
						else paste(rhs, "+", var)
					tkicursor(rhsEntry, "end")
					tkxview.moveto(rhsEntry, "1")
				}
			}
			else{
				function(){
					var <- getSelection(xBox)
					which <- tkcurselection(xBox$listbox)
					tkselection.clear(xBox$listbox, "0", "end")
					if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
					lhs <- tclvalue(lhsVariable)
					if (lhs == "" || tclvalue(tkselection.present(lhsEntry)) == "1"){
						tclvalue(lhsVariable) <- var
						tkselection.clear(lhsEntry)
						tkfocus(rhsEntry)
					}
					else {
						tkfocus(rhsEntry)
						rhs <- tclvalue(rhsVariable)
						rhs.chars <- rev(strsplit(rhs, "")[[1]])
						check.char <- if (length(rhs.chars) > 0){
								if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
									rhs.chars[1] else rhs.chars[2]
							}
							else ""
						tclvalue(rhsVariable) <- if (rhs == "" ||
								is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
								paste(rhs, var, sep="")
							else paste(rhs, "+", var)
					}
					tkicursor(rhsEntry, "end")
					tkxview.moveto(rhsEntry, "1")
				}
			}
		tkbind(xBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
		onPlus <- function(){
			rhs <- tclvalue(rhsVariable)
			var <- getSelection(xBox)
			tkselection.clear(xBox$listbox, "0", "end")										
			if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
			if (length(var) > 1){
				if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
				if (length(var) > 1) var <- paste(var, collapse=" + ")
			}
			tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onTimes <- function(){
			rhs <- tclvalue(rhsVariable)
			var <- getSelection(xBox)
			tkselection.clear(xBox$listbox, "0", "end")						
			if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
			if (length(var) > 1){
				if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
				var <- trim.blanks(var)
				if (length(var) > 1) var <- paste(var, collapse="*")
				tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
			}
			else tclvalue(rhsVariable) <- paste(rhs, if (!check) "*", sep="")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onColon <- function(){
			rhs <- tclvalue(rhsVariable)
			var <- getSelection(xBox)
			tkselection.clear(xBox$listbox, "0", "end")						
			if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
			if (length(var) > 1){
				if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
				var <- trim.blanks(var)
				if (length(var) > 1) var <- paste(var, collapse=":")
				tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
			}
			else tclvalue(rhsVariable) <- paste(rhs, if (!check) ":", sep="")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onSlash <- function(){
			rhs <- tclvalue(rhsVariable)
			if (!checkAddOperator(rhs)) return()
			tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onIn <- function(){
			rhs <- tclvalue(rhsVariable)
			if (!checkAddOperator(rhs)) return()
			tclvalue(rhsVariable) <- paste(rhs, "%in% ")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onMinus <- function(){
			rhs <- tclvalue(rhsVariable)
			if (!checkAddOperator(rhs)) return()
			tclvalue(rhsVariable) <- paste(rhs, "- ")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onPower <- function(){
			rhs <- tclvalue(rhsVariable)
			if (!checkAddOperator(rhs)) return()
			tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onLeftParen <- function(){
			tkfocus(rhsEntry)
			rhs <- tclvalue(rhsVariable)
			tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		onRightParen <- function(){
			rhs <- tclvalue(rhsVariable)
			if (!checkAddOperator(rhs)) return()
			tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
			tkicursor(rhsEntry, "end")
			tkxview.moveto(rhsEntry, "1")
		}
		outerOperatorsFrame <- tkframe(frame)
		operatorsFrame <- tkframe(outerOperatorsFrame)
		plusButton <- buttonRcmdr(operatorsFrame, text="+", width="3", command=onPlus)
		timesButton <- buttonRcmdr(operatorsFrame, text="*", width="3", command=onTimes)
		colonButton <- buttonRcmdr(operatorsFrame, text=":", width="3", command=onColon)
		slashButton <- buttonRcmdr(operatorsFrame, text="/", width="3", command=onSlash)
		inButton <- buttonRcmdr(operatorsFrame, text="%in%", width="5", command=onIn)
		minusButton <- buttonRcmdr(operatorsFrame, text="-", width="3", command=onMinus)
		powerButton <- buttonRcmdr(operatorsFrame, text="^", width="3", command=onPower)
		leftParenButton <- buttonRcmdr(operatorsFrame, text="(", width="3", command=onLeftParen)
		rightParenButton <- buttonRcmdr(operatorsFrame, text=")", width="3", command=onRightParen)
		
		tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, minusButton,
			powerButton, leftParenButton, rightParenButton, sticky="w")
		formulaFrame <- tkframe(frame)
		if (hasLhs){
			tkgrid(labelRcmdr(outerOperatorsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Model Formula:     "), fg="blue"), operatorsFrame)
			lhsVariable <- if (currentModel) tclVar(currentFields$lhs) else tclVar("")
			rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
			rhsEntry <- ttkentry(formulaFrame, width="50", textvariable=rhsVariable)
			rhsXscroll <- ttkscrollbar(formulaFrame,
				orient="horizontal", command=function(...) tkxview(rhsEntry, ...))
			tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
			lhsEntry <- ttkentry(formulaFrame, width="10", textvariable=lhsVariable)
			lhsScroll <- ttkscrollbar(formulaFrame,
				orient="horizontal", command=function(...) tkxview(lhsEntry, ...))
			tkconfigure(lhsEntry, xscrollcommand=function(...) tkset(lhsScroll, ...))
			tkgrid(labelRcmdr(formulaFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Objective variable")), lhsEntry, labelRcmdr(formulaFrame, text=gettext(domain="R-RcmdrPlugin.EZR","~ Explanatory variables")), rhsEntry, sticky="w")
			tkgrid(lhsScroll, labelRcmdr(formulaFrame, text=""), rhsXscroll, sticky="w")
			tkgrid.configure(lhsScroll, sticky="ew")
		}
		else{
			rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
			rhsEntry <- ttkentry(formulaFrame, width="50", textvariable=rhsVariable)
			rhsXscroll <- ttkscrollbar(formulaFrame,
				orient="horizontal", command=function(...) tkxview(rhs, ...))
			tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
			tkgrid(labelRcmdr(formulaFrame, text="   ~ "), rhsEntry, sticky="w")
			tkgrid(labelRcmdr(formulaFrame, text=""), rhsXscroll, sticky="w")
		}
		tkgrid.configure(rhsXscroll, sticky="ew")
	})

	
modelFormulaCox <- defmacro(frame=top, hasLhs=TRUE, expr={  # from RcmdrPlugin.SurvivalT
    checkAddOperator <- function(rhs){
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        if (length(rhs.chars) < 1) return(FALSE)
        check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
                rhs.chars[1] else rhs.chars[2]
        !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
        }
    .variables <- Variables()
    word <- paste("\\[", gettext(domain="R-RcmdrPlugin.EZR","factor"), "\\]", sep="")
    variables <- paste(.variables,
        ifelse(is.element(.variables, Factors()), paste("[", gettext(domain="R-RcmdrPlugin.EZR","factor"), "]", sep=""), ""))
    xBox <- variableListBox(frame, variables, title=gettext(domain="R-RcmdrPlugin.EZR","Variables (double-click to formula)"), listHeight=8)
    onDoubleClick <- if (!hasLhs){
        function(){
            var <- getSelection(xBox)
            if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
            tkfocus(rhsEntry)
            rhs <- tclvalue(rhsVariable)
            rhs.chars <- rev(strsplit(rhs, "")[[1]])
            check.char <- if (length(rhs.chars) > 0){
                if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
                    rhs.chars[1] else rhs.chars[2]
                }
                else ""
            tclvalue(rhsVariable) <- if (rhs == "" ||
                is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
                    paste(rhs, var, sep="")
                else paste(rhs, "+", var)
            tkicursor(rhsEntry, "end")
            tkxview.moveto(rhsEntry, "1")
            }
        }
    else{
        function(){
            var <- getSelection(xBox)
            if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
            lhs <- tclvalue(SurvivalTimeVariable)
            lhs2 <- tclvalue(StatusVariable)            
            if (lhs == "") tclvalue(SurvivalTimeVariable) <- var
            else 
            {
            if (lhs2 == "") tclvalue(StatusVariable) <- var
            else {
                tkfocus(rhsEntry)
                rhs <- tclvalue(rhsVariable)
                rhs.chars <- rev(strsplit(rhs, "")[[1]])
                check.char <- if (length(rhs.chars) > 0){
                    if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
                        rhs.chars[1] else rhs.chars[2]
                    }
                    else ""
                tclvalue(rhsVariable) <- if (rhs == "" ||
                    is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
                        paste(rhs, var, sep="")
                    else paste(rhs, "+", var)
                }
            }    
            tkicursor(rhsEntry, "end")
            tkxview.moveto(rhsEntry, "1")
            }
        }
    tkbind(xBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
    onPlus <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "+ ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onTimes <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "*", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onColon <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ":", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onSlash <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onIn <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "%in% ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onMinus <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "- ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onPower <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onLeftParen <- function(){
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onRightParen <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    outerOperatorsFrame <- tkframe(frame)
    operatorsFrame <- tkframe(outerOperatorsFrame)
    plusButton <- buttonRcmdr(operatorsFrame, text="+", width="3", command=onPlus)
    timesButton <- buttonRcmdr(operatorsFrame, text="*", width="3", command=onTimes)
    colonButton <- buttonRcmdr(operatorsFrame, text=":", width="3", command=onColon)
    slashButton <- buttonRcmdr(operatorsFrame, text="/", width="3", command=onSlash)
    inButton <- buttonRcmdr(operatorsFrame, text="%in%", width="5", command=onIn)
    minusButton <- buttonRcmdr(operatorsFrame, text="-", width="3", command=onMinus)
    powerButton <- buttonRcmdr(operatorsFrame, text="^", width="3", command=onPower)
    leftParenButton <- buttonRcmdr(operatorsFrame, text="(", width="3", command=onLeftParen)
    rightParenButton <- buttonRcmdr(operatorsFrame, text=")", width="3", command=onRightParen)

    tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, minusButton,
        powerButton, leftParenButton, rightParenButton, sticky="w")
    formulaFrame <- tkframe(frame)
    if (hasLhs){
        tkgrid(labelRcmdr(outerOperatorsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Model Formula:     "), fg="blue"), operatorsFrame)
        SurvivalTimeVariable <- if (currentModel) tclVar(currentFields$SurvivalTimeVariable) else tclVar("")
        StatusVariable <- if (currentModel) tclVar(currentFields$StatusVariable) else tclVar("")
        rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
        rhsEntry <- ttkentry(formulaFrame, width="50", textvariable=rhsVariable)
        rhsXscroll <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(rhs, ...))
        tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
        lhsEntry <- ttkentry(formulaFrame, width="10", textvariable=SurvivalTimeVariable)
        lhsScroll <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(lhsEntry, ...))
        tkconfigure(lhsEntry, xscrollcommand=function(...) tkset(lhsScroll, ...))
        lhsEntry2 <- ttkentry(formulaFrame, width="10", textvariable=StatusVariable)
        lhsScroll2 <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(lhsEntry2, ...))
        tkconfigure(lhsEntry2, xscrollcommand=function(...) tkset(lhsScroll2, ...))
        tkgrid(labelRcmdr(formulaFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Time")), lhsEntry, labelRcmdr(formulaFrame, text=gettext(domain="R-RcmdrPlugin.EZR",", Event")), lhsEntry2, labelRcmdr(formulaFrame, text=gettext(domain="R-RcmdrPlugin.EZR","~ Explanatory variables")), rhsEntry, sticky="w")
        tkgrid(labelRcmdr(formulaFrame, text=""), lhsScroll, labelRcmdr(formulaFrame, text=""), lhsScroll2, labelRcmdr(formulaFrame, text=""), rhsXscroll, sticky="w")
        tkgrid.configure(lhsScroll, sticky="ew")
        }
    else{
        rhsVariable <- tclVar("")
        rhsEntry <- ttkentry(formulaFrame, width="50", textvariable=rhsVariable)
        rhsXscroll <- ttkscrollbar(formulaFrame,
            orient="horizontal", command=function(...) tkxview(rhs, ...))
        tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
        tkgrid(labelRcmdr(formulaFrame, text="   ~ "), rhsEntry, sticky="w")
        tkgrid(labelRcmdr(formulaFrame, text=""), rhsXscroll, sticky="w")
        }
    tkgrid.configure(rhsXscroll, sticky="ew")
    })
    
		
listCoxModels <- function(envir=.GlobalEnv, ...) {  # from RcmdrPlugin.SurvivalT
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
#        function(.x) "coxph" == (class(eval(parse(text=.x), envir=envir))[1]))]
        function(.x) "coxph" == (class(get(.x, envir=envir))[1]))]
    }		

	
listLMModels <- function(envir=.GlobalEnv, ...) {
	objects <- ls(envir=envir, ...)
	if (length(objects) == 0) NULL
	else objects[sapply(objects,
				function(.x) "lm" == (class(get(.x, envir=envir))[1]))]
}


StatMedLoadDataSet <- function() {
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Load data set"), "#####", sep=""))
#	file <- tclvalue(tkgetOpenFile(filetypes=
#							gettext(domain="R-RcmdrPlugin.EZR",'{"R Data Files" {".RData" ".rda" ".Rda" ".RDA"}} {"All Files" {"*"}}')))	
	file <- tclvalue(tkgetOpenFile(filetypes=
							gettextRcmdr('{"All Files" {"*"}} {"R Data Files" {".RData" ".rda" ".Rda" ".RDA"}}')))
	if (file == "") return()
	setBusyCursor()
	on.exit(setIdleCursor())
	command <- paste('load("', file,'")', sep="")
	dsname <- justDoIt(command)
	logger(command)
	if (class(dsname)[1] !=  "try-error") activeDataSet(dsname)
	tkfocus(CommanderWindow())

}
	
	
StatMedReadDataSet <- function() {
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Read Text Data From File, Clipboard, or URL"))
	optionsFrame <- tkframe(top)
	dsname <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","Dataset"))
	entryDsname <- ttkentry(optionsFrame, width="20", textvariable=dsname)
	radioButtons(optionsFrame, "location", buttons=c("local", "clipboard", "url"), 
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("Local file system", "Clipboard", "Internet URL")), title=gettext(domain="R-RcmdrPlugin.EZR","Location of Data File"))
	headerVariable <- tclVar("1")
	headerCheckBox <- tkcheckbutton(optionsFrame, variable=headerVariable)
	fillVariable <- tclVar("1")
	fillCheckBox <- tkcheckbutton(optionsFrame, variable=fillVariable)
	blankVariable <- tclVar("1")
	blankCheckBox <- tkcheckbutton(optionsFrame, variable=blankVariable)
	##   clipboardVariable <- tclVar("0")
	##   clipboardCheckBox <- tkcheckbutton(optionsFrame, variable=clipboardVariable)
	radioButtons(optionsFrame, "delimiter", buttons=c("whitespace", "commas", "tabs"), initialValue="commas",
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("White space", "Commas", "Tabs")), title=gettext(domain="R-RcmdrPlugin.EZR","Field Separator"))
	otherButton <- ttkradiobutton(delimiterFrame, variable=delimiterVariable, value="other")
	otherVariable <- tclVar("")
	otherEntry <- ttkentry(delimiterFrame, width="4", textvariable=otherVariable)
	radioButtons(optionsFrame, "decimal", buttons=c("period", "comma"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("Period [.]", "Comma [,]")), title=gettext(domain="R-RcmdrPlugin.EZR","Decimal-Point Character"))
	missingVariable <- tclVar("NA")
	missingEntry <- ttkentry(optionsFrame, width="8", textvariable=missingVariable)
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Read Data From Text File"), "#####", sep=""))
		closeDialog()
		dsnameValue <- trim.blanks(tclvalue(dsname))
		if (dsnameValue == ""){
			errorCondition(recall=StatMedReadDataSet,
				message=gettext(domain="R-RcmdrPlugin.EZR","You must enter a name for the data set."))
			return()
		}
		if (!is.valid.name(dsnameValue)){
			errorCondition(recall=StatMedReadDataSet,
				message=paste('"', dsnameValue, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				StatMedReadDataSet()
				return()
			}
		}
		##        clip <- tclvalue(clipboardVariable) == "1"
		location <- tclvalue(locationVariable)
		file <- if (location == "clipboard") "clipboard" 
			else if (location == "local") tclvalue(tkgetOpenFile(filetypes=
							gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"Text Files" {".txt" ".TXT" ".dat" ".DAT" ".csv" ".CSV"}}')))
			else {
				initializeDialog(subdialog, title=gettext(domain="R-RcmdrPlugin.EZR","Internet URL"))
				onOKsub <- function(){
					closeDialog(subdialog)
				}
				urlFrame <- tkframe(subdialog)
				urlVar <- tclVar("")
				url <- ttkentry(urlFrame, font=getRcmdr("logFont"), width="30", textvariable=urlVar)
				urlXscroll <- ttkscrollbar(urlFrame,
					orient="horizontal", command=function(...) tkxview(url, ...))
				tkconfigure(url, xscrollcommand=function(...) tkset(urlXscroll, ...))
				subOKCancelHelp()
				tkgrid(url, sticky="w")
				tkgrid(urlXscroll, sticky="ew")
				tkgrid(urlFrame, sticky="nw")
				tkgrid(subButtonsFrame, sticky="w")
				dialogSuffix(subdialog, rows=2, columns=1, focus=url, onOK=onOKsub, force.wait=TRUE)
				tclvalue(urlVar)
			}
		if (file == "") {
			if (getRcmdr("grab.focus")) tkgrab.release(top)
			tkdestroy(top)
			return()
		}
		head <- tclvalue(headerVariable) == "1"
		fill <- tclvalue(fillVariable)
		if (fill == 0){
			fill <- ""
		}else{
			fill <- ", fill=TRUE"
		}
		delimiter <- tclvalue(delimiterVariable)
		del <- if (delimiter == "whitespace") ""
			else if (delimiter == "commas") ","
			else if (delimiter == "tabs") "\\t"
			else tclvalue(otherVariable)
		blank <- tclvalue(blankVariable)
		miss <- tclvalue(missingVariable)
		if (blank == 1) {
			miss <- paste('c("", "', miss, '")', sep="")			
		} else {
			miss <- paste('"', miss, '"', sep="")
		}
		dec <- if (tclvalue(decimalVariable) == "period") "." else ","
		command <- paste('read.table("', file,'", header=', head,
			', sep="', del, '", na.strings=', miss, ', dec="', dec, '"', fill, ', quote="\\"", strip.white=TRUE)', sep="")
		logger(paste(dsnameValue, " <- ", command, sep=""))
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error"){
# 			assign(dsnameValue, result, envir=.GlobalEnv)
# 			logger(paste(dsnameValue, "<-", command))
		    doItAndPrint(paste(dsnameValue, "<-", command))
			activeDataSet(dsnameValue)
		}
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="read.table")
	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for data set:")), entryDsname, sticky="w")
	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Variable names in file:")), headerCheckBox, sticky="w")
	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Adjust for different column numbers:")), fillCheckBox, sticky="w")
	##    tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Read data from clipboard:")), clipboardCheckBox, sticky="w")
	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Apply NA for blank cells in character variables:")), blankCheckBox, sticky="w")
	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Characters indicating NA cells:")), missingEntry, sticky="w")
	tkgrid(locationFrame, sticky="w")
	tkgrid(labelRcmdr(delimiterFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Other")), otherButton,
		labelRcmdr(delimiterFrame, text=gettext(domain="R-RcmdrPlugin.EZR","  Specify:")), otherEntry, sticky="w")
	tkgrid(delimiterFrame, sticky="w", columnspan=2)
	tkgrid(decimalFrame, sticky="w")
	tkgrid(optionsFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=5, columns=1)
}


StatMedImportSPSS <- function() {
	Library("foreign")
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Import SPSS Data Set"))
	dsname <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","Dataset"))
	entryDsname <- ttkentry(top, width="20", textvariable=dsname)
	asFactor <- tclVar("1")
	asFactorCheckBox <- tkcheckbutton(top, variable=asFactor)
	maxLevels <- tclVar("Inf")
	entryMaxLevels <- ttkentry(top, width="5", textvariable=maxLevels)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Import SPSS Data Set"), "#####", sep=""))
		closeDialog()
		dsnameValue <- trim.blanks(tclvalue(dsname))
		if (dsnameValue == ""){
			errorCondition(recall=StatMedImportSPSS,
				message=gettext(domain="R-RcmdrPlugin.EZR","You must enter the name of a data set."))
			return()
		}
		if (!is.valid.name(dsnameValue)){
			errorCondition(recall=StatMedImportSPSS,
				message=paste('"', dsnameValue, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				importSPSS()
				return()
			}
		}
		file <- tclvalue(tkgetOpenFile(
				filetypes=gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"SPSS files" {".sav" ".SAV" ".por" ".POR"}}')))
		if (file == "") {
			tkfocus(CommanderWindow())
			return()
		}
		factor <- tclvalue(asFactor) == "1"
		levels <- as.numeric(tclvalue(maxLevels))
		command <- paste('read.spss("', file,'", use.value.labels=', factor,
			", max.value.labels=", levels, ", to.data.frame=TRUE)", sep="")
		logger(paste(dsnameValue, " <- ", command, sep=""))
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error"){
# 			assign(dsnameValue, result, envir=.GlobalEnv)
# 			logger(paste(dsnameValue, "<-", command))
		    doItAndPrint(paste(dsnameValue, "<-", command))
			activeDataSet(dsnameValue)
		}
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="read.spss")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for data set:")), entryDsname, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Convert value labels\nto factor levels"), justify="left"),
		asFactorCheckBox, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Maximum number\nof value labels\nfor factor conversion"), justify="left"),
		entryMaxLevels, sticky="w")
	tkgrid(buttonsFrame, columnspan="2", sticky="w")
	tkgrid.configure(entryDsname, sticky="w")
	tkgrid.configure(asFactorCheckBox, sticky="w")
	tkgrid.configure(entryMaxLevels, sticky="w")
	dialogSuffix(rows=4, columns=2, focus=entryDsname)
}


StatMedImportMinitab <- function() {
	Library("foreign")
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Import Minitab Data Set"))
	dsname <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","Dataset"))
	entryDsname <- ttkentry(top, width="20", textvariable=dsname)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Import Minitab Data Set"), "#####", sep=""))
		closeDialog()
		dsnameValue <- trim.blanks(tclvalue(dsname))
		if (dsnameValue == ""){
			errorCondition(recall=StatMedImportMinitab,
				message=gettext(domain="R-RcmdrPlugin.EZR","You must enter the name of a data set."))
			return()
		}
		if (!is.valid.name(dsnameValue)){
			errorCondition(recall=StatMedImportMinitab,
				message=paste('"', dsnameValue, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				importMinitab()
				return()
			}
		}
		file <- tclvalue(tkgetOpenFile(
				filetypes=gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"Minitab portable files" {".mtp" ".MTP"}}')))
		if (file == "") {
			tkfocus(CommanderWindow())
			return()
		}
		command <- paste('read.mtp("', file,'")', sep="")
		datalist <- justDoIt(command)
		lengths <- sapply(datalist, length)
		datalist <- datalist[lengths != 0]
		lengths <- lengths[lengths != 0]
		if (!all(lengths == length(datalist[[1]]))){
			Message(message=
					paste(gettext(domain="R-RcmdrPlugin.EZR","Minitab data set contains elements of unequal length.\nData set cannot be converted.")),
				type="error")
			tkdestroy(top)
			tkfocus(CommanderWindow())
			return()
		}
# 		assign(dsnameValue, as.data.frame(datalist), envir=.GlobalEnv)
# 		logger(paste(dsnameValue, " <- as.data.frame(", command, ")", sep=""))
		doItAndPrint(paste(dsnameValue, " <- as.data.frame(", command, ")", sep=""))
		activeDataSet(dsnameValue)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="read.mtp")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for data set:")), entryDsname, sticky="e")
	tkgrid(buttonsFrame, columnspan="2", sticky="w")
	tkgrid.configure(entryDsname, sticky="w")
	dialogSuffix(rows=2, columns=2, focus=entryDsname)
}


StatMedImportSTATA <- function() {
	Library("foreign")
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Import Stata Data Set"))
	dsname <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","Dataset"))
	entryDsname <- ttkentry(top, width="20", textvariable=dsname)
	asFactor <- tclVar("1")
	asFactorCheckBox <- tkcheckbutton(top, variable=asFactor)
	asDate <- tclVar("1")
	asDateCheckBox <- tkcheckbutton(top, variable=asDate)
	asMissingType <- tclVar("1")
	asMissingTypeCheckBox <- tkcheckbutton(top, variable=asMissingType)
	asConvertUnderscore <- tclVar("1")
	asConvertUnderscoreCheckBox <- tkcheckbutton(top, variable=asConvertUnderscore)
	asWarnMissingLabels <- tclVar("1")
	asWarnMissingLabelsCheckBox <- tkcheckbutton(top, variable=asWarnMissingLabels)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Import Stata Data Set"), "#####", sep=""))	
		closeDialog()
		dsnameValue <- trim.blanks(tclvalue(dsname))
		if (dsnameValue == ""){
			errorCondition(recall=StatMedImportSTATA,
				message=gettext(domain="R-RcmdrPlugin.EZR","You must enter the name of a data set."))
			return()
		}
		if (!is.valid.name(dsnameValue)){
			errorCondition(recall=StatMedImportSTATA,
				message=paste('"', dsnameValue, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				importSTATA()
				return()
			}
		}
		file <- tclvalue(tkgetOpenFile(
				filetypes=gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"Stata datasets" {".dta" ".DTA"}}')))
		if (file == "") {
			tkfocus(CommanderWindow())
			return()
		}
		convert.date <- tclvalue(asDate) == "1"
		factor <- tclvalue(asFactor) == "1"
		missingtype <- tclvalue(asMissingType) == "1"
		convertunderscore <- tclvalue(asConvertUnderscore) == "1"
		warnmissinglabels <- tclvalue(asWarnMissingLabels) == "1"
		command <- paste('read.dta("', file,'", convert.dates=', convert.date,
			", convert.factors=", factor, ", missing.type=", missingtype,
			", convert.underscore=", convertunderscore, ", warn.missing.labels=TRUE)", sep="")
		logger(paste(dsnameValue, " <- ", command, sep=""))
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error"){
# 			assign(dsnameValue, result, envir=.GlobalEnv)
# 			logger(paste(dsnameValue, "<-", command))
		    doItAndPrint(paste(dsnameValue, "<-", command))
			activeDataSet(dsnameValue)
		}
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="read.dta")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for data set:")), entryDsname, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Convert value labels\nto factor levels"), justify="left"),
		asFactorCheckBox, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Convert dates to R format"), justify="left"),
		asDateCheckBox, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Multiple missing types (>=Stata 8)"), justify="left"),
		asMissingTypeCheckBox, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Convert underscore to period"), justify="left"),
		asConvertUnderscoreCheckBox, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Warn on missing labels"), justify="left"),
		asWarnMissingLabelsCheckBox, sticky="w")
	tkgrid(buttonsFrame, columnspan="2", sticky="w")
	tkgrid.configure(entryDsname, sticky="w")
	tkgrid.configure(asFactorCheckBox, sticky="w")
	tkgrid.configure(asDateCheckBox, sticky="w")
	tkgrid.configure(asMissingTypeCheckBox, sticky="w")
	tkgrid.configure(asWarnMissingLabelsCheckBox, sticky="w")
	dialogSuffix(rows=4, columns=2, focus=entryDsname)
}


StatMedLoadWorkspace <- function() {
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Load work space file"), "#####", sep=""))
	file <- tclvalue(tkgetOpenFile(filetypes=
				gettext(domain="R-RcmdrPlugin.EZR",'{"R Data Files" {".RData"}} {"All Files" {"*"}}')))
	if (file == "") return()
	command <- paste('load("', file,'")', sep="")
	dsname <- justDoIt(command)
	logger(command)
	tkfocus(CommanderWindow())
}
	
	
trim.col.na <- function(dat){
# Remove variables with only missing values (occurs sometimes with modified Excel file)
	colsup <- NULL
	for (i in 1:ncol(dat))
	{
		if (length(dat[is.na(dat[,i])==T,i]) ==length(dat[,i]))
			colsup <- c(colsup,i)
	}
	if (length(colsup) > 0)
		dat <- dat[,-colsup]
	dat
}


if(.Platform$OS.type == 'windows')
StatMedImportRODBCtable <- function(){
	# load the RODBC package and stops the program if not available
	Library("RODBC")
	#if(!require(RODBC))
	#	stop("This function requires the RODBC package.")
	# close all databases in case of error
	on.exit(odbcCloseAll())
# Enter the name of data set, by default : Dataset
	initializeDialog(title = gettext(domain="R-RcmdrPlugin.EZR","Import from Excel, Access or dBase data set"))
	dsname <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","Dataset"))
	entryDsname <- ttkentry(top, width = "35", textvariable = dsname)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Import from Excel, Access or dBase data set"), "#####", sep=""))
		closeDialog()
		dsnameValue <- trim.blanks(tclvalue(dsname))
		if(dsnameValue == ""){
			errorCondition(recall = StatMedImportRODBCtable,
				message = gettext(domain="R-RcmdrPlugin.EZR","You must enter the name of a data set."))
			return()
		}
		if(!is.valid.name(dsnameValue)){
			errorCondition(recall = StatMedImportRODBCtable,
				message = paste('"', dsnameValue, '" ',
					gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep = ""))
			return()
		}
		if(is.element(dsnameValue, listDataSets())){
			if("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				importRODBCtable()
				return()
			}
		}
		File <- tclvalue(tkgetOpenFile(filetypes = gettext(domain="R-RcmdrPlugin.EZR",
					'{"All Files" {"*"}} {"MS Access database" {*.mdb ".MDB"}} {"MS Access 2007 database" {*.accdb ".ACCDB"}} {"dBase-like file" {*.dbf ".DBF"}} {"MS Excel file" {*.xls ".XLS" *.xlsx ".XLSX"}}'
					)))
		if(File == ""){
			tkfocus(CommanderWindow())
			return()
		}
		sop <- match(".", rev(strsplit(File, NULL)[[1]]))[1]
		ext <- tolower(substring(File, nchar(File) - sop + 2, nchar(File)))
		channel <- switch(EXPR = ext,
			xls = odbcConnectExcel(File),
			xlsx = odbcConnectExcel2007(File),
			mdb = odbcConnectAccess(File),
			accdb = odbcConnectAccess2007(File),
			dbf = odbcConnectDbase(File))
		# For Excel and Access cases, need to select a particular sheet or table
		if(ext != "dbf"){
			tabdat <- sqlTables(channel)
			names(tabdat) <- tolower(names(tabdat))
			if(ext == "mdb" || ext == "accdb")
				tabdat <- tabdat[tabdat$table_type == "TABLE", 3]
			if(ext == "xls" || ext == "xlsx"){
				tabname <- tabdat$table_name
				tabdat <- ifelse(tabdat$table_type =="TABLE",
					substring(tabname, 2, nchar(tabname) - 2),
					substring(tabname, 1, nchar(tabname) - 1))
			}
			# if there are several tables
			if(length(tabdat)>1)
				fil <- tk_select.list(sort(tabdat),
					title = gettext(domain="R-RcmdrPlugin.EZR","Select one table"))
			else
				fil <- tabdat
			if(fil == ""){
				errorCondition(message=gettext(domain="R-RcmdrPlugin.EZR","No table selected"))
				return()
			}
			if(ext == "xls" || ext == "xlsx")
				fil <- paste("[", fil, "$]", sep = "")
		}
		# dBase file
		else{
			sop <- match(".", rev(strsplit(File, NULL)[[1]]))[1]
			root <- tolower(substring(File, 1, nchar(File) - sop))
			revstr <- rev(strsplit(root, NULL)[[1]])
			sop <- if(is.na(match(c("/", "\\"), revstr)[1]))
					length(revstr) else match(c("/", "\\"), revstr)[1] - 1
			toor <- revstr[seq(sop)]
			fil <- paste(rev(toor), collapse = "")
		}
		# Retrieve the data
		dat <- sqlQuery(channel = channel, query = paste("select * from", fil))
		names(dat)<- trim.blanks(names(dat))
		dat <- trim.col.na(dat)
		odbcCloseAll()
		gassign(dsnameValue, as.data.frame(dat))
		command <- paste("sqlQuery(channel = ",channel,", select * from ", fil,")",
				sep = "")
		logger(paste(dsnameValue, " <- ", command, sep = ""))
		activeDataSet(dsnameValue)
		tkfocus(CommanderWindow())
	}  ## End of function onOK
	OKCancelHelp(helpSubject="odbcConnect")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name of data set:  ")),
		entryDsname, sticky="e")
	tkgrid(buttonsFrame, columnspan="2", sticky="w")
	tkgrid.configure(entryDsname, sticky="w")
	dialogSuffix(rows=2, columns=2, focus=entryDsname)
}


StatMedImportExcel <- function(){
    Library("XLConnect")
	Library("methods")
    initializeDialog(title = gettextRcmdr("Import Excel Data Set"))
    dsname <- tclVar(gettextRcmdr("Dataset"))
    entryDsname <- ttkentry(top, width = "35", textvariable = dsname)
    onOK <- function(){
        closeDialog()
	setBusyCursor()
	on.exit(setIdleCursor())
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if(dsnameValue == ""){
            errorCondition(recall = StatMedImportExcel,
                           message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if(!is.valid.name(dsnameValue)){
            errorCondition(recall = StatMedImportExcel,
                           message = paste('"', dsnameValue, '" ',
                                           gettextRcmdr("is not a valid name."), sep = ""))
            return()
        }
        if(is.element(dsnameValue, listDataSets())){
            if("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
                StatMedImportExcel()
                return()
            }
        }
#        File <- tclvalue(tkgetOpenFile(filetypes = gettextRcmdr(
#            '{"All Files" {"*"}} {{"MS Excel file" {".xls" ".XLS"}} "MS Excel 2007 file" {".xlsx" ".XLSX"}}'
#        ), parent=CommanderWindow()))
        File <- tclvalue(tkgetOpenFile(filetypes = gettextRcmdr(
            '{"All Files" {"*"}} {"MS Excel 2007 file" {".xlsx" ".XLSX"}} {"MS Excel file" {".xls" ".XLS"}}'
        ), parent=CommanderWindow()))
        if(File == ""){
            tkfocus(CommanderWindow())
            return()
        }
        command <- paste('loadWorkbook("', File, '")', sep="")
        doItAndPrint(paste(".Workbook <- ", command, sep=""))
        worksheets <- getSheets(.Workbook)
        if(length(worksheets)>1)
            worksheet <- tk_select.list(worksheets,
                                        title = gettextRcmdr("Select one table"))
        else
            worksheet <- worksheets
        if(worksheet == ""){
            errorCondition(message=gettextRcmdr("No table selected"))
            return()
        }
        command <- paste('readWorksheet(.Workbook, "', worksheet, '")', sep="")
        logger(paste(dsnameValue, " <- ", command, sep=""))
        result <- justDoIt(command)
        if (class(result)[1] !=  "try-error"){
            gassign(dsnameValue, result)
        }
        logger("remove(.Workbook)")
        justDoIt("remove(.Workbook, envir=.GlobalEnv)")
        if (class(result)[1] !=  "try-error"){
            factors <- sapply(get(dsnameValue, envir=.GlobalEnv), is.character)
            if (any(factors)){
                factors <- which(factors)
                command <- paste(dsnameValue, "[, c(", paste(factors, collapse=", "), 
                                 ")] <- lapply(", dsnameValue, "[, c(", 
                                 paste(factors, collapse=", "), "), drop=FALSE], as.factor)",
                                 sep="")
                doItAndPrint(command)
            }
            activeDataSet(dsnameValue)
        }
    }
    OKCancelHelp(helpSubject="readWorksheet")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name of data set:  ")),
           entryDsname, sticky="e")
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
    tkgrid.configure(entryDsname, sticky="w")
    dialogSuffix(focus=entryDsname)
}


StatMedCopyDataset <- function(){
    dataSets <- listDataSets()
    .activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR", "Copy data set"))
    dsname <- tclVar("NewDataset")
    dsnameFrame <- tkframe(top)
    entryDsname <- ttkentry(dsnameFrame, width="20", textvariable=dsname)
    dataSet1Box <- variableListBox(top, dataSets, title=gettext(domain="R-RcmdrPlugin.EZR","Original Data Set"),
                                   initialSelection=if (is.null(.activeDataSet)) NULL else which(.activeDataSet == dataSets) - 1)
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Copy data set"), "#####", sep=""))
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == getSelection(dataSet1Box)) {
            errorCondition(recall=StatMedCopyDataset,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must enter a different data set name."))
            return()
        }
        if (dsnameValue == "") {
            errorCondition(recall=StatMedCopyDataset,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall=StatMedCopyDataset,
                           message=paste('"', dsnameValue, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
                closeDialog()
                StatMedCopyDataset()
                return()
            }
        }
        name1 <- getSelection(dataSet1Box)
        if (length(name1) == 0){
            errorCondition(recall=StatMedCopyDataset,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must select a data set."))
            return()
        }
        command <- paste(dsnameValue, " <- ", name1, sep="")
        doItAndPrint(command)
        activeDataSet(dsnameValue)
        closeDialog()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp()
	tkgrid(labelRcmdr(dsnameFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Name for new data set:  ")), entryDsname)
    tkgrid(dsnameFrame, sticky="w", columnspan=2)
    tkgrid(getFrame(dataSet1Box), sticky="nw")
    tkgrid(buttonsFrame, sticky="w", columnspan=2)
    dialogSuffix()
}


StatMedRenameDataset <- function(){
    dataSets <- listDataSets()
    .activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR", "Rename data set"))
    dsname <- tclVar("NewName")
    dsnameFrame <- tkframe(top)
    entryDsname <- ttkentry(dsnameFrame, width="20", textvariable=dsname)
    dataSet1Box <- variableListBox(top, dataSets, title=gettext(domain="R-RcmdrPlugin.EZR","Original Data Set"),
                                   initialSelection=if (is.null(.activeDataSet)) NULL else which(.activeDataSet == dataSets) - 1)
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Rename data set"), "#####", sep=""))
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == getSelection(dataSet1Box)) {
            errorCondition(recall=StatMedRenameDataset,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must enter a different data set name."))
            return()
        }
        if (dsnameValue == "") {
            errorCondition(recall=StatMedRenameDataset,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall=StatMedRenameDataset,
                           message=paste('"', dsnameValue, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
                closeDialog()
                StatMedRenameDataset()
                return()
            }
        }
        name1 <- getSelection(dataSet1Box)
        if (length(name1) == 0){
            errorCondition(recall=StatMedRenameDataset,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must select a data set."))
            return()
        }
        command <- paste(dsnameValue, " <- ", name1, sep="")
        doItAndPrint(command)
        activeDataSet(dsnameValue)
        command <- paste("remove(", name1, ")", sep="")
        doItAndPrint(command)
        activeDataSet(dsnameValue)
        closeDialog()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp()
    tkgrid(labelRcmdr(dsnameFrame, text=gettext(domain="R-RcmdrPlugin.EZR","New name for the data set:  ")), entryDsname)
    tkgrid(dsnameFrame, sticky="w", columnspan=2)
    tkgrid(getFrame(dataSet1Box), sticky="nw")
    tkgrid(buttonsFrame, sticky="w", columnspan=2)
    dialogSuffix()
}


StatMedMergeDatasets <- function(){
    dataSets <- listDataSets()
    .activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR", "Merge data sets"))
    dsname <- tclVar("MergedDataset")
    dsnameFrame <- tkframe(top)
    entryDsname <- ttkentry(dsnameFrame, width="20", textvariable=dsname)
    dataSet1Box <- variableListBox(top, dataSets, title=gettext(domain="R-RcmdrPlugin.EZR","First Data Set (pick one)"),
                                   initialSelection=if (is.null(.activeDataSet)) NULL else which(.activeDataSet == dataSets) - 1)
    dataSet2Box <- variableListBox(top, dataSets, title=gettext(domain="R-RcmdrPlugin.EZR","Second Data Set (pick one)"))
    commonVar <- tclVar("0")
    commonFrame <- tkframe(top)
    commonButton <- ttkcheckbutton(commonFrame, variable=commonVar)    
    radioButtons(top, "direction", buttons=c("rows", "columns"), 
                 labels=gettext(domain="R-RcmdrPlugin.EZR",c("Merge rows", "Merge columns")), title=gettext(domain="R-RcmdrPlugin.EZR","Direction of Merge"))
    radioButtons(top, "columnmerge", buttons=c("rownumber", "columns"), 
                 labels=gettext(domain="R-RcmdrPlugin.EZR",c("Merge by row number", "Merge by specified columns")), title=gettext(domain="R-RcmdrPlugin.EZR","Matching method to merge columns"))
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Merge data sets"), "#####", sep=""))
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=StatMedMergeDatasets,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall=StatMedMergeDatasets,
                           message=paste('"', dsnameValue, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
                closeDialog()
                StatMedMergeDatasets()
                return()
            }
        }
        name1 <- getSelection(dataSet1Box)
        name2 <- getSelection(dataSet2Box)
        if (length(name1) == 0){
            errorCondition(recall=StatMedMergeDatasets,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must select a data set."))
            return()
        }
        if (length(name2) == 0){
            errorCondition(recall=StatMedMergeDatasets,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You must select a data set."))
            return()
        }
        if (name1 == name2){
            errorCondition(recall=StatMedMergeDatasets,
                           message=gettext(domain="R-RcmdrPlugin.EZR","You cannot merge a data set with itself."))
            return()
        }
        common <- if (tclvalue(commonVar) == "1") TRUE else FALSE
        direction <- tclvalue(directionVariable)
        columnmerge <- tclvalue(columnmergeVariable)
        if (direction == "rows"){
            command <- paste(dsnameValue, " <- mergeRows(", name1, ", ", name2,
                             ", common.only=", common, ")", sep="")
            doItAndPrint(command)	
			activeDataSet(dsnameValue)
        }
        else {		
			if (columnmerge == "columns"){
				command <- paste(dsnameValue, " <- NULL", sep="")
				doItAndPrint(command)
				initializeDialog(subdialog, title=gettext(domain="R-RcmdrPlugin.EZR","Columns to merge datasets"))
				onOKsub <- function() {
					column.name.1 <- getSelection(column1Box)
					column.name.2 <- getSelection(column2Box)
					if (length(column.name.1) == 0){
						errorCondition(recall=NULL,
                        message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables"))
						return()
					}
					if (length(column.name.2) == 0){
						errorCondition(recall=NULL,
                        message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables"))
						return()
					}
					closeDialog(subdialog)
					command <- paste(dsnameValue, " <- merge(", name1, ", ", name2, ", all=", !common, ', by.x="', column.name.1, '", by.y="', column.name.2, '")', sep="")
					doItAndPrint(command)
					activeDataSet(dsnameValue)
				}
				subOKCancelHelp()
				list1 <- listVariables(name1)
				list2 <- listVariables(name2)				
				column1Box <- variableListBox(subdialog, list1, title=gettext(domain="R-RcmdrPlugin.EZR","Column name for matching in dataset 1(pick one)"), listHeight=10)
				column2Box <- variableListBox(subdialog, list2, title=gettext(domain="R-RcmdrPlugin.EZR","Column name for matching in dataset 2(pick one)"), listHeight=10)
				tkgrid(getFrame(column1Box), getFrame(column2Box), sticky="nw")				
				tkgrid(subButtonsFrame, sticky="w", columnspan=2)
				dialogSuffix(subdialog, focus=subdialog, force.wait=TRUE)
			} else {
				command <- paste(dsnameValue, " <- merge(", name1, ", ", name2, ", all=", !common, ', by="row.names")', sep="")
				doItAndPrint(command)
				command <- paste("rownames(", dsnameValue, ") <- ", dsnameValue, "$Row.names", sep="")
				doItAndPrint(command)
				command <- paste(dsnameValue, "$Row.names <- NULL", sep="")
				doItAndPrint(command)
				activeDataSet(dsnameValue)
			}
        }
        closeDialog()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="mergeRows")
    tkgrid(labelRcmdr(dsnameFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Name for merged data set:  ")), entryDsname)
    tkgrid(dsnameFrame, sticky="w", columnspan=2)
    tkgrid(getFrame(dataSet1Box), getFrame(dataSet2Box), sticky="nw")
    tkgrid(commonButton, labelRcmdr(commonFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Merge only common\nrows or columns")), 
           sticky="nw")
    tkgrid(directionFrame, commonFrame, sticky="sw")
    tkgrid(columnmergeFrame, sticky="sw")
    tkgrid(buttonsFrame, sticky="w", columnspan=2)
    dialogSuffix()
}


StatMedSaveDataSet <- function() {
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Save active data set"), "#####", sep=""))
	if (activeDataSetP() == FALSE){
            logger(gettext(domain="R-RcmdrPlugin.EZR","There is no active data set."))
            return()					
	}
	file <- tclvalue(tkgetSaveFile(filetypes=
				gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"R Data Files" {".rda" ".Rda" ".RDA" ".RData"}}'),
			defaultextension="rda", initialfile=paste(activeDataSet(), "rda", sep=".")))
	if (file == "") return()
	command <- paste('save("', activeDataSet(), '", file="', file, '")', sep="")
	justDoIt(command)
	logger(command)
}


StatMedExportDataSet <- function() {
	if (activeDataSetP() == FALSE){
            logger(gettext(domain="R-RcmdrPlugin.EZR","There is no active data set."))
            return()					
	}
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Export active data set (Text)"), "#####", sep=""))
	dsname <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Export Active Data Set"))
	checkBoxes(frame="optionsFrame", boxes=c("colnames", "rownames", "quotes"),
		initialValues=c(1,0,1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Write variable names:", "Write row names:", "Quotes around character values:")))
	missingVariable <- tclVar("NA")
	missingEntry <- ttkentry(optionsFrame, width="8", textvariable=missingVariable)
	radioButtons(name="delimiter", buttons=c("spaces", "tabs", "commas"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Spaces", "Tabs", "Commas")),
		initialValue="commas", title=gettext(domain="R-RcmdrPlugin.EZR","Field Separator"))
	otherButton <- ttkradiobutton(delimiterFrame, variable=delimiterVariable, value="other")
	otherVariable <- tclVar("")
	otherEntry <- ttkentry(delimiterFrame, width="4", textvariable=otherVariable)
	onOK <- function(){
		closeDialog()
		col <- tclvalue(colnamesVariable) == 1
		row <- tclvalue(rownamesVariable) == 1
		quote <- tclvalue(quotesVariable) == 1
		delim <- tclvalue(delimiterVariable)
		missing <- tclvalue(missingVariable)
		sep <- if (delim == "tabs") "\\t"
			else if (delim == "spaces") " "
			else if (delim == "commas") ","
			else trim.blanks(tclvalue(otherVariable))
		saveFile <- tclvalue(tkgetSaveFile(filetypes=gettext(domain="R-RcmdrPlugin.EZR",'{"Text Files" {".txt" ".TXT" ".dat" ".DAT" ".csv" ".CSV"}} {"All Files" {"*"}}'),
				defaultextension="txt", initialfile=paste(dsname, ".txt", sep="")))
		if (saveFile == "") {
			tkfocus(CommanderWindow())
			return()
		}
		command <- paste("write.table(", dsname, ', "', saveFile, '", sep="', sep,
			'", col.names=', col, ", row.names=", row, ", quote=", quote,
			', na="', missing, '")', sep="")
		justDoIt(command)
		logger(command)
		Message(paste(gettext(domain="R-RcmdrPlugin.EZR","Active dataset exported to file"), saveFile), type="note")
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="write.table")
	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Missing values:")), missingEntry, sticky="w")
	tkgrid(optionsFrame, sticky="w")
	tkgrid(labelRcmdr(delimiterFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Other")), otherButton,
		labelRcmdr(delimiterFrame, text=gettext(domain="R-RcmdrPlugin.EZR","  Specify:")), otherEntry, sticky="w")
	tkgrid(delimiterFrame, stick="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=3, columns=1)
}


StatMedExportStata <- function() {
	Library("foreign")
	if (activeDataSetP() == FALSE){
            logger(gettext(domain="R-RcmdrPlugin.EZR","There is no active data set."))
            return()					
	}
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Export active data set (Stata)"), "#####", sep=""))
	dsname <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Export Active Data Set"))
	onOK <- function(){
		closeDialog()
		saveFile <- tclvalue(tkgetSaveFile(filetypes=gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"Stata datasets" {".dta" ".DTA"}}'),
				defaultextension="dta", initialfile=paste(dsname, ".dta", sep="")))
		if (saveFile == "") {
			tkfocus(CommanderWindow())
			return()
		}
		command <- paste("write.dta(", dsname, ', "', saveFile, '")', sep="")
		justDoIt(command)
		logger(command)
		Message(paste(gettext(domain="R-RcmdrPlugin.EZR","Active dataset exported to Stata file"), saveFile), type="note")
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="write.dta")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=3, columns=1)
}


StatMedVariableCheck <- function(){
  logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Show variables in active data set"), "#####", sep=""))
  command <- paste("str(", activeDataSet(), ")", sep="")
  doItAndPrint(command)
  invisible(NULL)
  }
  
  
StatMedSubsetDataSet <- function(){
	dataSet <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Subset Data Set"))
	allVariablesFrame <- tkframe(top)
	allVariables <- tclVar("1")
	allVariablesCheckBox <- tkcheckbutton(allVariablesFrame, variable=allVariables)
	variablesBox <- variableListBox(top, Variables(), selectmode="multiple",
		initialSelection=NULL, title=gettext(domain="R-RcmdrPlugin.EZR","Variables (select one or more)"), listHeight=10)
	subsetVariable <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<all cases>"))
	subsetFrame <- tkframe(top)
	subsetEntry <- ttkentry(subsetFrame, width="60", textvariable=subsetVariable)
	subsetScroll <- ttkscrollbar(subsetFrame, orient="horizontal",
		command=function(...) tkxview(subsetEntry, ...))
	tkconfigure(subsetEntry, xscrollcommand=function(...) tkset(subsetScroll, ...))
	newDataSetName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<same as active data set>"))
	justshowVariablesFrame <- tkframe(top)
	justshowVariables <- tclVar("0")
	justshowVariablesCheckBox <- tkcheckbutton(justshowVariablesFrame, variable=justshowVariables)
	dataSetNameFrame <- tkframe(top)
	dataSetNameEntry <- ttkentry(dataSetNameFrame, width="25", textvariable=newDataSetName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Create subset data set"), "#####", sep=""))
		justshow <- tclvalue(justshowVariables)
		newName <- trim.blanks(tclvalue(newDataSetName))
		if (newName == gettext(domain="R-RcmdrPlugin.EZR","<same as active data set>")) newName <- ActiveDataSet()
		if (!is.valid.name(newName)){
			errorCondition(recall=StatMedSubsetDataSet,
				message=paste('"', newName, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (justshow==0 & is.element(newName, listDataSets())) {
			if ("no" == tclvalue(checkReplace(newName, type=gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				closeDialog()
				StatMedSubsetDataSet()
				return()
			}
		}
		selectVars <- if (tclvalue(allVariables) == "1") ""
			else {
				x <- getSelection(variablesBox)
				if (0 > length(x)) {
					errorCondition(recall=StatMedSubsetDataSet,
						message=gettext(domain="R-RcmdrPlugin.EZR","No variables were selected."))
					return()
				}
				paste(", select=c(", paste(x, collapse=","), ")", sep="")
			}
		closeDialog()
		cases <- tclvalue(subsetVariable)
		selectCases <- if (cases == gettext(domain="R-RcmdrPlugin.EZR","<all cases>")) ""
			else paste(", subset=", cases, sep="")
		if (selectVars == "" && selectCases ==""){
			errorCondition(recall=StatMedSubsetDataSet,
				message=gettext(domain="R-RcmdrPlugin.EZR","New data set same as active data set."))
			return()
		}
		if (justshow==0){
			command <- paste(newName, " <- subset(", ActiveDataSet(), selectCases, selectVars, ")",
				sep="")
			logger(command)
			result <- justDoIt(command)
			if (class(result)[1] !=  "try-error") activeDataSet(newName)
		} else {
			command <- paste("subset(", ActiveDataSet(), selectCases, selectVars, ")", sep="")
			doItAndPrint(command)
		}
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="subset")
	tkgrid(labelRcmdr(allVariablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Include all variables")),
		allVariablesCheckBox, sticky="w")
	tkgrid(allVariablesFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","   OR"), fg="red"), sticky="w")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(labelRcmdr(subsetFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Condition to extract samples")), sticky="w")
	tkgrid(labelRcmdr(subsetFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Example 1: age>50 & Sex==0, Example 2: age<50 | Sex==1")), sticky="w")  
	tkgrid(subsetEntry, sticky="w")
	tkgrid(subsetScroll, sticky="ew")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(labelRcmdr(justshowVariablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","View data only (not create data set)")),
		justshowVariablesCheckBox, sticky="w")
	tkgrid(justshowVariablesFrame, sticky="w")
	tkgrid(labelRcmdr(dataSetNameFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Name for new data set")), sticky="w")
	tkgrid(dataSetNameEntry, sticky="w")
	tkgrid(dataSetNameFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}


StatMedRenameVariables <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Rename variables"))
	variableBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"),
		selectmode="multiple", initialSelection=NULL, listHeight=10)
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Rename variables"), "#####", sep=""))
		variables <- getSelection(variableBox)
		closeDialog()
		nvariables <- length(variables)
		if (nvariables < 1) {
			errorCondition(recall=renameVariables, message=gettext(domain="R-RcmdrPlugin.EZR","No variables selected."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		unordered.names <- names(get(.activeDataSet))
#        unordered.names <- names(eval(parse(text=.activeDataSet)))
		which.variables <- match(variables, unordered.names)
		initializeDialog(subdialog, title=gettext(domain="R-RcmdrPlugin.EZR","Variable Names"))
		newnames <- rep("", nvariables)
		onOKsub <- function() {
			closeDialog(subdialog)
			for (i in 1:nvariables){
				newnames[i] <- eval(parse(text=paste("tclvalue(newName", i, ")", sep="")))
			}
			if (any(newnames == "")){
				errorCondition(recall=renameVariables, message=gettext(domain="R-RcmdrPlugin.EZR","A variable name is empty."))
				return()
			}
			test.names <- newnames == make.names(newnames)
			if (!all(test.names)){
				errorCondition(recall=renameVariables,
					message=paste(gettext(domain="R-RcmdrPlugin.EZR","The following variable names are not valid:\n"),
						paste(newnames[!test.names], collapse=", ")))
				return()
			}
			all.names <- names(get(.activeDataSet))
#            all.names <- eval(parse(text=paste("names(", .activeDataSet, ")")))
			all.names[which.variables] <- newnames
			if (length(unique(all.names)) != length(all.names)){
				errorCondition(recall=renameVariables, message=gettext(domain="R-RcmdrPlugin.EZR","Variable names are not unique"))
				return()
			}
			command <- paste("names(", .activeDataSet, ")[c(", paste(which.variables, collapse=","),
				")] <- c(", paste('"', newnames, '"', collapse=",", sep=""), ")", sep="")
			result <- justDoIt(command)
			logger(command)
			if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
			tkfocus(CommanderWindow())
		}
		subOKCancelHelp()
		tkgrid(labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Old Name"), fg="blue"),
			labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","New name"), fg="blue"), sticky="w")
		for (i in 1:nvariables){
			valVar <- paste("newName", i, sep="")
			assign(valVar, tclVar(""))
			assign(paste("entry", i, sep=""), ttkentry(subdialog, width="20",
#                textvariable=eval(parse(text=valVar))))
					textvariable=get(valVar)))
			tkgrid(labelRcmdr(subdialog, text=variables[i]), get(paste("entry", i, sep="")), sticky="w")
#            tkgrid(labelRcmdr(subdialog, text=variables[i]), eval(parse(text=paste("entry", i, sep=""))), sticky="w")
		}
		tkgrid(subButtonsFrame, sticky="w", columnspan=2)
		dialogSuffix(subdialog, rows=nvariables+2, columns=2, focus=entry1, onOK=onOKsub, force.wait=TRUE)
	}
	OKCancelHelp(helpSubject="names")
	tkgrid(getFrame(variableBox), sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=2, columns=1)
}


StatMedDeleteVariable <- function(){
	dataSet <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Delete variables from data set"))
	variablesBox <- variableListBox(top, Variables(),
		title=gettext(domain="R-RcmdrPlugin.EZR","Variable(s) to delete (pick one or more)"), selectmode="multiple",
		initialSelection=NULL, listHeight=15)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Delete variables from data set"), "#####", sep=""))
		variables <- getSelection(variablesBox)
		closeDialog()
		if (length(variables) == 0) {
			errorCondition(recall=StatMedDeleteVariable, message=gettext(domain="R-RcmdrPlugin.EZR","You must select one or more variables."))
			return()
		}
		if (length(variables) == 1){
			response <- tclvalue(RcmdrTkmessageBox(message=sprintf(gettext(domain="R-RcmdrPlugin.EZR","Delete %s?\nPlease confirm."), variables), icon="warning", type="okcancel", default="cancel"))
			if (response == "cancel") {
				onCancel()
				return()
			}
		}
		else{
			response <- tclvalue(RcmdrTkmessageBox(message=
						sprintf(gettext(domain="R-RcmdrPlugin.EZR","Delete %d variables?\nPlease confirm."), length(variables)),
					icon="warning", type="okcancel", default="cancel"))
			if (response == "cancel") {
				onCancel()
				return()
			}
		}
		for (variable in variables){
			eval(parse(text=paste(dataSet, "$", variable, "<- NULL", sep="")), envir=.GlobalEnv)
			logger(paste(dataSet, "$", variable, " <- NULL", sep=""))
		}
		activeDataSet(dataSet, flushModel=FALSE)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="NULL")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=2, columns=1)
}


StatMedStack <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Stack variables to long format data set"))
	variableBox <- variableListBox(top, Variables(), selectmode="multiple",
		title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick two or more)"), listHeight=10)
	factorName <- tclVar("")
	factorNameField <- ttkentry(top, width="20", textvariable=factorName)
	variableName <- tclVar("")
	variableNameField <- ttkentry(top, width="20", textvariable=variableName)
	datasetName <- tclVar("")
	datasetNameField <- ttkentry(top, width="20", textvariable=datasetName)
	checkBoxes(frame="checkboxFrame", boxes=c("othervar"), initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Include other variables in new data set")))	
	#    subsetBox(model=TRUE)
    StatMedSubsetBox()
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Stack variables to long format data set"), "#####", sep=""))
		variables <- getSelection(variableBox)
		facname <- tclvalue(factorName)
		varname <- tclvalue(variableName)
		dsname <- tclvalue(datasetName)
		othervar <- tclvalue(othervarVariable)
		closeDialog()
		if (length(variables) < 2) {
			errorCondition(recall=StatMedStack,
				message=gettext(domain="R-RcmdrPlugin.EZR","You must select at least two variables."))
			return()
		}
		if (!is.valid.name(facname)){
			errorCondition(recall=StatMedStack,
				message=paste('"', facname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (!is.valid.name(varname)){
			errorCondition(recall=StatMedStack,
				message=paste('"', varname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (!is.valid.name(dsname)){
			errorCondition(recall=Stack,
				message=paste('"', dsname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(dsname, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsname, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				Stack()
				return()
			}
		}
#		command <- paste(dsname, " <- stack(", activeDataSet(), "[, c(",
#			paste(paste('"', variables, '"', sep=""), collapse=","), ")])", sep="")
#		logger(command)
#		result <- justDoIt(command)
#		command <- paste("names(", dsname, ') <- c("', varname, '", "', facname, '")',
#			sep="")
#		logger(command)
#		justDoIt(command)

		dataSet <- ActiveDataSet()
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			doItAndPrint(paste("TempDF <- ", dataSet))
		}
        else {
			doItAndPrint(paste("TempDF <- subset(", dataSet, ",", subset, ")") )
        }

		nvar <- length(variables)
#		RepeatedData <- variables[1]
#		RepeatedData2 <- paste('"', variables[1], '"', sep="")
#		for(i in 2:nvar){
#			RepeatedData <- paste(RepeatedData, ", ", variables[i], sep="")
#			RepeatedData2 <- paste(RepeatedData2, ', "', variables[i], '"', sep="")
#			}
		logger("#Convert to long format")
		doItAndPrint("n <- length(TempDF[,1])")
#		doItAndPrint("TempDF$TempIDforReshaping <- c(1:n)")

		if(othervar==0){
			command <- paste("TempDF <- data.frame(", variables[1], "=TempDF$", variables[1], sep="")
			for (i in 2:nvar){
				command <- paste(command, ", ", variables[i], "=TempDF$", variables[i], sep="")
			}
			command <- paste(command, ")", sep="")
			doItAndPrint(command)
		}
		
#		command <- paste('TempDF2 <- reshape(TempDF, idvar="TempIDforReshaping", varying=list(c("', variables[1], sep="")
		command <- paste('TempDF2 <- reshape(TempDF, varying=list(c("', variables[1], sep="")
		for (i in 2:nvar){
			command <- paste(command, '", "', variables[i], sep="")
		}
		command <- paste(command, '")), v.names="', varname, '", timevar="', facname, '", direction="long")', sep="")
		doItAndPrint(command)
		command <- paste('RepeatNumber <- c("', variables[1], sep="")
		for (i in 2:nvar){
			command <- paste(command, '", "', variables[i], sep="")
		}
		command <- paste(command, '")', sep="")
		doItAndPrint(command)
		doItAndPrint(paste("TempDF2$", facname, " <- RepeatNumber[TempDF2$", facname, "]", sep=""))
#		doItAndPrint("TempDF2$TempIDforReshaping <- NULL")
		result <- doItAndPrint(paste(dsname, " <- TempDF2", sep=""))
		if (class(result)[1] !=  "try-error") activeDataSet(dsname)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="stack")
	tkgrid(getFrame(variableBox), sticky="nw", columnspan=2)
	tkgrid(labelRcmdr(top, text=""))
	tkgrid(labelRcmdr(top,
			text=gettext(domain="R-RcmdrPlugin.EZR","Name for stacked data set:")), datasetNameField, sticky="w")
	tkgrid(labelRcmdr(top,
			text=gettext(domain="R-RcmdrPlugin.EZR","Name for stacked variable data in new data set:")), variableNameField, sticky="w")
	tkgrid(labelRcmdr(top,
			text=gettext(domain="R-RcmdrPlugin.EZR","Name for factor to identify stacked variables in new data set:")), factorNameField, sticky="w")
	 tkgrid(checkboxFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=5, columns=2, preventGrabFocus=TRUE)
}


StatMedSort <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Sort rows"))
	variablesBox <- variableListBox(top, Variables(), initialSelection=NULL, title=gettext(domain="R-RcmdrPlugin.EZR","Variable for sorting"), listHeight=10)
    optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="decreasing", buttons=gettext(domain="R-RcmdrPlugin.EZR",c("Ascending", "Descending")), values=c("FALSE", "TRUE"),
    labels=gettext(domain="R-RcmdrPlugin.EZR",c("Ascending", "Descending")), title=gettext(domain="R-RcmdrPlugin.EZR","Sorting order"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Sort rows"), "#####", sep=""))
		dataSet <- activeDataSet()
		x <- getSelection(variablesBox)
		if (length(x) == 0) {
            errorCondition(recall=StatMedSort, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            return()			
		}
		closeDialog()
        decreasing <- tclvalue(decreasingVariable)
		command <- paste("TempList <- order(", dataSet, "$", x, ", decreasing=", decreasing, ")", sep="")
		doItAndPrint(command)
		command <- paste(dataSet, " <- ", dataSet, "[TempList,]", sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="order")
	tkgrid(getFrame(variablesBox), sticky="nw")
    tkgrid(decreasingFrame, labelRcmdr(optionsFrame, text="    "), sticky="nw")
    tkgrid(optionsFrame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}


StatMedCountMissing <- function(){
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Count missing observations of specified variables"))
	variableBox <- variableListBox(top, Variables(), selectmode="multiple",
		title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"), listHeight=15)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Count missing observations of specified variables"), "#####", sep=""))
		variables <- getSelection(variableBox)
		closeDialog()
		if (length(variables) == 0) {
			errorCondition(recall=StatMedCountMissing, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		for (name in variables){
			command <- paste("sum(is.na(", .activeDataSet, "$", name, ")) ###", name, gettext(domain="R-RcmdrPlugin.EZR",": Number of missing observations"), sep="")
			doItAndPrint(command)
		}
	}
	OKCancelHelp(helpSubject="is.na")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
	tkgrid(getFrame(variableBox), sticky="nw")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=4, columns=2, preventGrabFocus=TRUE)
}


StatMedFilterNA <- function(){
	dataSet <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Remove rows with missing data in specified variables"))
	variablesBox <- variableListBox(top, Variables(), selectmode="multiple", initialSelection=NULL,
		title=gettext(domain="R-RcmdrPlugin.EZR","Variables to remove rows with missing data (pick one or more)"), listHeight=15)
	newDataSetName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<same as active data set>"))
	dataSetNameFrame <- tkframe(top)
	dataSetNameEntry <- ttkentry(dataSetNameFrame, width="25", textvariable=newDataSetName)
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Remove rows with missing data in specified variables"), "#####", sep=""))
		x <- getSelection(variablesBox)
		closeDialog()
		newName <- trim.blanks(tclvalue(newDataSetName))
		.activeDataSet <- ActiveDataSet()
		if (newName == gettext(domain="R-RcmdrPlugin.EZR","<same as active data set>")) newName <- .activeDataSet
		if (!is.valid.name(newName)){
			errorCondition(recall=StatMedFilterNA,
				message=paste('"', newName, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(newName, listDataSets())) {
			if ("no" == tclvalue(checkReplace(newName, gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				filterNA()
				return()
			}
		}	
		if (length(x) == 0) {
				errorCondition(recall=StatMedFilterNA, message=gettext(domain="R-RcmdrPlugin.EZR","No variables were selected."))
				return()
		}
		command <- paste(newName, " <- ", .activeDataSet, "[complete.cases(", .activeDataSet, "$", x[1], sep="")
		if (length(x)>1){
			for (i in 2:length(x)){
				command <- paste(command, ", ", .activeDataSet, "$", x[i], sep="")
			}
		}
		command <- paste(command, "),]", sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(newName)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="complete.cases")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(labelRcmdr(dataSetNameFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Name for new data set")), sticky="w")
	tkgrid(dataSetNameEntry, sticky="w")
	tkgrid(dataSetNameFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

  
StatMedNAgroup <- function(){		
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Convert missing observations to a group"))
  	dataSet <- activeDataSet()
	variablesBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable including missing data"), listHeight=15)
	newVariableName <- tclVar("")
	newVariableNameEntry <- ttkentry(top, width="20", textvariable=newVariableName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Convert missing observations to a group"), "#####", sep=""))
		var <- trim.blanks(getSelection(variablesBox))
        	if (length(var) == 0){
  	          errorCondition(recall=StatMedNAgroup, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            	return()
        }		
		newVar <- trim.blanks(tclvalue(newVariableName))
		if (!is.valid.name(newVar)){
			errorCondition(recall=StatMedNAgroup,
				message=paste('"', newVar, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (newVar == var){
			errorCondition(recall=StatMedNAgroup, message=gettext(domain="R-RcmdrPlugin.EZR","New variable name must be different from the original name."))
			return()
		}
		closeDialog()
		command <-  paste("if(sum(is.na(", dataSet, "$", var, "))>0) ", dataSet,"$",newVar, " <- as.factor(ifelse(is.na(", dataSet, "$", var, '), "NA", as.character(', dataSet, "$", var, ")))", sep="")
		result <- doItAndPrint(command)
		command <-  paste("if(sum(is.na(", dataSet, "$", var, '))==0) cat(gettext(domain="R-RcmdrPlugin.EZR","There was no missing data."), "\n")', sep="")
		result <- doItAndPrint(command)
		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
		doItAndPrint(paste("if(sum(is.na(", dataSet, "$", var, '))>0) cat(gettext(domain="R-RcmdrPlugin.EZR","New variable"), "', newVar, '",  gettext(domain="R-RcmdrPlugin.EZR","was made."), "\n")', sep="") )
		doItAndPrint(paste("if(sum(is.na(", dataSet, "$", var, "))>0) table(", dataSet, "$", newVar, ")", sep="") )
		tkfocus(CommanderWindow())
	}
    OKCancelHelp(helpSubject="is.na")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","New variable name")), newVariableNameEntry, sticky="w")
	tkgrid.configure(newVariableNameEntry, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=3, columns=2)
	}

	
StatMedNewVar <- function(){		
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Bin numeric variable with specified threshold"))
  	dataSet <- activeDataSet()
	variablesBox <- variableListBox(top, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Select one numeric variable"), listHeight=15)
	newVariableName <- tclVar("")
	newVariableNameEntry <- ttkentry(top, width="20", textvariable=newVariableName)
	threshold <- tclVar("")
	thresholdEntry <- ttkentry(top, width="20", textvariable=threshold)
	radioButtons(name="grouping", buttons=c("equalgreater", "greater"), values=c(">=", ">"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c(">= (equal to or greater than)", "> (greater than)")), title=gettext(domain="R-RcmdrPlugin.EZR","Threshold"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Bin numeric variable with specified threshold"), "#####", sep=""))
		var <- trim.blanks(getSelection(variablesBox))
        	if (length(var) == 0){
  	          errorCondition(recall=StatMedNewVar, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            	return()
            	}
		newVar <- trim.blanks(tclvalue(newVariableName))
		if (!is.valid.name(newVar)){
			errorCondition(recall=StatMedNewVar,
				message=paste('"', newVar, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		threshold <- tclvalue(threshold)
        	if (length(threshold) == 0){
  	          errorCondition(recall=StatMedNewVar, message=gettext(domain="R-RcmdrPlugin.EZR","Input threshold to bin a numeric variable."))
            	return()
            	}
		grouping <- as.character(tclvalue(groupingVariable))
		closeDialog()
		command <-  paste(dataSet,"$",newVar, " <- ifelse(", dataSet, "$", var, grouping, threshold, ", 1 , 0)", sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
#		logger(paste("#", gettext(domain="R-RcmdrPlugin.EZR","New variable"), " ", newVar, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), "(", threshold, gettext(domain="R-RcmdrPlugin.EZR","<=:1, >:0"), sep="") )
		logger(paste("#", gettext(domain="R-RcmdrPlugin.EZR","New variable"), " ", newVar, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), sep="") )
		doItAndPrint(paste("table(", dataSet, "$", newVar, ", exclude=NULL)", sep="") )
		tkfocus(CommanderWindow())
	}
        OKCancelHelp(helpSubject="ifelse")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","New variable name")), newVariableNameEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Threshold to bin a numeric variable.")), thresholdEntry, sticky="w")
	tkgrid.configure(newVariableNameEntry, sticky="w")
	tkgrid.configure(thresholdEntry, sticky="w")
	tkgrid(groupingFrame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=3, columns=2)
	}

	
StatMedNewVar2 <- function(){		
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Bin numeric variable to more than 2 groups with specified thresholds"))
  	dataSet <- activeDataSet()
	variablesBox <- variableListBox(top, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Select one numeric variable"), listHeight=15)

	thresholdFrame <- tkframe(top)
	newVariableNameFrame <- tkframe(thresholdFrame)
	newVariableNameVariable <- tclVar("")	
	newVariableNameField <- ttkentry(thresholdFrame, width="20", textvariable=newVariableNameVariable)
	threshold1Frame <- tkframe(thresholdFrame)
	threshold1Variable <- tclVar("")	
	threshold1Field <- ttkentry(thresholdFrame, width="10", textvariable=threshold1Variable)
	threshold2Frame <- tkframe(thresholdFrame)
	threshold2Variable <- tclVar("")	
	threshold2Field <- ttkentry(thresholdFrame, width="10", textvariable=threshold2Variable)
	threshold3Frame <- tkframe(thresholdFrame)
	threshold3Variable <- tclVar("")	
	threshold3Field <- ttkentry(thresholdFrame, width="10", textvariable=threshold3Variable)
	threshold4Frame <- tkframe(thresholdFrame)
	threshold4Variable <- tclVar("")	
	threshold4Field <- ttkentry(thresholdFrame, width="10", textvariable=threshold4Variable)
	threshold5Frame <- tkframe(thresholdFrame)
	threshold5Variable <- tclVar("")	
	threshold5Field <- ttkentry(thresholdFrame, width="10", textvariable=threshold5Variable)
	levelname1Frame <- tkframe(thresholdFrame)
	levelname1Variable <- tclVar("<no group>")	
	levelname1Field <- ttkentry(thresholdFrame, width="20", textvariable=levelname1Variable)
	levelname2Frame <- tkframe(thresholdFrame)
	levelname2Variable <- tclVar("<no group>")	
	levelname2Field <- ttkentry(thresholdFrame, width="20", textvariable=levelname2Variable)
	levelname3Frame <- tkframe(thresholdFrame)
	levelname3Variable <- tclVar("<no group>")	
	levelname3Field <- ttkentry(thresholdFrame, width="20", textvariable=levelname3Variable)
	levelname4Frame <- tkframe(thresholdFrame)
	levelname4Variable <- tclVar("<no group>")	
	levelname4Field <- ttkentry(thresholdFrame, width="20", textvariable=levelname4Variable)
	levelname5Frame <- tkframe(thresholdFrame)
	levelname5Variable <- tclVar("<no group>")	
	levelname5Field <- ttkentry(thresholdFrame, width="20", textvariable=levelname5Variable)
	levelname6Frame <- tkframe(thresholdFrame)
	levelname6Variable <- tclVar("<no group>")	
	levelname6Field <- ttkentry(thresholdFrame, width="20", textvariable=levelname6Variable)
	radioButtons(name="grouping", buttons=c("equalgreater", "greater"), values=c(">=", ">"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c(">= (equal to or greater than)", "> (greater than)")), title=gettext(domain="R-RcmdrPlugin.EZR","Threshold"))

	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Bin numeric variable to more than 2 groups with specified thresholds"), "#####", sep=""))
		var <- trim.blanks(getSelection(variablesBox))
        	if (length(var) == 0){
  	          errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            	return()
            	}
		newVar <- trim.blanks(tclvalue(newVariableNameVariable))
		if (!is.valid.name(newVar)){
			errorCondition(recall=StatMedNewVar2,
				message=paste('"', newVar, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		threshold1 <- tclvalue(threshold1Variable)
		threshold2 <- tclvalue(threshold2Variable)
		threshold3 <- tclvalue(threshold3Variable)
		threshold4 <- tclvalue(threshold4Variable)
		threshold5 <- tclvalue(threshold5Variable)
		levelname1 <- tclvalue(levelname1Variable)
		levelname2 <- tclvalue(levelname2Variable)
		levelname3 <- tclvalue(levelname3Variable)
		levelname4 <- tclvalue(levelname4Variable)
		levelname5 <- tclvalue(levelname5Variable)
		levelname6 <- tclvalue(levelname6Variable)
		grouping <- as.character(tclvalue(groupingVariable))
		if(grouping==">="){
			right <- ", right=FALSE)"
		} else {
			right <- ", right=TRUE)"
		}
		levels <- 0
		breaks <- ", breaks=c(-Inf, "
		labels <- ", labels=c("
		if (levelname1 == "<no group>"){
  	        errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","Input at least two groups."))
            return()
        } else {
			levels <- levels + 1
			labels <- paste(labels, '"', levelname1, '"', sep="")
		}
		if (levelname2 != "<no group>"){
			if (length(threshold1) == 0){
  	          errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","Input threshold to bin a numeric variable."))
            return()
            }
			levels <- levels + 1
			breaks <- paste(breaks, threshold1, sep="")
			labels <- paste(labels, ', "', levelname2, '"', sep="")
		}
		if (levelname3 != "<no group>"){
			if (length(threshold2) == 0){
  	          errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","Input threshold to bin a numeric variable."))
            return()
            }
			levels <- levels + 1
			breaks <- paste(breaks, ", ", threshold2, sep="")
			labels <- paste(labels, ', "', levelname3, '"', sep="")
		}
		if (levelname4 != "<no group>"){
			if (length(threshold3) == 0){
  	          errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","Input threshold to bin a numeric variable."))
            return()
            }
			levels <- levels + 1
			breaks <- paste(breaks, ", ", threshold3, sep="")
			labels <- paste(labels, ', "', levelname4, '"', sep="")
		}
		if (levelname5 != "<no group>"){
			if (length(threshold4) == 0){
  	          errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","Input threshold to bin a numeric variable."))
            return()
            }
			levels <- levels + 1
			breaks <- paste(breaks, ", ", threshold4, sep="")
			labels <- paste(labels, ', "', levelname5, '"', sep="")
		}
		if (levelname6 != "<no group>"){
			if (length(threshold5) == 0){
  	          errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","Input threshold to bin a numeric variable."))
            return()
            }
			levels <- levels + 1
			breaks <- paste(breaks, ", ", threshold5, sep="")
			labels <- paste(labels, ', "', levelname6, '"', sep="")
		}
		if (levels < 2){
  	          errorCondition(recall=StatMedNewVar2, message=gettext(domain="R-RcmdrPlugin.EZR","Input at least two groups."))
            return()
        }
		breaks <- paste(breaks, ", Inf)", sep="")
		labels <- paste(labels, ")", sep="")
		closeDialog()
		command <-  paste(dataSet,"$",newVar, " <- cut(", dataSet, "$", var, breaks, labels, right, sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
		logger(paste("#", gettext(domain="R-RcmdrPlugin.EZR","New variable"), " ", newVar, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), sep="") )
		doItAndPrint(paste("table(", dataSet, "$", newVar, ", exclude=NULL)", sep="") )
		tkfocus(CommanderWindow())
	}
        OKCancelHelp(helpSubject="ifelse")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(labelRcmdr(newVariableNameFrame, text=gettext(domain="R-RcmdrPlugin.EZR","New variable name:")), newVariableNameField, sticky = "w")
	tkgrid(newVariableNameFrame, labelRcmdr(thresholdFrame, text="  "), sticky="w")

	tkgrid(labelRcmdr(thresholdFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Input thresholds and level names."), fg="blue"), sticky="w")

	tkgrid(labelRcmdr(levelname1Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Level group name"), " 1:", sep="")), levelname1Field, sticky = "w")
	tkgrid(levelname1Frame, labelRcmdr(thresholdFrame, text="  "), sticky="w")
	tkgrid(labelRcmdr(levelname2Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Level group name"), " 2:", sep="")), levelname2Field, sticky = "w")
	tkgrid(labelRcmdr(threshold1Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Threshold"), " 1:", sep="")), threshold1Field, sticky = "w")
	tkgrid(levelname2Frame, labelRcmdr(thresholdFrame, text="  "), threshold1Frame, sticky="w")
	tkgrid(labelRcmdr(levelname3Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Level group name"), " 3:", sep="")), levelname3Field, sticky = "w")
	tkgrid(labelRcmdr(threshold2Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Threshold"), " 2:", sep="")), threshold2Field, sticky = "w")
	tkgrid(levelname3Frame, labelRcmdr(thresholdFrame, text="  "), threshold2Frame, sticky="w")
	tkgrid(labelRcmdr(levelname4Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Level group name"), " 4:", sep="")), levelname4Field, sticky = "w")
	tkgrid(labelRcmdr(threshold3Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Threshold"), " 3:", sep="")), threshold3Field, sticky = "w")
	tkgrid(levelname4Frame, labelRcmdr(thresholdFrame, text="  "), threshold3Frame, sticky="w")
	tkgrid(labelRcmdr(levelname5Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Level group name"), " 5:", sep="")), levelname5Field, sticky = "w")
	tkgrid(labelRcmdr(threshold4Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Threshold"), " 4:", sep="")), threshold4Field, sticky = "w")
	tkgrid(levelname5Frame, labelRcmdr(thresholdFrame, text="  "), threshold4Frame, sticky="w")
	tkgrid(labelRcmdr(levelname6Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Level group name"), " 6:", sep="")), levelname6Field, sticky = "w")
	tkgrid(labelRcmdr(threshold5Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Threshold"), " 5:", sep="")), threshold5Field, sticky = "w")
	tkgrid(levelname6Frame, labelRcmdr(thresholdFrame, text="  "), threshold5Frame, sticky="w")
	tkgrid(thresholdFrame, sticky="w")
	tkgrid(groupingFrame, sticky="nw")

	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=3, columns=2)
	}

	
StatMedCompute <- function(){
	onDoubleClick <-function(){
		var <- trim.blanks(getSelection(variablesBox))
		word <- paste("\\[", gettext(domain="R-RcmdrPlugin.EZR","factor"), "\\]", sep="")
		if (length(grep(word, var)) == 1)
			var <- trim.blanks(sub(word, "",  var))
		tkfocus(compute)
		expr <- tclvalue(computeVar)
		tclvalue(computeVar) <- if (expr == "") var
			else paste(expr, var, sep=if (rev(strsplit(expr, "")[[1]])[1] =="(" ) "" else " ")
		tkicursor(compute, "end")
		tkxview.moveto(compute, "1")
	}
	dataSet <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Create new variable"))
	.variables <- Variables()
	variables <- paste(.variables, ifelse(is.element(.variables, Factors()), gettext(domain="R-RcmdrPlugin.EZR","[factor]"), ""))
	variablesBox <- variableListBox(top, variables, title=gettext(domain="R-RcmdrPlugin.EZR","Current variables (double-click to expression)"), listHeight=15)
	tkbind(variablesBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
	variablesFrame <- tkframe(top)
	newVariableName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","variable"))
	newVariable <- ttkentry(variablesFrame, width="20", textvariable=newVariableName)
	computeFrame <- tkframe(top)
	computeVar <- tclVar("")
	compute <- ttkentry(computeFrame, font=getRcmdr("logFont"), width="60", textvariable=computeVar)
	computeXscroll <- ttkscrollbar(computeFrame,
		orient="horizontal", command=function(...) tkxview(compute, ...))
	tkconfigure(compute, xscrollcommand=function(...) tkset(computeXscroll, ...))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Create new variable"), "#####", sep=""))
		closeDialog()
		newVar <- trim.blanks(tclvalue(newVariableName))
		if (!is.valid.name(newVar)){
			errorCondition(recall=StatMedCompute,
				message=paste('"', newVar, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		express <- tclvalue(computeVar)
		check.empty <- gsub(";", "", gsub(" ", "", express))
		if ("" == check.empty) {
			errorCondition(recall=StatMedCompute,
				message=gettext(domain="R-RcmdrPlugin.EZR","No expression specified."))
			return()
		}
		if (is.element(newVar, Variables())) {
			if ("no" == tclvalue(checkReplace(newVar, gettext(domain="R-RcmdrPlugin.EZR","Variable")))){
				StatMedCompute()
				return()
			}
		}
		command <-  paste(dataSet,"$",newVar, " <- with(", ActiveDataSet(),
			", ", express, ")", sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
		logger(paste("#", gettext(domain="R-RcmdrPlugin.EZR","New variable"), " ", newVar, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), sep="") )
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="ifelse")
	tkgrid(getFrame(variablesBox), sticky="nw", columnspan=2)
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","New variable name")), sticky="w")
	tkgrid(newVariable, labelRcmdr(variablesFrame, text="     "), sticky="w")
	tkgrid(labelRcmdr(computeFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Expression to compute")), sticky="w")
	tkgrid(compute, sticky="w")
	tkgrid(computeXscroll, sticky="ew")
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Example 1: ifelse(age > 50 & Sex == 0, 1, 0)")), sticky="w")  
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Example 2: ifelse(age < 50 | Sex == 1, 1, 0)")), sticky="w")  
	tkgrid(variablesFrame, sticky="nw")
	tkgrid(computeFrame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=3, columns=2, focus=compute)
}

	
StatMedLog <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Logarithmic transformation"))
	variableBox <- variableListBox(top, Numeric(), selectmode="multiple",
		title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"), listHeight=15)
	radioButtons(name="base", buttons=c("common", "natural", "binary"), values=c("10", "exp(1)", "2"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("Common logarithm (base=10)", "Natural logarithm (base=e)", "Binary logarithm (base=2)")), title=gettext(domain="R-RcmdrPlugin.EZR","Base of logarithmic transformation"))
	logName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<same as variables>"))
	logNameField <- ttkentry(top, width="20", textvariable=logName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Logarithmic transformation"), "#####", sep=""))
		variables <- getSelection(variableBox)
		closeDialog()
		if (length(variables) == 0) {
			errorCondition(recall=StatMedLog, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		logname <- trim.blanks(tclvalue(logName))
		base <- as.character(tclvalue(baseVariable))
		.activeDataSet <- ActiveDataSet()
		for (name in variables){
			lname <- if (logname == gettext(domain="R-RcmdrPlugin.EZR","<same as variables>")) name
				else if (length(variables) == 1) logname
				else paste(logname, name, sep="")
			if (!is.valid.name(lname)){
				errorCondition(recall=StatMedLog,
					message=paste('"', lname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
				return()
			}
			if (is.element(lname, Variables())) {
				if ("no" == tclvalue(checkReplace(lname))){
					StatMedLog()
					return()
				}
			}
			command <- paste(.activeDataSet, "$", lname, " <- log(", .activeDataSet, "$", name, ", base=", base, ")", sep="")
			result <- justDoIt(command)
			logger(command)
			if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
		logger(paste("#", gettext(domain="R-RcmdrPlugin.EZR","New variable"), " ", lname, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), sep="") )
		}
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="log")
	tkgrid(getFrame(variableBox), baseFrame, sticky="nw")
	tkgrid(labelRcmdr(top,
			text=gettext(domain="R-RcmdrPlugin.EZR","New variable name or prefix for multiple variables:")),
		logNameField, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=4, columns=2, preventGrabFocus=TRUE)
}

	
StatMedNumericToFactor <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Convert Numeric Variables to Factors"))
	variableBox <- variableListBox(top, Numeric(), selectmode="multiple",
		title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"), listHeight=15)
	radioButtons(name="levels", buttons=c("names", "numbers"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("Supply level names", "Use numbers")), title=gettext(domain="R-RcmdrPlugin.EZR","Factor Levels"))
	factorName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<same as variables>"))
	factorNameField <- ttkentry(top, width="20", textvariable=factorName)
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Convert numeric variables to factors"), "#####", sep=""))
		variables <- getSelection(variableBox)
		closeDialog()
		if (length(variables) == 0) {
			errorCondition(recall=StatMedNumericToFactor, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		facname <- trim.blanks(tclvalue(factorName))
		.activeDataSet <- ActiveDataSet()
		cmd <- paste("apply(", .activeDataSet, "[c(", paste(
				paste('"', variables, '"', sep=""),
				collapse=","), ")], 2, function(x) sort(unique(x)))", sep="")
		levs <- eval(parse(text=cmd), envir=.GlobalEnv)
		sameLevels <- (length(variables) == 1) ||
			((is.matrix(levs)) && (all(0 == apply(levs, 1, var))))
		for (name in variables){
			fname <- if (facname == gettext(domain="R-RcmdrPlugin.EZR","<same as variables>")) name
				else if (length(variables) == 1) facname
				else paste(facname, name, sep="")
			if (!is.valid.name(fname)){
				errorCondition(recall=StatMedNumericToFactor,
					message=paste('"', fname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
				return()
			}
			if (is.element(fname, Variables())) {
				if ("no" == tclvalue(checkReplace(fname))){
					StatMedNumericToFactor()
					return()
				}
			}
			levelsType <- tclvalue(levelsVariable)
			env <- environment()
			if (((name == variables[1]) || (!sameLevels)) && (levelsType == "names")){
				values <- sort(unique(eval(parse(text=paste(.activeDataSet, "$", name, sep="")),
							envir=.GlobalEnv)))
				nvalues <- length(values)
				if (nvalues > 30) {
					errorCondition(recall=StatMedNumericToFactor,
						message=sprintf(gettext(domain="R-RcmdrPlugin.EZR","Number of levels (%d) too large."), nvalues))
					return()
				}
				initializeDialog(subdialog,
					title=paste(gettext(domain="R-RcmdrPlugin.EZR","Level Names for"),
						if(sameLevels && length(variables) > 1) "Factors" else fname))
				names <- rep("", nvalues)
				onOKsub <- function() {
					closeDialog(subdialog)
					for (i in 1:nvalues){
						names[i] <- eval(parse(text=paste("tclvalue(levelName", i, ")", sep="")))
					}
					if (length(unique(names)) != nvalues){
						errorCondition(recall=StatMedNumericToFactor,
							message=gettext(domain="R-RcmdrPlugin.EZR","Levels names are not unique."))
						return()
					}
					if (any(names == "")){
						errorCondition(recall=StatMedNumericToFactor,
							message=gettext(domain="R-RcmdrPlugin.EZR","A level name is empty."))
						return()
					}
					assign("labels", paste(paste("'", names, "'", sep=""), collapse=","),
						envir=env)
				}
				subOKCancelHelp()
				tkgrid(labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Numeric value")), labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Level name")), sticky="w")
				for (i in 1:nvalues){
					valVar <- paste("levelName", i, sep="")
					assign(valVar, tclVar(""))
					assign(paste("entry", i, sep=""), ttkentry(subdialog, width="20",
							textvariable=get(valVar)))
#                        textvariable=eval(parse(text=valVar))))
					tkgrid(labelRcmdr(subdialog, text=values[i]), get(paste("entry", i, sep="")), sticky="w")
#                    tkgrid(labelRcmdr(subdialog, text=values[i]), eval(parse(text=paste("entry", i, sep=""))), sticky="w")
				}
				tkgrid(subButtonsFrame, sticky="w", columnspan=2)
				dialogSuffix(subdialog, rows=nvalues+2, columns=2, focus=entry1, onOK=onOKsub, force.wait=TRUE)
			}
			if (levelsType == "names"){
				if (!exists("labels", mode="character")) return()
				command <- paste("factor(", .activeDataSet, "$", name,
					", labels=c(", labels, "))", sep="")
				result <- justDoIt(paste(.activeDataSet, "$", fname, " <- ", command, sep=""))
				logger(paste(.activeDataSet,"$", fname," <- ", command, sep=""))
				if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet)
				tkfocus(CommanderWindow())
			}
			else{
				command <- paste("as.factor(", .activeDataSet, "$", name, ")", sep="")
				result <- justDoIt(paste(.activeDataSet, "$", fname, " <- ", command, sep=""))
				logger(paste(.activeDataSet, "$", fname," <- ", command, sep=""))
				if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
				tkfocus(CommanderWindow())
			}
		}
	}
	OKCancelHelp(helpSubject="factor")
	tkgrid(getFrame(variableBox), levelsFrame, sticky="nw")
	tkgrid(labelRcmdr(top,
			text=gettext(domain="R-RcmdrPlugin.EZR","New variable name or prefix for multiple variables:")),
		factorNameField, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	tkgrid.configure(numbersButton, sticky="w")
	tkgrid.configure(namesButton, sticky="w")
	dialogSuffix(rows=4, columns=2, preventGrabFocus=TRUE)
}

		
StatMedBinVariable <- function(){
# Author: Dan Putler (revision by J. Fox, 2 Feb 05)
#    if (!checkActiveDataSet()) return()
#    if (!checkNumeric()) return()
	env <- environment()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Bin a Numeric Variable"))
	variableFrame <- tkframe(top)
	variableBox <- variableListBox(variableFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable to bin (pick one)"), listHeight=15)
	newVariableFrame <- tkframe(variableFrame)
	newVariableName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","variable"))
	newVariable <- ttkentry(newVariableFrame, width="18", textvariable=newVariableName)
	binsFrame <- tkframe(top)
	binsVariable <- tclVar("3")
	slider <- tkscale(binsFrame, from=2, to=20, showvalue=TRUE, variable=binsVariable,
		resolution=1, orient="horizontal")
	optionsFrame <- tkframe(top)
	radioButtons(optionsFrame, name="levels", buttons=c("specify", "numbers", "ranges"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("Specify names", "Numbers", "Ranges")), title=gettext(domain="R-RcmdrPlugin.EZR","Level Names"))
	radioButtons(optionsFrame, name="method", buttons=c("intervals", "proportions", "natural"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("Equal-width bins", "Equal-count bins", "Natural breaks\n(from K-means clustering)")),
		title=gettext(domain="R-RcmdrPlugin.EZR","Binning Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Bin a Numeric Variable"), "#####", sep=""))
		levels <- tclvalue(levelsVariable)
		bins <- as.numeric(tclvalue(binsVariable))
		varName <- getSelection(variableBox)
		closeDialog()
		if (length(varName) == 0){
			errorCondition(recall=StatMedBinVariable, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		newVar <- tclvalue(newVariableName)
		if (is.element(newVar, Variables())) {
			if ("no" == tclvalue(checkReplace(newVar))){
				binVariable()
				return()
			}
		}
		if (!is.valid.name(newVar)){
			errorCondition(message=paste('"', newVar, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""),
				recall=StatMedBinVariable)
			return()
		}
		method <- tclvalue(methodVariable)
		if (levels == "specify"){
			initializeDialog(subdialog, title=gettext(domain="R-RcmdrPlugin.EZR","Bin Names"))
			onOKsub <- function() {
				closeDialog(subdialog)
				level <- character(bins)
				for (i in 1:bins){
					level[i] <- eval(parse(text=paste("tclvalue(levelName", i, ")", sep="")))
				}
				if (length(unique(level)) != length(level)){
					errorCondition(window=subdialog, message=gettext(domain="R-RcmdrPlugin.EZR","Level names must be unique."),
						recall=onOK)
					return()
				}
				assign("levelNames", level, envir=env)
			}
			subOKCancelHelp()
			tkgrid(labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Bin"), fg="blue"),
				labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Name"), fg="blue"), sticky="w")
			for (i in 1:bins){
				valVar <- paste("levelName", i, sep="")
				assign(valVar, tclVar(i))
				assign(paste("entry", i, sep=""), ttkentry(subdialog, width="20",
						textvariable=get(valVar)))
#                    textvariable=eval(parse(text=valVar))))
				tkgrid(labelRcmdr(subdialog, text=as.character(i)), get(paste("entry", i, sep="")), sticky="w")
#                tkgrid(labelRcmdr(subdialog, text=as.character(i)), eval(parse(text=paste("entry", i, sep=""))), sticky="w")
			}
			tkgrid(subButtonsFrame, sticky="w", columnspan=2)
			dialogSuffix(subdialog, focus=entry1, rows=bins+1, columns=2, bindReturn=FALSE, force.wait=TRUE)
		}
		labels <- if (levels == "numbers") "FALSE"
			else if (levels == "ranges") "NULL"
			else {
				if (!exists("levelNames")){
					onCancel()
					binVariable()
					return()
				}
				paste("c('", paste(levelNames,  collapse="','"), "')", sep="")
			}
		.activeDataSet <- ActiveDataSet()
		command <- paste(.activeDataSet,"$",newVar, " <- ",
			"bin.var(", .activeDataSet,"$", varName, ", bins=", bins,
			", method=", "'", method, "', labels=", labels, ")", sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
		logger(paste("#", gettext(domain="R-RcmdrPlugin.EZR","New variable"), " ", newVar, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), sep="") )
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="bin.var")
	tkgrid(labelRcmdr(newVariableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","New variable name"), fg="blue"), sticky="w")
	tkgrid(newVariable, sticky="w")
	tkgrid(getFrame(variableBox), labelRcmdr(variableFrame, text="    "), newVariableFrame, sticky="nw")
	tkgrid(variableFrame, sticky="w")
	tkgrid(labelRcmdr(binsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of bins:")), slider, sticky="s")
	tkgrid(binsFrame, sticky="w")
	tkgrid(levelsFrame, labelRcmdr(optionsFrame, text="    "), methodFrame, sticky="nw")
	tkgrid(optionsFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

	
StatMedFactorToNumeric <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Convert factors or character variables of numeric data to numeric variables"))
	variableBox <- variableListBox(top, Variables(), selectmode="multiple",
		title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"), listHeight=15)
	factorName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<same as variables>"))
	factorNameField <- ttkentry(top, width="20", textvariable=factorName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Convert factors or character variables of numeric data to numeric variables"), "#####", sep=""))
		variables <- getSelection(variableBox)
		closeDialog()
		if (length(variables) == 0) {
			errorCondition(recall=StatMedFactorToNumeric, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		facname <- trim.blanks(tclvalue(factorName))
		.activeDataSet <- ActiveDataSet()
		cmd <- paste("apply(", .activeDataSet, "[c(", paste(
				paste('"', variables, '"', sep=""),
				collapse=","), ")], 2, function(x) sort(unique(x)))", sep="")
		for (name in variables){
			fname <- if (facname == gettext(domain="R-RcmdrPlugin.EZR","<same as variables>")) name
				else if (length(variables) == 1) facname
				else paste(facname, name, sep="")
			if (!is.valid.name(fname)){
				errorCondition(recall=StatMedFactorToNumeric,
					message=paste('"', fname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
				return()
			}
			if (is.element(fname, Variables())) {
				if ("no" == tclvalue(checkReplace(fname))){
					numericToFactor()
					return()
				}
			}	
			command <- paste("as.numeric(as.character(", .activeDataSet, "$", name, "))", sep="")
			result <- justDoIt(paste(.activeDataSet, "$", fname, " <- ", command, sep=""))
			logger(paste(.activeDataSet, "$", fname," <- ", command, sep=""))
			if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
			tkfocus(CommanderWindow())
		}
	}
	OKCancelHelp(helpSubject="numeric")
	tkgrid(getFrame(variableBox), sticky="nw")
	tkgrid(labelRcmdr(top,
			text=gettext(domain="R-RcmdrPlugin.EZR","New variable name or prefix for multiple variables:")),
		factorNameField, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=4, columns=2, preventGrabFocus=TRUE)
}


StatMedReorderFactor <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Reorder Factor Levels"))
	variableBox <- variableListBox(top, Factors(), title=gettext(domain="R-RcmdrPlugin.EZR","Factor (pick one)"), listHeight=15)
	orderedFrame <- tkframe(top)
	orderedVariable <- tclVar("0")
	orderedCheckBox <- tkcheckbutton(orderedFrame, variable=orderedVariable)
	factorName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<same as original>"))
	factorNameField <- ttkentry(top, width="20", textvariable=factorName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Reorder Factor Levels"), "#####", sep=""))
		variable <- getSelection(variableBox)
		closeDialog()
		if (length(variable) == 0) {
			errorCondition(recall=StatMedReorderFactor, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		name <- trim.blanks(tclvalue(factorName))
		if (name == gettext(domain="R-RcmdrPlugin.EZR","<same as original>")) name <- variable
		if (!is.valid.name(name)){
			errorCondition(recall=StatMedReorderFactor,
				message=paste('"', name, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(name, Variables())) {
			if ("no" == tclvalue(checkReplace(name))){
				reorderFactor()
				return()
			}
		}
		.activeDataSet <- ActiveDataSet()
		old.levels <- eval(parse(text=paste("levels(", .activeDataSet, "$", variable, ")",
					sep="")), envir=.GlobalEnv)
		nvalues <- length(old.levels)
		ordered <- tclvalue(orderedVariable)
		if (nvalues > 30) {
			errorCondition(recall=StatMedReorderFactor,
				message=sprintf(gettext(domain="R-RcmdrPlugin.EZR","Number of levels (%d) too large."), nvalues))
			return()
		}
		initializeDialog(subdialog, title=gettext(domain="R-RcmdrPlugin.EZR","Reorder Levels"))
		order <- 1:nvalues
		onOKsub <- function() {
			closeDialog(subdialog)
			opt <- options(warn=-1)
			for (i in 1:nvalues){
				order[i] <- as.numeric(eval(parse(text=paste("tclvalue(levelOrder", i, ")", sep=""))))
			}
			options(opt)
			if (any(sort(order) != 1:nvalues) || any(is.na(order))){
				errorCondition(recall=StatMedReorderFactor,
					message=paste(gettext(domain="R-RcmdrPlugin.EZR","Order of levels must include all integers from 1 to "), nvalues, sep=""))
				return()
			}
			levels <- old.levels[order(order)]
			ordered <- if (ordered == "1") ", ordered=TRUE" else ""
			command <- paste("factor(", .activeDataSet, "$", variable,
				", levels=c(", paste(paste("'", levels, "'", sep=""), collapse=","), ")",
				ordered, ")", sep="")
			result <- justDoIt(paste(.activeDataSet, "$", name, " <- ", command, sep=""))
			logger(paste(.activeDataSet,"$", name," <- ", command, sep=""))
			if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
		}
		subOKCancelHelp()
		tkgrid(labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Old Levels"), fg="blue"),
			labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","New order"), fg="blue"), sticky="w")
		for (i in 1:nvalues){
			valVar <- paste("levelOrder", i, sep="")
			assign(valVar, tclVar(i))
			assign(paste("entry", i, sep=""), ttkentry(subdialog, width="2",
					textvariable=get(valVar)))
#                textvariable=eval(parse(text=valVar))))
			tkgrid(labelRcmdr(subdialog, text=old.levels[i]), get(paste("entry", i, sep="")), sticky="w")
#            tkgrid(labelRcmdr(subdialog, text=old.levels[i]), eval(parse(text=paste("entry", i, sep=""))), sticky="w")
		}
		tkgrid(subButtonsFrame, sticky="w", columnspan=2)
		dialogSuffix(subdialog, focus=entry1, rows=nvalues+1, columns=2, force.wait=TRUE)
	}
	OKCancelHelp(helpSubject="factor")
	tkgrid(getFrame(variableBox), sticky="nw")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Name for factor")), sticky="w")
	tkgrid(factorNameField, sticky="w")
	tkgrid(labelRcmdr(orderedFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Make ordered factor")), orderedCheckBox, sticky="w")
	tkgrid(orderedFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=5, columns=1, preventGrabFocus=TRUE)
}


StatMedDropUnusedFactorLevels <- function(){
    dataSet <- activeDataSet()
    initializeDialog(title=gettextRcmdr("Drop Unused Factor Levels"))
    allfactorsVariable <- tclVar("0")
    allFrame <- tkframe(top)
    allfactorsCheckBox <- ttkcheckbutton(allFrame, variable = allfactorsVariable)
    variablesBox <- variableListBox(top, Factors(),
        title=gettextRcmdr("Factors(s) to drop levels (pick one or more)"), selectmode="multiple",
        initialSelection=NULL)
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Drop unused factor levels"), "#####", sep=""))
        all <- tclvalue(allfactorsVariable)
        variables <- getSelection(variablesBox)
        closeDialog()
        if (all == 0 && length(variables) == 0) {
            errorCondition(recall=StatMedDropUnusedFactorLevels, message=gettextRcmdr("You must select one or more variables."))
            return()
        }
        response <- tclvalue(RcmdrTkmessageBox(message=gettextRcmdr("Drop unused factor levels\nPlease confirm."), 
            icon="warning", type="okcancel", default="cancel"))
        if (response == "cancel") {
            onCancel()
            return()
        }
        if (all == 1) command <- paste(dataSet, " <- droplevels(", dataSet, ")", sep="")
        else{
				command <- ""
				for (variable in variables){
                command <- paste(command, dataSet, "$", variable, " <- droplevels(", dataSet, "$", variable, ")\n", sep="")
            }
        }
        doItAndPrint(command)
        activeDataSet(dataSet, flushModel=FALSE, flushDialogMemory=FALSE)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="droplevels")
    tkgrid(allfactorsCheckBox, labelRcmdr(allFrame, text=gettextRcmdr("all factors")), sticky="w")
    tkgrid(allFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("OR"), fg="red"), sticky="w")
    tkgrid(getFrame(variablesBox), sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}


StatMedRecodeDialog <- function () {
  processRecode <- function(recode) {
    parts <- strsplit(recode, "=")[[1]]
    if (length(grep(",", parts[1])) > 0) 
      paste("c(", parts[1], ") = ", parts[2], sep = "")
    else paste(parts, collapse = "=")
  }
  dataSet <- activeDataSet()
  defaults <- list (initial.asFactor = 1, initial.variables = NULL, initial.name = gettextRcmdr ("variable"),
                    initial.recode.directives="")
  dialog.values <- getDialog ("StatMedRecodeDialog", defaults)
  initializeDialog(title = gettextRcmdr("Recode Variables"))
  variablesBox <- variableListBox(top, Variables(), selectmode = "multiple", 
                                  title = gettextRcmdr("Variables to recode (pick one or more)"),
                                  initialSelection = varPosn (dialog.values$initial.variables, "all"))
  variablesFrame <- tkframe(top)
  newVariableName <- tclVar(dialog.values$initial.name)
  newVariable <- ttkentry(variablesFrame, width = "20", textvariable = newVariableName)
  recodesFrame <- tkframe(top)
  recodes <- tktext(recodesFrame, bg = "white", font = getRcmdr("logFont"), 
                    height = "5", width = "40", wrap = "none")
  recodesXscroll <- ttkscrollbar(recodesFrame, orient = "horizontal", 
                                 command = function(...) tkxview(recodes, ...))
  recodesYscroll <- ttkscrollbar(recodesFrame, command = function(...) tkyview(recodes, 
                                                                               ...))
  tkconfigure(recodes, xscrollcommand = function(...) tkset(recodesXscroll, 
                                                            ...))
  tkconfigure(recodes, yscrollcommand = function(...) tkset(recodesYscroll, 
                                                            ...))
  tkinsert(recodes, "1.0", dialog.values$initial.recode.directives)
  asFactorFrame <- tkframe(top)
  asFactorVariable <- tclVar(dialog.values$initial.asFactor)
  asFactorCheckBox <- ttkcheckbutton(asFactorFrame, variable = asFactorVariable)
  onOK <- function() {
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Recode variables"), "#####", sep=""))
    asFactor <- tclvalue(asFactorVariable) == "1"
    save.recodes <- trim.blanks(tclvalue(tkget(recodes, "1.0", "end")))
    recode.directives <- gsub("\n", "; ", save.recodes)
    check.empty <- gsub(";", "", gsub(" ", "", recode.directives))
    if ("" == check.empty) {
      errorCondition(recall = StatMedRecodeDialog, message = gettextRcmdr("No recode directives specified."))
      return()
    }
    if (0 != length(grep("'", recode.directives))) {
      errorCondition(recall = StatMedRecodeDialog, message = gettextRcmdr("Use only double-quotes (\" \") in recode directives"))
      return()
    }
    recode.directives <- strsplit(recode.directives, ";")[[1]]
    recode.directives <- paste(sapply(recode.directives, 
                                      processRecode), collapse = ";")
    recode.directives <- sub(" *; *$", "", recode.directives)
    variables <- getSelection(variablesBox)
    closeDialog()
    if (length(variables) == 0) {
      errorCondition(recall = StatMedRecodeDialog, message = gettextRcmdr("You must select a variable."))
      return()
    }
    multiple <- if (length(variables) > 1) 
      TRUE
    else FALSE
    name <- trim.blanks(tclvalue(newVariableName))
    #        save.recodes <- gsub("; ", "\\\n", trim.blanks(recode.directives))  
    putDialog ("StatMedRecodeDialog", list (initial.asFactor = asFactor, initial.variables = variables,
                                     initial.name = name, initial.recode.directives=save.recodes))
#    command <- paste(dataSet, " <- within(", dataSet, ", {", sep="")
    nvar <- length(variables)
    for (i in 1:nvar) {
      variable <- variables[nvar - i + 1]
      newVar <- if (multiple) 
        paste(name, variable, sep = "")
      else name
      if (!is.valid.name(newVar)) {
        errorCondition(recall = StatMedRecodeDialog, message = paste("\"", 
                                                              newVar, "\" ", gettextRcmdr("is not a valid name."), 
                                                              sep = ""))
        return()
      }
      if (is.element(newVar, Variables())) {
        if ("no" == tclvalue(checkReplace(newVar))) {
          StatMedRecodeDialog()
          return()
        }
      }
#      command <- paste(command, "\n  ", newVar, " <- Recode(", variable, ", '", 
#                       recode.directives, "', as.factor.result=", asFactor, 
#                       ")", sep = "")  
      command <- paste(dataSet, "$", newVar, " <- Recode(", dataSet, "$", variable, ", '", 
                       recode.directives, "', as.factor.result=", asFactor, 
                       ")\n", sep = "")  
    }
#    command <- paste(command, "\n})", sep="")
    result <- doItAndPrint(command)
    if (class(result)[1] != "try-error")
      activeDataSet(dataSet, flushModel = FALSE, flushDialogMemory = FALSE)
    #     else{
    #       if (getRcmdr("use.markdown")) removeLastRmdBlock()
    #       if (getRcmdr("use.knitr")) removeLastRnwBlock()
    #    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "StatMedRecodeDialog", reset = "StatMedRecodeDialog", apply = "StatMedRecodeDialog")
  tkgrid(getFrame(variablesBox), sticky = "nw")
  tkgrid(labelRcmdr(variablesFrame, text = ""))
  tkgrid(labelRcmdr(variablesFrame, text = gettextRcmdr("New variable name or prefix for multiple recodes: ")), 
         newVariable, sticky = "w")
  tkgrid(asFactorCheckBox, labelRcmdr(asFactorFrame, text = gettextRcmdr("Make (each) new variable a factor")), 
         sticky = "w")
  tkgrid(labelRcmdr(asFactorFrame, text = ""))
  tkgrid(labelRcmdr(recodesFrame, text = gettextRcmdr("Enter recode directives"), 
                    fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
  tkgrid(recodes, recodesYscroll, sticky = "nw")
  tkgrid(recodesXscroll)
  tkgrid(variablesFrame, sticky = "w")
  tkgrid(asFactorFrame, sticky = "w")
  tkgrid(recodesFrame, sticky = "w")
  tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
  tkgrid.configure(recodesXscroll, sticky = "ew")
  tkgrid.configure(recodesYscroll, sticky = "ns")
  dialogSuffix(bindReturn = FALSE)
}


StatMedSetContrasts <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Define contrasts for a factor"))
	variableBox <- variableListBox(top, Factors(), title=gettext(domain="R-RcmdrPlugin.EZR","Factor (pick one)"), listHeight=15)
	radioButtons(name="contrasts", buttons=c("treatment", "sum", "helmert", "poly", "specify"),
		values=c("contr.Treatment", "contr.Sum", "contr.helmert", "contr.poly", "specify"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("Treatment (dummy) contrasts", "Sum (deviation) contrasts", "Helmert contrasts",
				"Polynomial contrasts", "Other (specify)")), title=gettext(domain="R-RcmdrPlugin.EZR","Contrasts"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Define contrasts for a factor"), "#####", sep=""))
		variable <- getSelection(variableBox)
		closeDialog()
		if (length(variable) == 0) {
			errorCondition(recall=StaMedSetContrasts, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		contrasts <- tclvalue(contrastsVariable)
		if (contrasts != "specify"){
			command <- paste("contrasts(", ActiveDataSet(), "$", variable, ') <- "', contrasts, '"', sep="")
			result <- justDoIt(command)
			logger(command)
			if (class(result)[1] !=  "try-error") activeDataSet(ActiveDataSet())
			tkfocus(CommanderWindow())
		}
		else{
			initializeDialog(subdialog, title=gettext(domain="R-RcmdrPlugin.EZR","Specify Contrasts"))
			tkgrid(labelRcmdr(subdialog, text=gettext(domain="R-RcmdrPlugin.EZR","Enter Contrast Coefficients"), fg="blue"), sticky="w")
			env <- environment()
			tableFrame <- tkframe(subdialog)
			row.names <- eval(parse(text=paste("levels(", ActiveDataSet(), "$", variable, ")")))
			row.names <- substring(paste(abbreviate(row.names, 12), "            "), 1, 12)
			nrows <- length(row.names)
			ncols <- nrows - 1
			make.col.names <- paste("labelRcmdr(tableFrame, text='", gettext(domain="R-RcmdrPlugin.EZR","Contrast Name:"), "')", sep="")
			for (j in 1:ncols) {
				varname <- paste(".col.", j, sep="")
				assign(varname, tclVar(paste(".", j, sep="")), envir=env)
				make.col.names <- paste(make.col.names, ", ",
					"ttkentry(tableFrame, width='12', textvariable=", varname, ")", sep="")
			}
			eval(parse(text=paste("tkgrid(", make.col.names, ", sticky='w')", sep="")), envir=env)
			for (i in 1:nrows){
				make.row <- paste("labelRcmdr(tableFrame, text='", row.names[i], "')")
				for (j in 1:ncols){
					varname <- paste(".tab.", i, ".", j, sep="")
					assign(varname, tclVar("0"), envir=env)
					make.row <- paste(make.row, ", ", "ttkentry(tableFrame, width='5', textvariable=",
						varname, ")", sep="")
				}
				eval(parse(text=paste("tkgrid(", make.row, ", sticky='w')", sep="")), envir=env)
			}
			tkgrid(tableFrame, sticky="w")
			onOKsub <- function(){
				closeDialog(subdialog)
				cell <- 0
				values <- rep(NA, nrows*ncols)
				for (j in 1:ncols){
					for (i in 1:nrows){
						cell <- cell + 1
						varname <- paste(".tab.", i, ".", j, sep="")
						values[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep=""))))
					}
				}
				values <- na.omit(values)
				if (length(values) != nrows*ncols){
					errorCondition(subdialog, recall=StatMedSetContrasts,
						message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",
								"Number of valid entries in contrast matrix(%d)\nnot equal to number of levels (%d) * number of contrasts (%d)."), length(values), nrows, ncols))
					return()
				}
				if (qr(matrix(values, nrows, ncols))$rank < ncols) {
					errorCondition(subdialog, recall=StatMedSetContrasts, message=gettext(domain="R-RcmdrPlugin.EZR","Contrast matrix is not of full column rank"))
					return()
				}
				contrast.names <- rep("", ncols)
				for (j in 1:ncols){
					varname <- paste(".col.", j, sep="")
					contrast.names[j] <- eval(parse(text=paste("tclvalue(", varname,")", sep="")))
				}
				if (length(unique(contrast.names)) < ncols) {
					errorCondition(subdialog, recall=StatMedSetContrasts, message=gettext(domain="R-RcmdrPlugin.EZR","Contrast names must be unique"))
					return()
				}
				command <- paste("matrix(c(", paste(values, collapse=","), "), ", nrows, ", ", ncols,
					")", sep="")
# 				assign(".Contrasts", justDoIt(command), envir=.GlobalEnv)
# 				logger(paste(".Contrasts <- ", command, sep=""))
				doItAndPrint(paste(".Contrasts <- ", command, sep=""))
				command <- paste("colnames(.Contrasts) <- c(",
					paste("'", contrast.names, "'", sep="", collapse=", "), ")", sep="")
				justDoIt(command)
				logger(command)
				command <- paste("contrasts(", ActiveDataSet(), "$", variable, ") <- .Contrasts", sep="")
				result <- justDoIt(command)
				logger(command)
				justDoIt("remove(.Contrasts, envir=.GlobalEnv)")
				logger("remove(.Contrasts)")
				if (class(result)[1] !=  "try-error") activeDataSet(ActiveDataSet(), flushModel=FALSE)
				tkfocus(CommanderWindow())
			}
			subOKCancelHelp(helpSubject="contrasts")
			tkgrid(tableFrame, sticky="w")
			tkgrid(labelRcmdr(subdialog, text=""))
			tkgrid(subButtonsFrame, sticky="w")
			dialogSuffix(subdialog, rows=5, columns=1, focus=subdialog, force.wait=TRUE)
		}
	}
	OKCancelHelp(helpSubject="contrasts")
	tkgrid(getFrame(variableBox), sticky="nw")
	tkgrid(contrastsFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

  
StatMedDummy <- function(){		
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Create dummy variables"))
  	dataSet <- activeDataSet()
	variablesBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Select one variable to make dummy variables"), listHeight=15)
	newVariableName <- tclVar(".Dummy.")
	newVariableNameEntry <- ttkentry(top, width="20", textvariable=newVariableName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Create dummy variables"), "#####", sep=""))
		var <- trim.blanks(getSelection(variablesBox))
        	if (length(var) == 0){
  	          errorCondition(recall=StatMedDummy, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            	return()
        }		
		newVar <- trim.blanks(tclvalue(newVariableName))
		if (newVar == ""){
			errorCondition(recall=StatMedDummy, message=gettext(domain="R-RcmdrPlugin.EZR","Define characters to indentify dummy variables"))
			return()
		}
		closeDialog()
		groups <- eval(parse(text=paste("levels(factor(", ActiveDataSet(), "$", var, "))", sep="")))
		ngroups <- length(groups)
		for (i in 1:ngroups) {
			newvarname <- paste(var, newVar, groups[i], sep="")
			for(j in 1:nchar(newvarname)){
				char <- substring(newvarname, j, j)
				substring(newvarname, j, j) <- ifelse(char=="/" | char=="*" | char=="-" | char=="+" | char==" " | char=="(" | char==")", ".", char)
			}
			if (!is.valid.name(newvarname)){
				errorCondition(recall=StatMedDummy,
					message=paste('"', newvarname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
				return()
			}
			command <- paste(ActiveDataSet(), "$", newvarname, " <- ifelse(", ActiveDataSet(), "$", var, '=="', groups[i], '", 1, 0)', sep="")
			result <- justDoIt(command)
			logger(command)
			logger(paste("###", gettext(domain="R-RcmdrPlugin.EZR","Dummy variable"), " ", newvarname, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), sep=""))
			doItAndPrint(paste("table(", ActiveDataSet(), "$", newvarname, ", exclude=NULL)", sep="") )
		}
		if (class(result)[1] !=  "try-error") activeDataSet(ActiveDataSet(), flushModel=FALSE)			
		logger(gettext(domain="R-RcmdrPlugin.EZR","Input all dummy variables except for the referece group into the model."))
		tkfocus(CommanderWindow())
	}
    OKCancelHelp()
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Characters to identify dummy variables")), newVariableNameEntry, sticky="w")
	tkgrid.configure(newVariableNameEntry, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=3, columns=2)
	}


StatMedDatediff <- function(){		
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Compute difference between two date variables"))
  	dataSet <- activeDataSet()
	startBox <- variableListBox(top, title=gettext(domain="R-RcmdrPlugin.EZR","Select start date"), listHeight=12)
	stopBox <- variableListBox(top, title=gettext(domain="R-RcmdrPlugin.EZR","Select end date"), listHeight=12)
	newVariableName <- tclVar("")
	newVariableNameEntry <- ttkentry(top, width="20", textvariable=newVariableName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Compute difference between two date variables"), "#####", sep=""))
		start <- trim.blanks(getSelection(startBox))
		stop <- trim.blanks(getSelection(stopBox))
		if (length(start) == 0 || length(stop) == 0){
  	          errorCondition(recall=StatMedDatediff, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            	return()
            	}
		newVar <- trim.blanks(tclvalue(newVariableName))
		if (!is.valid.name(newVar)){
			errorCondition(recall=StatMedDatediff,
				message=paste('"', newVar, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}		
        ymd <- as.character(tclvalue(ymdVariable))
		switch(ymd, 
			"a" = ymd <- ', "%Y-%m-%d"',
			"b" = ymd <- ', "%Y/%m/%d"',
			"c" = ymd <- ', "%y-%m-%d"',
			"d" = ymd <- ', "%y/%m/%d"',
			"e" = ymd <- ', "%m-%d-%Y"',
			"f" = ymd <- ', "%m/%d/%Y"',
			"g" = ymd <- ', "%m-%d-%y"',
			"h" = ymd <- ', "%m/%d/%Y"',
		)
		closeDialog()
				command <-  paste(dataSet,"$",newVar, " <- with(", dataSet, ", as.numeric(as.Date(", stop, ymd,") - as.Date(", start, ymd,")))", sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
		logger(paste("#", gettext(domain="R-RcmdrPlugin.EZR","New variable"), " ", newVar, " ", gettext(domain="R-RcmdrPlugin.EZR","was made."), sep="") )
		tkfocus(CommanderWindow())
	}
    OKCancelHelp(helpSubject="as.Date")
	radioButtons(name="ymd",
        buttons=c("A", "B", "C", "D", "E", "F", "G", "H"),
		values=c("a", "b", "c", "d", "e", "f", "g", "h"), initialValue="a",
        labels=c("1999-12-31", "1999/12/31", "99-12-31", "99/12/31", "12-31-1999", "12/31/1999", "12-31-99", "12/31/99"), 
		title=gettext(domain="R-RcmdrPlugin.EZR","Select format"))
    tkgrid(ymdFrame, sticky="w")
	tkgrid(getFrame(startBox), getFrame(stopBox), sticky="nw")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","New variable name"), fg="blue"), newVariableNameEntry, sticky="w")
	tkgrid.configure(newVariableNameEntry, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=3, columns=2)
	}


StatMedRenewDataSet <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Renew active data set"))
	checkBoxes(frame="chrtofac", boxes=c("chrtofac"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR","Convert all character variables to factors"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Renew active data set"), "#####", sep=""))
		dataSet <- activeDataSet()
		chrtofac <- tclvalue(chrtofacVariable)
		closeDialog()
		if(chrtofac==1){
			doItAndPrint(paste(dataSet, " <- ChrToFactor(", dataSet, ")", sep="")) 
		}
		activeDataSet(dataSet, flushModel=FALSE)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="factor")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Renew currently active data set"), fg="blue"), sticky="w")
	tkgrid(chrtofac, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}


StatMedChrToFactor <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Convert all character variables to factors"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Convert all character variables to factors"), "#####", sep=""))
		dataSet <- activeDataSet()
		closeDialog()
		doItAndPrint(paste(dataSet, " <- ChrToFactor(", dataSet, ")", sep="")) 
		activeDataSet(dataSet, flushModel=FALSE)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="factor")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Convert all character variables to factors."), fg="blue"), sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}


StatMedGraphOptions <- function(){
defaults <- list(window.size="Medium", window.type="width=7, height=7", lwd="1", las="1", family="sans", cex="1")
dialog.values <- getDialog("StatMedGraphOptions", defaults)
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Graph settings"))
    optionsFrame <- tkframe(top)	
    radioButtons(optionsFrame, name="window.size", buttons=gettext(domain="R-RcmdrPlugin.EZR",c("Small", "Medium", "Large")), values=c("Small", "Medium", "Large"), initialValue=dialog.values$window.size, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Small", "Medium", "Large")), title=gettext(domain="R-RcmdrPlugin.EZR","Graph size"))
    radioButtons(optionsFrame, name="window.type", buttons=gettext(domain="R-RcmdrPlugin.EZR",c("Square", "Horizontal", "Vertical")), values=c("width=7, height=7", "width=9, height=6", "width=6, height=9"), initialValue=dialog.values$window.type, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Square", "Horizontal rectangle", "Vertical rectangle")), title=gettext(domain="R-RcmdrPlugin.EZR","Graph shape"))
    radioButtons(optionsFrame, name="lwd", buttons=gettext(domain="R-RcmdrPlugin.EZR",c("Thin", "Medium", "Thick")), values=c("1", "2", "3"), initialValue=dialog.values$lwd, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Thin", "Medium", "Thick")), title=gettext(domain="R-RcmdrPlugin.EZR","Line width"))
	radioButtons(optionsFrame, name="las", buttons=gettext(domain="R-RcmdrPlugin.EZR",c("ParallelAxis", "Horizontal", "PerpendicularAxis", "Vertical")), values=c("0", "1", "2", "3"), initialValue=dialog.values$las, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Parallel to axis", "Horizontal", "Perpendicular to axis", "Vertical")), title=gettext(domain="R-RcmdrPlugin.EZR","Axis label style"))
	radioButtons(optionsFrame, name="family", buttons=c("standard", "sans", "serif", "mono"), values=c("", "sans", "serif", "mono"), initialValue=dialog.values$family, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Standard", "Sans", "Serif", "Mono")), title=gettext(domain="R-RcmdrPlugin.EZR","Font"))
	radioButtons(optionsFrame, name="cex", buttons=gettext(domain="R-RcmdrPlugin.EZR",c("Small", "Medium", "Large")), values=c("1", "1.25", "1.5"), initialValue=dialog.values$cex, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Small", "Medium", "Large")), title=gettext(domain="R-RcmdrPlugin.EZR","Font size"))
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Graph settings"), "#####", sep=""))
#		doItAndPrint('unlockBinding("lwd", as.environment("package:RcmdrPlugin.EZR"))')
#		doItAndPrint('unlockBinding("lwd", as.environment("package:Rcmdr"))')
#		justDoIt('unlockBinding("las", as.environment("package:RcmdrPlugin.EZR"))')
#		justDoIt('unlockBinding("font", as.environment("package:RcmdrPlugin.EZR"))')
#		justDoIt('unlockBinding("cex", as.environment("package:RcmdrPlugin.EZR"))')
#		justDoIt('unlockBinding("window.type", as.environment("package:RcmdrPlugin.EZR"))')
#		justDoIt('unlockBinding("par.option", as.environment("package:RcmdrPlugin.EZR"))')
        size <- tclvalue(window.sizeVariable)
		type <- tclvalue(window.typeVariable)
		lwd <- tclvalue(lwdVariable)
        las <- tclvalue(lasVariable)
        font <- tclvalue(familyVariable)
        cex <- tclvalue(cexVariable)
		putDialog("StatMedGraphOptions", list(window.size=size, window.type=type, lwd=lwd, las=las, family=font, cex=cex))
        closeDialog()
		if (size=="Medium"){
			switch(type,
#				"width=7, height=7" = window.type <<- "width=10.5, height=10.5",
#				"width=9, height=6" = window.type <<- "width=13.5, height=9",
#				"width=6, height=9" = window.type <<- "width=9, height=13.5"
#				"width=7, height=7" = assign("window.type", "width=10.5, height=10.5", envir=.GlobalEnv),
#				"width=9, height=6" = assign("window.type", "width=13.5, height=9", envir=.GlobalEnv),
#				"width=6, height=9" = assign("window.type", "width=9, height=13.5", envir=.GlobalEnv)
				"width=7, height=7" = justDoIt('window.type <- "width=7, height=7"'),
				"width=9, height=6" = justDoIt('window.type <- "width=9, height=6"'),
				"width=6, height=9" = justDoIt('window.type <- "width=6, height=9"')
				)	
		}
		if (size=="Large"){
			switch(type,
#				"width=7, height=7" = window.type <<- "width=10.5, height=10.5",
#				"width=9, height=6" = window.type <<- "width=13.5, height=9",
#				"width=6, height=9" = window.type <<- "width=9, height=13.5"
#				"width=7, height=7" = assign("window.type", "width=10.5, height=10.5", envir=.GlobalEnv),
#				"width=9, height=6" = assign("window.type", "width=13.5, height=9", envir=.GlobalEnv),
#				"width=6, height=9" = assign("window.type", "width=9, height=13.5", envir=.GlobalEnv)
				"width=7, height=7" = justDoIt('window.type <- "width=10.5, height=10.5"'),
				"width=9, height=6" = justDoIt('window.type <- "width=13.5, height=9"'),
				"width=6, height=9" = justDoIt('window.type <- "width=9, height=13.5"')
				)	
		}
		if (size=="Small"){
			switch(type,
#				"width=7, height=7" = window.type <<- "width=5, height=5",
#				"width=9, height=6" = window.type <<- "width=6, height=4",
#				"width=6, height=9" = window.type <<- "width=4, height=6"
#				"width=7, height=7" = assign("window.type", "width=5, height=5", envir=.GlobalEnv),
#				"width=9, height=6" = assign("window.type", "width=6, height=4", envir=.GlobalEnv),
#				"width=6, height=9" = assign("window.type", "width=4, height=6", envir=.GlobalEnv)
				"width=7, height=7" = justDoIt('window.type <- "width=5, height=5"'),
				"width=9, height=6" = justDoIt('window.type <- "width=6, height=4"'),
				"width=6, height=9" = justDoIt('window.type <- "width=4, height=6"')
				)	
		}		
		window.type <- get("window.type", envir=.GlobalEnv)
		#		doItAndPrint(paste("windows(", window.type, ")", sep=""))
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", window.type, ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", window.type, ")", sep=""))} else {doItAndPrint(paste("x11(", window.type, ")", sep=""))}
		
#		par.option <<- paste("lwd=", lwd, ", las=", las, ', family="', font, '", cex=', cex, sep="")
#		assign("par.option", paste("lwd=", lwd, ", las=", las, ', family="', font, '", cex=', cex, sep=""), envir=.GlobalEnv)
		par.option <- paste("lwd=", lwd, ", las=", las, ', family="', font, '", cex=', cex, ", mgp=c(3.0,1,0)", sep="")
		justDoIt(paste("par.option <- '", par.option, "'", sep=""))
		par.option <- get("par.option", envir=.GlobalEnv)

		par.lwd <- paste("lwd=", lwd, sep="")
		justDoIt(paste("par.lwd <- '", par.lwd, "'", sep=""))
		justDoIt(paste("par.cex <- '", cex, "'", sep=""))
		
		doItAndPrint(paste("par(", par.option, ")", sep=""))
		doItAndPrint('plot(sin, xlim=c(0,10), main="Sample")')
		tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="par", apply="StatMedGraphOptions", reset="StatMedGraphOptions")
    tkgrid(window.sizeFrame, labelRcmdr(optionsFrame, text="    "), window.typeFrame, labelRcmdr(optionsFrame, text="    "), lwdFrame, labelRcmdr(optionsFrame, text="    "), lasFrame, labelRcmdr(optionsFrame, text="    "), familyFrame, labelRcmdr(optionsFrame, text="    "), cexFrame, sticky="w")
	tkgrid(optionsFrame, sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
    }

	
StatMedChangePalette <- function(){
defaults <- list(palette="Standard")
dialog.values <- getDialog("StatMedChangePalette", defaults)
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Graph colors"))
    optionsFrame <- tkframe(top)	
    radioButtons(optionsFrame, name="palette", buttons=c("Standard", "Gray4", "Gray8", "Heat", "Cold"), values=c("Standard", "Gray4", "Gray8", "Heat", "Cold"), initialValue=dialog.values$palette,
       labels=gettext(domain="R-RcmdrPlugin.EZR",c("Standard color", "Gray (4 levels)", "Gray (8 levels)", "Heat", "Cold")), title=gettext(domain="R-RcmdrPlugin.EZR","Colors"))
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Graph colors"), "#####", sep=""))	
        palette.type <- tclvalue(paletteVariable)
putDialog("StatMedChangePalette", list(palette=palette.type))
        closeDialog()
		switch (palette.type,
			"Standard"=doItAndPrint('palette("default")'),
			"Gray4"=doItAndPrint("palette(gray(rep(c(0, 0.3, 0.6, 0.9),2)))"),
			"Gray8"=doItAndPrint("palette(gray(seq(0, 1, length=8)))"),
			"Heat"=doItAndPrint("palette(heat.colors(8))"),
			"Cold"=doItAndPrint("palette(cm.colors(8))")
		)
		
		if(getRversion() < '3.0.0') {
#        doItAndPrint(paste("windows(", window.type, "); par(", par.option, ")", sep=""))
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		}
		doItAndPrint('plot(0,1, type="n", yaxt="n", ylab="", xlim=c(0,9), ylim=c(0,1), xlab="Color number")')
		doItAndPrint("for (i in 1:8) {rect(i-0.5, 0.05, i+0.5, 0.95, col = i)}")
		doItAndPrint("axis(1, at=1:8)")
		tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="palette", apply="StatMedChangePalette", reset="StatMedChangePalette")
    tkgrid(paletteFrame, labelRcmdr(optionsFrame, text="    "), sticky="w")
	tkgrid(optionsFrame, sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
    }
	
		
StatMedSetPalette <- function() {
    cval <- function(x,y) -sum((x-y)^2)
    contrasting <- function(x)
        optim(rep(127, 3),cval,lower=0,upper=255,method="L-BFGS-B",y=x)$par
    # the following local function from Thomas Lumley via r-help
    convert <- function (color){
        rgb <- col2rgb(color)/255
        L <- c(0.2, 0.6, 0) %*% rgb
        ifelse(L >= 0.2, "#000060", "#FFFFA0")
        }
    env <- environment()
    pal <- palette()
    pickColor <- function(initialcolor, parent){
        tclvalue(.Tcl(paste("tk_chooseColor", .Tcl.args(title = "Select a Color",
            initialcolor=initialcolor, parent=parent))))
        }
	Library("tcltk")
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Graph detailed colors"))
    hexcolor <- colorConverter(toXYZ = function(hex,...) {
        rgb <- t(col2rgb(hex))/255
        colorspaces$sRGB$toXYZ(rgb,...) },
        fromXYZ = function(xyz,...) {
            rgb <- colorspaces$sRGB$fromXYZ(xyz,..)
            rgb <- round(rgb,5)
            if (min(rgb) < 0 || max(rgb) > 1) as.character(NA)
            else rgb(rgb[1],rgb[2],rgb[3])},
            white = "D65", name = "#rrggbb")
    cols <- t(col2rgb(pal))
    hex <- convertColor(cols, from="sRGB", to=hexcolor, scale.in=255, scale.out=NULL)
    for (i in 1:8) assign(paste("hex", i, sep="."), hex[i], envir=env)
    paletteFrame <- tkframe(top)
    button1 <- tkbutton(paletteFrame, text=hex[1], bg = hex[1],
        fg=convert(hex[1]),
        command=function() {
            color <- pickColor(hex[1], parent=button1)
            fg <- convert(color)
            tkconfigure(button1, bg=color, fg=fg)
            assign("hex.1", color, envir=env)
            }
        )
    button2 <- tkbutton(paletteFrame, text=hex[2], bg = hex[2],
        fg=convert(hex[2]),
        command=function() {
            color <- pickColor(hex[2], parent=button2)
            fg <- convert(color)
            tkconfigure(button2, bg=color, fg=fg)
            assign("hex.2", color, envir=env)
            }
        )
     button3 <- tkbutton(paletteFrame, text=hex[3], bg = hex[3],
        fg=convert(hex[3]),
        command=function() {
            color <- pickColor(hex[3], parent=button3)
            fg <- convert(color)
            tkconfigure(button3, bg=color, fg=fg)
            assign("hex.3", color, envir=env)
            }
        )
     button4 <- tkbutton(paletteFrame, text=hex[4], bg = hex[4],
        fg=convert(hex[4]),
        command=function() {
            color <- pickColor(hex[4], parent=button4)
            fg <- convert(color)
            tkconfigure(button4, bg=color, fg=fg)
            assign("hex.4", color, envir=env)
            }
        )
     button5 <- tkbutton(paletteFrame, text=hex[5], bg = hex[5],
        fg=convert(hex[5]),
        command=function() {
            color <- pickColor(hex[5], parent=button5)
            fg <- convert(color)
            tkconfigure(button5, bg=color, fg=fg)
            assign("hex.5", color, envir=env)
            }
        )
     button6 <- tkbutton(paletteFrame, text=hex[6], bg = hex[6],
        fg=convert(hex[6]),
        command=function() {
            color <- pickColor(hex[6], parent=button6)
            fg <- convert(color)
            tkconfigure(button6, bg=color, fg=fg)
            assign("hex.6", color, envir=env)
            }
        )
     button7 <- tkbutton(paletteFrame, text=hex[7], bg = hex[7],
        fg=convert(hex[7]),
        command=function() {
            color <- pickColor(hex[7], parent=button7)
            fg <- convert(color)
            tkconfigure(button7, bg=color, fg=fg)
            assign("hex.7", color, envir=env)
            }
        )
     button8 <- tkbutton(paletteFrame, text=hex[8], bg = hex[8],
        fg=convert(hex[8]),
        command=function() {
            color <- pickColor(hex[8], parent=button8)
            fg <- convert(color)
            tkconfigure(button8, bg=color, fg=fg)
            assign("hex.8", color, envir=env)
            }
        )
     onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Graph detailed colors"), "#####", sep=""))
        closeDialog(top)
		if(getRversion() < '3.0.0') {
#        doItAndPrint(paste("windows(", window.type, "); par(", par.option, ")", sep=""))
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		}
        palette(c(hex.1, hex.2, hex.3, hex.4, hex.5, hex.6, hex.7, hex.8))
		logger(paste('palette(c("', hex.1, '", "', hex.2, '", "', hex.3, '", "', hex.4, '", "', hex.5, '", "', hex.6, '", "', hex.7, '", "', hex.8, '"))', sep=""))
		doItAndPrint('plot(0,1, type="n", yaxt="n", ylab="", xlim=c(0,9), ylim=c(0,1), xlab="Color number")')
		doItAndPrint("for (i in 1:8) {rect(i-0.5, 0.05, i+0.5, 0.95, col = i)}")
		doItAndPrint("axis(1, at=1:8)")
        Message(gettext(domain="R-RcmdrPlugin.EZR","Color palette reset.", type="note"))
        }
    OKCancelHelp(helpSubject="palette")
    tkgrid(button1, button2, button3, button4, button5, button6, button7, button8)
    tkgrid(paletteFrame)
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2)
    }

	
StatMedNumericalSummaries <- function(){
	Library("tcltk")
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Numerical Summaries"))
    xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"), listHeight=15)
    checkBoxes(frame="checkBoxFrame", boxes=c("graph", "mean", "sd"), initialValues=c("0", "1", "1"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show graph", "Mean", "Standard Deviation")))
    quantilesVariable <- tclVar("1")
    quantilesFrame <- tkframe(top)
    quantilesCheckBox <- tkcheckbutton(quantilesFrame, variable=quantilesVariable)
    quantiles <- tclVar("0, .25, .5, .75, 1")
    quantilesEntry <- ttkentry(quantilesFrame, width="20", textvariable=quantiles)
    StatMedGroupsBox(recall=StatMedNumericalSummaries, label=gettext(domain="R-RcmdrPlugin.EZR","Summarize by:"), initialLabel=gettext(domain="R-RcmdrPlugin.EZR","Summarize by groups"))
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Numerical summaries"), "#####", sep=""))
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=StatMedNumericalSummaries, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            return()
            }
		graph <- tclvalue(graphVariable)
        closeDialog()
        quants <- paste("c(", gsub(",+", ",", gsub(" ", ",", tclvalue(quantiles))), ")", sep="")
        .activeDataSet <- ActiveDataSet()
        vars <- if (length(x) == 1) paste('"', x, '"', sep="") 
            else paste("c(", paste('"', x, '"', collapse=", ", sep=""), ")", sep="")
        vars <- paste(.activeDataSet, "[,", vars, "]", sep="")
        stats <- paste("c(",
            paste(c('"mean"', '"sd"', '"quantiles"')
                [c(tclvalue(meanVariable), tclvalue(sdVariable), tclvalue(quantilesVariable)) == 1], 
                collapse=", "), ")", sep="")
        if (stats == "c()"){
             errorCondition(recall=StatMedNumericalSummaries, message=gettext(domain="R-RcmdrPlugin.EZR","No statistics selected."))
            return()
            }               
        command <- if (.groups != FALSE) {
            grps <- paste(.activeDataSet, "$", .groups, sep="")
            paste("res <- numSummary(", vars, ", groups=", grps, ", statistics=", stats, 
				", quantiles=", quants, ")", sep="")
            }
        else  paste("res <- numSummary(", vars, ", statistics=", stats, 
			", quantiles=", quants, ")", sep="")
        doItAndPrint(command) 
		doItAndPrint('colnames(res$table)[1:2] <- gettext(domain="R-RcmdrPlugin.EZR", colnames(res$table)[1:2])')
		if (graph==1){
			for (i in 1:length(x)){			
				if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
				if (.groups == FALSE){
					doItAndPrint(paste("dummyX <- rep(0, length(", .activeDataSet, "$", x[i], "))"))		
					doItAndPrint(paste("dot.plot(dummyX, ", .activeDataSet, "$", x[i], ', xlab="", ylab="', x[i], '")', sep=""))			
				} else {
					groupNames <- paste(.activeDataSet, "$", .groups, collapse="*")
					doItAndPrint(paste("dot.plot(", .activeDataSet, "$", .groups, ", ", .activeDataSet, "$", x[i], ', xlab="', .groups, '", ylab="', x[i], '")', sep=""))
				}
			}
		}
		doItAndPrint("res")
		doItAndPrint("remove(res)")
		tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="numSummary")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
    tkgrid(getFrame(xBox), sticky="nw")    
    tkgrid(checkBoxFrame, sticky="w")
    tkgrid(labelRcmdr(quantilesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Quantiles")), quantilesCheckBox,
        labelRcmdr(quantilesFrame, text=gettext(domain="R-RcmdrPlugin.EZR"," quantiles:")), quantilesEntry, sticky="w")
    tkgrid(quantilesFrame, sticky="w")
    tkgrid(groupsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
    }

	
StatMedQQPlot <- function () {
# this function modified by Martin Maechler
	requireNamespace("car")
	defaults <- list(initial.x = NULL, initial.identify = 0, initial.dist = "norm", initial.df = "",
			initial.chisqdf = "", initial.fdf1 = "", initial.fdf2 = "", initial.othername = "", 
			initial.otherparam = "")
	dialog.values <- getDialog("StatMedQQPlot", defaults)
	initializeDialog(title = gettextRcmdr("Quantile-Comparison (QQ) Plot"))
	xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), listHeight=15, 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	onOK <- function() {
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Quantile-comparison plot"), "#####", sep=""))
		x <- getSelection(xBox)
		initial.dist <-dist <- tclvalue(distVariable)
		identify <- tclvalue(identifyVariable)
		tdf <- tclvalue(tDfVariable)
		chisqdf <- tclvalue(chisqDfVariable)
		fdf1 <- tclvalue(FDf1Variable)
		fdf2 <- tclvalue(FDf2Variable)
		othername <- tclvalue(otherNameVariable)
		otherparam <- tclvalue(otherParamsVariable)
		putDialog ("StatMedQQPlot", list (initial.x = x, initial.dist = initial.dist,
						initial.identify = identify, initial.df = tdf, initial.chisqdf = chisqdf,
						initial.fdf1 = fdf1, initial.fdf2 = fdf2, initial.othername = othername, 
						initial.otherparam = otherparam))
		closeDialog()
		if (0 == length(x)) {
			errorCondition(recall = StatMedQQPlot, message = gettextRcmdr("You must select a variable."))
			return()
		}
		save <- options(warn = -1)
		on.exit(save)
		retryMe <- function(msg) {
			Message(message = msg, type = "error")
			QQPlot()
		}
		switch(dist, norm = {
					args <- "dist=\"norm\""
				}, t = {
					df <- tclvalue(tDfVariable)
					df.num <- as.numeric(df)
					if (is.na(df.num) || df.num < 1) {
						retryMe(gettextRcmdr("df for t must be a positive number."))
						return()
					}
					args <- paste("dist=\"t\", df=", df, sep = "")
				}, chisq = {
					df <- tclvalue(chisqDfVariable)
					df.num <- as.numeric(df)
					if (is.na(df.num) || df.num < 1) {
						retryMe(gettextRcmdr("df for chi-square must be a positive number."))
						return()
					}
					args <- paste("dist=\"chisq\", df=", df, sep = "")
				}, f = {
					df1 <- tclvalue(FDf1Variable)
					df2 <- tclvalue(FDf2Variable)
					df.num1 <- as.numeric(df1)
					df.num2 <- as.numeric(df2)
					if (is.na(df.num1) || df.num1 < 1 || is.na(df.num2) || 
							df.num2 < 1) {
						retryMe(gettextRcmdr("numerator and denominator \ndf for F must be positive numbers."))
						return()
					}
					args <- paste("dist=\"f\", df1=", df1, ", df2=", 
							df2, sep = "")
				}, {
					dist <- tclvalue(otherNameVariable)
					params <- tclvalue(otherParamsVariable)
					args <- paste("dist=\"", dist, "\", ", params, sep = "")
				})
		.activeDataSet <- ActiveDataSet()
		if ("1" == tclvalue(identifyVariable)) {
			RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
							gettextRcmdr(if (MacOSXP()) 
												"esc key to exit."
											else "right button to exit."), sep = ""), icon = "info", 
					type = "ok")
			idtext <- paste(", labels=rownames(", .activeDataSet, 
					"), id.method=\"identify\"", sep = "")
		}
		else idtext <- ""
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("qqPlot", "(", .activeDataSet, "$", 
				x, ", ", args, idtext, ")", sep = "")
		doItAndPrint(command)
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "qqPlot", apply="StatMedQQPlot", reset = "StatMedQQPlot")
	distFrame <- tkframe(top)
	distVariable <- tclVar(dialog.values$initial.dist)
	normalButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "norm")
	tButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "t")
	chisqButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "chisq")
	FButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "f")
	otherButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "other")
	tDfFrame <- tkframe(distFrame)
	tDfVariable <- tclVar(dialog.values$initial.df)
	tDfField <- ttkentry(tDfFrame, width = "6", textvariable = tDfVariable)
	chisqDfFrame <- tkframe(distFrame)
	chisqDfVariable <- tclVar(dialog.values$initial.chisqdf)
	chisqDfField <- ttkentry(chisqDfFrame, width = "6", textvariable = chisqDfVariable)
	FDfFrame <- tkframe(distFrame)
	FDf1Variable <- tclVar(dialog.values$initial.fdf1)
	FDf1Field <- ttkentry(FDfFrame, width = "6", textvariable = FDf1Variable)
	FDf2Variable <- tclVar(dialog.values$initial.fdf2)
	FDf2Field <- ttkentry(FDfFrame, width = "6", textvariable = FDf2Variable)
	otherParamsFrame <- tkframe(distFrame)
	otherParamsVariable <- tclVar(dialog.values$initial.otherparam)
	otherParamsField <- ttkentry(otherParamsFrame, width = "30", 
			textvariable = otherParamsVariable)
	otherNameVariable <- tclVar(dialog.values$initial.othername)
	otherNameField <- ttkentry(otherParamsFrame, width = "10", 
			textvariable = otherNameVariable)
	identifyVariable <- tclVar(dialog.values$initial.identify)
	identifyFrame <- tkframe(top)
	identifyCheckBox <- tkcheckbutton(identifyFrame, variable = identifyVariable)
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(identifyFrame, text = gettextRcmdr("Identify observations with mouse")), 
			identifyCheckBox, sticky = "w")
	tkgrid(identifyFrame, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Distribution"), 
					fg = "blue"), columnspan = 6, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Normal")), 
			normalButton, sticky = "w")
	tkgrid(labelRcmdr(tDfFrame, text = gettextRcmdr("df = ")), 
			tDfField, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = "t"), tButton, tDfFrame, 
			sticky = "w")
	tkgrid(labelRcmdr(chisqDfFrame, text = gettextRcmdr("df = ")), 
			chisqDfField, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Chi-square")), 
			chisqButton, chisqDfFrame, sticky = "w")
	tkgrid(labelRcmdr(FDfFrame, text = gettextRcmdr("Numerator df = ")), 
			FDf1Field, labelRcmdr(FDfFrame, text = gettextRcmdr("Denominator df = ")), 
			FDf2Field, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = "F"), FButton, FDfFrame, 
			sticky = "w")
	tkgrid(labelRcmdr(otherParamsFrame, text = gettextRcmdr("Specify: ")), 
			otherNameField, labelRcmdr(otherParamsFrame, text = gettextRcmdr("Parameters: ")), 
			otherParamsField, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Other")), 
			otherButton, otherParamsFrame, sticky = "w")
	tkgrid(distFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 5, columns = 1)
}


StatMedHistogram <- function(){
defaults <- list(x=NULL, group=NULL, color=0, bins="<auto>", scale="frequency", subset="")
dialog.values <- getDialog("StatMedHistogram", defaults)
currentFields$subset <- dialog.values$subset	#Valued of currentFields will be sent to subsetBox
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Histogram"))
    variablesFrame <- tkframe(top)
    xBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$x, "numeric"))
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    checkBoxes(frame="color", boxes=c("color"),initialValues=dialog.values$color,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Draw in color (when grouped)")))
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Histogram"), "#####", sep=""))
        x <- getSelection(xBox)
        group <- getSelection(groupBox)
		color <- tclvalue(colorVariable)
        bins <- tclvalue(binsVariable)
        scale <- tclvalue(scaleVariable)
		subset <- tclvalue(subsetVariable)
putDialog("StatMedHistogram", list(x=x, group=group, color=color, bins=tclvalue(binsVariable), scale=scale, subset=tclvalue(subsetVariable)))
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
        closeDialog()
        if (length(x) == 0){
            errorCondition(recall=StatMedHistogram, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable"))
            return()
            }
        opts <- options(warn=-1)
#        bins <- if (bins == gettext(domain="R-RcmdrPlugin.EZR","<auto>")) '"Sturges"' else as.numeric(bins)
        bins <- if (bins == gettext(domain="R-RcmdrPlugin.EZR","<auto>")) '"scott"' else as.numeric(bins)-1	#chabge default to Scott, bins <- bins - 1
        options(opts)
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		if (length(group)==0) {
			command <- paste("HistEZR(", subset1, ActiveDataSet(), subset2, "$", x, ', scale="',
            scale, '", breaks=', bins, ', xlab="', x, '", col="darkgray")', sep="")
			doItAndPrint(command)
        } else {
			groups <- eval(parse(text=paste("levels(factor(", subset1, ActiveDataSet(), subset2, "$", group, "))", sep="")))
			ngroup <- length(groups)
			if (color == 0){
				color <- NULL
			} else {
				color <- paste(", col=c(2:", ngroup+1, ")", sep="")
			}	
			doItAndPrint(paste("breaks <- hist(", subset1, ActiveDataSet(), subset2, "$", x, ", breaks='scott', right=FALSE, plot=FALSE)$breaks", sep=""))
			command <- paste("temp.y <- ", subset1, ActiveDataSet(), subset2, "[", subset1, ActiveDataSet(), subset2, "$", group, '=="', groups[1], '",]$', x, sep="")
			doItAndPrint(command)
			doItAndPrint("temp.h <- hist(temp.y, breaks=breaks, right=FALSE, plot=FALSE)$counts")
			if (ngroup >=2){
				for (i in 2:ngroup){
					command <- paste("temp.y <- ", subset1, ActiveDataSet(), subset2, "[", subset1, ActiveDataSet(), subset2, "$", group, '=="', groups[i], '",]$', x, sep="")
					doItAndPrint(command)
					doItAndPrint("temp.h <- rbind(temp.h, hist(temp.y, breaks=breaks, right=FALSE, plot=FALSE)$counts)")			
				}
				command <- paste("barplot(temp.h, beside=TRUE, space=c(0, 0.4), names.arg=breaks[-length(temp.h[1])], legend=levels(factor(", ActiveDataSet(), "$", group, ')), args.legend=list(title="', group, '", box.lty=0), axis.lty=1, axisnames=TRUE', color, ")", sep="")
				doItAndPrint(command)					
			}
		}
		activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Hist", apply="StatMedHistogram", reset="StatMedHistogram")
    radioButtons(name="scale", buttons=c("frequency", "percent", "density"), initialValue=dialog.values$scale,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Frequency counts", "Percentages", "Densities")), title=gettext(domain="R-RcmdrPlugin.EZR","Y axis (when not grouped)"))
    binsFrame <- tkframe(top)
    binsVariable <- tclVar(gettext(domain="R-RcmdrPlugin.EZR",dialog.values$bins))
    binsField <- ttkentry(binsFrame, width="6", textvariable=binsVariable)
	tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(color, sticky="w")
    tkgrid(labelRcmdr(binsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of sections (when not grouped)")), binsField, sticky="w")
    tkgrid(binsFrame, sticky="w")
    tkgrid(scaleFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(binsField, sticky="e")
    dialogSuffix(rows=4, columns=1)
    }

	
StatMedStemAndLeaf <- function(){
	defaults <- list(initial.x = NULL, initial.leafs.auto="1", initial.unit = 0,  initial.m = "auto", 
			initial.trim = 1, initial.depths = 1, initial.reverse = 1, initial.style = "Tukey") 
	dialog.values <- getDialog("StatMedStemAndLeaf", defaults)
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Stem and Leaf Display"), preventCrisp=TRUE)
    xBox <- variableListBox(top, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable (pick one)"), listHeight=10, initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    displayDigits <- tclVar("1")
    onDigits <- function(...){
        tclvalue(displayDigits) <- formatC(10^as.numeric(tclvalue(leafsDigitValue)),
            format="fg", big.mark=",")
        tclvalue(leafsAutoVariable) <- "0"
        }
	radioButtons(name = "parts", buttons = c("auto", "one", "two", 
					"five"), values = c("auto", "1", "2", "5"), labels = c(gettext(domain="R-RcmdrPlugin.EZR","Automatic"), 
					"   1", "   2", "   5"), title = gettext(domain="R-RcmdrPlugin.EZR","Parts Per Stem"), 
			initialValue = dialog.values$initial.m)
	radioButtons(name = "style", buttons = c("Tukey", "bare"), 
			labels = gettext(domain="R-RcmdrPlugin.EZR",c("Tukey", "Repeated stem digits")), 
			title = gettext(domain="R-RcmdrPlugin.EZR","Style of Divided Stems"), 
			initialValue = dialog.values$initial.style)
	checkBoxes(frame = "optionsFrame", boxes = c("trimOutliers", 
					"showDepths", "reverseNegative"), initialValues = c(dialog.values$initial.trim,
					dialog.values$initial.depths, dialog.values$initial.reverse),
			labels = gettext(domain="R-RcmdrPlugin.EZR",c("Trim outliers", "Show depths", 
							"Reverse negative leaves")))
#		radioButtons(name="parts", buttons=c("auto", "one", "two", "five"),
#        values=c("auto", "1", "2", "5"), labels=c(gettext(domain="R-RcmdrPlugin.EZR","Automatic"), "   1", "   2", "   5"),
#        title=gettext(domain="R-RcmdrPlugin.EZR","Parts Per Stem"))
#    radioButtons(name="style", buttons=c("Tukey", "bare"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Tukey", "Repeated stem digits")),
#        title=gettext(domain="R-RcmdrPlugin.EZR","Style of Divided Stems"))
#    checkBoxes(frame="optionsFrame", boxes=c("trimOutliers", "showDepths", "reverseNegative"),
#        initialValues=rep(1, 3), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Trim outliers", "Show depths", "Reverse negative leaves")))
	leafsFrame <- tkframe(top)
	leafsDigitValue <- tclVar(dialog.values$initial.unit) #tclVar("0")
	leafsDigitSlider <- tkscale(leafsFrame, from = -6, to = 6, 
			showvalue = FALSE, variable = leafsDigitValue, resolution = 1, 
			orient = "horizontal", command = onDigits)
	leafsDigitShow <- labelRcmdr(leafsFrame, textvariable = displayDigits, 
			width = 8, justify = "right")
	leafsAutoVariable <- tclVar("1") # tclVar(dialog.values$initial.leafs.auto)
	leafsDigitCheckBox <- tkcheckbutton(leafsFrame, variable = leafsAutoVariable)
#    leafsFrame <- tkframe(top)
#    leafsDigitValue <- tclVar("0")
#    leafsDigitSlider <- tkscale(leafsFrame, from=-6, to=6, showvalue=FALSE, variable=leafsDigitValue,
#        resolution=1, orient="horizontal", command=onDigits)
#    leafsDigitShow <- labelRcmdr(leafsFrame, textvariable=displayDigits, width=8, justify="right")
#    leafsAutoVariable <- tclVar("1")
#    leafsDigitCheckBox <- tkcheckbutton(leafsFrame, variable=leafsAutoVariable)
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Stem and Leaf Display"), "#####", sep=""))
		x <- getSelection(xBox)
		m <- tclvalue(partsVariable)
		style <- tclvalue (styleVariable)
		trim <- tclvalue (trimOutliersVariable)
		depths <- tclvalue (showDepthsVariable)
		reverse <- tclvalue (reverseNegativeVariable)
		unit <- if (tclvalue(leafsAutoVariable) == "1") 
					""
				else paste(", unit=", 10^as.numeric(tclvalue(leafsDigitValue)), 
							sep = "")
		putDialog ("StatMedStemAndLeaf", list(initial.x = x, initial.leafs.auto=tclvalue(leafsAutoVariable),
						initial.unit = as.numeric(tclvalue(leafsDigitValue)),  initial.m = m, 
						initial.trim = trim, initial.depths = depths, initial.reverse = reverse, 
						initial.style = style))
		closeDialog()
		if (length(x) == 0) {
			errorCondition(recall = StatMedStemAndLeaf, message = gettext(domain="R-RcmdrPlugin.EZR","You must select a variable"))
			return()
		}
		trim <- if (tclvalue(trimOutliersVariable) == "1") 
					""
				else ", trim.outliers=FALSE"
		depths <- if (tclvalue(showDepthsVariable) == "1") 
					""
				else ", depths=FALSE"
		reverse <- if (tclvalue(reverseNegativeVariable) == "1") 
					""
				else ", reverse.negative.leaves=FALSE"
		m <- if (tclvalue(partsVariable) == "auto") 
					""
				else paste(", m=", tclvalue(partsVariable), sep = "")
		style <- if (tclvalue(styleVariable) == "Tukey") 
					""
				else ", style=\"bare\""
		doItAndPrint("library(aplpack)")
		command <- paste("stem.leaf(", ActiveDataSet(), "$", 
				x, style, unit, m, trim, depths, reverse, ", na.rm=TRUE)", 
				sep = "")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
        }
	OKCancelHelp(helpSubject = "stem.leaf", apply = "StatMedStemAndLeaf", reset = "StatMedStemAndLeaf")
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(leafsFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Leafs Digit:  "), 
					fg = "blue"), labelRcmdr(leafsFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Automatic")), 
			leafsDigitCheckBox, labelRcmdr(leafsFrame, text = gettext(domain="R-RcmdrPlugin.EZR","  or set:"), 
					fg = "red"), leafsDigitShow, leafsDigitSlider, sticky = "w")
#    tkgrid(labelRcmdr(leafsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Leafs Digit:  "), fg="blue"),
#        labelRcmdr(leafsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Automatic")), leafsDigitCheckBox,
#        labelRcmdr(leafsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","  or set:"), fg="red"), leafsDigitShow, leafsDigitSlider, sticky="w")
	tkgrid(leafsFrame, sticky = "w")
	tkgrid(partsFrame, sticky = "w")
	tkgrid(styleFrame, sticky = "w")
	tkgrid(labelRcmdr(top, text = gettext(domain="R-RcmdrPlugin.EZR","Options"), fg = "blue"), 
			sticky = "w")
#    tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Options"), fg="blue"), sticky="w")
	tkgrid(optionsFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "w")
	tclvalue(leafsAutoVariable) <- dialog.values$initial.leafs.auto
#    tclvalue(leafsAutoVariable) <- "1"
	dialogSuffix(rows = 7, columns = 1, preventCrisp = TRUE)
}

	
StatMedBoxPlot <- function(){
defaults <- list(x=NULL, group=NULL, logy=0, whisker="90", subset = "")
dialog.values <- getDialog("StatMedBoxPlot", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Boxplot"))
	variablesFrame <- tkframe(top)
    xBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$x, "numeric"))
	checkBoxes(frame="logy", boxes=c("logy"),initialValues=dialog.values$logy,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Log y-axis")))
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
#    identifyVariable <- tclVar("0")
#    identifyFrame <- tkframe(top)
#    identifyCheckBox <- tkcheckbutton(identifyFrame, variable=identifyVariable)
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Boxplot"), "#####", sep=""))
        x <- getSelection(xBox)
		group <- getSelection(groupBox)
		logy <- tclvalue(logyVariable)
        whisker <- tclvalue(whiskerVariable)
		subset <- tclvalue(subsetVariable)
putDialog("StatMedBoxPlot", list(x=x, group=group, logy=logy, whisker=whisker, subset = tclvalue(subsetVariable)))
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}
        closeDialog()
        if (length(x) == 0){
            errorCondition(recall=StatMedBoxPlot, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable"))
            return()
            }
		if (logy==0){
			logy <- ""
		}
		else{
			logy <- ', log="y"'
		}
#        identifyPoints <- "1" == tclvalue(identifyVariable)
        .activeDataSet <- ActiveDataSet()
        var <- paste(subset1, .activeDataSet, subset2, "[complete.cases(", subset1, .activeDataSet, subset2, "$", x, "),]$", x, sep="")
        compgroup <- paste(subset1, .activeDataSet, subset2, "[complete.cases(", subset1, .activeDataSet, subset2, "$", x, "),]$", group, sep="")
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
        if (length(group) == 0) {
			if(whisker=="default"){
				command <- (paste("boxplot(", var, ', ylab="', x, '"', logy, ')', sep=""))
				doItAndPrint(command)
            } else {
				command <- (paste("boxdata <- boxplot(", var, ', ylab="', x, '"', logy, ', plot=FALSE)', sep=""))
				doItAndPrint(command)
				if(whisker=="90"){
					doItAndPrint(paste("boxdata$stats[1,] <- quantile(", var, ", .1, na.rm=TRUE)", sep=""))
					doItAndPrint(paste("boxdata$stats[5,] <- quantile(", var, ", .9, na.rm=TRUE)", sep=""))
					doItAndPrint(paste("boxdata.outliers <- ", var, "[", var , "<boxdata$stats[1,] | ", var, ">boxdata$stats[5,]]", sep=""))
					doItAndPrint("boxdata$out <- c(boxdata$out, boxdata.outliers)")
					doItAndPrint(paste("boxdata$group <- c(boxdata$group, rep(1, length(boxdata.outliers)))", sep=""))
				}
				if(whisker=="95"){
					doItAndPrint(paste("boxdata$stats[1,] <- quantile(", var, ", .05, na.rm=TRUE)", sep=""))
					doItAndPrint(paste("boxdata$stats[5,] <- quantile(", var, ", .95, na.rm=TRUE)", sep=""))
					doItAndPrint(paste("boxdata.outliers <- ", var, "[", var , "<boxdata$stats[1,] | ", var, ">boxdata$stats[5,]]", sep=""))
					doItAndPrint("boxdata$out <- c(boxdata$out, boxdata.outliers)")
					doItAndPrint(paste("boxdata$group <- c(boxdata$group, rep(1, length(boxdata.outliers)))", sep=""))
				}
				if(whisker=="maxmin"){
					doItAndPrint(paste("boxdata$stats[1,] <- min(", var, ", na.rm=TRUE)", sep=""))
					doItAndPrint(paste("boxdata$stats[5,] <- max(", var, ", na.rm=TRUE)", sep=""))
					doItAndPrint("boxdata$out <- NULL")
					doItAndPrint("boxdata$group <- NULL")
				}
				doItAndPrint(paste('bxp(boxdata, ylab="', x, '"', logy, ")", sep=""))
				doItAndPrint("remove(boxdata)")				
				doItAndPrint("remove(boxdata.outliers)")
			}
        } else {
			if(whisker=="default"){
				command <- (paste("boxplot(", x, "~ factor(", group, '), ylab="', x,
                '", xlab="', group,'"',
                ", data=", .activeDataSet, subset, logy, ")", sep=""))
				doItAndPrint(command)
			} else {
				command <- (paste("boxdata <- boxplot(", x, "~ factor(", group, '), ylab="', x,
                '", xlab="', group,'"',
                ", data=", .activeDataSet, subset, logy, ", plot=FALSE)", sep=""))
				doItAndPrint(command)
				groups <- eval(parse(text=paste("levels(factor(", subset1, .activeDataSet, subset2, "$", group, "))", sep="")))
				ngroup <- length(groups)
				doItAndPrint("boxdata$out <- NULL")
				doItAndPrint("boxdata$group <- NULL")
				for (i in 1:ngroup){
					if(whisker=="90"){
						doItAndPrint(paste("boxdata$stats[1,", i, "] <- quantile(", var, "[", compgroup, '=="', groups[i], '"], .1, na.rm=TRUE)', sep=""))
						doItAndPrint(paste("boxdata$stats[5,", i, "] <- quantile(", var, "[", compgroup, '=="', groups[i], '"], .9, na.rm=TRUE)', sep=""))
						doItAndPrint(paste("boxdata.outliers <- ", subset1, .activeDataSet, subset2, "[!is.na(", subset1, .activeDataSet, subset2, "$", x, ") & ", subset1, .activeDataSet, subset2, "$", group, '=="', groups[i], '",]$', x, "[", subset1, .activeDataSet, subset2, "[!is.na(", subset1, .activeDataSet, subset2, "$", x, ") & ", subset1, .activeDataSet, subset2, "$", group, '=="', groups[i], '",]$', x , "<boxdata$stats[1,", i, "] | ", subset1, .activeDataSet, subset2, "[!is.na(", subset1, .activeDataSet, subset2, "$", x, ") & ", subset1, .activeDataSet, subset2, "$", group, '=="', groups[i], '",]$', x , ">boxdata$stats[5,", i, "]]", sep=""))
						doItAndPrint("boxdata$out <- c(boxdata$out, boxdata.outliers)")
						doItAndPrint(paste("boxdata$group <- c(boxdata$group, rep(", i, ", length(boxdata.outliers)))", sep=""))
						}
					if(whisker=="95"){
						doItAndPrint(paste("boxdata$stats[1,", i, "] <- quantile(", var, "[", compgroup, '=="', groups[i], '"], .05, na.rm=TRUE)', sep=""))
						doItAndPrint(paste("boxdata$stats[5,", i, "] <- quantile(", var, "[", compgroup, '=="', groups[i], '"], .95, na.rm=TRUE)', sep=""))
						doItAndPrint(paste("boxdata.outliers <- ", subset1, .activeDataSet, subset2, "[!is.na(", subset1, .activeDataSet, subset2, "$", x, ") & ", subset1, .activeDataSet, subset2, "$", group, '=="', groups[i], '",]$', x, "[", subset1, .activeDataSet, subset2, "[!is.na(", subset1, .activeDataSet, subset2, "$", x, ") & ", subset1, .activeDataSet, subset2, "$", group, '=="', groups[i], '",]$', x , "<boxdata$stats[1,", i, "] | ", subset1, .activeDataSet, subset2, "[!is.na(", subset1, .activeDataSet, subset2, "$", x, ") & ", subset1, .activeDataSet, subset2, "$", group, '=="', groups[i], '",]$', x , ">boxdata$stats[5,", i, "]]", sep=""))
						doItAndPrint("boxdata$out <- c(boxdata$out, boxdata.outliers)")
						doItAndPrint(paste("boxdata$group <- c(boxdata$group, rep(", i, ", length(boxdata.outliers)))", sep=""))
					}
					if(whisker=="maxmin"){
						doItAndPrint(paste("boxdata$stats[1,", i, "] <- min(", var, "[", compgroup, '=="', groups[i], '"], na.rm=TRUE)', sep=""))
						doItAndPrint(paste("boxdata$stats[5,", i, "] <- max(", var, "[", compgroup, '=="', groups[i], '"], na.rm=TRUE)', sep=""))
						doItAndPrint("boxdata$out <- NULL")
						doItAndPrint("boxdata$group <- NULL")
					}
				}
				doItAndPrint(paste('bxp(boxdata, ylab="', x, '"', logy, ")", sep=""))
				doItAndPrint("remove(boxdata)")				
				doItAndPrint("remove(boxdata.outliers)")
			}
		}
#		if (identifyPoints) {
#			RcmdrTkmessageBox(title="Identify Points",
#					message=paste(gettext(domain="R-RcmdrPlugin.EZR","Use left mouse button to identify points,\n"),
#					gettext(domain="R-RcmdrPlugin.EZR",if (MacOSXP()) "esc key to exit." else "right button to exit."), sep=""),
#                    icon="info", type="ok")
#               doItAndPrint(paste("identify(rep(1, length(", var,
#                    ")), ", var, ", rownames(", .activeDataSet,"))", sep=""))
#       }
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="boxplot", apply="StatMedBoxPlot", reset="StatMedBoxPlot")
    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","    ")), getFrame(groupBox), sticky="w")
    tkgrid(variablesFrame, stick="w")
	tkgrid(logy, sticky="w")    
	radioButtons(name="whisker", buttons=c("ninety", "ninetyfive", "maxmin", "default"), values=c("90", "95", "maxmin", "default"), initialValue=dialog.values$whisker, labels=gettext(domain="R-RcmdrPlugin.EZR",c("10-90 percentiles", "5-95 percentiles", "Min-Max", "(1Q-1.5xIQR)-(3Q+1.5xIQR)")), 
		title=gettext(domain="R-RcmdrPlugin.EZR","Whisker range"))
    tkgrid(whiskerFrame, sticky="nw")	

#    tkgrid(labelRcmdr(identifyFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Identify outliers with mouse"), justify="left"),
#       identifyCheckBox, sticky="w")
#	tkgrid(identifyFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
    }
    
	
StatMedBarMeans <- function(){
defaults <- list(group1=NULL, group2=NULL, response=NULL, errorBars="bar.sds", subset="")
dialog.values <- getDialog("StatMedBarMeans", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE 
   initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Bar graph(Means)"))
    variablesFrame <- tkframe(top)
    group1Box <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable1(pick 0 or 1)"), listHeight=15, initialSelection=varPosn(dialog.values$group1, "all"))
    group2Box <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable2(pick 0 or 1)"), listHeight=15, initialSelection=varPosn(dialog.values$group2, "all"))
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Bar graph(Means)"), "#####", sep=""))
        group1 <- getSelection(group1Box)
        group2 <- getSelection(group2Box)
        response <- getSelection(responseBox)
        error.bars <- tclvalue(errorBarsVariable)
		subset <- tclvalue(subsetVariable)
putDialog("StatMedBarMeans", list(group1=group1, group2=group2, response=response, errorBars=error.bars, subset=tclvalue(subsetVariable)))
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
        closeDialog()
        if (length(response) == 0) {
            errorCondition(recall=StatMedBarMeans, message=gettext(domain="R-RcmdrPlugin.EZR","No response variable selected."))
            return()
          }
        dataSet <- ActiveDataSet()
		if (length(group1) == 0){
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint(paste("bar.sums <- sum(", subset1, dataSet, subset2, "$", response, ", na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.means <- mean(", subset1, dataSet, subset2, "$", response, ", na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.sds <- mean(", subset1, dataSet, subset2, "$", response, ", na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.n <- bar.sums/bar.means"))
			doItAndPrint(paste("bar.ses <- bar.sds/sqrt(bar.n)"))
			doItAndPrint("bar.sds <- ifelse(is.na(bar.sds), 0, bar.sds)")
			doItAndPrint("bar.ses <- ifelse(is.na(bar.ses), 0, bar.ses)")			

			if (error.bars == "none"){
				doItAndPrint(paste('barx <- barplot(bar.means, ylim=c(ifelse(min(bar.means)>0, 0, min(bar.means)*1.2), max(bar.means)*1.2), ylab="', response, '", axis.lty=1)',sep=""))
			}
			else{
				doItAndPrint(paste('barx <- barplot(bar.means, ylim=c(ifelse(min(bar.means)>0, 0, min(bar.means-', error.bars, ')*1.2), max(bar.means+', error.bars, ')*1.2), ylab="', response, '", axis.lty=1)',sep=""))
				doItAndPrint(paste("error.bar(barx, bar.means, ", error.bars, ")", sep=""))
			}
		}
		if (length(group1) == 1 && length(group2) == 0){
		    if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint(paste("bar.sums <- tapply(", subset1, dataSet, subset2, "$", response, ", factor(", subset1, dataSet, subset2, "$", group1, "), sum, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.means <- tapply(", subset1, dataSet, subset2, "$", response, ", factor(", subset1, dataSet, subset2, "$", group1, "), mean, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.sds <- tapply(", subset1, dataSet, subset2, "$", response, ", factor(", subset1, dataSet, subset2, "$", group1, "), sd, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.n <- bar.sums/bar.means"))
			doItAndPrint(paste("bar.ses <- bar.sds/sqrt(bar.n)"))
			doItAndPrint("bar.sds <- ifelse(is.na(bar.sds), 0, bar.sds)")
			doItAndPrint("bar.ses <- ifelse(is.na(bar.ses), 0, bar.ses)")			
			if (error.bars == "none"){
				doItAndPrint(paste('barx <- barplot(bar.means, ylim=c(ifelse(min(bar.means, na.rm=TRUE)>0, 0, min(bar.means, na.rm=TRUE)*1.2), max(bar.means, na.rm=TRUE)*1.2), xlab="', group1, '", ylab="', response, '", axis.lty=1)',sep=""))
			}
			else{
				doItAndPrint(paste('barx <- barplot(bar.means, ylim=c(ifelse(min(bar.means, na.rm=TRUE)>0, 0, min(bar.means-', error.bars, ', na.rm=TRUE)*1.2), max(bar.means+', error.bars, ', na.rm=TRUE)*1.2), xlab="', group1, '", ylab="', response, '", axis.lty=1)',sep=""))
				doItAndPrint(paste("error.bar(barx, bar.means, ", error.bars, ")", sep=""))
			}
		}
		if (length(group1) == 1 && length(group2) == 1){
			if (eval(parse(text=paste("min(table(", subset1, dataSet, subset2, "$", group1, ", ", subset1, dataSet, subset2, "$", group2, "))", sep="")))==0) {
				logger(gettext(domain="R-RcmdrPlugin.EZR","Graph not created when a group with 0 sample exists"))
			} else {			
			eval.bar.var <- eval(parse(text=paste("length(levels(factor(", subset1, dataSet, subset2, "$", group2, ")))", sep="")))
		    if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint(paste("bar.var <- length(levels(factor(", subset1, dataSet, subset2, "$", group2, ")))", sep=""))
			doItAndPrint(paste("bar.sums <- tapply(subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[1])$", response, ", subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[1])$", group1, ", sum, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.means <- tapply(subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[1])$", response, ", subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[1])$", group1, ", mean, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.sds <- tapply(subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[1])$", response, ", subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[1])$", group1, ", sd, na.rm=TRUE)", sep=""))
			
			if(eval.bar.var > 1){
				for (i in 2: eval.bar.var){
				doItAndPrint(paste("bar.sums <- c(bar.sums, tapply(subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[", i, "])$", response, ", subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[", i, "])$", group1, ", sum, na.rm=TRUE))", sep=""))
				doItAndPrint(paste("bar.means <- c(bar.means, tapply(subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[", i, "])$", response, ", subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[", i, "])$", group1, ", mean, na.rm=TRUE))", sep=""))
				doItAndPrint(paste("bar.sds <- c(bar.sds, tapply(subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[", i, "])$", response, ", subset(", subset1, dataSet, subset2, ", ", group2, "==levels(factor(", group2, "))[", i, "])$", group1, ", sd, na.rm=TRUE))", sep=""))
				doItAndPrint("bar.n <- bar.sums/bar.means")
				doItAndPrint("bar.ses <- bar.sds/sqrt(bar.n)")
				}
			}
			doItAndPrint(paste("bar.var2 <- length(levels(factor(", subset1, dataSet, subset2, "$", group1, ")))", sep=""))
			doItAndPrint("bar.means <- matrix(bar.means, bar.var2)")
			doItAndPrint("bar.sds <- matrix(bar.sds, bar.var2)")
			doItAndPrint("bar.ses <- matrix(bar.ses, bar.var2)")
			doItAndPrint("bar.sds <- ifelse(is.na(bar.sds), 0, bar.sds)")
			doItAndPrint("bar.ses <- ifelse(is.na(bar.ses), 0, bar.ses)")
			if (error.bars == "none"){
				doItAndPrint(paste('barx <- barplot(bar.means, beside=TRUE, ylim=c(ifelse(min(bar.means)>0, 0, min(bar.means)*1.2), max(bar.means)*1.2), xlab="', group2, '", ylab="', response, '", names.arg=levels(factor(', subset1, dataSet, subset2, "$", group2, ")), legend.text=levels(factor(", subset1, dataSet, subset2, "$", group1, ')), args.legend=list(title="', group1, '", box.lty=0), axis.lty=1)',  sep=""))
			}
			else{
				doItAndPrint(paste('barx <- barplot(bar.means, beside=TRUE, ylim=c(ifelse(min(bar.means)>0, 0, min(bar.means-', error.bars, ')*1.2), max(bar.means+', error.bars, ')*1.2), xlab="', group2, '", ylab="', response, '", names.arg=levels(factor(', subset1, dataSet, subset2, "$", group2, ")), legend.text=levels(factor(", subset1, dataSet, subset2, "$", group1, ')), args.legend=list(title="', group1, '", box.lty=0), axis.lty=1)',  sep=""))
				doItAndPrint(paste("error.bar(barx, bar.means, ", error.bars, ")",  sep=""))
			}
		}
		}
        activateMenus()
        tkfocus(CommanderWindow())
        }
    optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="errorBars", buttons=c("bar.ses", "bar.sds", "none"), values=c("bar.ses", "bar.sds", "none"),initialValue=dialog.values$errorBars, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Standard errors", "Standard deviations", "No error bars")),
        title=gettext(domain="R-RcmdrPlugin.EZR","Error Bars"))
#	errorBarsVariable <- tclVar("bar.sds")
#    seButton <- ttkradiobutton(optionsFrame, variable=errorBarsVariable, value="bar.ses")
#    sdButton <- ttkradiobutton(optionsFrame, variable=errorBarsVariable, value="bar.sds")
#    noneButton <- ttkradiobutton(optionsFrame, variable=errorBarsVariable, value="none")
    buttonsFrame <- tkframe(top)
    OKCancelHelp(helpSubject="barplot", apply="StatMedBarMeans", reset="StatMedBarMeans")
	tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(group1Box), getFrame(group2Box), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
#    tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Error Bars"), fg="blue"), sticky="w")
#    tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Standard errors")), seButton, sticky="w")
#    tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviations")), sdButton, sticky="w")
#    tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","No error bars")), noneButton, sticky="w")	
    tkgrid(errorBarsFrame, columnspan=2, sticky="w")
    tkgrid(optionsFrame, columnspan=2, sticky="w")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }
	
	
StatMedStripChart <- function(){
defaults <- list(group=NULL, response=NULL, logy=0, subset = "")
dialog.values <- getDialog("StatMedStripChart", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Strip Chart"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), listHeight=15, 
	initialSelection=varPosn(dialog.values$group, "all"))
	responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
	checkBoxes(frame="logy", boxes=c("logy"),initialValues=dialog.values$logy,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Log y-axis")))
    StatMedSubsetBox(model=TRUE)  
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Strip Chart"), "#####", sep=""))
		groups <- getSelection(groupBox)
		response <- getSelection(responseBox)
		logy <- tclvalue(logyVariable)
		.activeDataSet <- ActiveDataSet()
		subset <- tclvalue(subsetVariable)
putDialog("StatMedStripChart", list(group=groups, response=response, logy=logy, subset = tclvalue(subsetVariable)))
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				.subDataSet <- .activeDataSet
			} else {
				.subDataSet <- paste("subset(", .activeDataSet, ", ", subset, ")", sep="")
			}
		closeDialog()
		if (0 == length(response)) {
			errorCondition(recall=StatMedStripChart, message=gettext(domain="R-RcmdrPlugin.EZR","No response variable selected."))
			return()
		}
		if (logy==0){
			logy <- ""         
			logflag <- ""
      	   	}   
		else{
			logy <- ', log="y"'
			logflag <- ", log.flag=TRUE"
		}
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		if (length(groups) == 0){
			doItAndPrint(paste("dummyX <- rep(0, length(", .subDataSet, "$", response, "))"))		
			doItAndPrint(paste("dot.plot(dummyX, ", .subDataSet, "$", response, logflag, ', xlab="", ylab="', response, '")', sep=""))			
		} else {
			doItAndPrint(paste("dot.plot(", .subDataSet, "$", groups, ", ", .subDataSet, "$", response, logflag, ', xlab="', groups, '", ylab="', response, '")', sep=""))
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	buttonsFrame <- tkframe(top)
	OKCancelHelp(helpSubject="plot", apply="StatMedStripChart", reset="StatMedStripChart")
	tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(logy, sticky="w")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, columnspan=2, sticky="w")
	dialogSuffix(rows=3, columns=2)
}
	
	
StatMedOrderedChart <- function(){
defaults <- list(response=NULL, group=NULL, type="line", trend="FALSE", lowlim="<auto>", uplim="<auto>", logy=0, subset="")
dialog.values <- getDialog("StatMedOrderedChart", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Ordered Chart"))
    variablesFrame <- tkframe(top)
	groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Factors (pick zero or more)"), selectmode="multiple", listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
	responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
	optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="type", buttons=c("line", "box"), values=c("line", "box"), initialValue=dialog.values$type, 
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Line", "Box")), title=gettext(domain="R-RcmdrPlugin.EZR","Plot type"))
    radioButtons(optionsFrame, name="trend", buttons=c("inc", "dec"), values=c("FALSE", "TRUE"), initialValue=dialog.values$trend,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Increasing", "Decreasing")), title=gettext(domain="R-RcmdrPlugin.EZR","Order"))
	options2Frame <- tkframe(top)
	checkBoxes(frame="options2Frame", boxes=c("logy"),initialValues=dialog.values$logy,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Log y-axis")))
	options3Frame <- tkframe(top)
	lowlimFrame <- tkframe(options3Frame)
	lowlim <- tclVar(dialog.values$lowlim)
	lowlimField <- ttkentry(lowlimFrame, width="20", textvariable=lowlim)
	tkgrid(tklabel(lowlimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis lower limit")), lowlimField, sticky="w")
#	tkgrid(lowlimFrame, sticky="w")
	uplimFrame <- tkframe(options3Frame)
	uplim <- tclVar(dialog.values$uplim)
	uplimField <- ttkentry(uplimFrame, width="20", textvariable=uplim)
	tkgrid(tklabel(uplimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis upper limit")), uplimField, sticky="w")
#	tkgrid(uplimFrame, sticky="w")
	tkgrid(lowlimFrame, labelRcmdr(options3Frame, text="  "), uplimFrame, sticky="w")
     StatMedSubsetBox(model=TRUE)  
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Ordered Chart"), "#####", sep=""))
		groups <- getSelection(groupBox)
		response <- getSelection(responseBox)
		ylog <- tclvalue(logyVariable)
		type <- as.character(tclvalue(typeVariable))
		trend <- as.character(tclvalue(trendVariable))		
		lowlim <- tclvalue(lowlim)
		uplim <- tclvalue(uplim)
		.activeDataSet <- ActiveDataSet()
		subset <- tclvalue(subsetVariable)
putDialog("StatMedOrderedChart", list(response=response, group=getSelection(groupBox), type=type, trend=trend, lowlim=lowlim, uplim=uplim, logy=ylog, subset=tclvalue(subsetVariable)))
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				.subDataSet <- .activeDataSet
			} else {
				.subDataSet <- paste("subset(", .activeDataSet, ", ", subset, ")", sep="")
			}
		if (length(groups)==0) {
			groups <- "NULL"
		} else {
			groups <- paste(.subDataSet, "$", groups, sep="")
		}
		closeDialog()
		if (0 == length(response)) {
			errorCondition(recall=StatMedOrderedChart, message=gettext(domain="R-RcmdrPlugin.EZR","No response variable selected."))
			return()
		}
		if (lowlim=="<auto>") lowlim <- NULL
		if (uplim=="<auto>") uplim <- NULL
		
		if (ylog==0){
			ylog <- FALSE         
      	   	} else {
			ylog <- TRUE
		}
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		
		command <- paste("OrderedPlot(y=", .subDataSet, "$", response, ", group=", groups, ', type="', type, '", ylab="', response, '", ylog=', ylog, ", lowlim=", lowlim, ", uplim=", uplim, ', decreasing="', trend, '")', sep="") 
		doItAndPrint(command)
				activateMenus()
		tkfocus(CommanderWindow())
	}
	buttonsFrame <- tkframe(top)
	OKCancelHelp(helpSubject="plot", apply="StatMedOrderedChart", reset="StatMedOrderedChart")
	tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Grouping is valid only when line plot is selected."), fg="blue"), sticky="w")
	tkgrid(typeFrame, trendFrame, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(options2Frame, sticky="w")
    tkgrid(options3Frame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, columnspan=2, sticky="w")
	dialogSuffix(rows=3, columns=2)
}
	
	
StatMedScatterPlot <- function () {
#	require("car")
	defaults <- list(initial.x = NULL, initial.y = NULL, initial.jitterx = 0, initial.jittery = 0, 
			initial.logstringx = 0, initial.logstringy = 0, initial.log = 0, initial.box = 1, 
			initial.line = 1, initial.smooth = 0, initial.spread = 0, initial.span = 50,
			initial.subset = gettext ("<all valid cases>"), initial.ylab = gettext ("<auto>"), 
			initial.xlab = gettext(domain="R-RcmdrPlugin.EZR","<auto>"), initial.pch = gettext(domain="R-RcmdrPlugin.EZR","<auto>"), 
			initial.cexValue = 1, initial.cex.axisValue = 1, initial.cex.labValue = 1, initialGroup=NULL, initial.lines.by.group=1, subset="") 
	dialog.values <- getDialog("StatMedScatterPlot", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	initial.group <- dialog.values$initial.group
	.linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
	.groups <- if (is.null(initial.group)) FALSE else initial.group
	Library("tcltk")
	initializeDialog(title = gettext(domain="R-RcmdrPlugin.EZR","Scatterplot"))
	.numeric <- Numeric()
	variablesFrame <- tkframe(top)
	xBox <- variableListBox(variablesFrame, .numeric, title = gettext(domain="R-RcmdrPlugin.EZR","x-variable (pick one)"), listHeight=6,
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	yBox <- variableListBox(variablesFrame, .numeric, title = gettext(domain="R-RcmdrPlugin.EZR","y-variable (pick one)"), listHeight=6,
			initialSelection = varPosn (dialog.values$initial.y, "numeric"))
	optionsParFrame <- tkframe(top)
	checkBoxes(window = optionsParFrame, frame = "optionsFrame", 
			boxes = c("identify", "jitterX", "jitterY", "logX", "logY", 
					"boxplots", "lsLine", "smoothLine", "spread"), initialValues = c(dialog.values$initial.log, 
					dialog.values$initial.jitterx, dialog.values$initial.jittery, 
					dialog.values$initial.logstringx, dialog.values$initial.logstringy,
					dialog.values$initial.box, dialog.values$initial.line, dialog.values$initial.smooth,
					dialog.values$initial.spread),labels = gettext(domain="R-RcmdrPlugin.EZR",c("Identify points", 
							"Jitter x-variable", "Jitter y-variable", "Log x-axis", 
							"Log y-axis", "Marginal boxplots", "Least-squares line", 
							"Smooth line", "Show spread")), title = "Options")
	sliderValue <- tclVar(dialog.values$initial.span)
	slider <- tkscale(optionsFrame, from = 0, to = 100, showvalue = TRUE, 
			variable = sliderValue, resolution = 5, orient = "horizontal")
#	subsetBox(subset.expression = dialog.values$initial.subset)
    StatMedSubsetBox(model=TRUE)
	labelsFrame <- tkframe(top)
	xlabVar <- tclVar(dialog.values$initial.xlab)
	ylabVar <- tclVar(dialog.values$initial.ylab)
	xlabFrame <- tkframe(labelsFrame)
	xlabEntry <- ttkentry(xlabFrame, width = "25", textvariable = xlabVar)
	xlabScroll <- ttkscrollbar(xlabFrame, orient = "horizontal", 
			command = function(...) tkxview(xlabEntry, ...))
	tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll, 
						...))
	tkgrid(labelRcmdr(xlabFrame, text = gettext(domain="R-RcmdrPlugin.EZR","x-axis label"), 
					fg = "blue"), sticky = "w")
	tkgrid(xlabEntry, sticky = "w")
	tkgrid(xlabScroll, sticky = "ew")
	ylabFrame <- tkframe(labelsFrame)
	ylabEntry <- ttkentry(ylabFrame, width = "25", textvariable = ylabVar)
	ylabScroll <- ttkscrollbar(ylabFrame, orient = "horizontal", 
			command = function(...) tkxview(ylabEntry, ...))
	tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll, 
						...))
	tkgrid(labelRcmdr(ylabFrame, text = gettext(domain="R-RcmdrPlugin.EZR","y-axis label"), 
					fg = "blue"), sticky = "w")
	tkgrid(ylabEntry, sticky = "w")
	tkgrid(ylabScroll, sticky = "ew")
	tkgrid(xlabFrame, labelRcmdr(labelsFrame, text = "     "), 
			ylabFrame, sticky = "w")
	parFrame <- tkframe(optionsParFrame)
	pchVar <- tclVar(dialog.values$initial.pch)
	pchEntry <- ttkentry(parFrame, width = 25, textvariable = pchVar)
	cexValue <- tclVar(dialog.values$initial.cexValue)
	cex.axisValue <- tclVar(dialog.values$initial.cex.axisValue)
	cex.labValue <- tclVar(dialog.values$initial.cex.labValue)
	cexSlider <- tkscale(parFrame, from = 0.5, to = 2.5, showvalue = TRUE, 
			variable = cexValue, resolution = 0.1, orient = "horizontal")
	cex.axisSlider <- tkscale(parFrame, from = 0.5, to = 2.5, 
			showvalue = TRUE, variable = cex.axisValue, resolution = 0.1, 
			orient = "horizontal")
	cex.labSlider <- tkscale(parFrame, from = 0.5, to = 2.5, 
			showvalue = TRUE, variable = cex.labValue, resolution = 0.1, 
			orient = "horizontal")
	onOK <- function() {
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Scatterplot"), "#####", sep=""))
		x <- getSelection(xBox)
		y <- getSelection(yBox)
		jitter <- if ("1" == tclvalue(jitterXVariable) && "1" == 
						tclvalue(jitterYVariable)) 
					", jitter=list(x=1, y=1)"
				else if ("1" == tclvalue(jitterXVariable)) 
					", jitter=list(x=1)"
				else if ("1" == tclvalue(jitterYVariable)) 
					", jitter=list(y=1)"
				else ""
		logstring <- ""
		if ("1" == tclvalue(logXVariable)) 
			logstring <- paste(logstring, "x", sep = "")
		if ("1" == tclvalue(logYVariable)) 
			logstring <- paste(logstring, "y", sep = "")
		log <- tclvalue(identifyVariable)
		box <- tclvalue(boxplotsVariable)
		line <- tclvalue(lsLineVariable)
		smooth <-  tclvalue(smoothLineVariable)
		spread <- tclvalue(spreadVariable)
		span <- as.numeric(tclvalue(sliderValue))
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) 
					""
				else paste(", subset=", subset, sep = "")
		cex.axis <- as.numeric(tclvalue(cex.axisValue))
		cex <- as.numeric(tclvalue(cexValue))
		cex.lab <- as.numeric(tclvalue(cex.labValue))
		xlab <- trim.blanks(tclvalue(xlabVar))
		xlab <- if (xlab == gettext(domain="R-RcmdrPlugin.EZR","<auto>")) 
					""
				else paste(", xlab=\"", xlab, "\"", sep = "")
		ylab <- trim.blanks(tclvalue(ylabVar))
		ylab <- if (ylab == gettext(domain="R-RcmdrPlugin.EZR","<auto>")) 
					""
				else paste(", ylab=\"", ylab, "\"", sep = "")
		pch <- gsub(" ", ",", tclvalue(pchVar))
		putDialog ("StatMedScatterPlot", list (initial.x = x, initial.y = y, initial.jitterx = tclvalue(jitterXVariable),
						initial.jittery = tclvalue(jitterYVariable), initial.logstringx = tclvalue(logXVariable),
						initial.logstringy = tclvalue(logYVariable), initial.log = log, initial.box = box, 
						initial.line = line, initial.smooth = smooth, initial.spread = spread,
						initial.span = span, initial.subset = initial.subset, initial.xlab = tclvalue(xlabVar),
						initial.ylab = tclvalue(ylabVar), initial.cexValue = tclvalue(cexValue), 
						initial.cex.axisValue = tclvalue(cex.axisValue), initial.cex.labValue = tclvalue(cex.labValue), 
						initial.pch = pch, initial.group=if (.groups == FALSE) NULL else .groups,
						initial.lines.by.group=if (.linesByGroup) 1 else 0, subset=tclvalue(subsetVariable)))
		closeDialog()
		if ("" == pch) {
			errorCondition(recall = StatMedScatterPlot, message = gettext(domain="R-RcmdrPlugin.EZR","No plotting characters."))
			return()
		}
		pch <- if (trim.blanks(pch) == gettext(domain="R-RcmdrPlugin.EZR","<auto>")) 
					""
				else paste(", pch=c(", pch, ")", sep = "")
		if (length(x) == 0 || length(y) == 0) {
			errorCondition(recall = StatMedScatterPlot, message = gettext(domain="R-RcmdrPlugin.EZR","You must select two variables"))
			return()
		}
		if (x == y) {
			errorCondition(recall = StatMedScatterPlot, message = gettext(domain="R-RcmdrPlugin.EZR","x and y variables must be different"))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		log <- if (logstring != "") 
					paste(", log=\"", logstring, "\"", sep = "")
				else ""
		if ("1" == tclvalue(identifyVariable)) {
			RcmdrTkmessageBox(title = "Identify Points", message = paste(gettext(domain="R-RcmdrPlugin.EZR","Use left mouse button to identify points,\n"), 
							gettext(domain="R-RcmdrPlugin.EZR",if (MacOSXP()) 
												"esc key to exit."
											else "right button to exit."), sep = ""), icon = "info", 
					type = "ok")
			idtext <- ", id.method=\"identify\""
		}
		else idtext <- ""
		box <- if ("1" == tclvalue(boxplotsVariable)) 
					"'xy'"
				else "FALSE"
		line <- if ("1" == tclvalue(lsLineVariable)) 
					"lm"
				else "FALSE"
		smooth <- as.character("1" == tclvalue(smoothLineVariable))
		spread <- as.character("1" == tclvalue(spreadVariable))
		cex <- if (cex == 1) 
					""
				else paste(", cex=", cex, sep = "")
		cex.axis <- if (cex.axis == 1) 
					""
				else paste(", cex.axis=", cex.axis, sep = "")
		cex.lab <- if (cex.lab == 1) 
					""
				else paste(", cex.lab=", cex.lab, sep = "")
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		if (.groups == FALSE) {
			doItAndPrint(paste("scatterplot(", y, "~", x, log, 
							", reg.line=", line, ", smooth=", smooth, ", spread=", 
							spread, idtext, ", boxplots=", box, ", span=", 
							span/100, jitter, xlab, ylab, cex, cex.axis, 
							cex.lab, pch, ", data=", .activeDataSet, subset, 
							")", sep = ""))
		}
		else {
			doItAndPrint(paste("scatterplot(", y, "~", x, " | ", 
							.groups, log, ", reg.line=", line, ", smooth=", smooth, 
							", spread=", spread, idtext, ", boxplots=", box, 
							", span=", span/100, jitter, xlab, ylab, cex, 
							cex.axis, cex.lab, pch, ", by.groups=", .linesByGroup, 
							", data=", .activeDataSet, subset, ")", sep = ""))
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(scatterPlot, plotLinesByGroup = TRUE, initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
			initialLabel=if (is.null(initial.group)) gettext(domain="R-RcmdrPlugin.EZR","Plot by groups") else paste(gettext(domain="R-RcmdrPlugin.EZR","Plot by:"), initial.group))
	OKCancelHelp(helpSubject = "scatterplot", apply = "StatMedScatterPlot", reset = "StatMedScatterPlot")
	tkgrid(getFrame(xBox), getFrame(yBox), sticky = "nw")
	tkgrid(variablesFrame, sticky = "w")
	tkgrid(labelRcmdr(optionsFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Span for smooth")), 
			slider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Plotting Parameters"), 
					fg = "blue"), sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Plotting characters")), 
			pchEntry, stick = "w")
	tkgrid(labelRcmdr(parFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Point size")), 
			cexSlider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Axis text size")), 
			cex.axisSlider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Axis-labels text size")), 
			cex.labSlider, sticky = "w")
	tkgrid(optionsFrame, parFrame, sticky = "nw")
	tkgrid(optionsParFrame, sticky = "w")
	tkgrid(labelsFrame, sticky = "w")
	tkgrid(subsetFrame, sticky = "w")
	tkgrid(groupsFrame, sticky = "w")
	tkgrid(labelRcmdr(top, text = " "))
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 8, columns = 2)
}
		
	
StatMedScatterPlotMatrix <- function () {
#	require("car")
	defaults <- list(initial.variables = NULL, initial.line = 1, initial.smooth = 0, initial.spread = 0, 
			initial.span = 50, initial.diag = "density", initial.subset = gettext ("<all valid cases>"),
			initialGroup=NULL, initial.lines.by.group=1, subset="") 
	dialog.values <- getDialog("StatMedScatterPlotMatrix", defaults)
	currentFields$subset <- dialog.values$subset	
	currentModel <- TRUE
	initial.group <- dialog.values$initial.group
	.linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
	.groups <- if (is.null(initial.group)) FALSE else initial.group
	Library("tcltk")
	initializeDialog(title = gettext(domain="R-RcmdrPlugin.EZR","Scatterplot Matrix"))
	variablesBox <- variableListBox(top, Numeric(), title = gettext(domain="R-RcmdrPlugin.EZR","Select variables (three or more)"), 
			selectmode = "multiple", listHeight=10, initialSelection = varPosn (dialog.values$initial.variables, "numeric"))
	checkBoxes(frame = "optionsFrame", boxes = c("lsLine", "smoothLine", 
					"spread"), initialValues = c(dialog.values$initial.line, dialog.values$initial.smooth,
					dialog.values$initial.spread), labels = gettext(domain="R-RcmdrPlugin.EZR",c("Least-squares lines", 
							"Smooth lines", "Show spread")))
	sliderValue <- tclVar(dialog.values$initial.span)
	slider <- tkscale(optionsFrame, from = 0, to = 100, showvalue = TRUE, 
			variable = sliderValue, resolution = 5, orient = "horizontal")
	radioButtons(name = "diagonal", buttons = c("density", "histogram", 
					"boxplot", "oned", "qqplot", "none"), labels = gettext(domain="R-RcmdrPlugin.EZR",c("Density plots", 
							"Histograms", "Boxplots", "One-dimensional scatterplots", 
							"Normal QQ plots", "Nothing (empty)")), title = gettext(domain="R-RcmdrPlugin.EZR","On Diagonal"), 
			initialValue = dialog.values$initial.diag)
#	subsetBox(subset.expression = dialog.values$initial.subset)
    StatMedSubsetBox(model=TRUE)
	onOK <- function() {
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Scatterplot Matrix"), "#####", sep=""))
		variables <- getSelection(variablesBox)
		closeDialog()
		line <- if ("1" == tclvalue(lsLineVariable)) 
					"lm"
				else "FALSE"
		smooth <- as.character("1" == tclvalue(smoothLineVariable))
		spread <- as.character("1" == tclvalue(spreadVariable))
		span <- as.numeric(tclvalue(sliderValue))
		diag <- as.character(tclvalue(diagonalVariable))
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) ""
				else paste(", subset=", subset, sep="")
		.activeDataSet <- ActiveDataSet()
		putDialog("StatMedScatterPlotMatrix", list(initial.variables = variables, initial.line = tclvalue (lsLineVariable), 
						initial.smooth = tclvalue(smoothLineVariable),initial.spread = tclvalue (spreadVariable), 
						initial.span = span, initial.diag = diag, initial.subset = initial.subset, 
						initial.group=if (.groups == FALSE) NULL else .groups,
						initial.lines.by.group=if (.linesByGroup) 1 else 0, subset=tclvalue(subsetVariable)))
		if (length(variables) < 3) {
			errorCondition(recall = StatMedScatterPlotMatrix, message = gettext(domain="R-RcmdrPlugin.EZR","Fewer than 3 variable selected."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		if (.groups == FALSE) {
			command <- paste("scatterplotMatrix(~", paste(variables, 
							collapse = "+"), ", reg.line=", line, ", smooth=", 
					smooth, ", spread=", spread, ", span=", span/100, 
					", diagonal = '", diag, "', data=", .activeDataSet, 
					subset, ")", sep = "")
			logger(command)
			justDoIt(command)
		}
		else {
			command <- paste("scatterplotMatrix(~", paste(variables, 
							collapse = "+"), " | ", .groups, ", reg.line=", 
					line, ", smooth=", smooth, ", spread=", spread, 
					", span=", span/100, ", diagonal= '", diag, "', by.groups=", 
					.linesByGroup, ", data=", .activeDataSet, subset, 
					")", sep = "")
			logger(command)
			justDoIt(command)
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(scatterPlot, plotLinesByGroup = TRUE, initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
			initialLabel=if (is.null(initial.group)) gettext(domain="R-RcmdrPlugin.EZR","Plot by groups") else paste(gettext(domain="R-RcmdrPlugin.EZR","Plot by:"), initial.group))
	OKCancelHelp(helpSubject = "scatterplotMatrix", apply = "StatMedScatterPlotMatrix", reset = "StatMedScatterPlotMatrix")
	tkgrid(getFrame(variablesBox), sticky = "nw")
	tkgrid(labelRcmdr(optionsFrame, text = gettext(domain="R-RcmdrPlugin.EZR","Span for smooth")), 
			slider, sticky = "w")
	tkgrid(optionsFrame, sticky = "w")
	tkgrid(diagonalFrame, sticky = "w")
	tkgrid(subsetFrame, sticky = "w")
	tkgrid(groupsFrame, sticky = "w")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 6, columns = 2)
}


StatMedPlotMeans <- function(){
defaults <- list(group=NULL, response=NULL, errorBars="sd", confidence="0.95", graph="narrow", line="color", subset = "")
dialog.values <- getDialog("StatMedPlotMeans", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Plot Means"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Factors (pick one or two)"), selectmode="multiple", listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Plot Means"), "#####", sep=""))
        groups <- getSelection(groupBox)
        response <- getSelection(responseBox)
		graph <- as.character(tclvalue(graphVariable))
        error.bars <- tclvalue(errorBarsVariable)
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
		line <- tclvalue(lineVariable)
		if (line=="color") line <- ", lty=1, lwd=1, "
		if (line=="type") line <- ", col=1, lwd=1, "
		if (line=="width") line <- ", col=1, lty=1, "
		
putDialog("StatMedPlotMeans", list(group=groups, response=response, errorBars=error.bars, confidence=as.character(tclvalue(levelVariable)), graph=graph, line=tclvalue(lineVariable), subset = tclvalue(subsetVariable)))

        closeDialog()
        if (0 == length(groups)) {
            errorCondition(recall=StatMedPlotMeans, message=gettext(domain="R-RcmdrPlugin.EZR","No factors selected."))
            return()
            }
        if (2 < length(groups)) {
            errorCondition(recall=StatMedPlotMeans, message=gettext(domain="R-RcmdrPlugin.EZR","More than two factors selected."))
            return()
            }
        if (0 == length(response)) {
            errorCondition(recall=StatMedPlotMeans, message=gettext(domain="R-RcmdrPlugin.EZR","No response variable selected."))
            return()
            }
        .activeDataSet <- ActiveDataSet()
        level <- if (error.bars == "conf.int") paste(", level=", tclvalue(levelVariable), sep="") else ""
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		if (length(groups) == 1) doItAndPrint(paste("StatMedplotMeans(", subset1, .activeDataSet, subset2, "$", response,
            ", factor(", subset1, .activeDataSet, subset2, "$", groups[1],
            '), error.bars="', error.bars, '"', level, ', xlab="', groups[1], '", ylab="', response, '")', sep=""))
        else {
			if (graph == "narrow"){
            if (eval(parse(text=paste("length(levels(", subset1, .activeDataSet, subset2, "$", groups[1],
                ")) < length(levels(", subset1, .activeDataSet, subset2, "$", groups[2], "))", sep=""))))
                groups <- rev(groups)
            doItAndPrint(paste("StatMedplotMeans(", subset1, .activeDataSet, subset2, "$", response, ", as.factor(", subset1, .activeDataSet, subset2, "$", groups[1],
                "), as.factor(", subset1, .activeDataSet, subset2, "$", groups[2], '), error.bars="', error.bars, '"', level, ', xlab="', groups[1], '"', line, 'ylab="', response, '", legend.lab="', groups[2], '")', sep=""))
            }
			else{
				doItAndPrint(paste("dummyfactor <- paste(as.factor(", subset1, .activeDataSet, subset2, "$", groups[1], '), " : ", as.factor(', .activeDataSet, "$", groups[2], '), sep="")', sep=""))
				doItAndPrint(paste('xlab <- paste("', groups[1], '", " : ", "', groups[2], '", sep="")'))
				doItAndPrint(paste("StatMedplotMeans(", subset1, .activeDataSet, subset2, "$", response,
				', as.factor(dummyfactor), error.bars="', error.bars, '"', line, level, ', xlab=xlab, ylab="', response, '")', sep=""))
			}
		}
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="plotMeans", apply="StatMedPlotMeans", reset="StatMedPlotMeans")

	tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")

    optionsFrame <- tkframe(top)
	radioButtons(optionsFrame, name="errorBars", buttons=c("se", "sd", "conf.int", "none"), values=c("se", "sd", "conf.int", "none"), initialValue=dialog.values$errorBars, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Standard errors", "Standard deviations", "Confidence intervals", "No error bars")), title=gettext(domain="R-RcmdrPlugin.EZR","Error Bars"))	
	levelFrame <- tkframe(optionsFrame)
    levelVariable <- tclVar(dialog.values$confidence)
    levelField <- ttkentry(levelFrame, width="6", textvariable=levelVariable)
    tkgrid(labelRcmdr(levelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","   Level of confidence:")), levelField, sticky="w")	
    tkgrid(errorBarsFrame, labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","   ")), levelFrame, sticky="w")
    tkgrid(optionsFrame, columnspan=2, sticky="w")

    options2Frame <- tkframe(top)
	radioButtons(options2Frame, name="graph", buttons=c("narrow", "wide"), values=c("narrow", "wide"), initialValue=dialog.values$graph,
    labels=gettext(domain="R-RcmdrPlugin.EZR",c("Narrow view", "Wide view")), title=gettext(domain="R-RcmdrPlugin.EZR","When two factors were picked:"))
    radioButtons(options2Frame, name="line", buttons=c("color", "type", "width"), values=c("color", "type", "width"), initialValue=dialog.values$line,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Color", "Line type", "Line width")), title=gettext(domain="R-RcmdrPlugin.EZR","Line discrimination"))
	tkgrid(graphFrame, labelRcmdr(options2Frame, text="   "), lineFrame, sticky="w")
	tkgrid(options2Frame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }
	

StatMedLinePlot <- function(){
defaults <- list(data=NULL, group=NULL, axisLabel="", log=0, multi=0, subset = "")
dialog.values <- getDialog("StatMedLinePlot", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Line graph(Repeated measures)"))
    variablesFrame <- tkframe(top)
    dataBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Repeatedly measured data (pick at least 2)"), selectmode="multiple", listHeight=15, initialSelection=varPosn(dialog.values$data, "numeric"))
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    axisLabelFrame <- tkframe(top)
    axisLabelVariable <- tclVar(dialog.values$axisLabel)
    axisLabelField <- ttkentry(axisLabelFrame, width="40", textvariable=axisLabelVariable)

	optionsFrame <- tkframe(top)
	checkBoxes(frame="optionsFrame", boxes=c("log", "multi"), initialValues=c(dialog.values$log, dialog.values$multi),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Log y-axis", "Show different groups in separate graphs")))	

#    logFrame <- tkframe(top)
#    logVariable <- tclVar("0")
#    logCheckBox <- tkcheckbutton(logFrame, variable=logVariable)
#    multiFrame <- tkframe(top)
#    multiVariable <- tclVar("0")
#    multiCheckBox <- tkcheckbutton(multiFrame, variable=multiVariable)
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Line graph(Repeated measures)"), "#####", sep=""))
		dataSet <- ActiveDataSet()
        data <- getSelection(dataBox)
        group <- getSelection(groupBox)
        axisLabel <- tclvalue(axisLabelVariable)
        logy <- tclvalue(logVariable)
        multi <- tclvalue(multiVariable)
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
putDialog("StatMedLinePlot", list(data=data, group=group, axisLabel=axisLabel, log=logy, multi=multi, subset = tclvalue(subsetVariable))
)
        closeDialog()
		ndata <- length(data)
        if (ndata < 2) {
            errorCondition(recall=StatMedLinePlot, message=gettext(domain="R-RcmdrPlugin.EZR","Pick at least 2 repeatedly measured data"))
            return()
        }
		command <- paste("alldata <- c(", subset1, dataSet, subset2, "$", data[1], sep="")
		command2 <- paste('xlabels <- c("', data[1], '"', sep="")
		for (i in 2:ndata){
			command <- paste(command, ", ", subset1, dataSet, subset2, "$", data[i], sep="")
			command2 <- paste(command2, ', "', data[i], '"', sep="")
		}
		command <- paste(command, ")", sep="")
		command2 <- paste(command2, ")", sep="")
		doItAndPrint(command)
		doItAndPrint(command2)
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		if (logy==0){
			logy <- ""
			doItAndPrint("ylimu <- max(alldata, na.rm=TRUE)")
			doItAndPrint("yliml <- ylimu - (ylimu - min(alldata, na.rm=TRUE))*1.2")
			doItAndPrint("ylimu <- ylimu*1.2")
			}
		else{
			logy <- ', log="y"'
			doItAndPrint("ylimu <- max(alldata, na.rm=TRUE)")
			doItAndPrint("yliml <- min(alldata, na.rm=TRUE)")
		}
        if (axisLabel == gettext(domain="R-RcmdrPlugin.EZR","<use y-variable names>")){
            axisLabel <- if (legend) ""
                else if(length(y) == 1) y
                else paste(paste("(", 1:length(y), ") ", y, sep=""), collapse=", ")
        }
		if (length(group) == 0){
			command <- paste("y <- rbind(", subset1, dataSet, subset2, "$", data[1], sep="")
			for (i in 2:ndata){
				command <- paste(command, ", ", subset1, dataSet, subset2, "$", data[i], sep="")
			}
			command <- paste(command, ")", sep="")
			doItAndPrint(command)	
			command <- paste('matplot(y, type="o", lty=1, pch=1, col=1, ylab="', axisLabel, '", ylim=c(yliml, ylimu), axes=FALSE', logy, ")", sep="")
			doItAndPrint(command)
			doItAndPrint("box()")
			doItAndPrint("axis(2)")
			doItAndPrint(paste("axis(1, at=1:", ndata, ", labels=xlabels)", sep=""))
		}		
		
		if (length(group) == 1){	
			groups <- eval(parse(text=paste("levels(factor(", subset1, dataSet, subset2, "$", group, "))", sep="")))
			ngroup <- length(groups)
			groupmembers <- paste('c("', groups[1], '"', sep="")
			if (ngroup >= 2){
				for (i in 2:ngroup){
					groupmembers <- paste(groupmembers, ', "', groups[i], '"', sep="")
				}
			}
			groupmembers <- paste(groupmembers, ')', sep="")
			command <- paste("y <- rbind(", subset1, dataSet, subset2, "[", subset1, dataSet, subset2, "$", group, '=="', groups[1], '",]$', data[1], sep="")
			for (i in 2:ndata){
				command <- paste(command, ", ", subset1, dataSet, subset2, "[", subset1, dataSet, subset2, "$", group, '=="', groups[1], '",]$', data[i], sep="")
			}
			command <- paste(command, ")", sep="")
			doItAndPrint(command)	
			command <- paste('matplot(y, type="o", lty=1, pch=1, col=1, ylab="', axisLabel, '", ylim=c(yliml, ylimu), axes=FALSE', logy, ")", sep="")
			doItAndPrint(command)
			doItAndPrint("box()")
			doItAndPrint("axis(2)")
			doItAndPrint(paste("axis(1, at=1:", ndata, ", labels=xlabels)", sep=""))			
			
			if (ngroup >= 2){
				if (multi == 1){
					doItAndPrint(paste('legend("topright", "', group, "=", groups[1], '", box.lty=0)', sep=""))
				} else {
					command <- paste('legend("topright", ', groupmembers, ", col=1:", ngroup, ", lty=1:", ngroup, ", lwd=1:", ngroup, ', title="', group, '", box.lty=0)', sep="")
					doItAndPrint(command)	
				}
				for (j in 2:ngroup){			
					command <- paste("y <- rbind(", subset1, dataSet, subset2, "[", subset1, dataSet, subset2, "$", group, '=="', groups[j], '",]$', data[1], sep="")
					for (i in 2:ndata){
						command <- paste(command, ", ", subset1, dataSet, subset2, "[", subset1, dataSet, subset2, "$", group, '=="', groups[j], '",]$', data[i], sep="")
					}
					command <- paste(command, ")", sep="")
					doItAndPrint(command)
					if (multi == 1){
						if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
						command <- paste('matplot(y, type="o", lty=1, pch=1, col=1, ylab="', axisLabel, '", ylim=c(yliml, ylimu), axes=FALSE', logy, ")", sep="")
						doItAndPrint(command)
						doItAndPrint("box()")
						doItAndPrint("axis(2)")
						doItAndPrint(paste("axis(1, at=1:", ndata, ", labels=xlabels)", sep=""))
						doItAndPrint(paste('legend("topright", "', group, "=", groups[j], '", box.lty=0)', sep=""))
					} 
					else {
						command <- paste('matplot(y, type="o", lty=', j, ", pch=", j, ", lwd=", j, ", col=", j, ', ylab="', axisLabel, '", ylim=c(yliml, ylimu), axes=FALSE', logy, ", add=TRUE)", sep="")
						doItAndPrint(command)											
					}
				}
			}		
		}
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="matplot", apply="StatMedLinePlot", reset="StatMedLinePlot")
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
    tkgrid(getFrame(dataBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(labelRcmdr(axisLabelFrame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Label for y-axis"), ":  "), fg="blue"), axisLabelField, sticky="w")
#    tkgrid(axisLabelEntry, sticky="w")
    tkgrid(axisLabelFrame, sticky="w")
#    tkgrid(labelRcmdr(logFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Log y-axis")),
#        logCheckBox, sticky="w")
#    tkgrid(logFrame, sticky="w")
#    tkgrid(labelRcmdr(multiFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Show different groups in separate graphs")),
#        multiCheckBox, sticky="w")
 #   tkgrid(multiFrame, sticky="w")
	tkgrid(optionsFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, stick="w")
    dialogSuffix(rows=4, columns=1)
    }

	
StatMedMeanCI <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a mean"))
	variableFrame <- tkframe(top)
	mean <- tclVar("")
	meanEntry <- ttkentry(variableFrame, width="20", textvariable=mean)
	sd <- tclVar("")
	sdEntry <- ttkentry(variableFrame, width="20", textvariable=sd)
	variable2Frame <- tkframe(top)
	sample <- tclVar("")
	sampleEntry <- ttkentry(variable2Frame, width="20", textvariable=sample)
	CI <- tclVar("95")
	CIEntry <- ttkentry(variable2Frame, width="20", textvariable=CI)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a mean"), "#####", sep=""))
	mean <- tclvalue(mean)
	sd <- tclvalue(sd)
	sample <- tclvalue(sample)
	CI <- tclvalue(CI)
	closeDialog()
	if (length(mean) == 0 || length(sd) == 0 || length(sample) == 0){
			errorCondition(recall=StatMedMeanCI, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
	}
	doItAndPrint(paste("se <- ", sd, "/ sqrt(", sample, ")", sep=""))
	doItAndPrint(paste("CIL <- ", mean, " - qt((100+", CI, ")/200, ", sample, "-1)*se", sep=""))
	doItAndPrint(paste("CIH <- ", mean, " + qt((100+", CI, ")/200, ", sample, "-1)*se", sep=""))	
	doItAndPrint(paste('cat("', CI, '", gettext(domain="R-RcmdrPlugin.EZR","%CI"), " ", round(CIL,3), "-", round(CIH,3), "\n", sep="")'))	
	tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="epi.tests")
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Mean")), meanEntry, sticky="w")
	tkgrid.configure(meanEntry, sticky="w")
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation")), sdEntry, sticky="w")
	tkgrid.configure(sdEntry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size")), sampleEntry, sticky="w")
	tkgrid.configure(sampleEntry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval")), CIEntry, sticky="w")
	tkgrid.configure(CIEntry, sticky="w")
	tkgrid(variableFrame, sticky="nw")
	tkgrid(variable2Frame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}
	
	
StatMedSG <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Smirnov-Grubbs test for outliers"))
	variableBox <- variableListBox(top, Numeric(), selectmode="multiple",
		title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"), listHeight=15)
	radioButtons(name="remove", buttons=c("no", "yes"), values=c("0", "1"),
		labels=gettext(domain="R-RcmdrPlugin.EZR",c("No", "Yes")), title=gettext(domain="R-RcmdrPlugin.EZR","Create a variable converting outliers to NA"))
	newName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","<same as variables>"))
	newNameField <- ttkentry(top, width="20", textvariable=newName)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Smirnov-Grubbs test for outliers"), "#####", sep=""))
		variables <- getSelection(variableBox)
		closeDialog()
		if (length(variables) == 0) {
			errorCondition(recall=StatMedSG, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		newname <- trim.blanks(tclvalue(newName))
		remove <- tclvalue(removeVariable)
		.activeDataSet <- ActiveDataSet()
		if(remove==1){
			for (name in variables){
				nname <- if (newname == gettext(domain="R-RcmdrPlugin.EZR","<same as variables>")) name
					else if (length(variables) == 1) newname
					else paste(newname, name, sep="")
				if (!is.valid.name(nname)){
					errorCondition(recall=StatMedSG,
						message=paste('"', nname, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
					return()
				}
				if (is.element(nname, Variables())) {
					if ("no" == tclvalue(checkReplace(nname))){
						StatMedSG()
						return()
					}
				}
				command <- paste("(", .activeDataSet, "$", nname, " <- RemoveOutlier(", .activeDataSet, "$", name, ", return=1))", sep="")
				result <- doItAndPrint(command)
				if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
			}
		tkfocus(CommanderWindow())
		} else {
			for (name in variables){
				command <- paste("RemoveOutlier(", .activeDataSet, "$", name, ", return=0)", sep="")
				doItAndPrint(command)
			}
		}
	}
	OKCancelHelp()
	tkgrid(getFrame(variableBox), removeFrame, sticky="nw")
	tkgrid(labelRcmdr(top,
			text=gettext(domain="R-RcmdrPlugin.EZR","New variable name or prefix for multiple variables:")),
		newNameField, sticky="w")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=4, columns=2, preventGrabFocus=TRUE)
}

		
StatMedSingleSampleTTest <- function(){
defaults <- list(x=NULL, mu="0.0", confidence="0.95", alternative="two.sided", subset = "")
dialog.values <- getDialog("StatMedSingleSampleTTest", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Single-Sample t-Test"))
    xBox <- variableListBox(top, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$x, "numeric"))
	StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Single-Sample t-Test"), "#####", sep=""))
        x <- getSelection(xBox)
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceVariable)
        mu <- tclvalue(muVariable)
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
putDialog("StatMedSingleSampleTTest", list(x=x, mu=mu, confidence=tclvalue(confidenceVariable), alternative=alternative, subset = tclvalue(subsetVariable)))	
        if (length(x) == 0){
            errorCondition(recall=StatMedSingleSampleTTest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            return()
            }
		closeDialog()
        doItAndPrint(paste("(res <- t.test(", subset1, ActiveDataSet(), subset2, "$", x,
            ", alternative='", alternative, "', mu=", mu, ", conf.level=", level,
            "))", sep=""))
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "mean"), " = ", res$estimate, ", ", gettext(domain="R-RcmdrPlugin.EZR", "95% CI"), " ", res$conf.int[1], "-", res$conf.int[2], ", ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="t.test", apply="StatMedSingleSampleTTest", reset="StatMedSingleSampleTTest")
    radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),initialValue=dialog.values$alternative, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Population mean != mu0", "Population mean < mu0", "Population mean > mu0")),
        title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))
    rightFrame <- tkframe(top)	
    confidenceFrame <- tkframe(rightFrame)
    confidenceVariable <- tclVar(dialog.values$confidence)
    confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceVariable)
    muFrame <- tkframe(rightFrame)
    muVariable <- tclVar(dialog.values$mu)
    muField <- ttkentry(muFrame, width="6", textvariable=muVariable)

#    confidenceFrame <- tkframe(rightFrame)
#    confidenceLevel <- tclVar(".95")
#    confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
#    muFrame <- tkframe(rightFrame)
#    muVariable <- tclVar("0.0")
#    muField <- ttkentry(muFrame, width="8", textvariable=muVariable)
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(labelRcmdr(rightFrame, text=""), sticky="w")
    tkgrid(labelRcmdr(muFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Null hypothesis: mu = ")), muField, sticky="w")
    tkgrid(muFrame, sticky="w")
    tkgrid(labelRcmdr(confidenceFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence Level: ")), confidenceField, sticky="w")
    tkgrid(confidenceFrame, sticky="w")
    tkgrid(alternativeFrame, sticky="nw")
    tkgrid(rightFrame, sticky="nw")
    tkgrid.configure(confidenceField, sticky="e")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }

	
StatMedKS <- function(){
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Kolmogorov-Smimov test for normal distribution"))
    variablesFrame <- tkframe(top)
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable (pick one)"), listHeight=15)
	StatMedSubsetBox()
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Kolmogorov-Smimov test for normal distribution"), "#####", sep=""))
        response <- getSelection(responseBox)
        if (length(response) == 0) {
            errorCondition(recall=StatMedKS, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
        subset <- tclvalue(subsetVariable)
        .activeDataSet <- ActiveDataSet()
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset <- .activeDataSet
		} else {
			subset <- paste("subset(", .activeDataSet, ", ", subset, ")", sep="")
		}
        closeDialog()
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("hist2(", subset, "$", response, ', freq=F, main="", xlab="', response, '", ylab="", col="darkgray")', sep="")
		doItAndPrint(command)
		command <- paste("curve(dnorm(x, mean=mean(", subset, "$", response, "[!is.na(", subset, "$", response, ")]), sd=sd(", subset, "$", response, "[!is.na(", subset, "$", response, ")])), add=T)", sep="") 
		doItAndPrint(command)
		doItAndPrint(paste("skewness.kurtosis(", subset, "$", response, ")", sep=""))
        doItAndPrint(paste("ks.test(", subset, "$", response, ', "pnorm", mean=mean(', subset, "$", response, "[!is.na(", subset, "$", response, ")]), sd=sd(", subset, "$", response, "[!is.na(", subset, "$", response, ")]))", sep=""))	
		n <- eval(parse(text=paste("length(", subset, "$", response, ")", sep="")))	
		logger(paste(gettext(domain="R-RcmdrPlugin.EZR","# Shapiro-Wilk test can be performed only when the sample size is less than 5000. (Sample size ="), " ", n, ")", sep=""))
		if(n <= 5000){
			doItAndPrint(paste("shapiro.test(", subset, "$", response, ")", sep=""))
		}
		tkfocus(CommanderWindow())
        tkdestroy(top)
        }
    OKCancelHelp(helpSubject="ks.test")
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=5, columns=1)
    }

	
StatMedFTest <- function(){
defaults <- list(group=NULL, response=NULL, confidence="0.95", alternative="two.sided", subset = "")
dialog.values <- getDialog("StatMedFTest", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Two-variances F-test"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Groups (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
	StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Two-variances F-test"), "#####", sep=""))
        group <- getSelection(groupBox)
        response <- getSelection(responseBox)
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceVariable)
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
			subset <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
			subset <- paste(", subset=", subset, sep="")
		}
putDialog("StatMedFTest", list(group=group, response=response, confidence=level, alternative=alternative, subset = tclvalue(subsetVariable)))	
        if (length(group) == 0) {
            errorCondition(recall=StatMedFTest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups variable."))
            return()
            }
        if (length(response) == 0) {
            errorCondition(recall=StatMedFTest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
            closeDialog()
        .activeDataSet <- ActiveDataSet()
        doItAndPrint(paste("tapply(", subset1, .activeDataSet, subset2, "$", response, ", ", subset1,
            .activeDataSet, subset2, "$", group, ",  var, na.rm=TRUE)", sep=""))
        doItAndPrint(paste("(res <- var.test(", response, " ~ ", group,
            ", alternative='", alternative, "', conf.level=", level,
            ", data=", .activeDataSet, subset, "))", sep=""))
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "F test"), " ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")
		tkfocus(CommanderWindow())
        tkdestroy(top)
        }
    OKCancelHelp(helpSubject="var.test", apply="StatMedFTest", reset="StatMedFTest")
    radioButtons(name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Difference < 0", "Difference > 0")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))

    confidenceFrame <- tkframe(top)
    confidenceVariable <- tclVar(dialog.values$confidence)
    confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceVariable)
		
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
#    groupsLabel(groupsBox=groupBox)
    tkgrid(labelRcmdr(confidenceFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence Level:  "), fg="blue"), confidenceField, sticky="w")
    tkgrid(alternativeFrame, sticky="w")
    tkgrid(confidenceFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=5, columns=1)
    }

	
StatMedBartlett <- function(){
defaults <- list(group=NULL, response=NULL, subset = "")
dialog.values <- getDialog("StatMedBartlett", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Bartlett's test"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Groups (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
	StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Bartlett's test"), "#####", sep=""))
        group <- getSelection(groupBox)
        response <- getSelection(responseBox)
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
			subset <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
			subset <- paste(", subset=", subset, sep="")
		}
putDialog("StatMedBartlett", list(group=group, response=response, subset = tclvalue(subsetVariable)))	
        if (length(group) == 0) {
            errorCondition(recall=StatMedBartlett, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups variable."))
            return()
            }
        if (length(response) == 0) {
            errorCondition(recall=StatMedBartlett, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        doItAndPrint(paste("tapply(", subset1, .activeDataSet, subset2, "$", response, ", ", subset1,
            .activeDataSet, subset2, "$", group, ",  var, na.rm=TRUE)", sep=""))
        doItAndPrint(paste("(res <- bartlett.test(", response, " ~ ", group,
            ", data=", .activeDataSet, subset, "))", sep=""))
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "Bartlett test"), " ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")
        tkfocus(CommanderWindow())
        tkdestroy(top)
        }
    OKCancelHelp(helpSubject="bartlett.test", apply="StatMedBartlett", reset="StatMedBartlett")
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
#    groupsLabel(groupsBox=groupBox)
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=5, columns=1)
    }

	
StatMedTtest <- function(){
defaults <- list(group=NULL, response=NULL, confidence="0.95", alternative="two.sided", variances="TRUE", graph="bar", subset = "")
dialog.values <- getDialog("StatMedTtest", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Two-sample t-test"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(),selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variables with two levels (pick at least one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
 #change to "multiple" to accept analyses for multiple factors
 #change to "Variables()" to accept numeric variabels as grouping variable
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Two-sample t-test"), "#####", sep=""))	
	    group <- getSelection(groupBox)
        response <- getSelection(responseBox)
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceVariable)
        variances <- as.character(tclvalue(variancesVariable))
        graph <- as.character(tclvalue(graphVariable))
		    subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}			
putDialog("StatMedTtest", list(group=group, response=response, confidence=tclvalue(confidenceVariable), alternative=alternative, variances=variances, graph=graph, subset = tclvalue(subsetVariable)))
        if (length(group) == 0) {
            errorCondition(recall=StatMedTtest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups variable."))
            return()
            }
        if (length(response) == 0) {
            errorCondition(recall=StatMedTtest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
        closeDialog()
	nvar = length(group)
	doItAndPrint("group.names <- NULL")
	doItAndPrint("group.means <- NULL")
	doItAndPrint("group.sds <- NULL")
	doItAndPrint("group.p <- NULL")
	for (i in 1:nvar) {
	        doItAndPrint(paste("(res <- t.test(", response, "~factor(", group[i],
        	    "), alternative='", alternative, "', conf.level=", level,
          	  ", var.equal=", variances,
            	", data=", ActiveDataSet(), subset, "))", sep=""))
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			if (graph == "box"){
			    command <- (paste("boxplot(", response, "~ factor(", group[i], '), ylab="', response,
                '", xlab="', group[i], '"',
                ", data=", ActiveDataSet(), subset, ")", sep=""))
            logger(command)
            justDoIt(command)
			}
			if (graph == "point"){
				command <- paste("StatMedplotMeans(", subset1, ActiveDataSet(), subset2, "$", response,
				", factor(", subset1, ActiveDataSet(), subset2, "$", group[i],
				'), ylab="', response, '", xlab="', group[i], '", error.bars="sd", level=0.95)', sep="")
				logger(command)
				justDoIt(command)
			}
			doItAndPrint(paste("bar.means <- tapply(", subset1, ActiveDataSet(), subset2, "$", response, ", factor(", subset1, ActiveDataSet(), subset2, "$", group[i], "), mean, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.sds <- tapply(", subset1, ActiveDataSet(), subset2, "$", response, ", factor(", subset1, ActiveDataSet(), subset2, "$", group[i], "), sd, na.rm=TRUE)", sep=""))
			if (graph == "bar"){
				doItAndPrint("bar.sds <- ifelse(is.na(bar.sds), 0, bar.sds)")
				doItAndPrint(paste('barx <- barplot(bar.means, ylim=c(ifelse(min(bar.means, na.rm=TRUE)>0, 0, min(bar.means-bar.sds, na.rm=TRUE)*1.2), max(bar.means+bar.sds, na.rm=TRUE)*1.2), xlab="', group[i], '", ylab="', response, '", axis.lty=1)',sep=""))
				doItAndPrint(paste("error.bar(barx, bar.means, bar.sds)", sep=""))
			}
			group.levels <- eval(parse(text=paste("levels(factor(", subset1, ActiveDataSet(), subset2, "$", group[i], "))", sep="")))
			for (j in 1:2){
				doItAndPrint(paste('group.names <- c(group.names, "', group[i], "=", group.levels[j], '")', sep=""))
				doItAndPrint(paste("group.means <- c(group.means, bar.means[", j, "])", sep=""))
				doItAndPrint(paste("group.sds <- c(group.sds, bar.sds[", j, "])", sep=""))
				if (j == 1){
					doItAndPrint("group.p <- c(group.p, signif(res$p.value,digits=3))")
				} else {
					doItAndPrint('group.p <- c(group.p, "")')	
				}
			}
			doItAndPrint("remove(res)")	
	}
		doItAndPrint("summary.ttest <- data.frame(mean=group.means, sd=group.sds, p.value=group.p)")
		doItAndPrint("rownames(summary.ttest) <- group.names")
		doItAndPrint('colnames(summary.ttest) <- gettext(domain="R-RcmdrPlugin.EZR",colnames(summary.ttest))')
		doItAndPrint("summary.ttest")	
#		doItAndPrint("remove(summary.ttest)")				
	tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="t.test", apply="StatMedTtest", reset="StatMedTtest")
    optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Difference < 0", "Difference > 0")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))
		
    confidenceFrame <- tkframe(optionsFrame)
    confidenceVariable <- tclVar(dialog.values$confidence)
    confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceVariable)

    radioButtons(optionsFrame, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue=dialog.values$variances,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes (t-test)", "No (Welch test)")), title=gettext(domain="R-RcmdrPlugin.EZR","Assume equal variances?"))
    radioButtons(optionsFrame, name="graph", buttons=c("box", "bar", "point"), values=c("box", "bar", "point"), initialValue=dialog.values$graph,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("BoxGraph", "BarGraph", "LinePlot")), title=gettext(domain="R-RcmdrPlugin.EZR","Graphs"))
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
	tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(labelRcmdr(confidenceFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence Level"), fg="blue"))
    tkgrid(confidenceField, sticky="w")
    groupsLabel(groupsBox=groupBox)
    tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), confidenceFrame, labelRcmdr(optionsFrame, text="    "), 
        variancesFrame, labelRcmdr(optionsFrame, text="    "), graphFrame, sticky="nw")
    tkgrid(optionsFrame, sticky="nw")
#    tkgrid(confidenceFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
    }


StatMedPairedTtest <- function(){
defaults <- list(xBox=NULL, yBox=NULL, confidence="0.95", alternative="two.sided", subset = "")
dialog.values <- getDialog("StatMedPairedTtest", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Paired t-test"))
    .numeric <- Numeric()
    variablesFrame <- tkframe(top)
    xBox <- variableListBox(variablesFrame, .numeric, title=gettext(domain="R-RcmdrPlugin.EZR","First variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$xBox, "numeric"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettext(domain="R-RcmdrPlugin.EZR","Second variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$yBox, "numeric"))
	StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Paired t-test"), "#####", sep=""))
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceVariable)
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}		
putDialog("StatMedPairedTtest", list(xBox=x, yBox=y, confidence=level, alternative=alternative, subset = tclvalue(subsetVariable)))
        if (length(x) == 0 | length(y) == 0){
            errorCondition(recall=StatMedPairedTtest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            return()
            }
        if (x == y){
            errorCondition(recall=StatMedPairedTtest, message=gettext(domain="R-RcmdrPlugin.EZR","Variables must be different."))
            return()
            }		
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        doItAndPrint(paste("(res <- t.test(", subset1, .activeDataSet, subset2, "$", x, ", ",
            subset1, .activeDataSet, subset2, "$", y,
            ", alternative='", alternative, "', conf.level=", level,
            ", paired=TRUE))", sep=""))

        doItAndPrint(paste("mean1 <- mean(", subset1, .activeDataSet, subset2, "$", x, ", na.rm=TRUE)", sep=""))
        doItAndPrint(paste("mean2 <- mean(", subset1, .activeDataSet, subset2, "$", y, ", na.rm=TRUE)", sep=""))
        doItAndPrint(paste("sd1 <- sd(", subset1, .activeDataSet, subset2, "$", x, ", na.rm=TRUE)", sep=""))
        doItAndPrint(paste("sd2 <- sd(", subset1, .activeDataSet, subset2, "$", y, ", na.rm=TRUE)", sep=""))
		doItAndPrint('summary.ttest <- data.frame(mean=c(mean1, mean2), sd=c(sd1, sd2), p.value=c(signif(res$p.value, digit=3),""))')
		doItAndPrint(paste('rownames(summary.ttest) <- c("', x, '", "', y, '")', sep=""))
		doItAndPrint('colnames(summary.ttest) <- gettext(domain="R-RcmdrPlugin.EZR",colnames(summary.ttest))')
		doItAndPrint("summary.ttest")	
						
		tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="t.test", apply="StatMedPairedTtest", reset="StatMedPairedTtest")
    radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Difference < 0", "Difference > 0")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))

    confidenceFrame <- tkframe(top)
    confidenceVariable <- tclVar(dialog.values$confidence)
    confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceVariable)

    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(labelRcmdr(confidenceFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence Level"), fg="blue"))
    tkgrid(confidenceField, sticky="w")
    tkgrid(alternativeFrame, sticky="nw")
    tkgrid(confidenceFrame, sticky="nw")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }

	
StatMedANOVA <- function(){
		Library("multcomp")
		Library("abind")
defaults <- list(group=NULL, response=NULL, variances="TRUE", pairwise=0, dunnett=0, bonferroni=0, holm=0, actmodel=0, graph="bar", subset = "")
dialog.values <- getDialog("StatMedANOVA", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

		initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","One-way ANOVA"))
		UpdateModelNumber()
		modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
		modelFrame <- tkframe(top)
		model <- ttkentry(modelFrame, width="20", textvariable=modelName)
		variablesFrame <- tkframe(top)
		groupBox <- variableListBox(variablesFrame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variables (pick at least one)"), listHeight=12, initialSelection=varPosn(dialog.values$group, "all"))
		responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=12, initialSelection=varPosn(dialog.values$response, "numeric"))

#tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison not performed when more than one grouping variables are picked."), fg="blue"), sticky="w")
		optionsFrame <- tkframe(top)
checkBoxes(frame="optionsFrame", boxes=c("bonferroni", "holm", "pairwise", "dunnett"), initialValues=c(dialog.values$bonferroni, dialog.values$holm, dialog.values$pairwise, dialog.values$dunnett),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Pairwise comparison (Bonferroni)", "Pairwise comparison (Holm)","Pairwise comparison (Tukey)", "Pairwise comparison (Dunnett)")))
#tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","The first group in alphabetical will be treated as the reference group."), fg="blue"), sticky="w")
		options2Frame <- tkframe(top)
checkBoxes(frame="options2Frame", boxes="actmodel", initialValues=dialog.values$actmodel,labels=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses"))

#		pairwiseVariable <- tclVar("0")
#		pairwiseCheckBox <- tkcheckbutton(optionsFrame, variable=pairwiseVariable)
#		dunnettVariable <- tclVar("0")
#		dunnettCheckBox <- tkcheckbutton(optionsFrame, variable=dunnettVariable)
#		bonferroniVariable <- tclVar("0")
#		bonferroniCheckBox <- tkcheckbutton(optionsFrame, variable=bonferroniVariable)
#		holmVariable <- tclVar("0")
#		holmCheckBox <- tkcheckbutton(optionsFrame, variable=holmVariable)
#		actmodelVariable <- tclVar("0")
#		actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)

		StatMedSubsetBox(model=TRUE)
		onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","One-way ANOVA"), "#####", sep=""))
			modelValue <- trim.blanks(tclvalue(modelName))
			group <- getSelection(groupBox)			
			response <- getSelection(responseBox)
			variances <- as.character(tclvalue(variancesVariable))
			graph <- as.character(tclvalue(graphVariable))
			pairwise <- tclvalue(pairwiseVariable)
			dunnett <- tclvalue(dunnettVariable)
			bonferroni <- tclvalue(bonferroniVariable)
			holm <- tclvalue(holmVariable)
			actmodel <- tclvalue(actmodelVariable)
		    subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}
putDialog("StatMedANOVA", list(group=group, response=response, variances=variances, pairwise=pairwise, dunnett=dunnett, bonferroni=bonferroni, holm=holm, actmodel=actmodel, graph=graph, subset=tclvalue(subsetVariable)))
			if (!is.valid.name(modelValue)){
				UpdateModelNumber(-1)
				errorCondition(recall=StatMedANOVA, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue))
				return()
			}
			if (is.element(modelValue, listAOVModels())) {
				if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
					UpdateModelNumber(-1)
					tkdestroy(top)
					oneWayAnova()
					return()
				}
			}
			closeDialog()
			if (length(group) == 0){
				errorCondition(recall=StatMedANOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups factor."))
				return()
			}
			if (length(response) == 0){
				errorCondition(recall=StatMedANOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
				return()
			}
			.activeDataSet <- ActiveDataSet()
			
			nvar = length(group)
			doItAndPrint("group.names <- NULL")
			doItAndPrint("group.means <- NULL")
			doItAndPrint("group.sds <- NULL")
			doItAndPrint("group.p <- NULL")
			for (i in 1:nvar) {
				if(variances=="TRUE"){
					command <- paste(modelValue, " <- aov(", response, " ~ factor(", group[i], "), data=", .activeDataSet, subset, ", na.action=na.omit)", sep="")
					justDoIt(command)
					logger(command)
				}
#		    	assign(modelValue, justDoIt(command), envir=.GlobalEnv)			
#				doItAndPrint(paste("numSummary(", subset1, .activeDataSet, subset2, "$", response, " , groups=", subset1, .activeDataSet, subset2, "$", group[i], ', statistics=c("mean", "sd"))', sep=""))					
				if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}

				#bar.means and bar.sds are required to show summary.anova even for "box" or "point"	
				doItAndPrint(paste("bar.means <- tapply(", subset1, ActiveDataSet(), subset2, "$", response, ", factor(", subset1, ActiveDataSet(), subset2, "$", group[i], "), mean, na.rm=TRUE)", sep=""))
				doItAndPrint(paste("bar.sds <- tapply(", subset1, ActiveDataSet(), subset2, "$", response, ", factor(", subset1, ActiveDataSet(), subset2, "$", group[i], "), sd, na.rm=TRUE)", sep=""))

				if (graph == "box"){
					command <- (paste("boxplot(", response, "~ factor(", group[i], '), ylab="', response,
					'", xlab="', group[i], '"',
					", data=", ActiveDataSet(), subset, ")", sep=""))
				logger(command)
				justDoIt(command)
				}
				if (graph == "point"){
				command <- paste("StatMedplotMeans(", subset1, ActiveDataSet(), subset2, "$", response,
				", factor(", subset1, ActiveDataSet(), subset2, "$", group[i],
				'), ylab="', response, '", xlab="', group[i], '", error.bars="sd", level=0.95)', sep="")
				logger(command)
				justDoIt(command)
				}
				if (graph == "bar"){
				doItAndPrint(
'error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}')
				doItAndPrint(paste('barx <- barplot(bar.means, ylim=c(ifelse(min(bar.means, na.rm=TRUE)>0, 0, min(bar.means-bar.sds, na.rm=TRUE)*1.2), max(bar.means+bar.sds, na.rm=TRUE)*1.2), xlab="', group[i], '", ylab="', response, '", axis.lty=1)',sep=""))
				doItAndPrint(paste("error.bar(barx, bar.means, bar.sds)", sep=""))
				}
				group.levels <- eval(parse(text=paste("levels(factor(", subset1, ActiveDataSet(), subset2, "$", group[i], "))", sep="")))
				for (j in 1:length(group.levels)){
				doItAndPrint(paste('group.names <- c(group.names, "', group[i], "=", group.levels[j], '")', sep=""))
				doItAndPrint(paste("group.means <- c(group.means, bar.means[", j, "])", sep=""))
				doItAndPrint(paste("group.sds <- c(group.sds, bar.sds[", j, "])", sep=""))
				if (j == 1 & variances=="TRUE"){
					doItAndPrint(paste("res <- summary(lm(", response, " ~ factor(", group[i], "), data=", .activeDataSet, subset, "))", sep=""))
					doItAndPrint('group.p <- c(group.p, signif(pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail=FALSE), digits=3))')
					doItAndPrint("remove(res)")
				} else if(j == 1 & variances=="FALSE"){
					doItAndPrint(paste("res <- oneway.test(", response, " ~ factor(", group[i], "), data=", .activeDataSet, subset, ", var.equal=FALSE)", sep=""))
					doItAndPrint('group.p <- c(group.p, signif(res$p.value, digits=3))')				
				} else {
					doItAndPrint('group.p <- c(group.p, "")')	
				}
				}	
				if(variances=="TRUE") doItAndPrint(paste("summary(", modelValue, ")", sep=""))				
			}
			doItAndPrint("summary.anova <- data.frame(mean=group.means, sd=group.sds, p.value=group.p)")
			doItAndPrint("rownames(summary.anova) <- group.names")
			doItAndPrint('colnames(summary.anova) <- gettext(domain="R-RcmdrPlugin.EZR",colnames(summary.anova))')
			doItAndPrint("summary.anova")	
#			doItAndPrint("remove(summary.anova)")				
			if (bonferroni == 1 && nvar == 1 && variances=="TRUE"){
					dataSet=ActiveDataSet()
					doItAndPrint(paste("pairwise.t.test(", subset1, dataSet, subset2, "$", response, ", ", subset1, dataSet, subset2, "$", group, ", var.equal=", variances, ', p.adj="bonferroni")', sep=""))
			}
			if (holm == 1 && nvar == 1 && variances=="TRUE"){
					dataSet=ActiveDataSet()
					doItAndPrint(paste("pairwise.t.test(", subset1, dataSet, subset2, "$", response, ", ", subset1, dataSet, subset2, "$", group, ", var.equal=", variances, ', p.adj="holm")', sep=""))
			}
			if (pairwise == 1 && nvar == 1 && variances=="TRUE") {
				if (eval(parse(text=paste("length(levels(factor(", subset1, .activeDataSet, subset2, "$", group, "))) < 3"))))
					Message(message=gettext(domain="R-RcmdrPlugin.EZR","Factor has fewer than 3 levels; pairwise comparisons omitted."),
						type="warning")
				# the following lines modified by Richard Heiberger and subsequently by J. Fox
				else {
#					command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", group, ' = "Tukey"))', sep="")
					command <- paste("TukeyHSD(", modelValue, ', "factor(', group, ')")', sep="")
					doItAndPrint(command)
					if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
					command <- paste("plot(TukeyHSD(", modelValue, ', "factor(', group, ')"))', sep="")
					doItAndPrint(command)
#					doItAndPrint("confint(.Pairs) # confidence intervals")
#					doItAndPrint("cld(.Pairs) # compact letter display")
#					justDoIt("old.oma <- par(oma=c(0,5,0,0))")
#					logger("old.oma <- par(oma=c(0,5,0,0))")
#					justDoIt("plot(confint(.Pairs))")
#					logger("plot(confint(.Pairs))")
#					justDoIt("par(old.oma)")
#					logger("par(old.oma)")
#					logger("remove(.Pairs)")
#					remove(.Pairs, envir=.GlobalEnv)
				}
			}
			if (dunnett == 1 && nvar == 1 && variances=="TRUE"){
					doItAndPrint(paste("group.factor <- factor(", subset1, .activeDataSet, subset2, "$", group, ")", sep=""))
					command <- paste("res <- aov(", response, " ~ group.factor, data=", .activeDataSet, subset, ")", sep="")
					doItAndPrint(command)
					command <- 'summary(glht(res, linfct=mcp(group.factor="Dunnett")))'
					doItAndPrint(command)
			}
			if (actmodel==1) activeModel(modelValue)
			tkfocus(CommanderWindow())
		}
		OKCancelHelp(helpSubject="anova", model=TRUE, apply="StatMedANOVA", reset="StatMedANOVA")
		tkgrid(labelRcmdr(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model: ")), model, sticky="w")
		tkgrid(modelFrame, sticky="w", columnspan=2)
		tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
		tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
		tkgrid(variablesFrame, sticky="w")	
		
		options0Frame <- tkframe(top)
		radioButtons(options0Frame, name="graph", buttons=c("box", "bar", "point"), values=c("box", "bar", "point"), initialValue=dialog.values$graph,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("BoxGraph", "BarGraph", "LinePlot")), title=gettext(domain="R-RcmdrPlugin.EZR","Graphs"))
		radioButtons(options0Frame, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue=dialog.values$variances,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes (ANOVA)", "No (Welch test)")), title=gettext(domain="R-RcmdrPlugin.EZR","Assume equal variances?"))
		tkgrid(graphFrame, labelRcmdr(options0Frame, text="    "), variancesFrame, sticky="nw")
		tkgrid(options0Frame, sticky="nw")
	#		tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison not performed when more than one grouping variables are picked."), fg="blue"), sticky="w")
#		tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Bonferroni)")), bonferroniCheckBox, sticky="w")
#		tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Holm)")), holmCheckBox, sticky="w")
#		tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Tukey)")), pairwiseCheckBox, sticky="w")
#		tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Dunnett)")), dunnettCheckBox, sticky="w")
#		tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","The first group in alphabetical will be treated as the reference group."), fg="blue"), sticky="w")
#		tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")

tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison and active model keeping not performed for Welch test."), fg="blue"), sticky="w")
tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison not performed when more than one grouping variables are picked."), fg="blue"), sticky="w")
		tkgrid(optionsFrame, sticky="w", columnspan=2)
tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","The first group in alphabetical will be treated as the reference group."), fg="blue"), sticky="w")
		tkgrid(options2Frame, sticky="w", columnspan=2)
		tkgrid(subsetFrame, sticky="w")
		tkgrid(buttonsFrame, columnspan=2, sticky="w")
		dialogSuffix(rows=4, columns=2)
	}	


StatMedRepANOVA <- function(){
defaults <- list(group=NULL, data=NULL, line="color", bonferroni=0, holm=0, actmodel=0, subset = "")
dialog.values <- getDialog("StatMedRepANOVA", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Repeated-measures ANOVA"))
	UpdateModelNumber()
	modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width="20", textvariable=modelName)
    variablesFrame <- tkframe(top)
    dataBox <- variableListBox(variablesFrame, Numeric(),selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Repeatedly measured data (pick at least 2)"), listHeight=15, initialSelection=varPosn(dialog.values$data, "numeric"))
    groupBox <- variableListBox(variablesFrame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (pick 0, 1, or more)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
	optionsFrame <- tkframe(top)	
checkBoxes(frame="optionsFrame", boxes=c("bonferroni", "holm", "actmodel"), initialValues=c(dialog.values$bonferroni, dialog.values$holm, dialog.values$actmodel),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Pairwise comparison (Bonferroni)", "Pairwise comparison (Holm)", "Keep results as active model for further analyses")))	
#	bonferroniVariable <- tclVar("0")
#	bonferroniCheckBox <- tkcheckbutton(optionsFrame, variable=bonferroniVariable)
#	holmVariable <- tclVar("0")
#	holmCheckBox <- tkcheckbutton(optionsFrame, variable=holmVariable)
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Repeated-measures ANOVA"), "#####", sep=""))
		modelValue <- trim.blanks(tclvalue(modelName))
        data <- getSelection(dataBox)
		group <- getSelection(groupBox)
		bonferroni <- tclvalue(bonferroniVariable)
		holm <- tclvalue(holmVariable)
		actmodel <- tclvalue(actmodelVariable)
		dataSet <- ActiveDataSet()
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			doItAndPrint(paste("TempDF <- ", dataSet))
		}
        else {
			doItAndPrint(paste("TempDF <- subset(", dataSet, ",", subset, ")") )
        }
		line <- tclvalue(lineVariable)
		if (line=="color") line <- ", lty=1, lwd=1"
		if (line=="type") line <- ", col=1, lwd=1"
		if (line=="width") line <- ", col=1, lty=1"
putDialog("StatMedRepANOVA", list(group=group, data=data, line=tclvalue(lineVariable), bonferroni=bonferroni, holm=holm, actmodel=actmodel, subset = tclvalue(subsetVariable)))
		if (!is.valid.name(modelValue)){
			UpdateModelNumber(-1)
			errorCondition(recall=StatMedRepANOVA, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue))
			return()
		}
		if (is.element(modelValue, listLMModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
				UpdateModelNumber(-1)
				tkdestroy(top)
				StatMedRepANOVA()
				return()
			}
		}
        if (length(data) < 2) {
            errorCondition(recall=StatMedRepANOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
		closeDialog()
		nvar <- length(data)
		RepeatedData <- data[1]
		RepeatedData2 <- paste('"', data[1], '"', sep="")
		for(i in 2:nvar){
			RepeatedData <- paste(RepeatedData, ", ", data[i], sep="")
			RepeatedData2 <- paste(RepeatedData2, ', "', data[i], '"', sep="")
			}
		nvar2 <- length(group)
		if (nvar2 >= 1){
			for(i in 1:nvar2){
				doItAndPrint(paste("TempDF$Factor", i, ".", group[i], " <- factor(TempDF$", group[i], ")", sep=""))
				doItAndPrint(paste("contrasts(TempDF$Factor", i, ".", group[i], ') <- "contr.Sum"', sep=""))
			}
		}
		if(nvar2 == 0){
			factors <- "1"
		}		
		if(nvar2 == 1){
			factors <- paste("Factor1.", group[1], sep="")
		}
		if(nvar2 >= 2){
			factors <- paste("Factor1.", group[1], sep="")
			for(i in 2:nvar2){
				factors <- paste(factors, "*Factor", i, ".", group[i], sep="")
			}
		}

		logger("#Convert to long format to draw graph")
		doItAndPrint("n <- length(TempDF[,1])")
		doItAndPrint("TempDF$TempID <- c(1:n)")
		command <- "TempDF2 <- data.frame(TempID=TempDF$TempID"
		for (i in 1:nvar){
			command <- paste(command, ", ", data[i], "=TempDF$", data[i], sep="")
		}
		if (length(group) == 0){
		}
		else{
			for (i in 1:length(group)){
				command <- paste(command, ", ", group[i], "=TempDF$", group[i], sep="")
			}
		}
		command <- paste(command, ")", sep="")
		doItAndPrint(command)
		doItAndPrint("TempDF2 <- na.omit(TempDF2)")		#delete rows with NA
		command <- paste('TempDF3 <- reshape(TempDF2, idvar="TempID", varying=list(c("', data[1], sep="")
		for (i in 2:nvar){
			command <- paste(command, '", "', data[i], sep="")
		}
		command <- paste(command, '")), v.names="data", direction="long")', sep="")
		doItAndPrint(command)
		command <- paste('RepeatNumber <- c("', data[1], sep="")
		for (i in 2:nvar){
			command <- paste(command, '", "', data[i], sep="")
		}
		command <- paste(command, '")', sep="")
		doItAndPrint(command)
		doItAndPrint("nvar <- length(TempDF3$time)")
		doItAndPrint("for (i in 1:nvar){TempDF3$time2[i] <- RepeatNumber[TempDF3$time[i]]}")
		if (length(group) == 0){
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint('StatMedplotMeans(TempDF3$data, factor(TempDF3$time2), error.bars="sd", xlab="", ylab="")')
		}
		if (length(group) == 1){
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}	
			doItAndPrint(paste("StatMedplotMeans(TempDF3$data, factor(TempDF3$time2), factor(TempDF3$", group[1], '), error.bars="sd", xlab="", ylab="", legend.lab="', group[1], '", ', line, ")", sep=""))
		}
		if (length(group) == 2){
			doItAndPrint(paste("for (i in 1:length(levels(factor(TempDF3$", group[1], ")))){windows(); par(", par.option, "); StatMedplotMeans(TempDF3$data[TempDF3$", group[1], "==levels(factor(TempDF3$", group[1], "))[i]], factor(TempDF3$time2[TempDF3$", group[1], "==levels(factor(TempDF3$", group[1], "))[i]]), factor(TempDF3$", group[2], "[TempDF3$", group[1], "==levels(factor(TempDF3$", group[1], "))[i]]), error.bars=", '"sd", xlab="", ylab="", legend.lab="', group[2], '", main=paste("', group[1], '", " : ", levels(factor(TempDF3$', group[1], "))[i]), ", line, ")}", sep=""))
		}
		command <- paste("lm(cbind(", RepeatedData, ") ~ ", factors, ", data=TempDF, na.action=na.omit)", sep="")
# 		logger(paste(modelValue, " <- ", command, sep = ""))
# 		assign(modelValue, justDoIt(command), envir = .GlobalEnv)
		doItAndPrint(paste(modelValue, " <- ", command, sep = ""))
		doItAndPrint(paste("time <- factor(c(", RepeatedData2, "))", sep=""))
		doItAndPrint("time <- data.frame(Time = time)")
		doItAndPrint(paste("res <- Anova(", modelValue, ', idata=time, idesign=~Time, type="III")', sep=""))
		if (actmodel==1) activeModel(modelValue)
		doItAndPrint("summary(res, multivariate=FALSE)")
		
		if (bonferroni == 1 && length(group) == 0){
			command <- paste("pairwise.pairedt.test(with(TempDF, cbind(", RepeatedData, ')), group=NULL, "', dataSet, '", p.adjust.method="bonferroni")', sep="")
			doItAndPrint(command)
		}
		if (bonferroni == 1 && length(group) == 1){
			command <- paste("pairwise.pairedt.test(with(TempDF, cbind(", RepeatedData, ")), TempDF$", factors, ', "', dataSet, '", p.adjust.method="bonferroni")', sep="")
			doItAndPrint(command)
		}
		if (holm == 1 && length(group) == 0){
			command <- paste("pairwise.pairedt.test(with(TempDF, cbind(", RepeatedData, ')), group=NULL, "', dataSet, '", p.adjust.method="holm")', sep="")
			doItAndPrint(command)
		}
		if (holm == 1 && length(group) == 1){
			command <- paste("pairwise.pairedt.test(with(TempDF, cbind(", RepeatedData, ")), TempDF$", factors, ', "', dataSet, '", p.adjust.method="holm")', sep="")
			doItAndPrint(command)
		}		
		doItAndPrint("remove(res)")
		tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Anova", model=TRUE, apply="StatMedRepANOVA", reset="StatMedRepANOVA")
	tkgrid(labelRcmdr(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model: ")), model, sticky="w")
	tkgrid(modelFrame, sticky="w", columnspan=2)
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
	tkgrid(getFrame(dataBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Graph not created when 3 or more grouping variables are picked.")), sticky="w")  
    radioButtons(name="line", buttons=c("color", "type", "width"), values=c("color", "type", "width"), initialValue=dialog.values$line,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Color", "Line type", "Line width")), title=gettext(domain="R-RcmdrPlugin.EZR","Line discrimination"))
	tkgrid(lineFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison not performed when more than one grouping variables are picked."), fg="blue"), sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Bonferroni)")), bonferroniCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Holm)")), holmCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
    }


StatMedMultiANOVA <- function(){
defaults <- list(group=NULL, data=NULL, interaction=1, actmodel=0, subset = "")
dialog.values <- getDialog("StatMedMultiANOVA", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Multi-way ANOVA"))
	UpdateModelNumber()
	modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width="20", textvariable=modelName)
    variablesFrame <- tkframe(top)
    dataBox <- variableListBox(variablesFrame, Numeric(),title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$data, "numeric"))
    groupBox <- variableListBox(variablesFrame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Factors (pick one or more)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    optionsFrame <- tkframe(top)
#    checkBoxes(window=optionsFrame, frame="interaction", boxes=c("interaction"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Include interaction term (when less than 4 grouping variables are picked)")))

	optionsFrame <- tkframe(top)	
checkBoxes(frame="optionsFrame", boxes="interaction", initialValues=dialog.values$interaction,labels=gettext(domain="R-RcmdrPlugin.EZR","Include interaction term (when less than 4 grouping variables are picked)"))	
	options2Frame <- tkframe(top)	
checkBoxes(frame="options2Frame", boxes="actmodel", initialValues=dialog.values$actmodel,labels=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses"))	

#	interactionVariable <- tclVar("1")
#	interactionCheckBox <- tkcheckbutton(optionsFrame, variable=interactionVariable)
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
	StatMedSubsetBox(model=TRUE)
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Multi-way ANOVA"), "#####", sep=""))
		modelValue <- trim.blanks(tclvalue(modelName))
        data <- getSelection(dataBox)
		group <- getSelection(groupBox)
		dataSet <- ActiveDataSet()
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			doItAndPrint(paste("TempDF <- ", dataSet))
		} else {
			doItAndPrint(paste("TempDF <- subset(", dataSet, ",", subset, ")") )
        }
		interaction <- tclvalue(interactionVariable)
		actmodel <- tclvalue(actmodelVariable)
		
	putDialog("StatMedMultiANOVA", list(group=group, data=data, interaction=interaction, actmodel=actmodel, subset = tclvalue(subsetVariable)))
		if (!is.valid.name(modelValue)){
			UpdateModelNumber(-1)
			errorCondition(recall=StatMedMultiANOVA, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue))
			return()
		}
		if (is.element(modelValue, listLMModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
				UpdateModelNumber(-1)
				tkdestroy(top)
				StatMedMultiANOVA()
				return()
			}
		}
        if (length(data) == 0) {
            errorCondition(recall=StatMedMultiANOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
		if (length(group) == 0) {
            errorCondition(recall=StatMedMultiANOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select at least one factor."))
            return()
            }
		closeDialog()
		nvar <- length(group)
		if (nvar <=3 & interaction==1){
			mark <- "*"
		} else {
			mark <- "+"
		}
		if (nvar >= 1){
			for(i in 1:nvar){
				doItAndPrint(paste("TempDF$Factor", i, ".", group[i], " <- factor(TempDF$", group[i], ")", sep=""))
				doItAndPrint(paste("contrasts(TempDF$Factor", i, ".", group[i], ') <- "contr.Sum"', sep=""))
			}
		}
		if(nvar == 1){
			factors <- paste(" + Factor1.", group[1], sep="")
		}
		if(nvar >= 2){
			factors <- paste(" + Factor1.", group[1], sep="")
			for (i in 2:nvar){
				factors <- paste(factors, mark, "Factor", i, ".", group[i], sep="")
			}
		}

		
		if (nvar == 1){
		    if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint(paste("bar.means <- tapply(TempDF$", data, ", TempDF$", group[1], ", mean, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.sds <- tapply(TempDF$", data, ", TempDF$", group[1], ", sd, na.rm=TRUE)", sep=""))
			doItAndPrint("bar.sds <- ifelse(is.na(bar.sds), 0, bar.sds)")
			doItAndPrint(paste('barx <- barplot(bar.means, ylim=c(ifelse(min(bar.means)>0, 0, min(bar.means-bar.sds)*1.2), max(bar.means+bar.sds)*1.2), xlab="', group[1], '", ylab="', data, '", axis.lty=1)',sep=""))
			doItAndPrint(paste("error.bar(barx, bar.means, bar.sds)", sep=""))
			}
		if (nvar == 2){
			if (eval(parse(text=paste("min(table(TempDF$", group[1], ", TempDF$", group[2], "))", sep="")))==0) {
				logger(gettext(domain="R-RcmdrPlugin.EZR","Graph not created when a group with 0 sample exists"))
			} else {			
			eval.bar.var <- eval(parse(text=paste("length(levels(factor(TempDF$", group[2], ")))", sep="")))
		    if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint(paste("bar.var <- length(levels(factor(TempDF$", group[2], ")))", sep=""))
			doItAndPrint(paste("bar.sums <- tapply(subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[1])$", data, ", subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[1])$", group[1], ", sum, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.means <- tapply(subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[1])$", data, ", subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[1])$", group[1], ", mean, na.rm=TRUE)", sep=""))
			doItAndPrint(paste("bar.sds <- tapply(subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[1])$", data, ", subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[1])$", group[1], ", sd, na.rm=TRUE)", sep=""))			
			if(eval.bar.var > 1){
				for (i in 2: eval.bar.var){
				doItAndPrint(paste("bar.sums <- c(bar.sums, tapply(subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[", i, "])$", data, ", subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[", i, "])$", group[1], ", sum, na.rm=TRUE))", sep=""))
				doItAndPrint(paste("bar.means <- c(bar.means, tapply(subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[", i, "])$", data, ", subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[", i, "])$", group[1], ", mean, na.rm=TRUE))", sep=""))
				doItAndPrint(paste("bar.sds <- c(bar.sds, tapply(subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[", i, "])$", data, ", subset(TempDF, ", group[2], "==levels(factor(", group[2], "))[", i, "])$", group[1], ", sd, na.rm=TRUE))", sep=""))
				doItAndPrint("bar.n <- bar.sums/bar.means")
				doItAndPrint("bar.ses <- bar.sds/sqrt(bar.n)")
				}
			}
			doItAndPrint(paste("bar.var2 <- length(levels(factor(TempDF$", group[1], ")))", sep=""))
			doItAndPrint("bar.means <- matrix(bar.means, bar.var2)")
			doItAndPrint("bar.sds <- matrix(bar.sds, bar.var2)")
			doItAndPrint("bar.ses <- matrix(bar.ses, bar.var2)")
			doItAndPrint("bar.sds <- ifelse(is.na(bar.sds), 0, bar.sds)")
			doItAndPrint("bar.ses <- ifelse(is.na(bar.ses), 0, bar.ses)")
			doItAndPrint(paste('barx <- barplot(bar.means, beside=TRUE, ylim=c(ifelse(min(bar.means)>0, 0, min(bar.means-bar.sds)*1.2), max(bar.means+bar.sds)*1.2), xlab="', group[2], '", ylab="', data, '", names.arg=levels(factor(TempDF$', group[2], ")), legend.text=levels(factor(TempDF$", group[1], ')), args.legend=list(title="', group[1], '", box.lty=0), axis.lty=1)',  sep=""))
			doItAndPrint("error.bar(barx, bar.means, bar.sds)")
			}		
		}
		groups.list <- paste(paste(group, "=TempDF$", group, sep=""), collapse=", ")
		doItAndPrint(paste("tapply(TempDF$", data, ", list(", groups.list,
					"), mean, na.rm=TRUE) # means", sep=""))
		doItAndPrint(paste("tapply(TempDF$", data, ", list(", groups.list,
					"), sd, na.rm=TRUE) # std. deviations", sep=""))
		doItAndPrint(paste("tapply(TempDF$", data, ", list(", groups.list,
					"), function(x) sum(!is.na(x))) # counts", sep=""))
		command <- paste("lm(", data, " ~ 1", factors, ", data=TempDF, na.action=na.omit)", sep="")
# 		logger(paste(modelValue, " <- ", command, sep = ""))
# 		assign(modelValue, justDoIt(command), envir = .GlobalEnv)
		doItAndPrint(paste(modelValue, " <- ", command, sep = ""))		
#		doItAndPrint("library(car)")
		doItAndPrint(paste("Anova(", modelValue, ', type="III")', sep=""))
		if (actmodel==1) activeModel(modelValue)
		tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Anova", model=TRUE, apply="StatMedMultiANOVA", reset="StatMedMultiANOVA")
	tkgrid(labelRcmdr(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model: ")), model, sticky="w")
	tkgrid(modelFrame, sticky="w", columnspan=2)
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
    tkgrid(getFrame(dataBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
#	tkgrid(labelRcmdr(optionsFrame, text=""), interaction, sticky="w")
	tkgrid(optionsFrame, sticky="nw")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Include interaction term (when less than 4 grouping variables are picked)")), interactionCheckBox, sticky="w")
 	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Graph not created when 3 or more grouping variables are picked.")), sticky="w")  
	tkgrid(options2Frame, sticky="nw")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}
	

StatMedANCOVA <- function(){
defaults <- list(group=NULL, data=NULL, covariate=NULL, actmodel=0, subset = "")
dialog.values <- getDialog("StatMedANCOVA", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","ANCOVA"))
	UpdateModelNumber()
	modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width="20", textvariable=modelName)
    variablesFrame <- tkframe(top)
    dataBox <- variableListBox(variablesFrame, Numeric(),title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=10, initialSelection=varPosn(dialog.values$data, "numeric"))
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (pick one)"), listHeight=10, initialSelection=varPosn(dialog.values$group, "all"))
    variables2Frame <- tkframe(top)
    covariateBox <- variableListBox(variables2Frame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Numeric variable for adjustment (pick one)"), listHeight=10, initialSelection=varPosn(dialog.values$covariate, "numeric"))	
    optionsFrame <- tkframe(top)
	checkBoxes(frame="optionsFrame", boxes="actmodel", initialValues=dialog.values$actmodel,labels=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses"))	
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","ANCOVA"), "#####", sep=""))
		modelValue <- trim.blanks(tclvalue(modelName))
        data <- getSelection(dataBox)
		group <- getSelection(groupBox)
		covariate <- getSelection(covariateBox)
		dataSet <- ActiveDataSet()
		actmodel <- tclvalue(actmodelVariable)
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			doItAndPrint(paste("TempDF <- ", dataSet))
		} else {
			doItAndPrint(paste("TempDF <- subset(", dataSet, ",", subset, ")") )
        }
putDialog("StatMedANCOVA", list(group=group, data=data, covariate=covariate, actmodel=actmodel, subset = tclvalue(subsetVariable)))
		if (!is.valid.name(modelValue)){
			UpdateModelNumber(-1)
			errorCondition(recall=StatMedANCOVA, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue))
			return()
		}
		if (is.element(modelValue, listLMModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
				UpdateModelNumber(-1)
				tkdestroy(top)
				StatMedMultiANOVA()
				return()
			}
		}
        if (length(data) == 0) {
            errorCondition(recall=StatMedANCOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
		if (length(group) == 0) {
            errorCondition(recall=StatMedANCOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select at least one factor."))
            return()
            }
		if (length(covariate) == 0) {
            errorCondition(recall=StatMedANCOVA, message=gettext(domain="R-RcmdrPlugin.EZR","You must select one numeric variable for adjustment."))
            return()
            }
		closeDialog()
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}		
		command <- paste("scatterplot(", data, " ~ ", covariate, " | factor(", group, "), reg.line=lm, smooth=FALSE, spread=FALSE, by.groups=TRUE, data=TempDF)", sep="")
		doItAndPrint(command)
#		doItAndPrint("library(car)")		
		interaction <- eval(parse(text=paste("signif(Anova(lm(", data, " ~ 1 + factor(", group, ") * ", covariate, ', data=TempDF, na.action=na.omit), type="III")$Pr[4], digits=3)', sep="")))
		doItAndPrint(paste('cat(gettext(domain="R-RcmdrPlugin.EZR","P value for interaction between grouping variable and covariate is"), ', " ", interaction, ', "\n")', sep=""))
		if(interaction < 0.05){
			logger(gettextRcmdr("ANCOVA not performed due to significant interaction between grouping variable and covariate."))
		} else {
			command <- paste(modelValue, " <- lm(", data, " ~ 1 + factor(", group, ") + ", covariate, ", data=TempDF, na.action=na.omit)", sep="")
			# 		logger(paste(modelValue, " <- ", command, sep = ""))
			# 		assign(modelValue, justDoIt(command), envir = .GlobalEnv)
			doItAndPrint(paste(modelValue, " <- ", command, sep = ""))
			doItAndPrint(paste("Anova(", modelValue, ', type="III")', sep=""))
			if (actmodel==1) activeModel(modelValue)
		}
		tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Anova", model=TRUE, apply="StatMedANCOVA", reset="StatMedANCOVA")
	tkgrid(labelRcmdr(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model: ")), model, sticky="w")
	tkgrid(modelFrame, sticky="w", columnspan=2)
    tkgrid(getFrame(dataBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(getFrame(covariateBox), labelRcmdr(variables2Frame, text="    "), sticky="nw")
    tkgrid(variables2Frame, sticky="nw")
	tkgrid(optionsFrame, sticky="nw")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}
	
	
StatMedCorrelation <- function(){
defaults <- list(x=NULL, alternative="two.sided", subset = "")
dialog.values <- getDialog("StatMedCorrelation", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Test for Pearson's correlation"))
  xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick two)"), listHeight=15, initialSelection=varPosn(dialog.values$x, "numeric"))
  radioButtons(name="alternative", buttons=c("two.sided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative,
               labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Correlation < 0", "Correlation > 0")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis")) 
	StatMedSubsetBox(model=TRUE)		   
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Test for Pearson's correlation"), "#####", sep=""))
    alternative <- as.character(tclvalue(alternativeVariable))
    x <- getSelection(xBox)
		    subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}
putDialog("StatMedCorrelation", list(x=x, alternative=alternative, subset = tclvalue(subsetVariable)))
    if (2 > length(x)) {
      errorCondition(recall=StatMedCorrelation,
        message=gettext(domain="R-RcmdrPlugin.EZR","Fewer than 2 variables selected."))
      return()
    }
    if(2 < length(x)) {
      errorCondition(recall=StatMedCorrelation,
        message=gettext(domain="R-RcmdrPlugin.EZR","More than 2 variables selected."))
      return()
    }			
    closeDialog()
    .activeDataSet <- ActiveDataSet()
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
    command2 <- paste("scatterplot(", x[1], "~", x[2],
        ", reg.line=lm, smooth=FALSE, spread=FALSE, boxplots='xy', span=0.5, data=", .activeDataSet, subset, ")", sep="")
    doItAndPrint(command2)  
    command <- paste("(res <- cor.test(", subset1, .activeDataSet, subset2, "$", x[1], ", ", subset1, .activeDataSet, subset2, "$", x[2],
        ', alternative="', alternative, '", method="pearson"))', sep="")
    doItAndPrint(command) 
	
	doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "correlation coefficient"), " = ", signif(res$estimate, digits=3), ", ", gettext(domain="R-RcmdrPlugin.EZR", "95% CI"), " ", signif(res$conf.int[1],digits=3), "-", signif(res$conf.int[2],digits=3), ", ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
	doItAndPrint("remove(res)")
		
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cor.test", apply="StatMedCorrelation", reset="StatMedCorrelation")
  tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(alternativeFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame,columnspan=2,sticky="w")
  dialogSuffix(rows=4, columns=1)
}
	
	
StatMedLinearRegression <- function(){
defaults <- list(x=NULL, y=NULL, wald=0, actmodel=0, diagnosis=0, stepwise1=0, stepwise2=0, stepwise3=0, subset = "")
dialog.values <- getDialog("StatMedLinearRegression", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Linear regression"))
    variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, Variables(), selectmode="multiple",
        title=gettext(domain="R-RcmdrPlugin.EZR","Explanatory variables (pick one or more)"), listHeight=15, initialSelection=varPosn(dialog.values$x, "all"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettext(domain="R-RcmdrPlugin.EZR","Response variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$y, "numeric"))
    UpdateModelNumber()
    modelName <- tclVar(paste("RegModel.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
	optionsFrame <- tkframe(top)
	
checkBoxes(frame="optionsFrame", boxes=c("wald", "actmodel", "diagnosis", "stepwise1", "stepwise2", "stepwise3"), initialValues=c(dialog.values$wald, dialog.values$actmodel, dialog.values$diagnosis, dialog.values$stepwise1, dialog.values$stepwise2, dialog.values$stepwise3),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Wald test for overall p-value for factors with >2 levels", "Keep results as active model for further analyses", "Show basic diagnostic plots", "Stepwise selection based on AIC", "Stepwise selection based on BIC", "Stepwise selection based on p-value")))	

#	waldVariable <- tclVar("0")
#	waldCheckBox <- tkcheckbutton(optionsFrame, variable=waldVariable)
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
#	stepwise1Variable <- tclVar("0")
#	stepwise2Variable <- tclVar("0")
#	stepwise3Variable <- tclVar("0")
#	stepwise1CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise1Variable)
#	stepwise2CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise2Variable)
#	stepwise3CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise3Variable)
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Linear regression"), "#####", sep=""))
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
		wald <- tclvalue(waldVariable)
		actmodel <- tclvalue(actmodelVariable)
		diagnosis <- tclvalue(diagnosisVariable)
		stepwise1 <- tclvalue(stepwise1Variable)
		stepwise2 <- tclvalue(stepwise2Variable)
		stepwise3 <- tclvalue(stepwise3Variable)
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>") || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
            }
        else{
            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
            }
putDialog("StatMedLinearRegression", list(x=x, y=y, wald=wald, actmodel=actmodel, diagnosis=diagnosis, stepwise1=stepwise1, stepwise2=stepwise2, stepwise3=stepwise3, subset = tclvalue(subsetVariable)))
        if (0 == length(y)) {
            UpdateModelNumber(-1)
            errorCondition(recall=StatMedLinearRegression, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
        if (0 == length(x)) {
            UpdateModelNumber(-1)
            errorCondition(recall=StatMedLinearRegression, message=gettext(domain="R-RcmdrPlugin.EZR","No explanatory variables selected."))
            return()
            }
        if (is.element(y, x)) {
            UpdateModelNumber(-1)
            errorCondition(recall=StatMedLinearRegression, message=gettext(domain="R-RcmdrPlugin.EZR","Response and explanatory variables must be different."))
            return()
            }			
		Library("aod")
		modelValue <- trim.blanks(tclvalue(modelName))
        if (!is.valid.name(modelValue)){
            UpdateModelNumber(-1)
            errorCondition(recall=StatMedLinearRegression, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue))
            return()
            }
        if (is.element(modelValue, listLinearModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
                UpdateModelNumber(-1)
                linearRegressionModel()
                return()
                }
            }
        command <- paste("lm(", y, "~", paste(x, collapse="+"),
            ", data=", ActiveDataSet(), subset, ")", sep="")
# 		logger(paste(modelValue, " <- ", command, sep = ""))
# 		assign(modelValue, justDoIt(command), envir = .GlobalEnv)
		doItAndPrint(paste(modelValue, " <- ", command, sep = ""))
        doItAndPrint(paste("(res <- summary(", modelValue, "))", sep=""))
        doItAndPrint(paste("vif(", modelValue, ")", sep=""))
		logger("###variance inflation factors")
		doItAndPrint('colnames(res$coefficients) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(res$coefficients))')
		doItAndPrint("res$coefficients")
		doItAndPrint("multireg.table <- res$coefficients")
		doItAndPrint("remove(res)")
		if (wald==1) doItAndPrint(paste("waldtest(", modelValue, ")", sep=""))
		if (diagnosis==1){
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}			
			doItAndPrint("oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))")
			doItAndPrint(paste("plot(", modelValue, ")", sep=""))
			doItAndPrint("par(oldpar)")			
		}
		if (stepwise1 == 1 | stepwise2 == 1 | stepwise3 == 1){
			command <- paste("TempDF <- with(", ActiveDataSet(), ", ", ActiveDataSet(), "[complete.cases(", paste(x, collapse=","), "),])", sep="")
			doItAndPrint(command)
			command <- paste("lm(", y, "~", paste(x, collapse="+"), ", data=TempDF", subset, ")", sep="")
			doItAndPrint(paste(modelValue, " <- ", command, sep=""))
			}
		if (stepwise1 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="AIC")', sep=""))
			doItAndPrint("summary(res)")
			if (wald==1) doItAndPrint("waldtest(res)")
			doItAndPrint("remove(res)")
			}
		if (stepwise2 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="BIC")', sep=""))
			doItAndPrint("summary(res)")
			if (wald==1) doItAndPrint("waldtest(res)")
			doItAndPrint("remove(res)")
		}
		if (stepwise3 == 1){
			subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
				|| trim.blanks(subset) == ""){
				subset <- ""
			}
			else{
				subset <- paste(", subset='", trim.blanks(subset), "'", sep="")
			}
			doItAndPrint(paste('step.p.lm(', modelValue, ', "TempDF", wald=', wald, subset, ")", sep=""))
		}
		if (actmodel==1) activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="lm", model=TRUE, apply="StatMedLinearRegression", reset="StatMedLinearRegression")
    tkgrid(labelRcmdr(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
    tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Wald test for overall p-value for factors with >2 levels")), waldCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on AIC")), stepwise1CheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on BIC")), stepwise2CheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on p-value")), stepwise3CheckBox, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, stick="w")
    tkgrid.configure(helpButton, sticky="e")
    dialogSuffix(rows=4, columns=1)
    }

	
StatMedMannW <- function(){
defaults <- list(group=NULL, response=NULL, alternative="two.sided", test="default", subset = "")
dialog.values <- getDialog("StatMedMannW", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Mann-Whitney U test"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variables with two levels (pick at least one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
    optionsFrame <- tkframe(top)
    StatMedSubsetBox(model=TRUE)   
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Mann-Whitney U test"), "#####", sep=""))
        group <- getSelection(groupBox)
        response <- getSelection(responseBox)
        alternative <- as.character(tclvalue(alternativeVariable))
        test <- as.character(tclvalue(testVariable))
        subset <- tclvalue(subsetVariable)
#        subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) ""
#            else paste(", subset=", subset, sep="")
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}			
			
putDialog("StatMedMannW", list(group=group, response=response, alternative=alternative, test=test, subset = tclvalue(subsetVariable)))
        if (length(group) == 0) {
            errorCondition(recall=StatMedMannW, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups variable."))
            return()
            }
        if (length(response) == 0) {
            errorCondition(recall=StatMedMannW, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
		closeDialog()
        .activeDataSet <- ActiveDataSet()
        nvar = length(group)
#		doItAndPrint("p.value <- NA")
#		doItAndPrint("groups <- NA")
		
		doItAndPrint("group.names <- NULL")
		doItAndPrint("group.median <- NULL")
		doItAndPrint("group.min <- NULL")
		doItAndPrint("group.max <- NULL")
		doItAndPrint("group.1Q <- NULL")
		doItAndPrint("group.3Q <- NULL")
		doItAndPrint("group.p <- NULL")
	
		if(eval(parse(text=paste('"res" %in% objects()')))) doItAndPrint("remove(res)")
	
		for (i in 1:nvar) {
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			command <- (paste("boxplot(", response, "~ factor(", group[i], '), ylab="', response,
                '", xlab="', group[i], '"',
                ", data=", ActiveDataSet(), subset, ")", sep=""))
            logger(command)
            justDoIt(command)
			if (test == "default"){
				doItAndPrint(paste("(res <- wilcox.test(", response, " ~ factor(", group[i], '), alternative="', 
				alternative, '", data=', .activeDataSet, subset, "))", sep=""))
				}
			else {
				doItAndPrint(paste("(res <- wilcox.test(", response, " ~ factor(", group[i], "), alternative='", 
				alternative, "', exact=", test=="exact", 
				", correct=", test=="correct",", data=", .activeDataSet, subset, "))", sep=""))
			}
			
#			doItAndPrint(paste("p.value[", i, "] <- signif(res$p.value, digits=3)", sep=""))
#			doItAndPrint(paste("groups[", i, '] <- "', group[i], '"', sep=""))
		
			group.levels <- eval(parse(text=paste("levels(factor(", subset1, ActiveDataSet(), subset2, "$", group[i], "))", sep="")))
			if(length(group.levels)!=2) next
			for (j in 1:2){
				doItAndPrint(paste('group.names <- c(group.names, "', group[i], "=", group.levels[j], '")', sep=""))
				doItAndPrint(paste("group.min <- c(group.min, with(", subset1, ActiveDataSet(), subset2, ", min(", response, "[", group[i], "=='", group.levels[j], "'], na.rm=TRUE)))", sep=""))
				doItAndPrint(paste("group.1Q <- c(group.1Q, with(", subset1, ActiveDataSet(), subset2, ", quantile(", response, "[", group[i], "=='", group.levels[j], "'], 0.25, na.rm=TRUE)))", sep=""))
				doItAndPrint(paste("group.median <- c(group.median, with(", subset1, ActiveDataSet(), subset2, ", median(", response, "[", group[i], "=='", group.levels[j], "'], na.rm=TRUE)))", sep=""))
				doItAndPrint(paste("group.3Q <- c(group.3Q, with(", subset1, ActiveDataSet(), subset2, ", quantile(", response, "[", group[i], "=='", group.levels[j], "'], 0.75, na.rm=TRUE)))", sep=""))
				doItAndPrint(paste("group.max <- c(group.max, with(", subset1, ActiveDataSet(), subset2, ", max(", response, "[", group[i], "=='", group.levels[j], "'], na.rm=TRUE)))", sep=""))
				
				if (j == 1){
					doItAndPrint("group.p <- c(group.p, signif(res$p.value,digits=3))")
				} else {
					doItAndPrint('group.p <- c(group.p, "")')	
				}
			}
				doItAndPrint("remove(res)")
			}
#		doItAndPrint("mannwhitney.table <- data.frame(p.value=p.value)")
#		doItAndPrint('colnames(mannwhitney.table) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(mannwhitney.table))')
#		doItAndPrint("rownames(mannwhitney.table) <- groups")
#		doItAndPrint("mannwhitney.table")
		
		doItAndPrint("mannwhitney.table <- data.frame(Minimum=group.min, Q1=group.1Q, Median=group.median, Q3=group.3Q, Maximum=group.max, p.value=group.p)")
		doItAndPrint("rownames(mannwhitney.table) <- group.names")
		doItAndPrint('colnames(mannwhitney.table)[c(2,4)] <- c("25%", "75%")') 
		doItAndPrint('colnames(mannwhitney.table) <- gettext(domain="R-RcmdrPlugin.EZR",colnames(mannwhitney.table))')
		doItAndPrint("mannwhitney.table")			
		
		tkfocus(CommanderWindow())
		}
    OKCancelHelp(helpSubject="wilcox.test", apply="StatMedMannW", reset="StatMedMannW")
    radioButtons(optionsFrame, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Difference < 0", "Difference > 0")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))
    radioButtons(optionsFrame, name="test", buttons=c("default", "exact", "normal", "correct"), 
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Default", "Exact", "Normal approximation", "Normal approximation with\ncontinuity correction")), initialValue=dialog.values$test,
        title=gettext(domain="R-RcmdrPlugin.EZR","Type of Test"))
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
     tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), testFrame, sticky="nw")
    tkgrid(optionsFrame, sticky="nw")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
}    

	
StatMedWilSign <- function(){
defaults <- list(x=NULL, y=NULL, alternative="two.sided", test="default", subset = "")
dialog.values <- getDialog("StatMedWilSign", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Wilcoxon's signed rank test"))
    .numeric <- Numeric()
    variablesFrame <- tkframe(top)
    xBox <- variableListBox(variablesFrame, .numeric, title=gettext(domain="R-RcmdrPlugin.EZR","First variable (pick one)"), listHeight=12, initialSelection=varPosn(dialog.values$x, "numeric"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettext(domain="R-RcmdrPlugin.EZR","Second variable (pick one)"), listHeight=12, initialSelection=varPosn(dialog.values$y, "numeric"))
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Wilcoxon's signed rank test"), "#####", sep=""))
        x <- getSelection(xBox)
        y <- getSelection(yBox)
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
        alternative <- as.character(tclvalue(alternativeVariable))
        test <- as.character(tclvalue(testVariable))
putDialog("StatMedWilSign", list(x=x, y=y, alternative=alternative, test=test, subset = tclvalue(subsetVariable)))
		closeDialog()
        if (length(x) == 0 | length(y) == 0) {
            errorCondition(recall=StatMedWilSign, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            return()
            }
        if (x == y) {
            errorCondition(recall=StatMedWilSign, message=gettext(domain="R-RcmdrPlugin.EZR","The two variables must be different."))
            return()
            }
        .activeDataSet <- ActiveDataSet()
        doItAndPrint(paste("median(", subset1, .activeDataSet, subset2, "$", x, " - ", subset1, .activeDataSet, subset2, "$", y, 
            ", na.rm=TRUE) # median difference", sep=""))
        if (test == "default"){
             doItAndPrint(paste("(res <- wilcox.test(", subset1, .activeDataSet, subset2, "$", x, ", ", 
                subset1, .activeDataSet, subset2, "$", y,
                ", alternative='", alternative,
                "', paired=TRUE))", sep=""))           
            }
        else if (test == "exact"){
            doItAndPrint(paste("(res <- wilcox.test(", subset1, .activeDataSet, subset2, "$", x, ", ", 
                subset1, .activeDataSet, subset2, "$", y,
                ", alternative='", alternative,
                "', exact=TRUE, paired=TRUE))", sep=""))
                }
        else {
            doItAndPrint(paste("(res <- wilcox.test(", subset1, .activeDataSet, subset2, "$", x, ", ", 
                subset1, .activeDataSet, subset2, "$", y,
                ", alternative='", alternative, "', correct=", test=="correct",
                ", exact=FALSE, paired=TRUE))", sep=""))
                }
		command <- paste('cat(gettext(domain="R-RcmdrPlugin.EZR", "Wilcoxon', "'", 's signed rank test")', ', "', gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), '\n")', sep="")
		doItAndPrint(command)
		doItAndPrint("remove(res)")
								
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="wilcox.test", apply="StatMedWilSign", reset="StatMedWilSign")
    radioButtons(name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative, 
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Difference < 0", "Difference > 0")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))
    radioButtons(name="test", buttons=c("default", "exact", "normal", "correct"), 
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Default", "Exact", "Normal approximation", "Normal approximation with\ncontinuity correction")), initialValue=dialog.values$test,
        title=gettext(domain="R-RcmdrPlugin.EZR","Type of Test"))
	tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(alternativeFrame, sticky="nw")
    tkgrid(testFrame, sticky="nw")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }
    
	
StatMedKruWalli <- function(){
defaults <- list(group=NULL, response=NULL, steeldwass=0, steel=0, bonferroni=0, holm=0, subset = "")
dialog.values <- getDialog("StatMedKruWalli", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Kruskal-Wallis test"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Groups (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
	optionsFrame <- tkframe(top)
checkBoxes(frame="optionsFrame", boxes=c("bonferroni", "holm", "steeldwass", "steel"), initialValues=c(dialog.values$bonferroni, dialog.values$holm, dialog.values$steeldwass, dialog.values$steel),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Pairwise comparison (Bonferroni)", "Pairwise comparison (Holm)", "Pairwise comparison (Steel-Dwass)", "Pairwise comparison (Steel)")))	
#	steeldwassVariable <- tclVar("0")
#	steeldwassCheckBox <- tkcheckbutton(optionsFrame, variable=steeldwassVariable)
#	steelVariable <- tclVar("0")
#	steelCheckBox <- tkcheckbutton(optionsFrame, variable=steelVariable)
#	bonferroniVariable <- tclVar("0")
#	bonferroniCheckBox <- tkcheckbutton(optionsFrame, variable=bonferroniVariable)
#	holmVariable <- tclVar("0")
#	holmCheckBox <- tkcheckbutton(optionsFrame, variable=holmVariable)
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Kruskal-Wallis test"), "#####", sep=""))
        group <- getSelection(groupBox)
        response <- getSelection(responseBox)
		steeldwass <- tclvalue(steeldwassVariable)
		steel <- tclvalue(steelVariable)
		bonferroni <- tclvalue(bonferroniVariable)
		holm <- tclvalue(holmVariable)
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}			
putDialog("StatMedKruWalli", list(group=group, response=response, steeldwass=steeldwass, steel=steel, bonferroni=bonferroni, holm=holm, subset = tclvalue(subsetVariable)))			
        if (length(group) == 0) {
            errorCondition(recall=StatMedKruWalli, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups variable."))
            return()
            }
        closeDialog()
        if (length(response) == 0) {
            errorCondition(recall=StatMedKruWalli, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a response variable."))
            return()
            }
        .activeDataSet <- ActiveDataSet()
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- (paste("boxplot(", response, "~ factor(", group, '), ylab="', response,
                '", xlab="', group, '"',
                ", data=", ActiveDataSet(), subset, ")", sep=""))
        doItAndPrint(command)
        doItAndPrint(paste("tapply(", paste(subset1, .activeDataSet, subset2, "$", response, sep=""),
            ", ", paste(subset1, .activeDataSet, subset2, "$", group, sep=""), ", median, na.rm=TRUE)", sep=""))
        doItAndPrint(paste("(res <- kruskal.test(", response, " ~ factor(", group, "), data=",
            .activeDataSet, subset, "))", sep=""))
			
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "Kruskal-Wallis test"), " ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")
			
		if (bonferroni==1){
			doItAndPrint(paste("pairwise.kruskal.test(", subset1, .activeDataSet, subset2, "$", response, ", ", subset1, .activeDataSet, subset2, "$", group, ', data.name="', .activeDataSet, '", p.adjust.method="bonferroni")', sep=""))
		}
		if (holm==1){
			doItAndPrint(paste("pairwise.kruskal.test(", subset1, .activeDataSet, subset2, "$", response, ", ", subset1, .activeDataSet, subset2, "$", group, ', data.name="', .activeDataSet, '", p.adjust.method="holm")', sep=""))
		}
		if (steeldwass==1){
			command <- paste("Steel.Dwass(", subset1, .activeDataSet, subset2, "$", response, ", ", subset1, .activeDataSet, subset2, "$", group, ")", sep="")
			doItAndPrint(command)
		}		
		if (steel==1){
			Library("mvtnorm")
			command <- paste("Steel(", subset1, .activeDataSet, subset2, "$", response, ", ", subset1, .activeDataSet, subset2, "$", group, ")", sep="")
			doItAndPrint(command)
		}		
		tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="kruskal.test", apply="StatMedKruWalli", reset="StatMedKruWalli")
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison"), fg="blue"), sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Bonferroni)")), bonferroniCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Holm)")), holmCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Steel-Dwass)")), steeldwassCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Steel)")), steelCheckBox, sticky="w")
	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","The first group in alphabetical will be treated as the reference group."), fg="blue"), sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=2, columns=2)
    }
	
	
StatMedFriedman <- function(){
defaults <- list(response=NULL, bonferroni=0, holm=0, subset = "")
dialog.values <- getDialog("StatMedFriedman", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Friedman test"))
	responseBox <- variableListBox(top, Numeric(), selectmode="multiple", 
			title=gettext(domain="R-RcmdrPlugin.EZR","Repeated-Measures Variables (pick two or more)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
	optionsFrame <- tkframe(top)
checkBoxes(frame="optionsFrame", boxes=c("bonferroni", "holm"), initialValues=c(dialog.values$bonferroni, dialog.values$holm),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Pairwise comparison (Bonferroni)", "Pairwise comparison (Holm)")))	
	
#	bonferroniVariable <- tclVar("0")
#	bonferroniCheckBox <- tkcheckbutton(optionsFrame, variable=bonferroniVariable)
#	holmVariable <- tclVar("0")
#	holmCheckBox <- tkcheckbutton(optionsFrame, variable=holmVariable)
    StatMedSubsetBox(model=TRUE)   
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Friedman test"), "#####", sep=""))
		responses <- getSelection(responseBox)
		bonferroni <- tclvalue(bonferroniVariable)
		holm <- tclvalue(holmVariable)
		    subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
			}
putDialog("StatMedFriedman", list(response=responses, bonferroni=bonferroni, holm=holm, subset = tclvalue(subsetVariable)))	
			closeDialog()
		if (length(responses) < 2) {
			errorCondition(recall=StatMedFriedman, message=gettext(domain="R-RcmdrPlugin.EZR","You must select at least two variables."))
			return()
			}
		.activeDataSet <- ActiveDataSet()
		command <- paste('na.omit(with(', subset1, .activeDataSet, subset2, ', cbind(', paste(responses, collapse=", "), ')))', sep="")
#		logger(paste(".Responses <- ", command, sep=""))
#		assign(".Responses", justDoIt(command), envir=.GlobalEnv)
		doItAndPrint(paste(".Responses <- ", command, sep=""))
		doItAndPrint("apply(.Responses, 2, median)")
		doItAndPrint("(res <- friedman.test(.Responses))")
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "Friedman test"), " ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")
		if (bonferroni==1){
			doItAndPrint(paste('pairwise.friedman.test(.Responses, "', .activeDataSet, '", p.adjust.method="bonferroni")', sep=""))
		}
		if (holm==1){
			doItAndPrint(paste('pairwise.friedman.test(.Responses, "', .activeDataSet, '", p.adjust.method="holm")', sep=""))
		}		
		logger("remove(.Responses)")
		remove(.Responses, envir=.GlobalEnv)
		tkfocus(CommanderWindow())
		}
	OKCancelHelp(helpSubject="friedman.test", apply="StatMedFriedman", reset="StatMedFriedman")
	tkgrid(getFrame(responseBox), sticky="nw")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison"), fg="blue"), sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Bonferroni)")), bonferroniCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Holm)")), holmCheckBox, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=2, columns=1)
	}
	

StatMedJT <- function(){
    Library("clinfun")
defaults <- list(response=NULL, group=NULL, alternative="two.sided", subset = "")
dialog.values <- getDialog("StatMedJT", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Jonckheere-Terpstra test"))
    variablesFrame <- tkframe(top)
   .factors <- Variables()
    responseBox <- variableListBox(variablesFrame, Numeric(), title=gettext(domain="R-RcmdrPlugin.EZR","Response Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "numeric"))
    groupBox <- variableListBox(variablesFrame, .factors, title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Jonckheere-Terpstra test"), "#####", sep=""))
        response <- getSelection(responseBox)
        group <- getSelection(groupBox)
        alternative <- as.character(tclvalue(alternativeVariable))
		.activeDataSet <- ActiveDataSet()
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			.subDataSet <- .activeDataSet
		} else {
			.subDataSet <- paste("subset(", .activeDataSet, ", ", subset, ")", sep="")
		}
putDialog("StatMedJT", list(response=response, group=group, alternative=alternative, subset = tclvalue(subsetVariable)))
		closeDialog()
        if (length(response) == 0 || length(group) == 0){
            errorCondition(recall=StatMedJT, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            return()
            }
        if (response == group) {
            	errorCondition(recall=StatMedJT, message=gettext(domain="R-RcmdrPlugin.EZR","Objective variable and grouping variable must be different."))
 	           return()
        } 
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- (paste("boxplot(", response, "~ factor(", group, '), ylab="', response,
                '", xlab="', group, '"',
                ", data=", .subDataSet, ")", sep=""))
        logger(command)
        justDoIt(command)
		command <- paste("(res <- jonckheere.test(", .subDataSet, "$", response, ", as.ordered(", .subDataSet, "$", group, '), alternative="', alternative, '"))', sep="")
        doItAndPrint(command)
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "Jonckheere-Terpstra test"), " ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")
		
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="jonckheere.test", apply="StatMedJT", reset="StatMedJT")
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Trend will be evaluated among groups in alphabetical order."), fg="blue"), sticky="w")
    optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="alternative", buttons=c("two", "inc", "dec"), values=c("two.sided", "increasing", "decreasing"), initialValue=dialog.values$alternative, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Increasing", "Decreasing")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))
    tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), sticky="nw")
    tkgrid(optionsFrame, sticky="nw")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
}

		
StatMedSpearman <- function(){
defaults <- list(x=NULL, alternative="two.sided", method="spearman", subset = "")
dialog.values <- getDialog("StatMedSpearman", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Spearman's rank correlation test"))
  xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick two)"), listHeight=15, initialSelection=varPosn(dialog.values$x, "numeric"))
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name="alternative", buttons=c("two.sided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "Correlation < 0", "Correlation > 0")), title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))  
  radioButtons(optionsFrame, name="method", buttons=c("Spearman", "Kendall"), values=c("spearman", "kendall"), initialValue=dialog.values$method,
               labels=gettext(domain="R-RcmdrPlugin.EZR",c("Spearman", "Kendall")), title=gettext(domain="R-RcmdrPlugin.EZR","Method"))  
    StatMedSubsetBox(model=TRUE)   
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Spearman's rank correlation test"), "#####", sep=""))
	alternative <- as.character(tclvalue(alternativeVariable))
    method <- as.character(tclvalue(methodVariable))
    x <- getSelection(xBox)
		    subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}
putDialog("StatMedSpearman", list(x=x, alternative=alternative, method=method, subset = tclvalue(subsetVariable)))			
    if (2 > length(x)) {
      errorCondition(recall=StatMedSpearman,
        message=gettext(domain="R-RcmdrPlugin.EZR","Fewer than 2 variables selected."))
      return()
    }
    if(2 < length(x)) {
      errorCondition(recall=StatMedSpearman,
        message=gettext(domain="R-RcmdrPlugin.EZR","More than 2 variables selected."))
      return()
    }
    closeDialog()
    .activeDataSet <- ActiveDataSet()
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
    command2 <- paste("scatterplot(", x[1], "~", x[2],
        ", reg.line=lm, smooth=FALSE, spread=FALSE, boxplots='xy', span=0.5, data=", .activeDataSet, subset, ")", sep="")
    doItAndPrint(command2)  
    command <- paste("(res <- cor.test(", subset1, .activeDataSet, subset2, "$", x[1], ", ", subset1, .activeDataSet, subset2, "$", x[2],
        ', alternative="', alternative, '", method="', method, '"))', sep="")
    doItAndPrint(command)  
	
	command <- paste('cat(gettext(domain="R-RcmdrPlugin.EZR", "Spearman', "'", 's rank correlation coefficient")', ', signif(res$estimate, digits=3), gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n")', sep="")
	doItAndPrint(command)
	doItAndPrint("remove(res)")

    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cor.test", apply="StatMedSpearman", reset="StatMedSpearman")
  tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
  tkgrid(getFrame(xBox), sticky="nw")
     tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), methodFrame, sticky="nw")
    tkgrid(optionsFrame, sticky="nw")
	tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame,columnspan=2,sticky="w")
  dialogSuffix(rows=4, columns=1)
}


StatMedFrequency <- function(){
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Frequency Distributions"))
    xBox <- variableListBox(top, Variables(), selectmode="multiple",
        title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick one or more)"), listHeight=15)
    optionsFrame <- tkframe(top)
    goodnessOfFitVariable <- tclVar("0")
    goodnessOfFitCheckBox <- tkcheckbutton(optionsFrame, variable=goodnessOfFitVariable)
    options2Frame <- tkframe(top)
    shownaVariable <- tclVar("1")
    shownaCheckBox <- tkcheckbutton(options2Frame, variable=shownaVariable)
    options3Frame <- tkframe(top)
    percentVariable <- tclVar("0")
    percentCheckBox <- tkcheckbutton(options3Frame, variable=percentVariable)
    options4Frame <- tkframe(top)
    graphVariable <- tclVar("0")
    graphCheckBox <- tkcheckbutton(options4Frame, variable=graphVariable)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Frequency Distributions"), "#####", sep=""))
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=StatMedFrequency, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            return()
            }
        goodnessOfFit <- tclvalue(goodnessOfFitVariable)
        if (length(x) > 1 && goodnessOfFit == "1"){
            errorCondition(recall=StatMedFrequency, 
                message=gettext(domain="R-RcmdrPlugin.EZR","Goodness-of-fit test not available when more than one variable is selected."))
            return()
            }
        showna <- tclvalue(shownaVariable)
		if (showna == 0){
			showna <- ""
		} else {
			showna <- ", exclude=NULL"
		}
        percent <- tclvalue(percentVariable)
        graph <- tclvalue(graphVariable)
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        for (variable in x){
            command <- paste("table(", .activeDataSet, "$", variable, showna, ")", sep="")
            doItAndPrint(paste("(.Table <- ", command, gettext(domain="R-RcmdrPlugin.EZR",")  # counts for "), variable, sep=""))
#            assign(".Table", justDoIt(command), envir=.GlobalEnv)
#            doItAndPrint(paste(".Table  # counts for", variable))
            if (percent==1) doItAndPrint(paste("round(100*.Table/sum(.Table), 2)  # percentages for", " ", variable))
			if (graph==1) {
					if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
					command <- paste('barplot(.Table, xlab="', variable, '", ylab="Frequency", axis.lty=1)', sep="")
					doItAndPrint(command)
					}
            }
        env <- environment()
        if (goodnessOfFit == 1){
            initializeDialog(subwin, title=gettext(domain="R-RcmdrPlugin.EZR","Goodness-of-Fit Test"))
            hypothesisFrame <- tkframe(subwin)
            levs <- eval(parse(text=paste("levels(", .activeDataSet, "$", x, ")", sep="")))
            n.levs <- length(levs)
            assign(".entry.1", tclVar(paste("1/", n.levs, sep="")), envir=env)
            make.entries <- "labelRcmdr(hypothesisFrame, text='Hypothesized probabilities:   ')"
            make.lev.names <- "labelRcmdr(hypothesisFrame, text='Factor levels:')"
            for (i in 1:n.levs) {
                entry.varname <- paste(".entry.", i, sep="")
                assign(entry.varname, tclVar(paste("1/", n.levs, sep="")), envir=env)
                make.entries <- paste(make.entries, ", ", "ttkentry(hypothesisFrame, width='5', textvariable=", 
                        entry.varname, ")", sep="")
                make.lev.names <- paste(make.lev.names, ", labelRcmdr(hypothesisFrame, text='", levs[i], "')", sep="")
                }
            eval(parse(text=paste("tkgrid(", make.lev.names, ", sticky='w')", sep="")), envir=env)
            eval(parse(text=paste("tkgrid(", make.entries, ", stick='w')", sep="")), envir=env)
            tkgrid(hypothesisFrame, sticky="w")
            onOKsub <- function(){
                probs <- rep(NA, n.levs)
                for (i in 1:n.levs){
                    entry.varname <- paste(".entry.", i, sep="")
                    res <- try(
                        entry <- eval(parse(text=eval(parse(text=paste("tclvalue(", entry.varname,")", sep="")), envir=env))),
                        silent=TRUE)
                    if (class(res) == "try-error"){
                        errorCondition(subwin, message=gettext(domain="R-RcmdrPlugin.EZR","Invalid entry."))
                        return()
                        }
                    if (length(entry) == 0){
                        errorCondition(subwin, message=gettext(domain="R-RcmdrPlugin.EZR","Missing entry."))
                        return()
                        }
                    opts <- options(warn=-1)
                    probs[i] <- as.numeric(entry)
                    options(opts)
                    }
                probs <- na.omit(probs)
                if (length(probs) != n.levs){
                    errorCondition(subwin, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR","Number of valid entries (%d)\nnot equal to number levels (%d)."), length(probs), 
                        n.levs))
                    return()
                    }
                if (any(probs < 0)){
                    errorCondition(subwin, message=gettext(domain="R-RcmdrPlugin.EZR","Negative probabilities not allowed."))
                    return()
                    }
                if (abs(sum(probs) - 1) > 0.001){
                    Message(message=gettext(domain="R-RcmdrPlugin.EZR","Probabilities rescaled to sum to 1."), type="warning")
                    probs <- probs/sum(probs)
                    }
                closeDialog(subwin)
                command <- paste("c(", paste(probs, collapse=","), ")", sep="")
#                logger(paste(".Probs <-", command))
#                assign(".Probs", justDoIt(command), envir=.GlobalEnv)
                doItAndPrint(paste(".Probs <-", command))
                doItAndPrint("chisq.test(.Table, p=.Probs)")
                logger("remove(.Probs)")
                remove(.Probs, envir=.GlobalEnv)
                }
            subOKCancelHelp(subwin)
            tkgrid(subButtonsFrame, sticky="w")
            dialogSuffix(subwin, rows=2, columns=1, onOK=onOKsub, focus=subwin)
            }            
        logger("remove(.Table)") 
        remove(.Table, envir=.GlobalEnv)  
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="table")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
    tkgrid(getFrame(xBox), sticky="nw")    
    tkgrid(labelRcmdr(optionsFrame, 
        text=gettext(domain="R-RcmdrPlugin.EZR","Chi-square goodness-of-fit test (for one variable only)")), 
            goodnessOfFitCheckBox, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(labelRcmdr(options2Frame, 
        text=gettext(domain="R-RcmdrPlugin.EZR","Show missing data")), 
            shownaCheckBox, sticky="w")
    tkgrid(options2Frame, sticky="w")
    tkgrid(labelRcmdr(options3Frame, 
        text=gettext(domain="R-RcmdrPlugin.EZR","Show percent")), 
            percentCheckBox, sticky="w")
    tkgrid(options3Frame, sticky="w")
    tkgrid(labelRcmdr(options4Frame, 
        text=gettext(domain="R-RcmdrPlugin.EZR","Show graph")), 
            graphCheckBox, sticky="w")
    tkgrid(options4Frame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }


StatMedProbCI <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a proportion"))
	variableFrame <- tkframe(top)
	sample <- tclVar("")
	sampleEntry <- ttkentry(variableFrame, width="20", textvariable=sample)
	event <- tclVar("")
	eventEntry <- ttkentry(variableFrame, width="20", textvariable=event)
	CI <- tclVar("95")
	CIEntry <- ttkentry(variableFrame, width="20", textvariable=CI)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a proportion"), "#####", sep=""))
		sample <- tclvalue(sample)
		event <- tclvalue(event)
		CI <- tclvalue(CI)
		closeDialog()
	if (length(sample) == 0 || length(event) == 0){
			errorCondition(recall=StatMedProbCI, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		doItAndPrint(paste("prop.conf(", event, ", ", sample, ", ", CI, ")", sep=""))
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Total number of samples")), sampleEntry, sticky="w")
	tkgrid.configure(sampleEntry, sticky="w")
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of events")), eventEntry, sticky="w")
	tkgrid.configure(eventEntry, sticky="w")
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval")), CIEntry, sticky="w")
	tkgrid.configure(CIEntry, sticky="w")
	tkgrid(variableFrame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

		
StatMedProbSingle <- function(){
defaults <- list(x=NULL, chisq=0, exact=1, continuity="TRUE", alternative="two.sided", p0="0.5", confidence="0.95", subset="")
dialog.values <- getDialog("StatMedProbSingle", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","One sample proportion test"))
    xBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Binary variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$x, "all"))
	StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","One sample proportion test"), "#####", sep=""))
        x <- getSelection(xBox)
		chisq <- tclvalue(chisqTestVariable)
        exact <- tclvalue(exactTestVariable)
		continuity <- tclvalue(continuityVariable)
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceVariable)
        p0 <- tclvalue(p0Variable)
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
putDialog("StatMedProbSingle", list(x=x, chisq=chisq, exact=exact, continuity=continuity, alternative=alternative, p0=p0, confidence=level, subset=tclvalue(subsetVariable)))		
        if (length(x) == 0){
            errorCondition(recall=StatMedProbSingle, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            return()
            }		
		closeDialog()
		doItAndPrint(paste("(.Table <- table(", subset1, ActiveDataSet(), subset2, "$", x, "))", sep="")) 
		if(chisq==1){
			command <- paste("(res <- prop.test(.Table[2], .Table[1]+ .Table[2], p=", p0, ', alternative="', alternative, '", conf.level=', level, ", correct=", continuity, "))", sep="")
			doItAndPrint(command)
		}
		if(exact==1){
			command <- paste("(res <- binom.test(.Table[2], .Table[1]+ .Table[2], p=", p0, ', alternative="', alternative, '", conf.level=', level, "))", sep="")
			doItAndPrint(command)
		}
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "Single-Sample Proportion Test"), " ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")	

		doItAndPrint("remove(.Table)")
			tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="binom.test", apply="StatMedProbSingle", reset="StatMedProbSingle")
    radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"), initialValue=dialog.values$alternative, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Population proportion p!=p0", "Population proportion p<p0", "Population proportion p>p0")),
        title=gettext(domain="R-RcmdrPlugin.EZR","Alternative Hypothesis"))
    rightFrame <- tkframe(top)
	
    confidenceFrame <- tkframe(rightFrame)
    confidenceVariable <- tclVar(dialog.values$confidence)
    confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceVariable)
    p0Frame <- tkframe(rightFrame)
    p0Variable <- tclVar(dialog.values$p0)
    p0Field <- ttkentry(p0Frame, width="6", textvariable=p0Variable)
	
    tkgrid(getFrame(xBox), sticky="nw")
	analysisFrame <- tkframe(top)
	checkBoxes(window=analysisFrame, frame="testsFrame", boxes=c("chisqTest", "exactTest"), 
		initialValues=c(dialog.values$chisq, dialog.values$exact), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Chi-square test", "Exact test")))
    radioButtons(analysisFrame, name="continuity", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue=dialog.values$continuity,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes", "No")), title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction of chi-square test")) 
	tkgrid(testsFrame, labelRcmdr(analysisFrame, text="   "), continuityFrame, sticky="w")
	tkgrid(analysisFrame, sticky="w")
	tkgrid(labelRcmdr(rightFrame, text=""), sticky="w")
    tkgrid(labelRcmdr(p0Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Null hypothesis: p=p0: p0 =")), p0Field, sticky="w")
    tkgrid(p0Frame, sticky="w")
    tkgrid(labelRcmdr(confidenceFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence Level: ")), confidenceField, sticky="w")
    tkgrid(confidenceFrame, sticky="w")
    tkgrid(alternativeFrame, sticky="nw")
    tkgrid(rightFrame, sticky="nw")
    tkgrid.configure(confidenceField, sticky="e")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }


StatMedProbDiffCI <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a difference between two proportions"))
	variableFrame <- tkframe(top)
	sample1 <- tclVar("")
	sample1Entry <- ttkentry(variableFrame, width="20", textvariable=sample1)
	event1 <- tclVar("")
	event1Entry <- ttkentry(variableFrame, width="20", textvariable=event1)
	variable2Frame <- tkframe(top)
	sample2 <- tclVar("")
	sample2Entry <- ttkentry(variable2Frame, width="20", textvariable=sample2)
	event2 <- tclVar("")
	event2Entry <- ttkentry(variable2Frame, width="20", textvariable=event2)
	CI <- tclVar("95")
	CIEntry <- ttkentry(variable2Frame, width="20", textvariable=CI)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a difference between two proportions"), "#####", sep=""))
		sample1 <- tclvalue(sample1)
		event1 <- tclvalue(event1)
		sample2 <- tclvalue(sample2)
		event2 <- tclvalue(event2)
		CI <- tclvalue(CI)
		closeDialog()
	if (length(sample1) == 0 || length(event1) == 0 || length(sample2) == 0 || length(event2) == 0){
			errorCondition(recall=StatMedProbDiffCI, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		doItAndPrint(paste("prop.diff.conf(", event1, ", ", sample1, ", ", event2, ", ", sample2, ", ", CI, ")", sep=""))
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in group 1")), sample1Entry, sticky="w")
	tkgrid.configure(sample1Entry, sticky="w")
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of events in group 1")), event1Entry, sticky="w")
	tkgrid.configure(event1Entry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in group 2")), sample2Entry, sticky="w")
	tkgrid.configure(sample2Entry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of events in group 2")), event2Entry, sticky="w")
	tkgrid.configure(event2Entry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval")), CIEntry, sticky="w")
	tkgrid.configure(CIEntry, sticky="w")
	tkgrid(variableFrame, sticky="nw")
	tkgrid(variable2Frame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedProbRatioCI <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a ratio of two proportions"))
	variableFrame <- tkframe(top)
	sample1 <- tclVar("")
	sample1Entry <- ttkentry(variableFrame, width="20", textvariable=sample1)
	event1 <- tclVar("")
	event1Entry <- ttkentry(variableFrame, width="20", textvariable=event1)
	variable2Frame <- tkframe(top)
	sample2 <- tclVar("")
	sample2Entry <- ttkentry(variable2Frame, width="20", textvariable=sample2)
	event2 <- tclVar("")
	event2Entry <- ttkentry(variable2Frame, width="20", textvariable=event2)
	CI <- tclVar("95")
	CIEntry <- ttkentry(variable2Frame, width="20", textvariable=CI)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Confidence interval for a ratio of two proportions"), "#####", sep=""))
		sample1 <- tclvalue(sample1)
		event1 <- tclvalue(event1)
		sample2 <- tclvalue(sample2)
		event2 <- tclvalue(event2)
		CI <- tclvalue(CI)
		closeDialog()
	if (length(sample1) == 0 || length(event1) == 0 || length(sample2) == 0 || length(event2) == 0){
			errorCondition(recall=StatMedProbRatioCI, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		doItAndPrint(paste("prop.ratio.conf(", event1, ", ", sample1, ", ", event2, ", ", sample2, ", ", CI, ")", sep=""))
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in group 1")), sample1Entry, sticky="w")
	tkgrid.configure(sample1Entry, sticky="w")
	tkgrid(tklabel(variableFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of events in group 1")), event1Entry, sticky="w")
	tkgrid.configure(event1Entry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in group 2")), sample2Entry, sticky="w")
	tkgrid.configure(sample2Entry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of events in group 2")), event2Entry, sticky="w")
	tkgrid.configure(event2Entry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval")), CIEntry, sticky="w")
	tkgrid.configure(CIEntry, sticky="w")
	tkgrid(variableFrame, sticky="nw")
	tkgrid(variable2Frame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

	
StatMedBarGraph <- function(){
defaults <- list(variable=NULL, group=NULL, group2=NULL, color=0, beside=0, percent=0, subset="")
dialog.values <- getDialog("StatMedBarGraph", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Bar graph(Frequencies)"))
	variablesFrame <- tkframe(top)
    variableBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$variable, "all"))
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable1(pick 0 or 1)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    group2Box <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable2(pick 0 or 1)"), listHeight=15, initialSelection=varPosn(dialog.values$group2, "all"))
	
optionsFrame <- tkframe(top)
checkBoxes(frame="optionsFrame", boxes=c("color", "beside", "percent"), initialValues=c(dialog.values$color, dialog.values$beside, dialog.values$percent),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Draw in color", "Show groups side by side", "Compare proportion in each group")))	

#    checkBoxes(frame="color", boxes=c("color"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Draw in color")))
#    checkBoxes(frame="beside", boxes=c("beside"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show groups side by side")))
#    checkBoxes(frame="percent", boxes=c("percent"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Compare proportion in each group")))
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Bar graph(Frequencies)"), "#####", sep=""))
        variable <- getSelection(variableBox)
        group <- getSelection(groupBox)
        group2 <- getSelection(group2Box)
		color <- tclvalue(colorVariable)
		beside <- tclvalue(besideVariable)
		percent <- tclvalue(percentVariable)
		variablemembers <- eval(parse(text=paste("length(levels(factor(", ActiveDataSet(), "$", variable, ")))", sep="")))
		if (color == 0){
			color <- NULL
		} else {
			color <- paste(", col=c(2:", variablemembers+1, ")", sep="")
		}	
		if (beside == 0){
			beside <- NULL
		} else {
			beside <- ", beside=TRUE"
		}	
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
putDialog("StatMedBarGraph", list(variable=variable, group=group, group2=group2, color=tclvalue(colorVariable), beside=tclvalue(besideVariable), percent=percent, subset=tclvalue(subsetVariable)))
        closeDialog()
        if (length(variable) == 0){
            errorCondition(recall=StatMedBarGraph, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
            return()
            }
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		if (length(group) == 0){
				command <- paste("barplot(table(", subset1, ActiveDataSet(), subset2, "$", variable, '), xlab="',
				variable, '", ylab="Frequency"', color, ", axis.lty=1)", sep="")
			}
		else if (length(group2) == 0){
			if(percent == 0){
				command <- paste("barplot(table(", subset1, ActiveDataSet(), subset2, "$", variable, ",", subset1, ActiveDataSet(), subset2, "$", group, '), xlab="',
				group, '", ylab="Frequency"', color, beside, ", legend=levels(factor(", subset1, ActiveDataSet(), subset2, "$", variable, ')), args.legend=list(title="', variable, '", box.lty=0), axis.lty=1)', sep="")
			}
			else{
				command <- paste("barplot(prop.table(table(", subset1, ActiveDataSet(), subset2, "$", variable, ",", subset1, ActiveDataSet(), subset2, "$", group, '),2), xlab="',
				group, '", ylab="Frequency"', color, beside, ", legend=levels(factor(", subset1, ActiveDataSet(), subset2, "$", variable, ')), args.legend=list(title="', variable, '", box.lty=0), axis.lty=1)', sep="")
			}
		} else {
				command <- paste('BarplotFor3Factors(First="', variable, '", Second="', group, '", Third="', group2, '", data="', subset1, ActiveDataSet(), subset2, '", prop=', percent, ", col=", tclvalue(colorVariable), ")", sep="")				
		}
        logger(command)
        justDoIt(command)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="barplot", apply="StatMedBarGraph", reset="StatMedBarGraph")
	tkgrid(getFrame(variableBox), labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","    ")), getFrame(groupBox), labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","    ")), getFrame(group2Box), sticky="w")
	tkgrid(variablesFrame, sticky="w")
#    tkgrid(color, sticky="w")
#    tkgrid(beside, sticky="w")
#    tkgrid(percent, sticky="w")
	tkgrid(optionsFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Side by side graph not created when 2 grouping variables are picked.")), sticky="w")  
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
    }

	
StatMedPieChart <- function(){
defaults <- list(variable=NULL, color=0, subset="")
dialog.values <- getDialog("StatMedPieChart", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	Library("colorspace")
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Pie chart(Frequencies)"))
    variableBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$variable, "all"))
    checkBoxes(frame="color", boxes=c("color"),initialValues=dialog.values$color,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Draw in color")))
    StatMedSubsetBox(model=TRUE)   
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Pie chart(Frequencies)"), "#####", sep=""))
        variable <- getSelection(variableBox)
 		color <- tclvalue(colorVariable)
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
		}
putDialog("StatMedPieChart", list(variable=variable, color=color, subset=tclvalue(subsetVariable)))
		closeDialog()
        if (length(variable) == 0){
            errorCondition(recall=StatMedPieChart, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable"))
            return()
            }
		.activeDataSet <- ActiveDataSet()
 		variablemembers <- eval(parse(text=paste("length(levels(factor(", subset1, .activeDataSet, subset2, "$", variable, ")))", sep="")))
		if (color == 0){
			gray = 0.8
			color <- ", col=(gray(c(0.9"
			if(variablemembers >= 2){
				for (i in 2:variablemembers){
					color <- paste(color, ", ", gray, sep="")
					gray <- gray - 0.1
				}
			}
			color <- paste(color, ")))", sep="")
		} else {
			color <- paste(", col=c(2:", variablemembers+1, ")", sep="")
		}	
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- (paste("pie(table(", subset1, .activeDataSet, subset2, "$", variable, "), labels=levels(factor(",
            .activeDataSet, "$", variable, ')), main="', variable, '"', color, ", clockwise=TRUE)", sep=""))
        logger(command)
        justDoIt(command)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="pie", apply="StatMedPieChart", reset="StatMedPieChart")
    tkgrid(getFrame(variableBox), sticky="nw")
    tkgrid(color, sticky="w")
	tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=3, columns=1)
    }

	
StatMedEnterTable <- function(){
    env <- environment()
    Library("abind")
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Enter and analyze two-way table"))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir=env)
    setUpTable <- function(...){
        tkdestroy(get(".tableFrame", envir=env))
        assign(".tableFrame", tkframe(outerTableFrame), envir=env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        make.col.names <- "labelRcmdr(.tableFrame, text='')"
        for (j in 1:ncols) {
            col.varname <- paste(".colname.", j, sep="")
            assign(col.varname, tclVar(j), envir=env)
            make.col.names <- paste(make.col.names, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                    col.varname, ")", sep="")
            }
        eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
        for (i in 1:nrows){
            varname <- paste(".tab.", i, ".1", sep="")
            assign(varname, tclVar("") , envir=env)
            row.varname <- paste(".rowname.", i, sep="")
            assign(row.varname, tclVar(i), envir=env)
            make.row <- paste("ttkentry(.tableFrame, width='5', textvariable=",
                row.varname, ")", sep="")
            make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                varname, ")", sep="")
            for (j in 2:ncols){
                varname <- paste(".tab.", i, ".", j, sep="")
                assign(varname, tclVar(""), envir=env)
                make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                    varname, ")", sep="")
                }
            eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
            }
        tkgrid(get(".tableFrame", envir=env), sticky="w")
        }
    rowColFrame <- tkframe(top)
    rowsValue <- tclVar("2")
    rowsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=rowsValue,
        resolution=1, orient="horizontal", command=setUpTable)
    rowsShow <- labelRcmdr(rowColFrame, textvariable=rowsValue, width=2, justify="right")
    colsValue <- tclVar("2")
    colsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=colsValue,
        resolution=1, orient="horizontal", command=setUpTable)
    colsShow <- labelRcmdr(rowColFrame, textvariable=colsValue, width=2, justify="right")
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Enter and analyze two-way table"), "#####", sep=""))
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        cell <- 0
        counts <- rep(NA, nrows*ncols)
        row.names <- rep("", nrows)
        col.names <- rep("", ncols)
        for (i in 1:nrows) row.names[i] <-
            eval(parse(text=paste("tclvalue(", paste(".rowname.", i, sep=""),")", sep="")))
        for (j in 1:ncols) col.names[j] <-
            eval(parse(text=paste("tclvalue(", paste(".colname.", j, sep=""),")", sep="")))
        for (i in 1:nrows){
            for (j in 1:ncols){
                cell <- cell+1
                varname <- paste(".tab.", i, ".", j, sep="")
                counts[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep=""))))
                }
            }
        counts <- na.omit(counts)
        if (length(counts) != nrows*ncols){
            errorCondition(recall=StatMedEnterTable, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR","Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), length(counts), nrows, ncols))
            return()
            }
        if (length(unique(row.names)) != nrows){
            errorCondition(recall=StatMedEnterTable, message=gettext(domain="R-RcmdrPlugin.EZR","Row names are not unique."))
            return()
            }
        if (length(unique(col.names)) != ncols){
            errorCondition(recall=StatMedEnterTable, message=gettext(domain="R-RcmdrPlugin.EZR","Column names are not unique."))
            return()
            }
        percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherVariable)
        closeDialog()
        command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", ncols,
            ", byrow=TRUE)", sep="")
#        assign(".Table", justDoIt(command), envir=.GlobalEnv)
#        logger(paste(".Table <- ", command, sep=""))
        doItAndPrint(paste(".Table <- ", command, sep=""))
        command <- paste("c(",paste(paste("'", row.names, "'", sep=""), collapse=", "), ")", sep="")
        justDoIt(paste("rownames(.Table) <- ", command, sep=""))
        logger(paste("rownames(.Table) <- ", command, sep=""))
        command <- paste("c(",paste(paste("'", col.names, "'", sep=""), collapse=", "), ")", sep="")
        justDoIt(paste("colnames(.Table) <- ", command, sep=""))
        logger(paste("colnames(.Table) <- ", command, sep=""))
        doItAndPrint(".Table  # Counts")
        if (percents == "row") doItAndPrint(gettext(domain="R-RcmdrPlugin.EZR","rowPercents(.Table) # Row Percentages"))
        if (percents == "column") doItAndPrint(gettext(domain="R-RcmdrPlugin.EZR","colPercents(.Table) # Column Percentages"))
        if (percents == "total") doItAndPrint(gettext(domain="R-RcmdrPlugin.EZR","totPercents(.Table) # Percentage of Total"))
        if (chisq == 1) {
            command <- "chisq.test(.Table, correct=TRUE)"
#            logger(paste(".Test <- ", command, sep=""))
#            assign(".Test", justDoIt(command), envir=.GlobalEnv)
            doItAndPrint(paste(".Test <- ", command, sep=""))
            doItAndPrint(".Test")
            if (expected == 1) doItAndPrint(".Test$expected # Expected Counts")
            warnText <- NULL
            if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
                gettext(domain="R-RcmdrPlugin.EZR","expected frequencies are less than 1"))
            if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                gettext(domain="R-RcmdrPlugin.EZR"," expected frequencies are less than 5"), sep="")
            if (!is.null(warnText)) Message(message=warnText,
                type="warning")
            if (chisqComp == 1) {
                command <- "round(.Test$residuals^2, 2) # Chi-square Components"
                doItAndPrint(command)
                }
            logger("remove(.Test)")
            remove(.Test, envir=.GlobalEnv)
            }
        if (fisher == 1) doItAndPrint("fisher.test(.Table)")

		if (fisher == 0 & chisq == 1){
            doItAndPrint("res <- chisq.test(.Table, correct=TRUE)")
			} else {
			doItAndPrint("res <- fisher.test(.Table)")
		}		
		doItAndPrint("summary.table <- data.frame(cbind(.Table, p.value=signif(res$p.value, digits=3)))")
		doItAndPrint('summary.table$p.value[2:length(.Table[,1])] <- ""')
		if(fisher == 0 & chisq == 1){
			doItAndPrint('colnames(summary.table)[length(.Table[1,])+1] <- gettext(domain="R-RcmdrPlugin.EZR", "Chisq.p.value")')	
		} else {
			doItAndPrint('colnames(summary.table)[length(.Table[1,])+1] <- gettext(domain="R-RcmdrPlugin.EZR", "Fisher.p.value")')
		}	
		doItAndPrint("remove(res)")	
		doItAndPrint("summary.table")
        logger("remove(.Table)")
        remove(.Table, envir=.GlobalEnv)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="fisher.test")
    radioButtons(name="percents", buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"), values=c("row", "column", "total", "none"),
        initialValue="none", labels=gettext(domain="R-RcmdrPlugin.EZR",c("Row percentages", "Column percentages",  "Percentages of total", "No percentages")), title=gettext(domain="R-RcmdrPlugin.EZR","Compute Percentages"))
    checkBoxes(frame="testsFrame", boxes=c("chisq", "chisqComponents", "expFreq", "fisher"), initialValues=c("0", "0", "0", "1"),
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Chi-square test with continuity correction", "Components of chi-square statistic",
            "Print expected frequencies", "Fisher's exact test")))
    tkgrid(labelRcmdr(rowColFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of Rows:")), rowsSlider, rowsShow, sticky="w")
    tkgrid(labelRcmdr(rowColFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of Columns:")), colsSlider, colsShow, sticky="w")
    tkgrid(rowColFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Enter counts:"), fg="blue"), sticky="w")
    tkgrid(outerTableFrame, sticky="w")
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Hypothesis Tests"), fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=7, columns=2)
    }
	
	
StatMedTwoWayTable <- function(){
#    Library("abind")
defaults <- list(row=NULL, column=NULL, percents="none", chisq=0, chisqComp=0, expected=0, fisher=1, continuity="TRUE", bonferroni=0, holm=0, subset="")
dialog.values <- getDialog("StatMedTwoWayTable", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    Library("abind")
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Create two-way table and compare two proportions (Fisher's exact test)"))
    variablesFrame <- tkframe(top)
   .factors <- Variables()
    rowBox <- variableListBox(variablesFrame, .factors, selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Row variable (pick one or more)"), listHeight=10, initialSelection=varPosn(dialog.values$row, "all"))
    columnBox <- variableListBox(variablesFrame, .factors, title=gettext(domain="R-RcmdrPlugin.EZR","Column variable (pick one)"), listHeight=10, initialSelection=varPosn(dialog.values$column, "all"))
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Create two-way table and compare two proportions (Fisher's exact test)"), "#####", sep=""))
        row <- getSelection(rowBox)
        column <- getSelection(columnBox)
        percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
		continuity <- tclvalue(continuityVariable)
		bonferroni <- tclvalue(bonferroniVariable)
		holm <- tclvalue(holmVariable)
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
			subset1 <- ""
			subset2 <- ""
			subset <- ""
		} else {
			subset1 <- "subset("
			subset2 <- paste(", ", subset, ")", sep="")
			subset <- paste(", subset=", subset, sep="")
		}
putDialog("StatMedTwoWayTable", list(row=row, column=column, percents=percents, chisq=chisq, chisqComp=chisqComp, expected=expected, fisher=fisher, continuity=continuity, bonferroni=bonferroni, holm=holm, subset=tclvalue(subsetVariable)))
        if (length(row) == 0 || length(column) == 0){
            errorCondition(recall=StatMedTwoWayTable, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            return()
            }	
		closeDialog()
	nvar = length(row)
	doItAndPrint("Fisher.summary.table <- NULL")
	for (i in 1:nvar) {
        	if (row[i] == column) {
            	errorCondition(recall=StatMedTwoWayTable, message=gettext(domain="R-RcmdrPlugin.EZR","Row and column variables are the same."))
 	           return()
            	}
		command <- paste("xtabs(~", row[i], "+", column, ", data=", ActiveDataSet(),
            	subset, ")", sep="")
#       	 	logger(paste(".Table <- ", command, sep=""))
#        		assign(".Table", justDoIt(command), envir=.GlobalEnv)
				doItAndPrint(paste(".Table <- ", command, sep=""))
        		doItAndPrint(".Table")
        if (percents == "row") doItAndPrint(gettext(domain="R-RcmdrPlugin.EZR","rowPercents(.Table) # Row Percentages"))
        if (percents == "column") doItAndPrint(gettext(domain="R-RcmdrPlugin.EZR","colPercents(.Table) # Column Percentages"))
        if (percents == "total") doItAndPrint(gettext(domain="R-RcmdrPlugin.EZR","totPercents(.Table) # Percentage of Total"))
        if (chisq == 1) {
            command <- paste("chisq.test(.Table, correct=", continuity, ")", sep="")
#            logger(paste(".Test <- ", command, sep=""))
#            assign(".Test", justDoIt(command), envir=.GlobalEnv)
            doItAndPrint(paste(".Test <- ", command, sep=""))
            doItAndPrint(".Test")
            if (expected == 1) doItAndPrint(".Test$expected # Expected Counts")
            warnText <- NULL
            if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
                gettext(domain="R-RcmdrPlugin.EZR","expected frequencies are less than 1"))
            if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                gettext(domain="R-RcmdrPlugin.EZR"," expected frequencies are less than 5"), sep="")
            if (!is.null(warnText)) Message(message=warnText,
                type="warning")
            if (chisqComp == 1) {
                command <- "round(.Test$residuals^2, 2) # Chi-square Components"
                doItAndPrint(command)
                }
            logger("remove(.Test)")
            remove(.Test, envir=.GlobalEnv)
            }
        if (fisher == 1) doItAndPrint("fisher.test(.Table)")
		if (fisher == 0 & chisq==1){
            doItAndPrint(paste("res <- chisq.test(.Table, correct=", continuity, ")", sep=""))
			} else {
			doItAndPrint("res <- fisher.test(.Table)")
		}
			doItAndPrint("Fisher.summary.table <- rbind(Fisher.summary.table, summary.table.twoway(table=.Table, res=res))")		
		doItAndPrint("remove(res)")	
		}        
	doItAndPrint('colnames(Fisher.summary.table)[length(Fisher.summary.table)] <-  gettext(domain="R-RcmdrPlugin.EZR", colnames(Fisher.summary.table)[length(Fisher.summary.table)])')
	doItAndPrint("Fisher.summary.table")
#	doItAndPrint("remove(Fisher.summary.table)")				
#	logger("remove(.Table)")
	
	if (bonferroni == 1 && nvar == 1){
		doItAndPrint(paste(".Table <- xtabs(~", column, "+", row[i], ", data=", ActiveDataSet(), subset, ")", sep=""))
		if(chisq==1){
			doItAndPrint('pairwise.prop2.test(.Table, p.adj="bonferroni", test.function=chisq.test)')
		}
		if(fisher==1){
			doItAndPrint('pairwise.prop2.test(.Table, p.adj="bonferroni", test.function=fisher.test)')
		}
	}
	if (holm == 1 && nvar == 1){
		doItAndPrint(paste(".Table <- xtabs(~", column, "+", row[i], ", data=", ActiveDataSet(), subset, ")", sep=""))
		if(chisq==1){
			doItAndPrint('pairwise.prop2.test(.Table, p.adj="holm", test.function=chisq.test)')
		}
		if(fisher==1){
			doItAndPrint('pairwise.prop2.test(.Table, p.adj="holm", test.function=fisher.test)')
		}
	}
        remove(.Table, envir=.GlobalEnv)
		logger("remove(.Table)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="xtabs", apply="StatMedTwoWayTable", reset="StatMedTwoWayTable")
    radioButtons(name="percents",
        buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
        values=c("row", "column", "total", "none"), initialValue=dialog.values$percents,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettext(domain="R-RcmdrPlugin.EZR","Compute Percentages"))

analysisFrame <- tkframe(top)
checkBoxes(window=analysisFrame, frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), initialValues=c(dialog.values$chisq, dialog.values$chisqComp, dialog.values$expected, dialog.values$fisher),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Chi-square test", "Components of chi-square statistic","Print expected frequencies", "Fisher's exact test")))

#	checkBoxes(window=analysisFrame, frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), 
#		initialValues=c("0", "0", "0", "1"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Chi-square test", "Components of chi-square statistic",
#           "Print expected frequencies", "Fisher's exact test")))

optionsFrame <- tkframe(top)
tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison not performed when more than one grouping variables are picked."), fg="blue"), sticky="w")
checkBoxes(frame="optionsFrame", boxes=c("bonferroni", "holm"), initialValues=c(dialog.values$bonferroni, dialog.values$holm),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Pairwise comparison (Bonferroni)", "Pairwise comparison (Holm)")))	
#	bonferroniVariable <- tclVar("0")
#	bonferroniCheckBox <- tkcheckbutton(optionsFrame, variable=dialog.values$bonferroniVariable)
#	holmVariable <- tclVar("0")
#	holmCheckBox <- tkcheckbutton(optionsFrame, variable=dialog.values$holmVariable)

    radioButtons(analysisFrame, name="continuity",
        buttons=c("yes", "no"),
        values=c("TRUE", "FALSE"), initialValue=dialog.values$continuity,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes", "No")), title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction of chi-square test"))

#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Bonferroni)")), bonferroniCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison (Holm)")), holmCheckBox, sticky="w")
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
    tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Hypothesis Tests"), fg="blue"), sticky="w")
	tkgrid(testsFrame, labelRcmdr(analysisFrame, text="   "), continuityFrame, sticky="w")
	tkgrid(analysisFrame, sticky="nw")
	options2Frame <- tkframe(top)
	tkgrid(labelRcmdr(options2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Pairwise comparison not performed when more than one grouping variables are picked."), fg="blue"), sticky="w")
	tkgrid(options2Frame, sticky="nw")
	tkgrid(optionsFrame, sticky="nw")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
}

	
StatMedMcNemar <- function(){
defaults <- list(row=NULL, column=NULL, continuity="TRUE", subset="")
dialog.values <- getDialog("StatMedMcNemar", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

#    Library("abind")
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Compare proportions of two paired samples (McNemar test)"))
    variablesFrame <- tkframe(top)
   .factors <- Variables()
    rowBox <- variableListBox(variablesFrame, .factors, title=gettext(domain="R-RcmdrPlugin.EZR","Row variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$row, "all"))
    columnBox <- variableListBox(variablesFrame, .factors, title=gettext(domain="R-RcmdrPlugin.EZR","Column variable (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$column, "all"))
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Compare proportions of two paired samples (McNemar test)"), "#####", sep=""))
        row <- getSelection(rowBox)
        column <- getSelection(columnBox)
		continuity <- tclvalue(continuityVariable)
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) ""
            else paste(", subset=", subset, sep="")
putDialog("StatMedMcNemar", list(row=row, column=column, continuity=continuity, subset=tclvalue(subsetVariable)))
        if (length(row) == 0 || length(column) == 0){
            errorCondition(recall=StatMedMcNemar, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            return()
            }
		closeDialog()

        if (row == column) {
            	errorCondition(recall=StatMedMcNemar, message=gettext(domain="R-RcmdrPlugin.EZR","Row and column variables are the same."))
 	           return()
        }
		command <- paste("xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
            	subset, ")", sep="")
#       	logger(paste(".Table <- ", command, sep=""))
#        		assign(".Table", justDoIt(command), envir=.GlobalEnv)
       	doItAndPrint(paste(".Table <- ", command, sep=""))
        doItAndPrint(".Table")
		command <- paste("(res <- mcnemar.test(.Table, correct=", continuity, "))", sep="")
        doItAndPrint(command)
		command <- paste('cat(gettext(domain="R-RcmdrPlugin.EZR", "McNemar', "'", 's test")', ', "', gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), '\n")', sep="")
		doItAndPrint(command)
		doItAndPrint("remove(res)")
        remove(.Table, envir=.GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="mcnemar.test", apply="StatMedMcNemar", reset="StatMedMcNemar")
    tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    radioButtons(name="continuity",
        buttons=c("yes", "no"),
        values=c("TRUE", "FALSE"), initialValue=dialog.values$continuity,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes", "No")), title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction"))
    tkgrid(continuityFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
    }

	
StatMedCochranQ <- function(){
defaults <- list(response=NULL, subset="")
dialog.values <- getDialog("StatMedCochranQ", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Compare proportions of more than two paired samples (Cochran Q test)"))
	responseBox <- variableListBox(top, Variables(), selectmode="multiple", 
			title=gettext(domain="R-RcmdrPlugin.EZR","Pick 2 or more paired binary variables"), listHeight=15, initialSelection=varPosn(dialog.values$response, "all"))
    StatMedSubsetBox(model=TRUE)   
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Compare proportions of more than two paired samples (Cochran Q test)"), "#####", sep=""))
		responses <- getSelection(responseBox)
		    subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
			}
putDialog("StatMedCochranQ", list(response=responses, subset=tclvalue(subsetVariable)))
		closeDialog()
		if (length(responses) < 2) {
			errorCondition(recall=StatMedCochranQ, message=gettext(domain="R-RcmdrPlugin.EZR","You must select at least two variables."))
			return()
			}
		.activeDataSet <- ActiveDataSet()
		command <- paste(".Table <- cbind(", subset1, .activeDataSet, subset2, "$", responses[1], sep="")
		for (i in 2:length(responses)){
			command <- paste(command, ", ", subset1, .activeDataSet, subset2, "$", responses[i], sep="")
		}
		command <- paste(command, ")", sep="")
		doItAndPrint(command)
		doItAndPrint("(res <- Cochran.Q.test(.Table))")
		command <- paste('cat(gettext(domain="R-RcmdrPlugin.EZR", "Cochran', "'", 's Q test")', ', "', gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), '\n")', sep="")
		doItAndPrint(command)
		doItAndPrint("remove(res)")		
        doItAndPrint("remove(.Table)")
		tkfocus(CommanderWindow())
		}
	OKCancelHelp(apply="StatMedCochranQ", reset="StatMedCochranQ")
	tkgrid(getFrame(responseBox), sticky="nw")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=2, columns=1)
	}

	
StatMedPropTrend <- function(){
#    Library("abind")
defaults <- list(response=NULL, group=NULL, subset="")
dialog.values <- getDialog("StatMedPropTrend", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Cochran-Armitage test for trend in proportions"))
    variablesFrame <- tkframe(top)
   .factors <- Variables()
    responseBox <- variableListBox(variablesFrame, .factors, title=gettext(domain="R-RcmdrPlugin.EZR","Binary varibale(Ex. No response=0, Response=1) (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "all"))
    groupBox <- variableListBox(variablesFrame, .factors, title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Cochran-Armitage test for trend in proportions"), "#####", sep=""))
        response <- getSelection(responseBox)
        group <- getSelection(groupBox)
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) ""
            else paste(", subset=", subset, sep="")			
putDialog("StatMedPropTrend", list(response=response, group=group, subset=tclvalue(subsetVariable)))
        if (length(response) == 0 || length(group) == 0){
            errorCondition(recall=StatMedPropTrend, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            return()
            }
        closeDialog()
        if (response == group) {
            	errorCondition(recall=StatMedPropTrend, message=gettext(domain="R-RcmdrPlugin.EZR","Binary variable and grouping variable must be different."))
 	           return()
        }
		command <- paste("xtabs(~", group, "+", response, ", data=", ActiveDataSet(),
            	subset, ")", sep="")
#       	logger(paste(".Table <- ", command, sep=""))
#        		assign(".Table", justDoIt(command), envir=.GlobalEnv)
       	doItAndPrint(paste(".Table <- ", command, sep=""))
        doItAndPrint(".Table")		
		command <- "(res <- prop.trend.test(.Table[,1], .Table[,1]+.Table[,2]))"
        doItAndPrint(command)
		doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR", "Cochran-Armitage test for trend in proportions"), " ", gettext(domain="R-RcmdrPlugin.EZR", "p.value"), " = ", signif(res$p.value, digits=3), "\n", sep="")')
		doItAndPrint("remove(res)")
        remove(.Table, envir=.GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="prop.trend.test", apply="StatMedPropTrend", reset="StatMedPropTrend")
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(groupBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Trend will be evaluated among groups in alphabetical order."), fg="blue"), sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
    }


StatMedLogisticRegression <- function(){
defaults <- list(lhs = "", rhs = "", waldVariable = 0,  rocVariable = 0, diagnosisVariable = 0, actmodelVariable = 0, stepwise1Variable = 0, stepwise2Variable = 0, stepwise3Variable = 0, subset = "")
dialog.values <- getDialog("StatMedLogisticRegression", defaults)
currentFields$lhs <- dialog.values$lhs			#Values in currentFields will be sent to modelFormula
currentFields$rhs <- dialog.values$rhs
currentFields$subset <- dialog.values$subset	

    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Logistic regression"))
    .activeModel <- ActiveModel()
    currentModel <- if (!is.null(.activeModel)) #if current model exists, input to modelFormula
        class(get(.activeModel, envir=.GlobalEnv))[1] == "glm"
#        eval(parse(text=paste("class(", .activeModel, ")[1] == 'glm'", sep="")),
#            envir=.GlobalEnv)
        else FALSE
#    if (currentModel) {
#        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv), glm=TRUE)
#        currentFields <- formulaFields(eval(parse(text=.activeModel),
#            envir=.GlobalEnv), glm=TRUE)
#        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
#        }
		
	currentModel <- TRUE
	
    StatMedModelFormula()
    UpdateModelNumber()
    modelName <- tclVar(paste("GLM.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
	optionsFrame <- tkframe(top)
	
	checkBoxes(frame="checkboxFrame", boxes=c("wald", "actmodel", "roc", "diagnosis", "stepwise1", "stepwise2", "stepwise3"), initialValues=c(dialog.values$waldVariable, dialog.values$actmodelVariable, dialog.values$rocVariable, dialog.values$diagnosisVariable, dialog.values$stepwise1Variabl, dialog.values$stepwise2Variabl, dialog.values$stepwise3Variabl),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Wald test for overall p-value for factors with >2 levels", "Keep results as active model for further analyses", "Show ROC curve", "Show basic diagnostic plots", "Stepwise selection based on AIC", "Stepwise selection based on BIC", "Stepwise selection based on p-value")))	

#	waldVariable <- tclVar("0")
#	waldCheckBox <- tkcheckbutton(optionsFrame, variable=waldVariable)
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
#	stepwise1Variable <- tclVar("0")
#	stepwise1CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise1Variable)
#	stepwise2Variable <- tclVar("0")
#	stepwise2CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise2Variable)
#	stepwise3Variable <- tclVar("0")
#	stepwise3CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise3Variable)
    StatMedSubsetBox(model=TRUE)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Logistic regression"), "#####", sep=""))
        modelValue <- trim.blanks(tclvalue(modelName))
        formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
		wald <- tclvalue(waldVariable)
		actmodel <- tclvalue(actmodelVariable)
		roc <- tclvalue(rocVariable)
		diagnosis <- tclvalue(diagnosisVariable)
		stepwise1 <- tclvalue(stepwise1Variable)
        stepwise2 <- tclvalue(stepwise2Variable)
		stepwise3 <- tclvalue(stepwise3Variable)
		subset <- tclvalue(subsetVariable)
#input values into dialog memory	
putDialog("StatMedLogisticRegression", list(lhs = tclvalue(lhsVariable), rhs = tclvalue(rhsVariable), waldVariable = wald,  actmodelVariable = actmodel, rocVariable = roc, diagnosisVariable = diagnosis, stepwise1Variable = stepwise1, stepwise2Variable = stepwise2, stepwise3Variable = stepwise3, subset=tclvalue(subsetVariable)))
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedLogisticRegression, model=TRUE, message=gettext(domain="R-RcmdrPlugin.EZR","Left-hand side of model empty."))
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedLogisticRegression, model=TRUE, message=gettext(domain="R-RcmdrPlugin.EZR","Right-hand side of model empty."))
            return()
            }
		if (!is.valid.name(modelValue)){
            errorCondition(recall=StatMedLogisticRegression, model=TRUE, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue))
            return()
            }
        if (is.element(modelValue, listGeneralizedLinearModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
                UpdateModelNumber(-1)
                closeDialog()
                StatMedLogisticRegression()
                return()
                }
            }
			closeDialog()
        if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>") || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
            }
        else{
            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
            }
		Library("aod")
			command <- paste("glm(", formula, ", family=binomial(logit), data=", ActiveDataSet(), subset, ")", sep="")
#        logger(paste(modelValue, " <- ", command, sep=""))
#        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(paste(modelValue, " <- ", command, sep=""))
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))	
		x <- strsplit(tclvalue(rhsVariable), split="\\+")
		command <- paste("TempDF <- with(", ActiveDataSet(), ", ", ActiveDataSet(), "[complete.cases(", tclvalue(lhsVariable), ", ", paste(x[[1]], collapse=","), "),])", sep="")
		doItAndPrint(command)
		doItAndPrint(paste("GLM.null <- glm(", tclvalue(lhsVariable), "~1, family=binomial(logit), data=TempDF", subset, ")", sep=""))
		doItAndPrint(paste("anova(", modelValue, ', GLM.null, test="Chisq")', sep=""))
        doItAndPrint(paste("vif(", modelValue, ")", sep=""))
		logger("###variance inflation factors")
		doItAndPrint(paste("odds <- data.frame(exp( summary(", modelValue, ")$coef[,1:2] %*% rbind(c(1,1,1), 1.96*c(0,-1,1))))", sep=""))
		doItAndPrint(paste("odds <- cbind(odds, summary(", modelValue, ")$coefficients[,4])", sep=""))
		doItAndPrint("odds <- signif(odds, digits=3)")
		doItAndPrint('names(odds) <- gettext(domain="R-RcmdrPlugin.EZR",c("odds ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
		doItAndPrint("odds")
		if (wald==1) doItAndPrint(paste("waldtest(", modelValue, ")", sep=""))
		if (roc==1){
			Library("pROC")
			doItAndPrint(paste("ROC <- roc(", tclvalue(lhsVariable), " ~ ", modelValue, "$fitted.values, data=TempDF", subset, ', ci=TRUE, direction="auto")', sep=""))
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}			
			doItAndPrint("plot(ROC)")
			doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR","Area under the curve"), signif(ROC$auc[1], digits=3), gettext(domain="R-RcmdrPlugin.EZR","95% CI"), signif(ROC$ci[1], digits=3), "-", signif(ROC$ci[3], digits=3), "\n")')
			doItAndPrint("remove(ROC)")
		}
		if (diagnosis==1){
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}			
			doItAndPrint("oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))")
			doItAndPrint(paste("plot(", modelValue, ")", sep=""))
			doItAndPrint("par(oldpar)")			
		}
		if (stepwise1 == 1 | stepwise2 == 1 | stepwise3 == 1){
			command <- paste("glm(", formula, ", family=binomial(logit), data=TempDF", subset, ")", sep="")
			doItAndPrint(paste(modelValue, " <- ", command, sep=""))
			}
		if (stepwise1 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="AIC")', sep=""))
			doItAndPrint("odds <- data.frame(exp( summary(res)$coef[,1:2] %*% rbind(c(1,1,1), 1.96*c(0,-1,1))))")
			doItAndPrint(paste("odds <- cbind(odds, summary(res)$coefficients[,4])", sep=""))
			doItAndPrint("odds <- signif(odds, digits=3)")
			doItAndPrint('names(odds) <- c("odds ratio", "lower .95", "upper .95", "p.value")')
			doItAndPrint("summary(res)")
			doItAndPrint("odds")
			if (wald==1) doItAndPrint("waldtest(res)")
			doItAndPrint("remove(res)")
		}
		if (stepwise2 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="BIC")', sep=""))
			doItAndPrint("odds <- data.frame(exp( summary(res)$coef[,1:2] %*% rbind(c(1,1,1), 1.96*c(0,-1,1))))")
			doItAndPrint(paste("odds <- cbind(odds, summary(res)$coefficients[,4])", sep=""))
			doItAndPrint("odds <- signif(odds, digits=3)")
			doItAndPrint('names(odds) <- c("odds ratio", "lower .95", "upper .95", "p.value")')
			doItAndPrint("summary(res)")
			doItAndPrint("odds")
			if (wald==1) doItAndPrint("waldtest(res)")
			doItAndPrint("remove(res)")
		}		
		if (stepwise3 == 1){
			subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
				|| trim.blanks(subset) == ""){
				subset <- ""
			}
			else{
				subset <- paste(", subset='", trim.blanks(subset), "'", sep="")
			}
			doItAndPrint(paste('step.p.glm(', modelValue, ', "TempDF", wald=', wald, subset, ")", sep=""))
		}
#		doItAndPrint("remove(odds)")
		if (actmodel==1) activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="generalizedLinearModel", apply="StatMedLogisticRegression", reset="StatMedLogisticRegression")
    helpButton <- buttonRcmdr(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(labelRcmdr(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(getFrame(xBox), sticky="w")
    tkgrid(outerOperatorsFrame, sticky="w")
    tkgrid(formulaFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Stratifing variable: + strata(#####)")), sticky="w")  

 tkgrid(checkboxFrame, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
    }


StatMedKaplanMeier  <- function(){
defaults <- list(event = "", timetoevent = "", group = "", strata = "", test = 0,  line = "color", place = "topright", xscale = "", posthoc = "", censor = 1, ci = 0, separatestrata = 0, atrisk = 0, point = "<none>", xlim = "<auto>", ylim = "<auto>", xlabel = "<auto>", ylabel = "<auto>", subset = "")
dialog.values <- getDialog("StatMedKaplanMeier", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Kaplan-Meier survival curve and logrank test"))
    variablesFrame <- tkframe(top)
    eventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Status indicator (censor=0, event=1) (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$event, "all"))
    timetoeventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time-to-event variable (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$timetoevent, "all"))
    variables2Frame <- tkframe(top)
    groupBox <- variableListBox(variables2Frame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (pick 0, 1, or more)"), listHeight=8, initialSelection=varPosn(dialog.values$group, "all"))
    strataBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Stratifying variable (pick 0 or 1)"), listHeight=8, initialSelection=varPosn(dialog.values$strata, "all"))
    plotoptionFrame <- tkframe(top)
    radioButtons(plotoptionFrame, name="test", buttons=c("logrank", "wilcoxon"), values=c("0", "1"), initialValue=dialog.values$test,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("logrank", "Peto-Peto-Wilcoxon")), title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
    radioButtons(plotoptionFrame, name="line", buttons=c("color", "type", "width"), values=c("color", "type", "width"), initialValue=dialog.values$line,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Color", "Line type", "Line width")), title=gettext(domain="R-RcmdrPlugin.EZR","Line discrimination"))
    radioButtons(plotoptionFrame, name="place", buttons=c("topright", "bottom", "mouse"), values=c("topright", "bottom", "mouse"), initialValue=dialog.values$place,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Upper right", "Bottom", "Mouse click")), title=gettext(domain="R-RcmdrPlugin.EZR","Legend"))
    radioButtons(plotoptionFrame, name="xscale", buttons=c("day", "daytomonth", "daytoyear", "monthtoyear"), values=c("", "30.4375", "365.25", "12"), initialValue=dialog.values$xscale,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("As is", "Day to month", "Day to year", "Month to year")), title=gettext(domain="R-RcmdrPlugin.EZR","X axis"))
    radioButtons(plotoptionFrame, name="posthoc", buttons=c("No", "Bonferroni", "Holm"), values=c("", "bon", "holm"), initialValue=dialog.values$posthoc,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("No", "Bonferroni", "Holm")), title=gettext(domain="R-RcmdrPlugin.EZR","Post-hoc test\n(when only one grouping\nvariable picked)"))
		
	plotoption2Frame <- tkframe(top)
	checkBoxes(window=plotoption2Frame, frame="censor", boxes=c("censor"), initialValues=c(dialog.values$censor),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show censoring marks")), title=gettext(domain="R-RcmdrPlugin.EZR","Options"))	
	checkBoxes(window=plotoption2Frame, frame="ci", boxes=c("ci"), initialValues=c(dialog.values$ci),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show 95% confidence intervals")), title=gettext(domain="R-RcmdrPlugin.EZR"," "))	
	checkBoxes(window=plotoption2Frame, frame="separatestrata", boxes=c("separatestrata"), initialValues=c(dialog.values$separatestrata),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show each strata separately")), title=gettext(domain="R-RcmdrPlugin.EZR"," "))	
	checkBoxes(window=plotoption2Frame, frame="atrisk", boxes=c("atrisk"), initialValues=c(dialog.values$atrisk),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show number at risk")), title=gettext(domain="R-RcmdrPlugin.EZR"," "))	

	#    checkBoxes(window=plotoption2Frame, frame="censor", boxes=c("censor"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show censoring marks")))
#   checkBoxes(window=plotoption2Frame, frame="ci", boxes=c("ci"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show 95% confidence intervals")))
#    checkBoxes(window=plotoption2Frame, frame="separatestrata", boxes=c("separatestrata"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show each strata separately")))
#    checkBoxes(window=plotoption2Frame, frame="atrisk", boxes=c("atrisk"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show number at risk")))
	axisFrame <- tkframe(top)
	axis2Frame <- tkframe(top)
	
	pointFrame <- tkframe(axisFrame)
	pointVariable <- tclVar(dialog.values$point)	
	pointField <- ttkentry(pointFrame, width="20", textvariable=pointVariable)
	xlimFrame <- tkframe(axis2Frame)
	xlimVariable <- tclVar(dialog.values$xlim)
	xlimField <- ttkentry(axis2Frame, width="20", textvariable=xlimVariable)
	ylimFrame <- tkframe(axis2Frame)
	ylimVariable <- tclVar(dialog.values$ylim)
	ylimField <- ttkentry(axis2Frame, width="20", textvariable=ylimVariable)
	xlabelFrame <- tkframe(axis2Frame)
	xlabelVariable <- tclVar(dialog.values$xlabel)
	xlabelField <- ttkentry(axis2Frame, width="20", textvariable=xlabelVariable)
	ylabelFrame <- tkframe(axis2Frame)
	ylabelVariable <- tclVar(dialog.values$ylabel)
	ylabelField <- ttkentry(axis2Frame, width="20", textvariable=ylabelVariable)
	
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Kaplan-Meier survival curve and logrank test"), "#####", sep=""))
    event <- getSelection(eventBox)
    timetoevent <- getSelection(timetoeventBox)
    group <- getSelection(groupBox)
    strata <- getSelection(strataBox)
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
	  sub1 <- ""
	  sub2 <- ""
      subset <- ""
	  }
    else{
	  sub1 <- "subset("
	  sub2 <- paste(", ", subset, ")", sep="")
      subset <- paste(", subset=", subset, sep="")
	  }
	if (length(strata) == 0){
		strata2 <- ""
	}
	else{
		strata2 <- paste("+strata(", strata, ")", sep="")
	}
    censor <- tclvalue(censorVariable)
    ci <- tclvalue(ciVariable)
    separatestrata <- tclvalue(separatestrataVariable)
	if (length(strata) == 0) separatestrata <- 0
    atrisk <- tclvalue(atriskVariable)
	point <- tclvalue(pointVariable)
	if (point == "<none>") {
		point <- ""
	} else {
		point <- paste(", time=", point, sep="")
	}
	test <- as.character(tclvalue(testVariable))
	line <- tclvalue(lineVariable)	
	par.lwd <- get("par.lwd", envir=.GlobalEnv)
	if (line=="color") {line <- paste("col=1:32, lty=1, ", par.lwd, ", ", sep=""); line2 <- paste("col=1:32, lty=1, ", par.lwd, ", ", sep="")}
	if (line=="type") {line <- paste("col=1, lty=1:32, ", par.lwd, ", ", sep=""); line2 <- paste("col=1, lty=1:32, ", par.lwd, ", ", sep="")}
	if (line=="width") {line <- paste("col=1, lty=1, ", par.lwd, ":8, ", sep=""); line2 <- paste("col=1, lty=1, ", par.lwd, ":8, ", sep="")}
	par.cex <- get("par.cex", envir=.GlobalEnv)
	place <- tclvalue(placeVariable)
	if(place=="mouse"){
		place <- "locator(1)"
	}else if (place=="topright"){
		place <- '"topright"'
	}else{
		place <- '"bottom", horiz=TRUE'
	}
	xscale <- tclvalue(xscaleVariable)
	xscale2 <- ""
	if (xscale!=""){
		xscale2 <- paste(" * ", xscale, sep="")
		xscale <- paste(", xscale=", xscale, sep="")
	}
	posthoc <- tclvalue(posthocVariable)
	xlim <- tclvalue(xlimVariable)
	ylim <- tclvalue(ylimVariable)
	if (xlim == "<auto>") {
		xlim <- ""
	} else {
		xlim <- paste(", xlim=c(", xlim, ")", sep="")
	}
	if (ylim == "<auto>") {
		ylim <- ""
	} else {
		ylim <- paste(", ylim=c(", ylim, ")", sep="")
	}
	xlabel <- tclvalue(xlabelVariable)
	ylabel <- tclvalue(ylabelVariable)
	if (xlabel == "<auto>") {
		xlabel <- paste(', xlab="', timetoevent, '"', sep="")
	} else {
		xlabel <- paste(', xlab="', xlabel, '"', sep="")
	}
	if (ylabel == "<auto>") {
		ylabel <- ', ylab="Probability"'
	} else {
		ylabel <- paste(', ylab="', ylabel, '"', sep="")
	}
	if (ci==0){
		conf.int <- "FALSE"
	}else{
		conf.int <- "TRUE"	
		if (line==paste("col=1:32, lty=1, ", par.lwd, ", ", sep="")) line <- paste("col=rep(1:32, each=3), lty=1, ", par.lwd, ", ", sep="")
		if (line==paste("col=1, lty=1:32, ", par.lwd, ", ", sep="")) line <- paste("col=1, lty=rep(1:32, each=3), ", par.lwd, ", ", sep="")
		if (line==paste("col=1, lty=1, ", par.lwd, ":8, ", sep="")) line <- paste("col=1, lty=1, lwd=rep(", substring(par.lwd, nchar(par.lwd),nchar(par.lwd)), ":8, each=3), ", sep="")		
	}
	if (censor==0){
		censor <- ", mark.time=FALSE"
	}else{
		censor <- ", mark.time=TRUE"	
	}
	dataSet <- activeDataSet()	
putDialog("StatMedKaplanMeier", list(event = event, timetoevent = timetoevent, group = group, strata = strata, test = test, line = tclvalue(lineVariable), place = tclvalue(placeVariable), xscale = tclvalue(xscaleVariable), posthoc = posthoc, censor = tclvalue(censorVariable), ci = ci, separatestrata = separatestrata, atrisk = atrisk, point = tclvalue(pointVariable), xlim = tclvalue(xlimVariable), ylim = tclvalue(ylimVariable), xlabel = tclvalue(xlabelVariable), ylabel = tclvalue(ylabelVariable), subset = tclvalue(subsetVariable)))
    if (length(event) != 1) {
      errorCondition(recall=StatMedKaplanMeier, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one status indicator (censor=0, event=1)"))
      return()
    }
    if (length(timetoevent) != 1) {
      errorCondition(recall=StatMedKaplanMeier, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one time-to-event variable"))
      return()
    }
	closeDialog()
    Library("survival")
    nvar <- length(group)
    if (nvar == 0){
	command <- paste("km <- survfit(Surv(", timetoevent, ",", event, "==1)~1, data=", ActiveDataSet(), subset, ', na.action = na.omit, conf.type="log-log")', sep="")
	doItAndPrint(command)
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
	if (atrisk==0){
		doItAndPrint(paste('plot(km, bty="l", ', line, "conf.int=", conf.int, censor, xlim, ylim, xlabel, ylabel, xscale, ")", sep=""))
	} else {
		doItAndPrint('mar <- par("mar")')
        doItAndPrint("mar[1] <- mar[1] + 1 + 0.5")
		doItAndPrint("par(mar=mar)")
        doItAndPrint("opar <- par(mar = mar)")
        doItAndPrint("on.exit(par(opar))")
		doItAndPrint(paste('plot(km, bty="l", ', line, "conf.int=", conf.int, censor, xlim, ylim, xlabel, ylabel, xscale, ")", sep=""))	
        doItAndPrint("xticks <- axTicks(1)")
        doItAndPrint(paste("n.atrisk <- nrisk(km, xticks", xscale2, ")", sep=""))
        doItAndPrint("axis(1, at = xticks, labels = n.atrisk, line = 3, tick = FALSE)")
        doItAndPrint('title(xlab = "Number at risk", line = 3, adj = 0)')
	}
	doItAndPrint("summary(km)")
	doItAndPrint(paste("summary.km(survfit=km", point, xscale2, ")", sep=""))
    } else {
	for (i in 1:nvar) {
	command <- paste("km <- survfit(Surv(", timetoevent, ",", event, "==1)~", group[i], strata2, ", data=", ActiveDataSet(), subset, ', na.action = na.omit, conf.type="log-log")', sep="")
	doItAndPrint(command)
	doItAndPrint("summary(km)")
#	doItAndPrint('legend <- c("0", "1")') #to create a legend vector. "0", "1" are dummy.
	if (length(strata) == 0 || separatestrata == 1){
		strata3 <- ""
		doItAndPrint(paste('len <- nchar("', group[i], '")', sep=""))
#		doItAndPrint("nvar2 <- length(names(km$strata))")	
#		doItAndPrint("k <- 1")
#		doItAndPrint("for (j in 1:nvar2){legend[k] <- levels(factor(substring(names(km$strata), len+2)))[j]; k <- k+1}")
		doItAndPrint("legend <- substring(names(km$strata), len+2)")
		}else{
			#To remove groups with n=0 by interaction() suggested by Dr. Yoshida		
		doItAndPrint(paste("legend <- levels(factor(interaction(", sub1, dataSet, sub2, "$", strata, ", ", 	sub1, dataSet, sub2, "$", group[i], ', sep=":")))', sep=""))
		strata3 <- paste(strata, " : ", sep="")
	}
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
	if (separatestrata == 0){
		if (atrisk==0){
			doItAndPrint(paste('plot(km, bty="l", ', line, "conf.int=", conf.int, censor, xlim, ylim, xlabel, ylabel, xscale, ")", sep=""))
			doItAndPrint(paste("legend (", place, ", legend, ", line2, ' box.lty=0, title="', strata3, group[i], '")', sep=""))
		} else{	
			doItAndPrint('mar <- par("mar")')
			doItAndPrint("mar[1] <- mar[1] + length(km$strata) + 0.5")
			doItAndPrint("mar[2] <- mar[2] + 2")
			doItAndPrint("par(mar=mar)")
			doItAndPrint("opar <- par(mar = mar)")
			doItAndPrint("on.exit(par(opar))")
			doItAndPrint(paste('plot(km, bty="l", ', line, "conf.int=", conf.int, censor, xlim, ylim, xlabel, ylabel, xscale, ")", sep=""))
			doItAndPrint("xticks <- axTicks(1)")
			doItAndPrint(paste("n.atrisk <- nrisk(km, xticks", xscale2, ")", sep=""))
			doItAndPrint("for (i in 1:length(km$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}")			
#			doItAndPrint(paste('#for (i in 1:length(km$strata)){for (j in 1:(length(xticks)-1)) {axis(1, at=c(xticks[j]+(xticks[2]-xticks[1])/3, xticks[j+1]-+(xticks[2]-xticks[1])/3), labels=c(" ", " "), line=4.6+i, ', line2, "lwd.ticks=0, tick = TRUE)}}", sep=""))			
			doItAndPrint(paste("for (i in 1:length(km$strata)){mtext(legend[i], at=-(xticks[2]-xticks[1])/2, side=1, line=4+i, cex=", par.cex, ")}", sep=""))			
			doItAndPrint('title(xlab = "Number at risk", line = 3.5, adj = 0)')
			doItAndPrint(paste("legend (", place, ", legend, ", line2, ' box.lty=0, title="', strata3, group[i], '")', sep=""))
		}
	}else{
		if (subset == ""){
			stratas <- eval(parse(text=paste("levels(factor(", dataSet, "$", strata, "))", sep="")))		
		}else{
			stratas <- eval(parse(text=paste("levels(factor(subset(", dataSet, ", ", tclvalue(subsetVariable), ")$", strata, "))", sep="")))
		}
		nstrata <- length(stratas)
		doItAndPrint("strata.names <- NULL")
		doItAndPrint("strata.p <- NULL")
		for(j in 1: nstrata){
#			command <- paste("km <- survfit(Surv(", timetoevent, ",", event, ")~", group[i], strata2, ", data=", dataSet, "[", dataSet, "$", strata, '=="', stratas[j], '",]', subset, ', na.action = na.omit, conf.type="log-log")', sep="")
			command <- paste("km <- survfit(Surv(", timetoevent, ",", event, "==1)~", group[i], ", data=", dataSet, "[", dataSet, "$", strata, '=="', stratas[j], '",]', subset, ', na.action = na.omit, conf.type="log-log")', sep="")
			doItAndPrint(command)
			doItAndPrint(paste('len <- nchar("', group[i], '")', sep=""))
			doItAndPrint("nvar2 <- length(names(km$strata))")	
#			doItAndPrint("k <- 1; legend <- NULL")
#			doItAndPrint("for (j in 1:nvar2){legend[k] <- levels(factor(substring(names(km$strata), len+2)))[j]; k <- k+1}")
			doItAndPrint("legend <- substring(names(km$strata), len+2)")
			main <- paste(', main="', strata, "=", stratas[j], '"', sep="")
			if (atrisk==0){
				doItAndPrint(paste('plot(km, bty="l", ', line, "conf.int=", conf.int, censor, xlim, ylim, xlabel, ylabel, main, xscale, ")", sep=""))
				doItAndPrint(paste("legend (", place, ", legend, ", line2, 'box.lty=0, title="', strata3, group[i], '")', sep=""))
			}else{
				doItAndPrint('mar <- par("mar")')
				doItAndPrint("mar[1] <- mar[1] + length(km$strata) + 0.5")
				doItAndPrint("mar[2] <- mar[2] + 2")
				doItAndPrint("par(mar=mar)")
				doItAndPrint("opar <- par(mar = mar)")
				doItAndPrint("on.exit(par(opar))")
				doItAndPrint(paste('plot(km, bty="l", ', line, "conf.int=", conf.int, censor, xlim, ylim, xlabel, ylabel, main, xscale, ")", sep=""))
				doItAndPrint("xticks <- axTicks(1)")
				doItAndPrint(paste("n.atrisk <- nrisk(km, xticks", xscale2, ")", sep=""))
				doItAndPrint("for (i in 1:length(km$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}")
#				doItAndPrint(paste('#for (i in 1:length(km$strata)){for (j in 1:(length(xticks)-1)) {axis(1, at=c(xticks[j]+(xticks[2]-xticks[1])/3, xticks[j+1]-+(xticks[2]-xticks[1])/3), labels=c(" ", " "), line=4.6+i, ', line2, "lwd.ticks=0, tick = TRUE)}}", sep=""))			
				doItAndPrint(paste("for (i in 1:length(km$strata)){mtext(legend[i], at=-(xticks[2]-xticks[1])/2, side=1, line=4+i, cex=", par.cex, ")}", sep=""))			
				doItAndPrint('title(xlab = "Number at risk", line = 3.5, adj = 0)')
				doItAndPrint(paste("legend (", place, ", legend, ", line2, ' box.lty=0, title="', strata3, group[i], '")', sep=""))			
			}
			if (subset == ""){
				levs <- eval(parse(text=paste("length(levels(factor(", dataSet, "[", dataSet, "$", strata, '=="', stratas[j], '",]$', group[i], ")))", sep="")))			
			} else {
				levs <- eval(parse(text=paste("length(levels(factor(subset(", dataSet, ", ", tclvalue(subsetVariable), ")[subset(", dataSet, ", ", tclvalue(subsetVariable), ")$", strata, '=="', stratas[j], '",]$', group[i], ")))", sep="")))
			}
			if (levs < 2){
				doItAndPrint(paste('strata.names <- c(strata.names, "', stratas[j], '")', sep=""))
				doItAndPrint(paste("strata.p <- c(strata.p, NA)", sep=""))			
			}else{
				command2 <- paste("res <- survdiff(Surv(", timetoevent, ",", event, "==1)~", group[i], strata2, ", data=", dataSet, "[", dataSet, "$", strata, '=="', stratas[j], '",]', subset, ", rho=", test, ", na.action = na.omit)", sep="")
				doItAndPrint(command2)
				doItAndPrint(paste('strata.names <- c(strata.names, "', stratas[j], '")', sep=""))
				doItAndPrint(paste("strata.p <- c(strata.p, signif(pchisq(c(res$chisq), df=length(res$n)-1, lower.tail=FALSE),digits=3))", sep=""))
				doItAndPrint("remove(res)")
			}
			if (j < nstrata) 	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
		}
		doItAndPrint(paste("strata.data <- data.frame(", strata, "=strata.names, p.value=strata.p)", sep=""))
 		logger("p-value calculated in each strata")
		doItAndPrint("strata.data")
#		doItAndPrint("remove(strata.data)")		
		command <- paste("km <- survfit(Surv(", timetoevent, ",", event, "==1)~", group[i], strata2, ", data=", ActiveDataSet(), subset, ', na.action = na.omit, conf.type="log-log")', sep="")
		doItAndPrint(command)	#To create km.summary.table in all strata
		}
#	command <- paste("km <- survfit(Surv(", timetoevent, ",", event, "==1)~", group[i], strata2, ", data=", ActiveDataSet(), subset, ', na.action = na.omit, conf.type="log-log")', sep="")
#	doItAndPrint(command)
	command2 <- paste("(res <- survdiff(Surv(", timetoevent, ",", event, "==1)~", group[i], strata2, ", data=", dataSet, subset, ", rho=", test, ', na.action = na.omit))', sep="")
	doItAndPrint(command2)
	if (i == 1){
		doItAndPrint(paste("km.summary.table <- summary.km(survfit=km, survdiff=res", point, xscale2, ")", sep=""))
	} else {
		doItAndPrint(paste("km.summary.table <- rbind(km.summary.table, summary.km(survfit=km, survdiff=res", point, xscale2, "))", sep=""))
	}
	if(nvar==1 & posthoc!=""){
		if (length(strata)==0) {
			command <- paste("pairwise.logrank.test(", sub1, dataSet, sub2, "$", timetoevent, ", ", sub1, dataSet, sub2, "$", event, ", ", sub1, dataSet, sub2, "$", group[i], ', strata=NULL, "', dataSet, '", p.adjust.method="', posthoc, '", rho=', test, ")", sep="")
		} else{
			command <- paste("pairwise.logrank.test(", sub1, dataSet, sub2, "$", timetoevent, ", ", sub1, dataSet, sub2, "$", event, ", ", sub1, dataSet, sub2, "$", group[i], ", strata=", sub1, dataSet, sub2, "$", strata, ', "', dataSet, '", p.adjust.method="', posthoc, '", rho=', test, ")", sep="")
		}
		doItAndPrint(command)
	}
	doItAndPrint("remove(res)")	
    }
	doItAndPrint("km.summary.table")
#	doItAndPrint("remove(km.summary.table)")	
	}
	doItAndPrint("remove(km)")
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="survfit", model=TRUE, apply="StatMedKaplanMeier", reset="StatMedKaplanMeier")
    tkgrid(getFrame(timetoeventBox), labelRcmdr(variablesFrame, text="    "), getFrame(eventBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(variables2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
    tkgrid(getFrame(groupBox), labelRcmdr(variables2Frame, text="    "), getFrame(strataBox), sticky="nw")
    tkgrid(variables2Frame, sticky="nw")
	
	tkgrid(testFrame, labelRcmdr(plotoptionFrame, text="   "), lineFrame, labelRcmdr(plotoptionFrame, text="   "), placeFrame, labelRcmdr(plotoptionFrame, text="   "), xscaleFrame, labelRcmdr(plotoptionFrame, text="   "), posthocFrame, sticky="w")
	tkgrid(plotoptionFrame, sticky="nw")

	tkgrid(censor, labelRcmdr(plotoption2Frame, text=" "), ci, labelRcmdr(plotoption2Frame, text=" "), separatestrata, labelRcmdr(plotoption2Frame, text=" "), atrisk, sticky="w")
	tkgrid(plotoption2Frame, sticky="nw")
	
	#	tkgrid(plotoptionFrame, plotoption2Frame, sticky="nw")

#    tkgrid(plotoption2_1Frame, plotoption2_2Frame, sticky="w")
#	tkgrid(plotoption2Frame, sticky="w")

	tkgrid(labelRcmdr(pointFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Time point to show survival rate")), pointField,  sticky = "w")
  tkgrid(pointFrame, sticky="w")	
	tkgrid(labelRcmdr(xlimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimField, sticky = "w")
	tkgrid(labelRcmdr(ylimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimField, sticky = "w")
   tkgrid(xlimFrame, labelRcmdr(axis2Frame, text="  "), ylimFrame, sticky="w")
	tkgrid(labelRcmdr(xlabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis lavel")), xlabelField, sticky = "w")
	tkgrid(labelRcmdr(ylabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis lavel")), ylabelField, sticky = "w")
   tkgrid(xlabelFrame, labelRcmdr(axis2Frame, text="  "), ylabelFrame, sticky="w")

#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Time point to show survival rate")), pointEntry, sticky="w")
#	tkgrid.configure(pointEntry, sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimEntry, sticky="w")
#	tkgrid.configure(xlimEntry, sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimEntry, sticky="w")
#	tkgrid.configure(ylimEntry, sticky="w")
  tkgrid(axisFrame, sticky="w")
  tkgrid(axis2Frame, sticky="w")
  
  StatMedSubsetBox(model=TRUE)
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedLogrankTrend  <- function(){
defaults <- list(event = "", timetoevent = "", group = "", subset = "")
dialog.values <- getDialog("StatMedLogrankTrend", defaults)

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Logrank trend test"))
    variablesFrame <- tkframe(top)
    eventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Status indicator (censor=0, event=1) (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$event, "all"))
    timetoeventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time-to-event variable (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$timetoevent, "all"))
    variables2Frame <- tkframe(top)
    groupBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$group, "all"))
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Logrank trend test"), "#####", sep=""))
    event <- getSelection(eventBox)
    timetoevent <- getSelection(timetoeventBox)
    group <- getSelection(groupBox)
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
	  sub1 <- ""
	  sub2 <- ""
      subset <- ""
	  }
    else{
	  sub1 <- "subset("
	  sub2 <- paste(", ", subset, ")", sep="")
      subset <- paste(", subset=", subset, sep="")
	  }
	dataSet <- activeDataSet()	
	putDialog("StatMedLogrankTrend", list(event = event, timetoevent = timetoevent, group = group, subset = tclvalue(subsetVariable)))
    if (length(event) != 1) {
      errorCondition(recall=StatMedLogrankTrend, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one status indicator (censor=0, event=1)"))
      return()
    }
    if (length(timetoevent) != 1) {
      errorCondition(recall=StatMedLogrankTrend, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one time-to-event variable"))
      return()
    }
    closeDialog()
	Library("survival")
    nvar <- length(group)
	command <- paste("(res <- survdiff(Surv(", timetoevent, ",", event, "==1)~", group, ", data=", dataSet, subset, ', na.action = na.omit))', sep="")
	doItAndPrint(command)
	doItAndPrint("logrank.trend(res)")
	doItAndPrint("remove(res)")
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="survdiff", model=TRUE, apply="StatMedLogrankTrend", reset="StatMedLogrankTrend")
    tkgrid(getFrame(timetoeventBox), labelRcmdr(variablesFrame, text="    "), getFrame(eventBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(getFrame(groupBox), labelRcmdr(variables2Frame, text="    "), sticky="nw")
    tkgrid(variables2Frame, sticky="nw")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Trend will be evaluated among groups in alphabetical order"), fg="blue"), sticky="w")
  StatMedSubsetBox()
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedCoxRegression  <- function(){
# add the class coxph to the modelClasses,  from fncCoxMode() in RcmdrPlugin.SurvivalT
    xx <- getRcmdr("modelClasses")
    bolCoxphExists = FALSE
    for(ii in 1:length(xx)){if (xx[ii] == "coxph") bolCoxphExists = TRUE}
    if (bolCoxphExists == FALSE) putRcmdr("modelClasses", c(getRcmdr("modelClasses"), "coxph"))
defaults <- list(SurvivalTimeVariable = "", StatusVariable = "", rhs = "", waldVariable = 0,  prophazVariable = 0, martinVariable = 0, basecurveVariable = 0, actmodelVariable = 0, stepwise1Variable = 0, stepwise2Variable = 0, stepwise3Variable = 0, subset = "")
dialog.values <- getDialog("StatMedCoxRegression", defaults)
currentFields$SurvivalTimeVariable <- dialog.values$SurvivalTimeVariable	
currentFields$StatusVariable <- dialog.values$StatusVariable
currentFields$rhs <- dialog.values$rhs
currentFields$subset <- dialog.values$subset	

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Cox proportional hazard regression"))
  .activeModel <- ActiveModel()
  currentModel <- if (!is.null(.activeModel)) 
        class(get(.activeModel, envir=.GlobalEnv))[1] == "coxph"
#    eval(parse(text=paste("class(", .activeModel, ")[1] == 'coxph'", sep="")), 
#         envir=.GlobalEnv) 
    else FALSE

	currentModel <- TRUE
#  if(currentModel){
#    currentFields <- formulaFields(eval(parse(text=.activeModel), 
#     envir=.GlobalEnv))
#    if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
#  }

  UpdateModelNumber()
  modelName <- tclVar(paste("CoxModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="30", textvariable=modelName)
  	optionsFrame <- tkframe(top)
	
	checkBoxes(frame="checkboxFrame", boxes=c("wald", "prophaz", "martin", "basecurve", "actmodel", "stepwise1", "stepwise2", "stepwise3"), initialValues=c(dialog.values$waldVariable, dialog.values$prophazVariable, dialog.values$martinVariable, dialog.values$basecurveVariable, dialog.values$actmodelVariable, dialog.values$stepwise1Variabl, dialog.values$stepwise2Variabl, dialog.values$stepwise3Variabl),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Wald test for overall p-value for factors with >2 levels", "Test proportional hazards assumption","Plot martingale residuals", "Show baseline survival curve", "Keep results as active model for further analyses", "Stepwise selection based on AIC", "Stepwise selection based on BIC", "Stepwise selection based on p-value")))	

#	waldVariable <- dialog.values$waldVariable
#	waldCheckBox <- tkcheckbutton(optionsFrame, variable=waldVariable)
#	prophazVariable <- dialog.values$prophazVariable
#	prophazCheckBox <- tkcheckbutton(optionsFrame, variable=prophazVariable)
#	basecurveVariable <- dialog.values$basecurveVariable
#	basecurveCheckBox <- tkcheckbutton(optionsFrame, variable=basecurveVariable)
#	actmodelVariable <- dialog.values$actmodelVariable
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
#	stepwise1Variable <- dialog.values$stepwise1Variable
#	stepwise1CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise1Variable)
#	stepwise2Variable <- dialog.values$stepwise2Variable
#	stepwise2CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise2Variable)
#	stepwise3Variable <- dialog.values$stepwise3Variable
#	stepwise3CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise3Variable)

  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Cox proportional hazard regression"), "#####", sep=""))
#    XXX <- getSelection(timeBox)
    modelValue <- trim.blanks(tclvalue(modelName))
		wald <- tclvalue(waldVariable)
		prophaz <- tclvalue(prophazVariable)
		martin <- tclvalue(martinVariable)
		basecurve <- tclvalue(basecurveVariable)
		actmodel <- tclvalue(actmodelVariable)
		stepwise1 <- tclvalue(stepwise1Variable)
		stepwise2 <- tclvalue(stepwise2Variable)
		stepwise3 <- tclvalue(stepwise3Variable)
		subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
putDialog("StatMedCoxRegression", list(SurvivalTimeVariable = tclvalue(SurvivalTimeVariable), StatusVariable = tclvalue(StatusVariable), rhs = tclvalue(rhsVariable), waldVariable = wald,  prophazVariable = prophaz, martinVariable = martin, basecurveVariable = basecurve, actmodelVariable = actmodel, stepwise1Variable = stepwise1, stepwise2Variable = stepwise2, stepwise3Variable = stepwise3, subset=tclvalue(subsetVariable)))
    if (!is.valid.name(modelValue)){
      errorCondition(recall=StatMedCoxRegression, 
        message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue), model=TRUE)
      return()
    }
#    check.empty <- gsub(" ", "", tclvalue(lhsVariable))
#    if ("" == check.empty) {
#      errorCondition(recall=StatMedCoxRegression,
#        message=gettext(domain="R-RcmdrPlugin.EZR","Left-hand side of model empty."), model=TRUE) 
#      return()
#    }
        check.empty <- gsub(" ", "", tclvalue(SurvivalTimeVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedCoxRegression, message=gettext(domain="R-RcmdrPlugin.EZR","Survival time variable of model empty."), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(StatusVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedCoxRegression, message=gettext(domain="R-RcmdrPlugin.EZR","Status variable of model empty."), model=TRUE)
            return()
            }
    check.empty <- gsub(" ", "", tclvalue(rhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=StatMedCoxRegression,
        message=gettext(domain="R-RcmdrPlugin.EZR","Right-hand side of model empty."), model=TRUE)
      return()
    }
    if (is.element(modelValue, listCoxModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
        UpdateModelNumber(-1)
        StatMedCoxRegression()
        return()
      }
    }	
    closeDialog()	
    Library("survival")
	Library("aod")
#    formula <- paste("Surv(", XXX, ", ", tclvalue(lhsVariable), ") ~ ", tclvalue(rhsVariable), sep="")

     formula <- paste("Surv(", tclvalue(SurvivalTimeVariable), ", ", tclvalue(StatusVariable), "==1)~ ", tclvalue(rhsVariable), sep="")

    command <- paste("coxph(", formula,
      ", data=", ActiveDataSet(), subset, ', method="breslow")', sep="")
#    logger(paste(modelValue, " <- ", command, sep=""))
#    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
    doItAndPrint(paste("(res <- summary(", modelValue, "))", sep=""))
#	doItAndPrint(paste("res <- ", command, sep=""))
#	doItAndPrint("res <- summary(res)")
	if(eval(parse(text="length(res$coefficients[,1])"))==1){
		doItAndPrint("cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)")
		doItAndPrint("rownames(cox.table) <- rownames(res$coefficients)")
		doItAndPrint('colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	} else {
		doItAndPrint("cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)")
		doItAndPrint("cox.table <- data.frame(cox.table)")
		doItAndPrint('colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	}
#	doItAndPrint("cox.table <- signif(cox.table, digits=3)")
	doItAndPrint("cox.table")
	if (wald==1) doItAndPrint(paste("waldtest(", modelValue, ")", sep=""))
	if (martin==1){
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}			
		doItAndPrint(paste("scatter.smooth(residuals(", modelValue, ', type="martingale"))', sep=""))	
		doItAndPrint("abline(h=0, lty=3)")	
		}
	if (prophaz == 1){
			doItAndPrint(paste("print(cox.zph(", modelValue, "))", sep=""))
	}
	if (basecurve ==1){
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}			
		doItAndPrint(paste("plot(survfit(", modelValue, "))", sep=""))
	}
	
	if (stepwise1 == 1 | stepwise2 == 1 | stepwise3 == 1){
		x <- strsplit(tclvalue(rhsVariable), split="\\+")
		command <- paste("TempDF <- with(", ActiveDataSet(), ", ", ActiveDataSet(), "[complete.cases(", paste(x[[1]], collapse=","), "),])", sep="")
		doItAndPrint(command)
		command <- paste("coxph(", formula, ", data=TempDF", subset, ', method="breslow")', sep="")
		doItAndPrint(paste(modelValue, " <- ", command, sep=""))
		}
	if (stepwise1 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="AIC")', sep=""))
			doItAndPrint("summary(res)")
			doItAndPrint("res2 <- summary(res)")
			if(eval(parse(text="length(res2$coefficients[,1])"))==1){
				doItAndPrint("cox.table <- signif(cbind(t(res2$conf.int[,c(1,3,4)]), p.value=res2$coefficients[,5]), digits=4)")
				doItAndPrint("rownames(cox.table) <- rownames(res2$coefficients)")
				doItAndPrint('colnames(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')	
				doItAndPrint("cox.table")
			} else if(eval(parse(text="length(res2$coefficients[,1])"))>1){
				doItAndPrint("cox.table <- signif(cbind(res2$conf.int[,c(1,3,4)], res2$coefficients[,5]), digits=4)")
				doItAndPrint("cox.table <- data.frame(cox.table)")
				doItAndPrint('names(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
				doItAndPrint("cox.table")
			}
			if (wald==1) doItAndPrint("waldtest(res)")
			}
	if (stepwise2 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="BIC")', sep=""))
			doItAndPrint("summary(res)")
			doItAndPrint("res2 <- summary(res)")
			if(eval(parse(text="length(res2$coefficients[,1])"))==1){
				doItAndPrint("cox.table <- signif(cbind(t(res2$conf.int[,c(1,3,4)]), p.value=res2$coefficients[,5]), digits=4)")
				doItAndPrint("rownames(cox.table) <- rownames(res2$coefficients)")
				doItAndPrint('colnames(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')	
				doItAndPrint("cox.table")
			} else if(eval(parse(text="length(res2$coefficients[,1])"))>1){
				doItAndPrint("cox.table <- signif(cbind(res2$conf.int[,c(1,3,4)], res2$coefficients[,5]), digits=4)")
				doItAndPrint("cox.table <- data.frame(cox.table)")
				doItAndPrint('names(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
				doItAndPrint("cox.table")
			}
			if (wald==1) doItAndPrint("waldtest(res)")
	}
	if (stepwise3 == 1){
			subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
				|| trim.blanks(subset) == ""){
				subset <- ""
			}
			else{
				subset <- paste(", subset='", trim.blanks(subset), "'", sep="")
			}
			doItAndPrint(paste('step.p.cox(', modelValue, ', "TempDF", wald=', wald, subset, ")", sep=""))
	}
	doItAndPrint("remove(res)")
#	doItAndPrint("remove(cox.table)")
	if (actmodel==1) activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="coxph", model=TRUE, apply="StatMedCoxRegression", reset="StatMedCoxRegression")
  tkgrid(tklabel(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")

  modelFormulaCox()
  StatMedSubsetBox(model=TRUE)
  
 tkgrid(getFrame(xBox), sticky="w")
#  tkgrid(getFrame(xBox), getFrame(timeBox), sticky="w")
  tkgrid(outerOperatorsFrame, sticky="w")
  tkgrid(formulaFrame, sticky="w")
  	tkgrid(labelRcmdr(top, text=paste("            ", gettext(domain="R-RcmdrPlugin.EZR","Stratifing variable: + strata(#####)"), sep="")), sticky="w") 
	
 tkgrid(checkboxFrame, sticky="w")
	
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Wald test for overall p-value for factors with >2 levels")), waldCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Test proportional hazards assumption")), prophazCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Show baseline survival curve")), basecurveCheckBox, sticky="w")	
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on AIC")), stepwise1CheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on BIC")), stepwise2CheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on p-value")), stepwise3CheckBox, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
}


StatMedAdjustedSurvival  <- function(){
defaults <- list(event = "", timetoevent = "", group = "", adjust = "", line = "color", place = "topright", xscale = "", censor = 1, atrisk = 0, xlim = "<auto>", ylim = "<auto>", xlabel = "<auto>", ylabel = "<auto>", subset = "")
dialog.values <- getDialog("StatMedAdjustedSurvival", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Adjusted survival curve"))
    variablesFrame <- tkframe(top)
    eventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Status indicator (censor=0, event=1) (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$event, "all"))
    timetoeventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time-to-event variable (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$timetoevent, "all"))
    variables2Frame <- tkframe(top)
    groupBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), listHeight=8, initialSelection=varPosn(dialog.values$group, "all"))
    adjustBox <- variableListBox(variables2Frame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables for adjustment (pick at least one)"), listHeight=8, initialSelection=varPosn(dialog.values$adjust, "all"))
    plotoptionFrame <- tkframe(top)
    radioButtons(plotoptionFrame, name="line", buttons=c("color", "type", "width"), values=c("color", "type", "width"), initialValue=dialog.values$line,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Color", "Line type", "Line width")), title=gettext(domain="R-RcmdrPlugin.EZR","Line discrimination"))
    radioButtons(plotoptionFrame, name="place", buttons=c("topright", "bottom", "mouse"), values=c("topright", "bottom", "mouse"), initialValue=dialog.values$place, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Upper right", "Bottom", "Mouse click")), title=gettext(domain="R-RcmdrPlugin.EZR","Legend"))
    radioButtons(plotoptionFrame, name="xscale", buttons=c("day", "daytomonth", "daytoyear", "monthtoyear"), values=c("", "30.4375", "365.25", "12"), initialValue=dialog.values$xscale, labels=gettext(domain="R-RcmdrPlugin.EZR",c("As is", "Day to month", "Day to year", "Month to year")), title=gettext(domain="R-RcmdrPlugin.EZR","X axis"))
	plotoption2Frame <- tkframe(top)
    checkBoxes(window=plotoption2Frame, frame="censor", boxes=c("censor"),initialValues=dialog.values$censor,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show censoring marks")), title=gettext(domain="R-RcmdrPlugin.EZR","Options"))
    checkBoxes(window=plotoption2Frame, frame="atrisk", boxes=c("atrisk"),initialValues=dialog.values$atrisk,labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show number at risk")), title=gettext(domain="R-RcmdrPlugin.EZR"," "))
	axisFrame <- tkframe(top)
	xlimFrame <- tkframe(axisFrame)
	xlimVariable <- tclVar(dialog.values$xlim)
	xlimField <- ttkentry(axisFrame, width="20", textvariable=xlimVariable)
	ylimFrame <- tkframe(axisFrame)
	ylimVariable <- tclVar(dialog.values$ylim)
	ylimField <- ttkentry(axisFrame, width="20", textvariable=ylimVariable)
	xlabelFrame <- tkframe(axisFrame)
	xlabelVariable <- tclVar(dialog.values$xlabel)
	xlabelField <- ttkentry(axisFrame, width="20", textvariable=xlabelVariable)
	ylabelFrame <- tkframe(axisFrame)
	ylabelVariable <- tclVar(dialog.values$ylabel)
	ylabelField <- ttkentry(axisFrame, width="20", textvariable=ylabelVariable)

onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Adjusted survival curve"), "#####", sep=""))
    event <- getSelection(eventBox)
    timetoevent <- getSelection(timetoeventBox)
    group <- getSelection(groupBox)
    adjust <- getSelection(adjustBox)
	dataSet <- activeDataSet()
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
		subdataSet <- dataSet
		naexcludeSubdataSet <- paste("subset(", dataSet, ", ", sep="")
		}
    else{
		subdataSet <- paste("subset(", dataSet, ", ", subset, ")", sep="")
		naexcludeSubdataSet <- paste("subset(", dataSet, ", (", subset, ") & ", sep="")	
		}
	line <- tclvalue(lineVariable)	
	par.lwd <- get("par.lwd", envir=.GlobalEnv)
	if (line=="color") line <- paste("col=1:32, lty=1, ", par.lwd, ", ", sep="")
	if (line=="type") line <- paste("col=1, lty=1:32, ", par.lwd, ", ", sep="")
	if (line=="width") line <- paste("col=1, lty=1, ", par.lwd, ":8, ", sep="")
	par.cex <- get("par.cex", envir=.GlobalEnv)	
	if(length(group)==0){line <- paste("col=1, lty=1, ", par.lwd, ", ", sep="")}
	place <- tclvalue(placeVariable)
	if(place=="mouse"){
		place <- "locator(1)"
	}else if (place=="topright"){
		place <- '"topright"'
	}else{
		place <- '"bottom", horiz=TRUE'
	}
    censor <- tclvalue(censorVariable)
    atrisk <- tclvalue(atriskVariable)
	xscale <- tclvalue(xscaleVariable)
	xscale2 <- ""
	if (xscale!=""){
		xscale2 <- paste(" * ", xscale, sep="")
		xscale <- paste(", xscale=", xscale, sep="")
	}
	xlim <- tclvalue(xlimVariable)
	ylim <- tclvalue(ylimVariable)
	if (xlim == "<auto>") {
		xlim <- ""
	} else {
		xlim <- paste(", xlim=c(", xlim, ")", sep="")
	}
	if (ylim == "<auto>") {
		ylim <- ""
	} else {
		ylim <- paste(", ylim=c(", ylim, ")", sep="")
	}
	xlabel <- tclvalue(xlabelVariable)
	ylabel <- tclvalue(ylabelVariable)
	if (xlabel == "<auto>") {
		xlabel <- paste(', xlab="', timetoevent, '"', sep="")
	} else {
		xlabel <- paste(', xlab="', xlabel, '"', sep="")
	}
	if (ylabel == "<auto>") {
		ylabel <- ', ylab="Probability"'
	} else {
		ylabel <- paste(', ylab="', ylabel, '"', sep="")
	}
	if (censor==0){
		censor <- ", mark.time=FALSE"
	}else{
		censor <- ", mark.time=TRUE"	
	}
putDialog("StatMedAdjustedSurvival", list(event = event, timetoevent = timetoevent, group = group, adjust = adjust, line = tclvalue(lineVariable), place = tclvalue(placeVariable), xscale = tclvalue(xscaleVariable), censor = tclvalue(censorVariable), atrisk = atrisk, xlim = tclvalue(xlimVariable), ylim = tclvalue(ylimVariable), xlabel = tclvalue(xlabelVariable), ylabel = tclvalue(ylabelVariable), subset = tclvalue(subsetVariable)))
    if (length(event) != 1) {
      errorCondition(recall=StatMedAdjustedSurvival, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one status indicator (censor=0, event=1)"))
      return()
    }
    if (length(timetoevent) != 1) {
      errorCondition(recall=StatMedAdjustedSurvival, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one time-to-event variable"))
      return()
    }
    if (length(adjust) == 0) {
      errorCondition(recall=StatMedAdjustedSurvival, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick at least one variable for adjustment."))
      return()
    }
    closeDialog()
    Library("survival")
	factor <- adjust[1]
	naexcludeSubdataSet <- paste(naexcludeSubdataSet, "(is.na(", adjust[1], ")==F", sep="")
	if(length(adjust)>1){
		for (i in 2:length(adjust)){
			factor <- paste(factor, " + ", adjust[i], sep="")
			naexcludeSubdataSet <- paste(naexcludeSubdataSet, " & is.na(", adjust[i], ")==F", sep="")
			}
	}
	factor2 <- factor
	naexcludeSubdataSet <- paste(naexcludeSubdataSet, "))", sep="")
	if (length(group)==1) factor2 <- paste(factor, " + strata(", group, ")", sep="")	
	command <- paste("coxmodel <- coxph(Surv(", timetoevent, ", ", event, "==1)~ ", factor2, ", data=", subdataSet, ', method="breslow")', sep="")
	doItAndPrint(command)
	doItAndPrint("cox <- survfit(coxmodel)")
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
	if(length(group)==1){		
			check.type <- eval(parse(text=paste(subdataSet, "$", group, sep="")))
			if(is.integer(check.type) | is.numeric(check.type)){
				doItAndPrint(paste('len <- nchar("', group, '")', sep=""))
				doItAndPrint("group.levels <- substring(names(cox$strata[cox$strata>0]),len+2)")
			} else {
				doItAndPrint("group.levels <- names(cox$strata[cox$strata>0])")
			}
	}
	if(atrisk==1){
		if(length(group)==0){
			doItAndPrint('mar <- par("mar")')
			doItAndPrint("mar[1] <- mar[1] + 1 + 0.5")
			doItAndPrint("par(mar=mar)")
			doItAndPrint("opar <- par(mar = mar)")
			doItAndPrint("on.exit(par(opar))")
			command3 <- paste("plot(cox, ", line, 'bty="l"', censor, xlim, ylim, xlabel, ylabel, xscale, ")", sep="") 
			doItAndPrint(command3)
			doItAndPrint("xticks <- axTicks(1)")
			doItAndPrint(paste("n.atrisk <- nrisk(cox, xticks", xscale2, ")", sep=""))
			doItAndPrint("axis(1, at = xticks, labels = n.atrisk, line = 3, tick = FALSE)")
			doItAndPrint('title(xlab = "Number at risk", line = 3, adj = 0)')
		} else {
			doItAndPrint('mar <- par("mar")')
			doItAndPrint("mar[1] <- mar[1] + length(cox$strata) + 0.5")
			doItAndPrint("mar[2] <- mar[2] + 2")
			doItAndPrint("par(mar=mar)")
			doItAndPrint("opar <- par(mar = mar)")
			doItAndPrint("on.exit(par(opar))")
			command3 <- paste("plot(cox, ", line, 'bty="l"', censor, xlim, ylim, xlabel, ylabel, xscale, ")", sep="") 			
			doItAndPrint(command3)			
			doItAndPrint("xticks <- axTicks(1)")
			doItAndPrint(paste("n.atrisk <- nrisk(cox, xticks", xscale2, ")", sep=""))
			doItAndPrint("for (i in 1:length(cox$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}")
#			doItAndPrint(paste('#for (i in 1:length(cox$strata)){for (j in 1:(length(xticks)-1)) {axis(1, at=c(xticks[j]+(xticks[2]-xticks[1])/3, xticks[j+1]-+(xticks[2]-xticks[1])/3), labels=c(" ", " "), line=4.6+i, ', line2, "lwd.ticks=0, tick = TRUE)}}", sep=""))			
			doItAndPrint(paste("for (i in 1:length(cox$strata)){mtext(group.levels[i], at=-(xticks[2]-xticks[1])/2, side=1, line=4+i, cex=", par.cex, ")}", sep=""))			
			doItAndPrint('title(xlab = "Number at risk", line = 3.5, adj = 0)')
#			doItAndPrint(paste("legend (", place, ", legend, ", line, ' box.lty=0, title="', strata3, group[i], '")', sep=""))			
		}
	} else {
		command3 <- paste("plot(cox, ", line, 'bty="l"', censor, xlim, ylim, xlabel, ylabel, xscale, ")", sep="") 
		doItAndPrint(command3)
	}
	if(length(group)==1){	
			doItAndPrint(paste("legend(", place, ', group.levels, title="', group, '", ', line, "box.lty=0)", sep=""))
	}
	doItAndPrint(paste('title("Survival curve adjusted for ', factor, '")', sep=""))
	doItAndPrint("summary(cox)")
	doItAndPrint("remove(cox)")
	doItAndPrint("remove(coxmodel)")
    tkfocus(CommanderWindow())
	}
  OKCancelHelp(helpSubject="coxph", model=TRUE, apply="StatMedAdjustedSurvival", reset="StatMedAdjustedSurvival")
    tkgrid(getFrame(timetoeventBox), labelRcmdr(variablesFrame, text="    "), getFrame(eventBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(variables2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
    tkgrid(getFrame(groupBox), labelRcmdr(variables2Frame, text="    "), getFrame(adjustBox), sticky="nw")
    tkgrid(variables2Frame, sticky="nw")
	tkgrid(lineFrame, labelRcmdr(plotoptionFrame, text="   "), placeFrame, labelRcmdr(plotoptionFrame, text="   "), xscaleFrame, sticky="w")
	tkgrid(plotoptionFrame, sticky="nw")
	tkgrid(censor, labelRcmdr(plotoption2Frame, text="  "), atrisk, sticky="w")
	tkgrid(plotoption2Frame, sticky="nw")

	tkgrid(labelRcmdr(xlimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimField, sticky = "w")
#  tkgrid(xlimFrame, sticky="w")
	tkgrid(labelRcmdr(ylimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimField, sticky = "w")
#  tkgrid(ylimFrame, sticky="w")
  tkgrid(xlimFrame, labelRcmdr(axisFrame, text="  "), ylimFrame, sticky="w")
	tkgrid(labelRcmdr(xlabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis lavel")), xlabelField, sticky = "w")
	tkgrid(labelRcmdr(ylabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis lavel")), ylabelField, sticky = "w")
   tkgrid(xlabelFrame, labelRcmdr(axisFrame, text="  "), ylabelFrame, sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimEntry, sticky="w")
#	tkgrid.configure(xlimEntry, sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimEntry, sticky="w")
#	tkgrid.configure(ylimEntry, sticky="w")
    tkgrid(axisFrame, sticky="w")
  StatMedSubsetBox(model=TRUE)
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedCumInc  <- function(){
defaults <- list(event = "", timetoevent = "", group = "", line = "color", place = "topright", xscale = "", posthoc = "", censor = 1, atrisk = 0, point = "<none>", plotevent = "<all>", xlim = "<auto>", ylim = "<auto>", xlabel = "<auto>", ylabel = "<auto>", subset = "")
dialog.values <- getDialog("StatMedCumInc", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Cumulative incidence of competing events and Gray test"))
    variablesFrame <- tkframe(top)
    eventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Status indicator (censor=0, event=1,2,3...) (pick one)"), listHeight=7, initialSelection=varPosn(dialog.values$event, "all"))
    timetoeventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time-to-event variable (pick one)"), listHeight=7, initialSelection=varPosn(dialog.values$timetoevent, "all"))
    variables2Frame <- tkframe(top)
    groupBox <- variableListBox(variables2Frame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (pick 0, 1, or more)"), listHeight=6, initialSelection=varPosn(dialog.values$group, "all"))
    plotoptionFrame <- tkframe(top)
    radioButtons(plotoptionFrame, name="line", buttons=c("color", "type", "width"), values=c("color", "type", "width"), initialValue=dialog.values$line,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Color", "Line type", "Line width")), title=gettext(domain="R-RcmdrPlugin.EZR","Line discrimination"))
    radioButtons(plotoptionFrame, name="place", buttons=c("topright", "bottom", "mouse"), values=c("topright", "bottom", "mouse"), initialValue=dialog.values$place,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Upper right", "Bottom", "Mouse click")), title=gettext(domain="R-RcmdrPlugin.EZR","Legend"))
    radioButtons(plotoptionFrame, name="xscale", buttons=c("day", "daytomonth", "daytoyear", "monthtoyear"), values=c("", "30.4375", "365.25", "12"), initialValue=dialog.values$xscale,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("As is", "Day to month", "Day to year", "Month to year")), title=gettext(domain="R-RcmdrPlugin.EZR","X axis"))
    radioButtons(plotoptionFrame, name="posthoc", buttons=c("No", "Bonferroni", "Holm"), values=c("", "bon", "holm"), initialValue=dialog.values$posthoc,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("No", "Bonferroni", "Holm")), title=gettext(domain="R-RcmdrPlugin.EZR","Post-hoc test (one event to show,\none grouping variable)"))
	plotoption2Frame <- tkframe(top)
	checkBoxes(window=plotoption2Frame, frame="censor", boxes=c("censor"), initialValues=c(dialog.values$censor),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show censoring marks")), title=gettext(domain="R-RcmdrPlugin.EZR","Options"))	
	checkBoxes(window=plotoption2Frame, frame="atrisk", boxes=c("atrisk"), initialValues=c(dialog.values$atrisk),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show number at risk")), title=gettext(domain="R-RcmdrPlugin.EZR"," "))	
	
#    checkBoxes(window=plotoption2Frame, frame="censor", boxes=c("censor"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show censoring marks")))
#    checkBoxes(window=plotoption2Frame, frame="atrisk", boxes=c("atrisk"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show number at risk")))
	axisFrame <- tkframe(top)
	axis2Frame <- tkframe(top)
	ploteventFrame <- tkframe(axisFrame)
	ploteventVariable <- tclVar(dialog.values$plotevent)	
	ploteventField <- ttkentry(ploteventFrame, width="20", textvariable=ploteventVariable)
	pointFrame <- tkframe(axisFrame)
	pointVariable <- tclVar(dialog.values$point)	
	pointField <- ttkentry(pointFrame, width="20", textvariable=pointVariable)
	xlimFrame <- tkframe(axis2Frame)
	xlimVariable <- tclVar(dialog.values$xlim)
	xlimField <- ttkentry(axis2Frame, width="20", textvariable=xlimVariable)
	ylimFrame <- tkframe(axis2Frame)
	ylimVariable <- tclVar(dialog.values$ylim)
	ylimField <- ttkentry(axis2Frame, width="20", textvariable=ylimVariable)
	xlabelFrame <- tkframe(axis2Frame)
	xlabelVariable <- tclVar(dialog.values$xlabel)
	xlabelField <- ttkentry(axis2Frame, width="20", textvariable=xlabelVariable)
	ylabelFrame <- tkframe(axis2Frame)
	ylabelVariable <- tclVar(dialog.values$ylabel)
	ylabelField <- ttkentry(axis2Frame, width="20", textvariable=ylabelVariable)
#	point <- tclVar("<none>")
#	pointEntry <- ttkentry(axisFrame, width="20", textvariable=point)
#	plotevent <- tclVar("<all>")
#	ploteventEntry <- ttkentry(axisFrame, width="20", textvariable=plotevent)
#	xlim <- tclVar("<auto>")
#	xlimEntry <- ttkentry(axisFrame, width="20", textvariable=xlim)
#	ylim <- tclVar("<auto>")
#	ylimEntry <- ttkentry(axisFrame, width="20", textvariable=ylim)
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Cumulative incidence of competing events and Gray test"), "#####", sep=""))
    event <- getSelection(eventBox)
    timetoevent <- getSelection(timetoeventBox)
    group <- getSelection(groupBox)
#    strata <- getSelection(strataBox)
	dataSet <- activeDataSet()
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
		subdataSet <- dataSet
		subset <- ""
    }
    else{
		subdataSet <- paste("subset(", dataSet, ", ", subset, ")", sep="")
		subset <- paste(", subset=", subset, sep="")
    }
	line <- tclvalue(lineVariable)
	par.lwd <- get("par.lwd", envir=.GlobalEnv)
	if (line=="color") line <- paste("col=1:32, lty=1, ", par.lwd, sep="")
	if (line=="type") line <- paste("col=1, lty=1:32, ", par.lwd, sep="")
	if (line=="width") line <- paste("col=1, lty=1, ", par.lwd, ":8", sep="")
	par.cex <- get("par.cex", envir=.GlobalEnv)
	point <- tclvalue(pointVariable)
	if (point == "<none>") {
		point <- ""
	} else {
		point <- paste(", time=", point, sep="")
	}
	place <- tclvalue(placeVariable)
	if(place=="mouse"){
		place <- "locator(1)"
	}else if (place=="topright"){
		place <- '"topright"'
	}else{
		place <- '"bottom", horiz=TRUE'
	}
#    color <- tclvalue(colorVariable)
    censor <- tclvalue(censorVariable)
    atrisk <- tclvalue(atriskVariable)
	if (censor==0){
		censor <- ", mark.time=FALSE"
	}else{
		censor <- ", mark.time=TRUE"	
	}
	plotevent <- tclvalue(ploteventVariable)
	if (plotevent == "<all>" | plotevent == "") {
		plotline <- 0
	} else {
		plotevent <- round(as.numeric(plotevent))
		nevents <- eval(parse(text=paste("length(levels(factor(", subdataSet, "$", event, ")))", sep="")))
		if (plotevent < 1 | plotevent > nevents){
			plotline <- 0
		} else {
			plotline <- plotevent
		}
	}	
	xscale <- tclvalue(xscaleVariable)
	xscale2 <- ""
	if (xscale!=""){
		xscale2 <- paste(" * ", xscale, sep="")
		xscale <- paste(", xscale=", xscale, sep="")
	}
	posthoc <- tclvalue(posthocVariable)
	xlim <- tclvalue(xlimVariable)
	ylim <- tclvalue(ylimVariable)
	if (xlim == "<auto>") {
		xlim <- ""
	} else {
		xlim <- paste(", xlim=c(", xlim, ")", sep="")
	}
	if (ylim == "<auto>") {
		ylim <- ", ylim=c(0, 1)"
	} else {
		ylim <- paste(", ylim=c(", ylim, ")", sep="")
	}
	xlabel <- tclvalue(xlabelVariable)
	ylabel <- tclvalue(ylabelVariable)
	if (xlabel == "<auto>") {
		xlabel <- paste(', xlab="', timetoevent, '"', sep="")
	} else {
		xlabel <- paste(', xlab="', xlabel, '"', sep="")
	}
	if (ylabel == "<auto>") {
		ylabel <- ', ylab="Cumulative incidence"'
	} else {
		ylabel <- paste(', ylab="', ylabel, '"', sep="")
	}	
putDialog("StatMedCumInc", list(event = event, timetoevent = timetoevent, group = group, line = tclvalue(lineVariable), place = tclvalue(placeVariable), xscale = tclvalue(xscaleVariable), posthoc = posthoc, censor = tclvalue(censorVariable), atrisk = atrisk, point = tclvalue(pointVariable), plotevent = tclvalue(ploteventVariable), xlim = tclvalue(xlimVariable), ylim = tclvalue(ylimVariable), xlabel = tclvalue(xlabelVariable), ylabel = tclvalue(ylabelVariable), subset = tclvalue(subsetVariable)))
    if (length(event) != 1) {
      errorCondition(recall=StatMedCumInc, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one status indicator (censor=0, event=1,2,3...)"))
      return()
    }
    if (length(timetoevent) != 1) {
      errorCondition(recall=StatMedCumInc, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one time-to-event variable"))
      return()
    }
#	if (length(strata) ==0){
#		strata <- ""
#	}
#	else{
#		strata <- paste(strata, ", ", sep="")
#	}
    closeDialog()	
    Library("survival")
    Library("cmprsk")
#    library(survival)
#    library(cmprsk)
#	justDoIt(paste("attach(",dataSet,")"))

	if(eval(parse(text=paste("min(", dataSet, "$", event, ", na.rm=TRUE)", sep="")))>0){	#no censoring in the dataset
		doItAndPrint(paste("DummyEventForCI <- ", dataSet, "$", event, sep=""))
		#only subset data will be used in the Surv() function, and therefore, all data should be included in the dummy data
		doItAndPrint('DummyEventForCI <- factor(DummyEventForCI, levels=c("0", levels(as.factor(DummyEventForCI)))) #Required for Surv() with mstate option')
		logger("#Making the smallest level of event as 0 to avoid the event with the smallest")
		logger("#event number will be treated as censoring when there are no censoring in the dataset.")
	} else {
		doItAndPrint(paste("DummyEventForCI <- ", dataSet, "$", event, sep=""))
		doItAndPrint("DummyEventForCI <- as.factor(DummyEventForCI) #Required for Surv() with mstate option")			
	}

    nvar <- length(group)
	nevent <- eval(parse(text=paste("length(levels(factor(", subdataSet, "$", event, "[", subdataSet, "$", event, ">0])))", sep="")))
    if (nvar == 0){	
	if(nevent==1){
		command <- paste("ci <- survfit(Surv(", timetoevent, ", ", event, ">0)~1, data=", dataSet, subset, ")", sep="")
		#Error message appears when etype option is chosen and there is only single group with only 1 event type.
		doItAndPrint(command)
		plotline <- 0
		doItAndPrint("ci$surv <- 1-ci$surv")
		doItAndPrint("tempCI <- 1-ci$lower")
		doItAndPrint("ci$lower <- 1-ci$upper")
		doItAndPrint("ci$upper <- tempCI")
		doItAndPrint("summary(ci)")	#To show cumulative incidence, substract from 1, and the add 1 for plot().		
		doItAndPrint("ci$surv <- 1-ci$surv")
		doItAndPrint("tempCI <- 1-ci$lower")
		doItAndPrint("ci$lower <- 1-ci$upper")
		doItAndPrint("ci$upper <- tempCI")		
	} else {
#		command <- paste("ci <- survfit(Surv(", timetoevent, ", ", event, ">0)~1, data=", dataSet, subset, ", etype=", event, ")", sep="")
		command <- paste("ci <- survfit(Surv(", timetoevent, ', DummyEventForCI, type="mstate")~1, data=', dataSet, subset, ")", sep="")
		doItAndPrint(command)
		command <- paste("res <- with(", dataSet, ", cuminc(", timetoevent, ", ", event, ", cencode=0", subset, ", na.action = na.omit))", sep="")
		doItAndPrint(command)
		doItAndPrint("print.ci.summary(ci=ci, res=res)")
	}
	if(nevent>1){
	if(plotline==0){
		for (j in 1:nevent){
			if(j==1) {doItAndPrint(paste("ci.summary.table <- summary.ci(ci=ci, res=res, event=", j, point, xscale2, ")", sep=""))	
			} else {
			doItAndPrint(paste("ci.summary.table <- rbind(ci.summary.table, summary.ci(ci=ci, res=res, event=", j, point, xscale2, "))", sep=""))		
			}
		}
	} else {
		doItAndPrint(paste("ci.summary.table <- summary.ci(ci=ci, res=res, event=", plotline, point, xscale2, ")", sep=""))	
	}
	}
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
	doItAndPrint(paste("compevents <- levels(factor(", subdataSet, "$", event, "))", sep=""))
	doItAndPrint("nevents <- length(compevents)")
	doItAndPrint('if (compevents[1]=="0") {compevents <- compevents[2:nevents]; nevents <- nevents - 1}')
	if (plotline==0){
		if(eval(parse(text=paste("length(levels(factor(", subdataSet, "$", event, "[", subdataSet, "$", event, ">0])))", sep="")))==1){line <- paste("col=1, lty=1, ", par.lwd, sep="")}
		if (atrisk==0){
			doItAndPrint(paste('plot(ci, fun="event", bty="l", conf.int=FALSE, ', line, xlim, ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
		} else {
			doItAndPrint('mar <- par("mar")')
			doItAndPrint("mar[1] <- mar[1] + 1 + 0.5")
			doItAndPrint("par(mar=mar)")
			doItAndPrint("opar <- par(mar = mar)")
			doItAndPrint("on.exit(par(opar))")
			doItAndPrint(paste('plot(ci, fun="event", bty="l", conf.int=FALSE, ', line, xlim, ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
			doItAndPrint("xticks <- axTicks(1)")
			doItAndPrint(paste("n.atrisk <- nrisk(ci, xticks", xscale2, ")", sep=""))
			doItAndPrint("axis(1, at = xticks, labels = n.atrisk, line = 3, tick = FALSE)")
			doItAndPrint('title(xlab = "Number at risk", line = 3, adj = 0)')
		}
		doItAndPrint(paste("legend(", place, ", compevents, ", line, ', box.lty=0, title="Competing events")', sep=""))
	}else{
		if (atrisk==0){
			doItAndPrint(paste("plot(ci[", plotline, '], fun="event", bty="l", lty=1:32, conf.int=FALSE', xlim, ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
		} else {
			doItAndPrint('mar <- par("mar")')
			doItAndPrint("mar[1] <- mar[1] + 1 + 0.5")
			doItAndPrint("par(mar=mar)")
			doItAndPrint("opar <- par(mar = mar)")
			doItAndPrint("on.exit(par(opar))")
			doItAndPrint(paste("plot(ci[", plotline, '], fun="event", bty="l", lty=1:32, conf.int=FALSE', xlim, ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
			doItAndPrint("xticks <- axTicks(1)")
			doItAndPrint(paste("n.atrisk <- nrisk(ci, xticks", xscale2, ")", sep=""))
			doItAndPrint("axis(1, at = xticks, labels = n.atrisk, line = 3, tick = FALSE)")
			doItAndPrint('title(xlab = "Number at risk", line = 3, adj = 0)')
		}
	}
    } else {
	for (i in 1:nvar) {
	if(nevent==1){
		command <- paste("ci <- survfit(Surv(", timetoevent, ", ", event, ">0)~", group[i], ", data=", dataSet, subset, ")", sep="")
		#Error message appears when etype option is chosen and there is only single group with only 1 event type.
		doItAndPrint(command)
		plotline <- 0
		doItAndPrint("ci$surv <- 1-ci$surv")
		doItAndPrint("tempCI <- 1-ci$lower")
		doItAndPrint("ci$lower <- 1-ci$upper")
		doItAndPrint("ci$upper <- tempCI")
		doItAndPrint("summary(ci)")	#To show cumulative incidence, substract from 1, and the add 1 for plot().		
		doItAndPrint("ci$surv <- 1-ci$surv")
		doItAndPrint("tempCI <- 1-ci$lower")
		doItAndPrint("ci$lower <- 1-ci$upper")
		doItAndPrint("ci$upper <- tempCI")		
		command <- paste("res <- with(", dataSet, ", cuminc(", timetoevent, ", ", event, ", ", group[i], ", cencode=0", subset, ", na.action = na.omit))", sep="")
		doItAndPrint(command)
	} else {
#	command <- paste("ci <- survfit(Surv(", timetoevent, ", ", event, ">0)~", group[i], ", data=", dataSet, subset, ", etype=", event, ")", sep="")
	command <- paste("ci <- survfit(Surv(", timetoevent, ', DummyEventForCI, type="mstate")~', group[i], ", data=", dataSet, subset, ")", sep="")
	doItAndPrint(command)
	command <- paste("res <- with(", dataSet, ", cuminc(", timetoevent, ", ", event, ", ", group[i], ", cencode=0", subset, ", na.action = na.omit))", sep="")
	doItAndPrint(command)
	doItAndPrint("print.ci.summary(ci=ci, res=res)")
	}
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
	doItAndPrint(paste("compevents <- levels(factor(", subdataSet, "$", event, "))", sep=""))
	doItAndPrint("nevents <- length(compevents)")
	doItAndPrint('if (compevents[1]=="0") {compevents <- compevents[2:nevents]; nevents <- nevents - 1}')
	doItAndPrint(paste('len <- nchar("', group[i], '")', sep=""))
	doItAndPrint("groups <- substring(names(ci$strata), len+2)")
	doItAndPrint("ngroups <- length(groups)")
	if(plotline==0){
		doItAndPrint('k <- 1; l <- 1; legend <- ""')
		doItAndPrint('for(j in 1:nevents){for(i in 1:ngroups){legend[k] <- paste(groups[i], ":", compevents[j]); ifelse(k==1,legendline <- (i-1)*nevents+j, legendline <- c(legendline, (i-1)*nevents+j)); k <- k+1 }}')	
		}else{		
#		doItAndPrint(paste("legend <- levels(factor(", subdataSet, "$", group[i], "))", sep=""))
		doItAndPrint("legend <- groups")
	}
		if (line==paste("col=1, lty=1, ", par.lwd, ":8", sep="") & par.lwd!="lwd=1") doItAndPrint(paste("legendline <- legendline + ", as.integer(substring(par.lwd, nchar(par.lwd),nchar(par.lwd))) - 1, sep=""))

	if (plotline==0){
			if (line==paste("col=legendline, lty=1, ", par.lwd, sep="")) line <- paste("col=1:32, lty=1, ", par.lwd, sep="")		#line cvariable changed for legend again changed for plot()
			if (line==paste("col=1, lty=legendline, ", par.lwd, sep="")) line <- paste("col=1, lty=1:32, ", par.lwd, sep="")
			if (line=="col=1, lty=1, lwd=legendline") line <- paste("col=1, lty=1, ", par.lwd, ":8", sep="")
			if (atrisk==0){
				doItAndPrint(paste('plot(ci, fun="event", bty="l", ', line, xlim, ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
#				doItAndPrint(paste("legend (", place, ", legend, ", line, ', box.lty=0, title="', strata3, group[i], '")', sep=""))
			}else{
				doItAndPrint('mar <- par("mar")')
				doItAndPrint("mar[1] <- mar[1] + length(ci$strata) + 0.5")
				doItAndPrint("mar[2] <- mar[2] + 2")
				doItAndPrint("par(mar=mar)")
				doItAndPrint("opar <- par(mar = mar)")
				doItAndPrint("on.exit(par(opar))")
				doItAndPrint(paste('plot(ci, fun="event", bty="l", ', line, xlim, ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
				doItAndPrint("xticks <- axTicks(1)")
				doItAndPrint(paste("n.atrisk <- nrisk(ci, xticks", xscale2, ")", sep=""))
				doItAndPrint("for (i in 1:length(ci$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}")
				doItAndPrint(paste("for (i in 1:length(ci$strata)){mtext(groups[i], side=1, at=-(xticks[2]-xticks[1])/2, line=4+i, cex=", par.cex, ")}", sep=""))
				doItAndPrint('title(xlab = "Number at risk", line = 3.5, adj = 0)')
#				doItAndPrint(paste("legend (", place, ", legend, ", line, ', box.lty=0, title="', strata3, group[i], '")', sep=""))			
			}
			if (line==paste("col=1:32, lty=1, ", par.lwd, sep="")) line <- paste("col=legendline, lty=1, ", par.lwd, sep="")
			if (line==paste("col=1, lty=1:32, ", par.lwd, sep="")) line <- paste("col=1, lty=legendline, ", par.lwd, sep="")
			if (line==paste("col=1, lty=1, ", par.lwd, ":8", sep="")) line <- "col=1, lty=1, lwd=legendline"
			doItAndPrint(paste("legend(", place, ", legend, box.lty=0, ", line, ', title="', group[i], ' : Competing events")', sep=""))
	}else{
			if (atrisk==0){
				doItAndPrint(paste("plot(ci[,", plotline, '], fun="event", bty="l", ', line, xlim, 	ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
#				doItAndPrint(paste("legend (", place, ", legend, ", line, ', box.lty=0, title="', strata3, group[i], '")', sep=""))
			}else{
				doItAndPrint('mar <- par("mar")')
				doItAndPrint("mar[1] <- mar[1] + length(ci$strata) + 0.5")
				doItAndPrint("mar[2] <- mar[2] + 2")
				doItAndPrint("par(mar=mar)")
				doItAndPrint("opar <- par(mar = mar)")
				doItAndPrint("on.exit(par(opar))")
				doItAndPrint(paste("plot(ci[,", plotline, '], fun="event", bty="l", ', line, xlim, ylim, xlabel, ylabel, censor, xscale, ")", sep=""))
				doItAndPrint("xticks <- axTicks(1)")
				doItAndPrint(paste("n.atrisk <- nrisk(ci, xticks", xscale2, ")", sep=""))
				doItAndPrint("for (i in 1:length(ci$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}")
				doItAndPrint(paste("for (i in 1:length(ci$strata)){mtext(groups[i], side=1, at=-(xticks[2]-xticks[1])/2, line=4+i, cex=", par.cex, ")}", sep=""))
				doItAndPrint('title(xlab = "Number at risk", line = 3.5, adj = 0)')
#				doItAndPrint(paste("legend (", place, ", legend, ", line, ', box.lty=0, title="', strata3, group[i], '")', sep=""))			
			}
			doItAndPrint(paste("legend(", place, ", legend, box.lty=0, ", line, ', title="', group[i], '")', sep=""))
	}
	doItAndPrint("res$Tests")
	if(nevent>1){
	if(plotline==0){
		for (j in 1:nevent){
			if(i==1 & j==1) {doItAndPrint(paste("ci.summary.table <- summary.ci(ci=ci, res=res, event=", j, point, xscale2, ")", sep=""))	
			} else {
			doItAndPrint(paste("ci.summary.table <- rbind(ci.summary.table, summary.ci(ci=ci, res=res, event=", j, point, xscale2, "))", sep=""))	
			}
		}
	} else {
		if (i == 1){
			if(plotline>0) doItAndPrint(paste("ci.summary.table <- summary.ci(ci=ci, res=res, event=", plotline, point, xscale2, ")", sep=""))	
		} else {
			if(plotline>0) doItAndPrint(paste("ci.summary.table <- rbind(ci.summary.table, summary.ci(ci=ci, res=res, event=", plotline, point, xscale2, "))", sep=""))	
		}
	}
	}
   }
	if(nvar==1 && plotline>0 && posthoc!=""){
		command <- paste("pairwise.gray.test(", subdataSet, "$", timetoevent, ", ", subdataSet, "$", event, ", ", subdataSet, "$", group[i], ', "', dataSet, '", p.adjust.method="', posthoc, '", endpoint=', plotline, ")", sep="")
		doItAndPrint(command)
	}
	doItAndPrint("remove(res)")
    }
#	if(plotline>0){
		if(nevent>1) doItAndPrint("ci.summary.table")
#		doItAndPrint("remove(ci.summary.table)")	
#	}
	doItAndPrint("remove(ci)")
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cuminc", apply="StatMedCumInc", reset="StatMedCumInc")
    tkgrid(getFrame(timetoeventBox), labelRcmdr(variablesFrame, text="    "), getFrame(eventBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(variables2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables"), fg="blue"), sticky="w")
#    tkgrid(getFrame(groupBox), labelRcmdr(variables2Frame, text="    "), getFrame(strataBox), sticky="nw")
    tkgrid(getFrame(groupBox), labelRcmdr(variables2Frame, text="    "), sticky="nw")
	tkgrid(variables2Frame, sticky="nw")
	tkgrid(lineFrame, labelRcmdr(plotoptionFrame, text="   "), placeFrame, labelRcmdr(plotoptionFrame, text="   "), xscaleFrame, labelRcmdr(plotoptionFrame, text="   "), posthocFrame, sticky="w")
	tkgrid(plotoptionFrame, sticky="nw")
#    tkgrid(color, sticky="w")

	tkgrid(censor, labelRcmdr(plotoption2Frame, text="  "), atrisk, sticky="w")
	tkgrid(plotoption2Frame, sticky="nw")

#	tkgrid(labelRcmdr(plotoption2Frame, text=""), censor, atrisk, sticky="w")
#	tkgrid(plotoption2Frame, sticky="nw")

	tkgrid(labelRcmdr(ploteventFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Code of event to show cumulative incidence rate")), ploteventField,  sticky = "w")
  tkgrid(ploteventFrame, sticky="w")	
	tkgrid(labelRcmdr(pointFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Time point to show survival rate")), pointField,  sticky = "w")
  tkgrid(pointFrame, sticky="w")	
	tkgrid(labelRcmdr(xlimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimField, sticky = "w")
#  tkgrid(xlimFrame, sticky="w")
	tkgrid(labelRcmdr(ylimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimField, sticky = "w")
# tkgrid(ylimFrame, sticky="w")
tkgrid(xlimFrame, labelRcmdr(axis2Frame, text="  "), ylimFrame, sticky="w")
	tkgrid(labelRcmdr(xlabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis lavel")), xlabelField, sticky = "w")
	tkgrid(labelRcmdr(ylabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis lavel")), ylabelField, sticky = "w")
   tkgrid(xlabelFrame, labelRcmdr(axis2Frame, text="  "), ylabelFrame, sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Code of event to show cumulative incidence rate")), ploteventEntry, sticky="w")
#	tkgrid.configure(ploteventEntry, sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Time point to show cumulative incidence rate")), pointEntry, sticky="w")
#	tkgrid.configure(pointEntry, sticky="w")
#	tkgrid(labelRcmdr(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Cumulative incidence rate shown only when one event specified"), fg="blue"), sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimEntry, sticky="w")
#	tkgrid.configure(xlimEntry, sticky="w")
#	tkgrid(tklabel(axisFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimEntry, sticky="w")
#	tkgrid.configure(ylimEntry, sticky="w")
    tkgrid(axisFrame, sticky="w")
    tkgrid(axis2Frame, sticky="w")
  StatMedSubsetBox(model=TRUE)
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedStackCumInc  <- function(){
defaults <- list(event = "", timetoevent = "", group = "", atrisk = 0, xlim = "<auto>", ylim = "<auto>", xlabel = "<auto>", ylabel = "<auto>", subset = "")
dialog.values <- getDialog("StatMedStackCumInc", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Stacked cumulative incidences"))
    variablesFrame <- tkframe(top)
    eventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Status indicator (censor=0, event=1,2,3...) (pick one)"), listHeight=7, initialSelection=varPosn(dialog.values$event, "all"))
    timetoeventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time-to-event variable (pick one)"), listHeight=7, initialSelection=varPosn(dialog.values$timetoevent, "all"))
    variables2Frame <- tkframe(top)
    groupBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), listHeight=6, initialSelection=varPosn(dialog.values$group, "all"))
    plotoptionFrame <- tkframe(top)

	checkBoxes(frame="plotoptionFrame", boxes="atrisk", initialValues=dialog.values$atrisk,labels=gettext(domain="R-RcmdrPlugin.EZR","Show number at risk"), title=gettext(domain="R-RcmdrPlugin.EZR","Options"))	
	
#    checkBoxes(window=plotoptionFrame, frame="atrisk", boxes=c("atrisk"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show number at risk")))

	xlimFrame <- tkframe(plotoptionFrame)
	xlimVariable <- tclVar(dialog.values$xlim)
	xlimField <- ttkentry(plotoptionFrame, width="20", textvariable=xlimVariable)
	ylimFrame <- tkframe(plotoptionFrame)
	ylimVariable <- tclVar(dialog.values$ylim)
	ylimField <- ttkentry(plotoptionFrame, width="20", textvariable=ylimVariable)
	xlabelFrame <- tkframe(plotoptionFrame)
	xlabelVariable <- tclVar(dialog.values$xlabel)
	xlabelField <- ttkentry(plotoptionFrame, width="20", textvariable=xlabelVariable)
	ylabelFrame <- tkframe(plotoptionFrame)
	ylabelVariable <- tclVar(dialog.values$ylabel)
	ylabelField <- ttkentry(plotoptionFrame, width="20", textvariable=ylabelVariable)

#	xlim <- tclVar("<auto>")
#	xlimEntry <- ttkentry(plotoptionFrame, width="20", textvariable=xlim)
#	ylim <- tclVar("<auto>")
#	ylimEntry <- ttkentry(plotoptionFrame, width="20", textvariable=ylim)
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Stacked cumulative incidences"), "#####", sep=""))
    event <- getSelection(eventBox)
    timetoevent <- getSelection(timetoeventBox)
    group <- getSelection(groupBox)
	dataSet <- activeDataSet()
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
		subdataSet <- dataSet
		subset <- ""
    }
    else{
		subdataSet <- paste("subset(", dataSet, ", ", subset, ")", sep="")
		subset <- paste(", subset=", subset, sep="")
    }
    atrisk <- tclvalue(atriskVariable)
	xlim <- tclvalue(xlimVariable)
	ylim <- tclvalue(ylimVariable)
	if (xlim == "<auto>") {
		xlim <- ""
	} else {
		xlim <- paste(", xlim=c(", xlim, ")", sep="")
	}
	if (ylim == "<auto>") {
		ylim <- ", ylim=c(0, 1)"
	} else {
		ylim <- paste(", ylim=c(", ylim, ")", sep="")
	}
	xlabel <- tclvalue(xlabelVariable)
	ylabel <- tclvalue(ylabelVariable)
	if (xlabel == "<auto>") {
		xlabel <- paste(', xlab="', timetoevent, '"', sep="")
	} else {
		xlabel <- paste(', xlab="', xlabel, '"', sep="")
	}
	if (ylabel == "<auto>") {
		ylabel <- ', ylab="Probability"'
	} else {
		ylabel <- paste(', ylab="', ylabel, '"', sep="")
	}
putDialog("StatMedStackCumInc", list(event = event, timetoevent = timetoevent, group = group, atrisk = atrisk, xlim = tclvalue(xlimVariable), ylim = tclvalue(ylimVariable), xlabel = tclvalue(xlabelVariable), ylabel = tclvalue(ylabelVariable), subset = tclvalue(subsetVariable)))
    if (length(event) != 1) {
      errorCondition(recall=StatMedStackCumInc, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one status indicator (censor=0, event=1,2,3...)"))
      return()
    }
    if (length(timetoevent) != 1) {
      errorCondition(recall=StatMedStackCumInc, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one time-to-event variable"))
      return()
    }	
    closeDialog()
    Library("survival")
    Library("cmprsk")
#    library(survival)
#    library(cmprsk)
    nvar <- length(group)
    if (nvar == 0){
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
	doItAndPrint(paste("stackcuminc(", subdataSet, "$", timetoevent, ", ", subdataSet, "$", event, xlim, ylim,  xlabel, ylabel, ", atrisk=", atrisk, ")", sep=""))
   } else {
	groups <- eval(parse(text=paste("levels(factor(", subdataSet, "$", group, "))", sep="")))
	for (i in groups){
		sub2dataSet <- paste("subset(", subdataSet, ", ", group, "=='", i, "')", sep="")	
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", paste(substring(get("par.option", envir=.GlobalEnv), 1, nchar(get("par.option", envir=.GlobalEnv))-8), "2.5,1,0)", sep=""), ")", sep=""))}
		doItAndPrint(paste("stackcuminc(", sub2dataSet, "$", timetoevent, ", ", sub2dataSet, "$", event, xlim, ylim,  xlabel, ylabel, ", atrisk=", atrisk, ", main='", group, " = ", i, "')", sep=""))
	}
   }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cuminc", apply="StatMedStackCumInc", reset="StatMedStackCumInc")
    tkgrid(getFrame(timetoeventBox), labelRcmdr(variablesFrame, text="    "), getFrame(eventBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(getFrame(groupBox), labelRcmdr(variables2Frame, text="    "), sticky="nw")
	tkgrid(variables2Frame, sticky="nw")
	tkgrid(plotoptionFrame, sticky="nw")
#	tkgrid(atrisk, labelRcmdr(plotoptionFrame, text=""), sticky="w")
	tkgrid(plotoptionFrame, sticky="nw")

	tkgrid(labelRcmdr(xlimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimField, sticky = "w")
#  tkgrid(xlimFrame, sticky="w")
	tkgrid(labelRcmdr(ylimFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimField, sticky = "w")
#  tkgrid(ylimFrame, sticky="w")
tkgrid(xlimFrame, labelRcmdr(plotoptionFrame, text="  "), ylimFrame, sticky="w")

	tkgrid(labelRcmdr(xlabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis lavel")), xlabelField, sticky = "w")
	tkgrid(labelRcmdr(ylabelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis lavel")), ylabelField, sticky = "w")
   tkgrid(xlabelFrame, labelRcmdr(plotoptionFrame, text="  "), ylabelFrame, sticky="w")

  #	tkgrid(tklabel(plotoptionFrame, text=gettext(domain="R-RcmdrPlugin.EZR","X axis range(Min, Max) Ex: 0, 365")), xlimEntry, sticky="w")
#	tkgrid.configure(xlimEntry, sticky="w")
#	tkgrid(tklabel(plotoptionFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Y axis range(Min, Max) Ex: 0.8, 1.0")), ylimEntry, sticky="w")
#	tkgrid.configure(ylimEntry, sticky="w")
  StatMedSubsetBox(model=TRUE)
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedCrr  <- function(){
defaults <- list(event = "", timetoevent = "", group = "", fcode = 1, wald = 0, stepwise1 = 0, stepwise2 = 0, stepwise3 = 0, subset = "")
dialog.values <- getDialog("StatMedCrr", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Fine-Gray proportional hazard regression for competing events"))
    variablesFrame <- tkframe(top)
	fcodeFrame <- tkframe(top)
    eventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Status indicator (censor=0, event=1,2,3...) (pick one)"), listHeight=10, initialSelection=varPosn(dialog.values$event, "all"))
    timetoeventBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time-to-event variable (pick one)"), listHeight=10, initialSelection=varPosn(dialog.values$timetoevent, "all"))
    groupBox <- variableListBox(top, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Explanatory (non-character) variables (pick one or more)"), listHeight=10, initialSelection=varPosn(dialog.values$group, "all"))

	fcodeFrame <- tkframe(top)
	fcodeVariable <- tclVar(dialog.values$fcode)	
	fcodeField <- ttkentry(fcodeFrame, width="20", textvariable=fcodeVariable)
	
#	fcode <- tclVar("1")
#	fcodeEntry <- ttkentry(fcodeFrame, width="10", textvariable=fcode)

  	optionsFrame <- tkframe(top)	
	checkBoxes(frame="optionsFrame", boxes=c("wald", "stepwise1", "stepwise2", "stepwise3"), initialValues=c(dialog.values$wald, dialog.values$stepwise1, dialog.values$stepwise2, dialog.values$stepwise3),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Wald test for overall p-value for factors with >2 levels", "Stepwise selection based on AIC", "Stepwise selection based on BIC", "Stepwise selection based on p-value")))	

#	waldVariable <- tclVar("0")
#	waldCheckBox <- tkcheckbutton(optionsFrame, variable=waldVariable)
#	stepwise1Variable <- tclVar("0")
#	stepwise1CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise1Variable)
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Fine-Gray proportional hazard regression for competing events"), "#####", sep=""))
    event <- getSelection(eventBox)
    timetoevent <- getSelection(timetoeventBox)
    group <- getSelection(groupBox)
    fcode <- tclvalue(fcodeVariable)
	wald <- tclvalue(waldVariable)
	stepwise1 <- tclvalue(stepwise1Variable)
	stepwise2 <- tclvalue(stepwise2Variable)
	stepwise3 <- tclvalue(stepwise3Variable)
	subset <- tclvalue(subsetVariable)	
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
      subset <- ""
    }
    else{
      subset <- paste(", subset=", subset, sep="")
    }	
putDialog("StatMedCrr", list(event = event, timetoevent = timetoevent, group = group, fcode = fcode, wald = wald, stepwise1 = stepwise1, stepwise2 = stepwise2, stepwise3 = stepwise3, subset = tclvalue(subsetVariable)))
    if (length(event) != 1) {
      errorCondition(recall=StatMedCrr, )
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one status indicator (censor=0, event=1,2,3...)")
      return()
    }
    if (length(timetoevent) != 1) {
      errorCondition(recall=StatMedCrr, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick one time-to-event variable"))
      return()
    }
    if (length(group) == 0) {
      errorCondition(recall=StatMedCrr, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick at least one explanatory variable"))
      return()
    }
    if (length(fcode) == 0) {
      errorCondition(recall=StatMedCrr, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Specify one event of interest"))
      return()
    }    	
    closeDialog()
    Library("survival")
    Library("cmprsk")
		Library("aod")
#    library(survival)
#    library(cmprsk)
	dataSet <- activeDataSet()
    nvar <- length(group)
    command <- paste("cov.matrix <- cbind(", group[1], "=", dataSet, "$", group[1], sep="")
    if (nvar >= 2){
		for (i in 2:nvar) {
		command <- paste(command, ", ", group[i], "=", dataSet, "$", group[i], sep="")
		}
    }
    command <- paste(command, ")", sep="")
    doItAndPrint(command)
    command2 <- paste("crr <- with(", dataSet, ", crr(", timetoevent, ", ", event, ", cov.matrix, failcode=", fcode, ", cencode=0", subset, ", na.action = na.omit))", sep="")
    doItAndPrint(command2)
    doItAndPrint("summary(crr)")
	
	if(eval(parse(text="length(summary(crr)$coef[,1])"))==1){
		doItAndPrint("crr.table <- signif(cbind(t(summary(crr)$conf.int[,c(1,3,4)]), p.value=summary(crr)$coef[,5]), digits=4)")
		doItAndPrint(paste('rownames(crr.table) <- "', group[1], '"', sep=""))
		doItAndPrint('colnames(crr.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	} else {
		doItAndPrint("crr.table <- signif(cbind(summary(crr)$conf.int[,c(1,3,4)], summary(crr)$coef[,5]), digits=4)")
		doItAndPrint('colnames(crr.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	}
#	doItAndPrint("crr.table <- signif(crr.table, digits=3)")
	doItAndPrint("crr.table")
	if (wald==1) doItAndPrint("waldtest.crr(crr, rownames(crr.table))")

	if (stepwise1 == 1 | stepwise2 == 1 | stepwise3 == 1){
		command <- paste("TempDF <- with(", ActiveDataSet(), ", ", ActiveDataSet(), "[complete.cases(", paste(group, collapse=", "), "),])", sep="")
		doItAndPrint(command)
		command <- paste('cov <- c("', group[1], '"', sep="")
		if (nvar >= 2){
			for (i in 2:nvar) {
			command <- paste(command, ', "', group[i], '"', sep="")
			}
		}
		command <- paste(command, ')', sep="")
		doItAndPrint(command)
		subset <- tclvalue(subsetVariable)
	    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
			|| trim.blanks(subset) == ""){
			subset <- ""
		}else{
			subset <- paste(", subset='", trim.blanks(subset), "'", sep="")
		}
	}
	if (stepwise1 == 1){
		doItAndPrint(paste('step.AIC.crr(crr, cov, "TempDF", BIC=0, waldtest=', wald, subset, ")", sep=""))
	}
	if (stepwise2 == 1){
		doItAndPrint(paste('step.AIC.crr(crr, cov, "TempDF", BIC=1, waldtest=', wald, subset, ")", sep=""))
	}
	if (stepwise3 == 1){
		doItAndPrint(paste('step.p.crr(crr, cov, "TempDF", wald=', wald, subset, ")", sep=""))
	}
	
	doItAndPrint("remove(crr)")
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="crr", apply="StatMedCrr", reset="StatMedCrr")
    tkgrid(getFrame(timetoeventBox), labelRcmdr(variablesFrame, text="  "), getFrame(eventBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	
	tkgrid(labelRcmdr(fcodeFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Input code of event of interest"), fg="blue"), fcodeField, sticky = "w")
  tkgrid(fcodeFrame, sticky="w")

#	tkgrid(tklabel(fcodeFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Input code of event of interest"), fg="blue"), fcodeEntry, sticky="w")
#	tkgrid.configure(fcodeEntry, sticky="w")
#   tkgrid(fcodeFrame, sticky="w")
	tkgrid(getFrame(groupBox), sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Dummy variables required for factors of more than 2 groups"), fg="blue"), sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Wald test for overall p-value for factors with >2 levels")), waldCheckBox, sticky="w")
# 	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on p-value")), stepwise1CheckBox, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
	StatMedSubsetBox(model=TRUE)
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedCoxTD  <- function(){
# add the class coxph to the modelClasses
    xx <- getRcmdr("modelClasses")
    bolCoxphExists = FALSE
    for(ii in 1:length(xx)){if (xx[ii] == "coxph") bolCoxphExists = TRUE}
    if (bolCoxphExists == FALSE) putRcmdr("modelClasses", c(getRcmdr("modelClasses"), "coxph"))

defaults <- list(SurvivalTimeVariable = "", StatusVariable = "", rhs = "", waldVariable = 0,  prophazVariable = 0, basecurveVariable = 0, actmodelVariable = 0, stepwise1Variable = 0, stepwise2Variable = 0, stepwise3Variable = 0, subset = "", timedependentcovariate = NULL, timepositive = NULL, timenegative = NULL)
dialog.values <- getDialog("StatMedCoxTD", defaults)
currentFields$SurvivalTimeVariable <- dialog.values$SurvivalTimeVariable	
currentFields$StatusVariable <- dialog.values$StatusVariable
currentFields$rhs <- dialog.values$rhs
currentFields$subset <- dialog.values$subset	

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Cox proportional hazard regression with time-dependent covariate"))
  .activeModel <- ActiveModel()
  currentModel <- if (!is.null(.activeModel)) 
        class(get(.activeModel, envir=.GlobalEnv))[1] == "coxph"
#    eval(parse(text=paste("class(", .activeModel, ")[1] == 'coxph'", sep="")), 
#         envir=.GlobalEnv) 
    else FALSE
	
	currentModel <- TRUE
	#  if(currentModel){
#    currentFields <- formulaFields(eval(parse(text=.activeModel), 
#     envir=.GlobalEnv))
#    if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
#  }
  UpdateModelNumber()
  modelName <- tclVar(paste("CoxModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="30", textvariable=modelName)
    variablesFrame <- tkframe(top)
    timedependentcovariateBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time-dependent (TD) covariate (pick one)"), listHeight=7, initialSelection=varPosn(dialog.values$timedependentcovariate, "all"))
    timepositiveBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time when TD covariate changes from 0 to 1(pick one)"), listHeight=7, initialSelection=varPosn(dialog.values$timepositive, "all"))
    timenegativeBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Time when TD covariate changes from 1 to 0(pick one)"), listHeight=7, initialSelection=varPosn(dialog.values$timenegative, "all"))
	textFrame <- tkframe(top)
  	optionsFrame <- tkframe(top)
	
	checkBoxes(frame="checkboxFrame", boxes=c("wald", "prophaz", "basecurve", "actmodel", "stepwise1", "stepwise2", "stepwise3"), initialValues=c(dialog.values$waldVariable, dialog.values$prophazVariable, dialog.values$basecurveVariable, dialog.values$actmodelVariable, dialog.values$stepwise1Variabl, dialog.values$stepwise2Variabl, dialog.values$stepwise3Variabl),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Wald test for overall p-value for factors with >2 levels", "Test proportional hazards assumption","Show baseline survival curve", "Keep results as active model for further analyses", "Stepwise selection based on AIC", "Stepwise selection based on BIC", "Stepwise selection based on p-value")))	
	
#	waldVariable <- tclVar("0")
#	waldCheckBox <- tkcheckbutton(optionsFrame, variable=waldVariable)
#	prophazVariable <- tclVar("0")
#	prophazCheckBox <- tkcheckbutton(optionsFrame, variable=prophazVariable)
#	basecurveVariable <- tclVar("0")
#	basecurveCheckBox <- tkcheckbutton(optionsFrame, variable=basecurveVariable)
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
#	stepwise1Variable <- tclVar("0")
#	stepwise1CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise1Variable)
#	stepwise2Variable <- tclVar("0")
#	stepwise2CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise2Variable)
#	stepwise3Variable <- tclVar("0")
#	stepwise3CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise3Variable)
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Cox proportional hazard regression with time-dependent covariate"), "#####", sep=""))
#    XXX <- getSelection(timeBox)
    modelValue <- trim.blanks(tclvalue(modelName))
	timedependentcovariate <- getSelection(timedependentcovariateBox)
	timepositive <- getSelection(timepositiveBox)
	timenegative <- getSelection(timenegativeBox)
		wald <- tclvalue(waldVariable)
		prophaz <- tclvalue(prophazVariable)
		basecurve <- tclvalue(basecurveVariable)
		actmodel <- tclvalue(actmodelVariable)
		stepwise1 <- tclvalue(stepwise1Variable)
		stepwise2 <- tclvalue(stepwise2Variable)
		stepwise3 <- tclvalue(stepwise3Variable)
		subset <- tclvalue(subsetVariable)
#    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
#        || trim.blanks(subset) == ""){
#      subset <- ""
#      putRcmdr("modelWithSubset", FALSE)
#    }
#    else{
#      subset <- paste(", subset=", subset, sep="")
#      putRcmdr("modelWithSubset", TRUE)
#    }	
putDialog("StatMedCoxTD", list(SurvivalTimeVariable = tclvalue(SurvivalTimeVariable), StatusVariable = tclvalue(StatusVariable), rhs = tclvalue(rhsVariable), waldVariable = wald,  prophazVariable = prophaz, basecurveVariable = basecurve, actmodelVariable = actmodel, stepwise1Variable = stepwise1, stepwise2Variable = stepwise2, stepwise3Variable = stepwise3, subset=tclvalue(subsetVariable), timedependentcovariate = timedependentcovariate, timepositive = timepositive, timenegative = timenegative))
		if (length(timedependentcovariate) == 0 || length(timepositive) == 0){
  	        errorCondition(recall=StatMedCoxTD, message=gettext(domain="R-RcmdrPlugin.EZR","Pick all required variables"))
            return()
        }
		if (length(timenegative) == 0){
			timenegative <- tclvalue(SurvivalTimeVariable)
		}
    if (!is.valid.name(modelValue)){
      errorCondition(recall=StatMedCoxTD, 
        message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue), model=TRUE)
      return()
    }	
#    check.empty <- gsub(" ", "", tclvalue(lhsVariable))
#    if ("" == check.empty) {
#      errorCondition(recall=StatMedCoxRegression,
#        message=gettext(domain="R-RcmdrPlugin.EZR","Left-hand side of model empty."), model=TRUE) 
#      return()
#    }
     check.empty <- gsub(" ", "", tclvalue(SurvivalTimeVariable))
     if ("" == check.empty) {
            errorCondition(recall=StatMedCoxTD, message=gettext(domain="R-RcmdrPlugin.EZR","Survival time variable of model empty."), model=TRUE)
            return()
     }
     check.empty <- gsub(" ", "", tclvalue(StatusVariable))
     if ("" == check.empty) {
            errorCondition(recall=StatMedCoxTD, message=gettext(domain="R-RcmdrPlugin.EZR","Status variable of model empty."), model=TRUE)
            return()
     }
	covariates <- ""
    check.empty <- gsub(" ", "", tclvalue(rhsVariable))
    if ("" == check.empty) {
		covariates <- "covariate_td"
    }
	else {
		covariates <- paste("covariate_td + ", tclvalue(rhsVariable), sep="")
	}
	if (is.element(modelValue, listCoxModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
        UpdateModelNumber(-1)
        StatMedCoxTD()
        return()
      }
    }
    closeDialog()
    Library("survival")
		Library("aod")
	dataSet <- activeDataSet()
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
        || trim.blanks(subset) == ""){
	  doItAndPrint(paste("TempDF <- ", dataSet, sep=""))
    }
   else{
	  doItAndPrint(paste("TempDF <- subset(", dataSet, ", ",subset, ")", sep=""))
    }		
#		doItAndPrint(paste("attach(", dataSet, ")"))
		command <-  paste("TempTD <- stsplit(TempDF, TempDF$", tclvalue(SurvivalTimeVariable), ", TempDF$", tclvalue(StatusVariable), ", TempDF$", timepositive, ", TempDF$", timedependentcovariate, ", TempDF$", timenegative, ")", sep="")
		result <- doItAndPrint(command)		
#    library(survival)
#    formula <- paste("Surv(", XXX, ", ", tclvalue(lhsVariable), ") ~ ", tclvalue(rhsVariable), sep="")
#     formula <- paste("Surv(", tclvalue(SurvivalTimeVariable), ", ", tclvalue(StatusVariable), ")~ ", tclvalue(rhsVariable), sep="")
    formula <- paste("Surv(start_td, stop_td, endpoint_td==1) ~ ", covariates, sep="")	 
#    command <- paste("coxph(", formula,
#      ", data=TempTD", subset, ', method="breslow")', sep="")
    command <- paste("coxph(", formula,
      ', data=TempTD, method="breslow")', sep="")	 	  
#    logger(paste(modelValue, " <- ", command, sep=""))
#    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
    doItAndPrint(paste("(res <- summary(", modelValue, "))", sep=""))
#	doItAndPrint(paste("res <- ", command, sep=""))
#	doItAndPrint("res <- summary(res)")
	if(eval(parse(text="length(res$coefficients[,1])"))==1){
		doItAndPrint("cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)")
		doItAndPrint("rownames(cox.table) <- rownames(res$coefficients)")
		doItAndPrint('colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	} else {
		doItAndPrint("cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)")
		doItAndPrint("cox.table <- data.frame(cox.table)")
		doItAndPrint('colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	}
#	doItAndPrint("cox.table <- signif(cox.table, digits=3)")
	doItAndPrint("cox.table")
	if (wald==1) doItAndPrint(paste("waldtest(", modelValue, ")", sep=""))
	if (prophaz == 1){
			doItAndPrint(paste("print(cox.zph(", modelValue, "))", sep=""))
	}
	if (basecurve ==1){
			doItAndPrint(paste("plot(survfit(", modelValue, "))", sep=""))
	}

	if (stepwise1 == 1 | stepwise2 == 1 | stepwise3 == 1){
		x <- strsplit(tclvalue(rhsVariable), split="\\+")
		if (length(x[[1]]>0)){
			command <- paste("TempDF <- with(TempTD, TempTD[complete.cases(covariate_td, ", paste(x[[1]], collapse=","), "),])", sep="")
		} else{
			command <- ("TempDF <- with(TempTD, TempTD[complete.cases(covariate_td),])")		
		}	
		doItAndPrint(command)
		command <- paste("coxph(", formula, ", data=TempDF", subset, ', method="breslow")', sep="")
		doItAndPrint(paste(modelValue, " <- ", command, sep=""))
		}
	if (stepwise1 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="AIC")', sep=""))
			doItAndPrint("summary(res)")
			doItAndPrint("res2 <- summary(res)")
			if(eval(parse(text="length(res2$coefficients[,1])"))==1){
				doItAndPrint("cox.table <- signif(cbind(t(res2$conf.int[,c(1,3,4)]), p.value=res2$coefficients[,5]), digits=4)")
				doItAndPrint("rownames(cox.table) <- rownames(res2$coefficients)")
				doItAndPrint('colnames(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')	
				doItAndPrint("cox.table")
			} else if(eval(parse(text="length(res2$coefficients[,1])"))>1){
				doItAndPrint("cox.table <- signif(cbind(res2$conf.int[,c(1,3,4)], res2$coefficients[,5]), digits=4)")
				doItAndPrint("cox.table <- data.frame(cox.table)")
				doItAndPrint('names(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
				doItAndPrint("cox.table")
			}
			if (wald==1) doItAndPrint("waldtest(res)")
	}
	if (stepwise2 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="BIC")', sep=""))
			doItAndPrint("summary(res)")
			doItAndPrint("res2 <- summary(res)")
			if(eval(parse(text="length(res2$coefficients[,1])"))==1){
				doItAndPrint("cox.table <- signif(cbind(t(res2$conf.int[,c(1,3,4)]), p.value=res2$coefficients[,5]), digits=4)")
				doItAndPrint("rownames(cox.table) <- rownames(res2$coefficients)")
				doItAndPrint('colnames(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')	
				doItAndPrint("cox.table")
			} else if(eval(parse(text="length(res2$coefficients[,1])"))>1){
				doItAndPrint("cox.table <- signif(cbind(res2$conf.int[,c(1,3,4)], res2$coefficients[,5]), digits=4)")
				doItAndPrint("cox.table <- data.frame(cox.table)")
				doItAndPrint('names(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
				doItAndPrint("cox.table")
			}
			doItAndPrint("cox.table")
			if (wald==1) doItAndPrint("waldtest(res)")
	}
	if (stepwise3 == 1){
			subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")
				|| trim.blanks(subset) == ""){
				subset <- ""
			}
			else{
				subset <- paste(", subset='", trim.blanks(subset), "'", sep="")
			}
			doItAndPrint(paste('step.p.coxtd(', modelValue, ', "TempDF", wald=', wald, subset, ")", sep=""))
	}
	doItAndPrint("remove(res)")
	if (actmodel==1) activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="coxph", model=TRUE, apply="StatMedCoxTD", reset="StatMedCoxTD")
  tkgrid(tklabel(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  
  modelFormulaCox()
  StatMedSubsetBox(model=TRUE)
  
 tkgrid(getFrame(xBox), sticky="w")
  tkgrid(outerOperatorsFrame, sticky="w")
  tkgrid(formulaFrame, sticky="w")
  	tkgrid(labelRcmdr(textFrame, text=paste("            ", gettext(domain="R-RcmdrPlugin.EZR","Stratifing variable: + strata(#####)"), sep="")), sticky="w")  
	tkgrid(textFrame, sticky="w")
	tkgrid(getFrame(timedependentcovariateBox), labelRcmdr(variablesFrame, text="  "), getFrame(timepositiveBox), getFrame(timenegativeBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")

 tkgrid(checkboxFrame, sticky="w")

#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Wald test for overall p-value for factors with >2 levels")), waldCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Test proportional hazards assumption")), prophazCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Show baseline survival curve")), basecurveCheckBox, sticky="w")	
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on AIC")), stepwise1CheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on BIC")), stepwise2CheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on p-value")), stepwise3CheckBox, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
}


StatMedROC <- function(){
defaults <- list(response=NULL, predictor=NULL, threshold=1, direction="auto", best="youden", cost="1", prevalence="0.5", subset = "")
dialog.values <- getDialog("StatMedROC", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE

    Library("pROC")
    Library("methods")
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","ROC curve analysis for quantitative test"))
    variablesFrame <- tkframe(top)
    responseBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Response (encoded as 0 or 1) (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$response, "all"))
    predictorBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Predictor (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$predictor, "all"))
	optionsFrame <- tkframe(top)
 radioButtons(optionsFrame, name="direction", buttons=c("auto", "higher", "lower"), initialValue=dialog.values$direction, 
    values=c("auto", "<", ">"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Automatic", ">=threshold as positive", "<=threshold as positive")),title=gettext(domain="R-RcmdrPlugin.EZR","Direction for comparison"))
 radioButtons(optionsFrame, name="best", buttons=c("youden", "closest.topleft"), initialValue=dialog.values$best, 
    values=c("youden", "closest.topleft"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Maximum sum of sensitivity + specificity", "Closest to the top-left corner")),title=gettext(domain="R-RcmdrPlugin.EZR","Optimal threshold"))
	
	checkBoxFrame <- tkframe(top)
	checkBoxes(frame="checkBoxFrame", boxes="threshold", initialValues=dialog.values$threshold,labels=gettext(domain="R-RcmdrPlugin.EZR","Show optimal threshold in graph"))
	
#    checkBoxes(frame="threshold", boxes=c("thres"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show optimal threshold in graph")))
 
    costFrame <- tkframe(top)
    costVariable <- tclVar(dialog.values$cost)
    costField <- ttkentry(costFrame, width="6", textvariable=costVariable)
    prevalenceFrame <- tkframe(top)
    prevalenceVariable <- tclVar(dialog.values$prevalence)
    prevalenceField <- ttkentry(prevalenceFrame, width="6", textvariable=prevalenceVariable)
	
#	costFrame <- tkframe(top)
#    costVariable <- tclVar("1")
#    costField <- ttkentry(costFrame, width="8", textvariable=costVariable)
#    prevalenceFrame <- tkframe(top)
#    prevalenceVariable <- tclVar("0.5")
#    prevalenceField <- ttkentry(prevalenceFrame, width="8", textvariable=prevalenceVariable)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","ROC curve analysis for quantitative test"), "#####", sep=""))
        response <- getSelection(responseBox)
        predictor <- getSelection(predictorBox)
	direction <- tclvalue(directionVariable)
	best <- tclvalue(bestVariable)
	cost <- tclvalue(costVariable)
	prevalence <- tclvalue(prevalenceVariable)
			subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
			}
putDialog("StatMedROC", list(response=response, predictor=predictor, threshold=tclvalue(thresholdVariable), direction=direction, best=best, cost=cost, prevalence=prevalence, subset = tclvalue(subsetVariable)))
        if (length(response) == 0 || length(predictor) == 0){
            errorCondition(recall=StatMedROC, message=gettext(domain="R-RcmdrPlugin.EZR","You must select two variables."))
            return()
            }
        closeDialog()
   	if (tclvalue(thresholdVariable) == "1"){
		pt <- paste(', print.thres="best", print.thres.best.method="', best, '", print.thres.best.weights=c(', cost, ", ", prevalence, ")", sep="")
		cpt <- paste(', "best", best.method="', best, '", best.weights=c(', cost, ", ", prevalence, ")", sep="")
		}
	else{
		pt <- ""
		cpt <- ""
	}
	command <- paste("ROC <- roc(", response, "~", predictor, ", data=", subset1, ActiveDataSet(), subset2,
            	', ci=TRUE, direction="', direction, '")', sep="")
    doItAndPrint(command)
#	doItAndPrint("if(ROC$thresholds[1]==-Inf) {ROC$thresholds[1:(length(levels(factor(ROC$predictor))))] <- as.numeric(levels(factor(ROC$predictor)))}")
#	doItAndPrint("if(ROC$thresholds[1]==Inf) {ROC$thresholds[1:(length(levels(factor(ROC$predictor))))] <- rev(as.numeric(levels(factor(ROC$predictor))))}")	

	doItAndPrint("if(ROC$thresholds[1]==-Inf){thre <- c(unique(sort(ROC$predictor)), Inf)}")
	doItAndPrint("if(ROC$thresholds[1]==Inf){thre <- c(unique(sort(ROC$predictor, decreasing=TRUE)), -Inf)}")

	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
#	doItAndPrint('plot(ROC$thresholds, ROC$sensitivities, ylim=c(0,1), type="l", ylab="Sensitivity/Specificity", xlab="Threshold")')
	doItAndPrint('plot(thre, ROC$sensitivities, ylim=c(0,1), type="l", ylab="Sensitivity/Specificity", xlab="Threshold")')
	doItAndPrint("par(new=T)")
#	doItAndPrint('plot(ROC$thresholds,ROC$specificities, ylim=c(0,1), type="l", lty=2, ylab="", xlab="")')
	doItAndPrint('plot(thre, ROC$specificities, ylim=c(0,1), type="l", lty=2, ylab="", xlab="")')
	doItAndPrint('legend("bottom", horiz=TRUE, c("Sensitivity", "Specificity"), lty=1:2, box.lty=0)')			
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}

	doItAndPrint(paste("co <- coords(ROC", cpt, ")", sep=""))
	if(eval(parse(text="class(co)"))=="matrix"){
		doItAndPrint("if(ROC$thresholds[1]==-Inf){co[1,] <- min(ROC$predictor[ROC$predictor>co[1,]])}")	###Change to exact values
		doItAndPrint("if(ROC$thresholds[1]==Inf)co[1,] <- max(ROC$predictor[ROC$predictor<co[1,]])")	###Change to exact values
		doItAndPrint("plot(ROC, print.thres=co[1,], grid=TRUE)")
		} else {
		doItAndPrint("if(ROC$thresholds[1]==-Inf){co[1] <- min(ROC$predictor[ROC$predictor>co[1]])}")	###Change to exact values
		doItAndPrint("if(ROC$thresholds[1]==Inf)co[1] <- max(ROC$predictor[ROC$predictor<co[1]])")	###Change to exact values
		doItAndPrint("plot(ROC, print.thres=co[1], grid=TRUE)")
	}
#	doItAndPrint('coords(ROC, "all")')
	doItAndPrint("if(ROC$thresholds[1]==-Inf){coords(ROC, x=c(-Inf, unique(sort(ROC$predictor)), Inf))}")
	doItAndPrint("if(ROC$thresholds[1]==Inf){coords(ROC, x=c(Inf, unique(sort(ROC$predictor, decreasing=TRUE)), -Inf))}")

	if(eval(parse(text="ROC$direction"))==">"){
		logger(gettext(domain="R-RcmdrPlugin.EZR","### <= threshold is considered positive"))
	}else{
		logger(gettext(domain="R-RcmdrPlugin.EZR","### >= threshold is considered positive"))
	}
	doItAndPrint("ROC")
	doItAndPrint('cat(gettext(domain="R-RcmdrPlugin.EZR","Area under the curve"), signif(ROC$auc[1], digits=3), gettext(domain="R-RcmdrPlugin.EZR","95% CI"), signif(ROC$ci[1], digits=3), "-", signif(ROC$ci[3], digits=3), "\n")')
	doItAndPrint("remove(ROC)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="roc", apply="StatMedROC", reset="StatMedROC")
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), getFrame(predictorBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
#    tkgrid(thresholdFrame, sticky="w")
	tkgrid(checkBoxFrame, sticky="w")
#    tkgrid(directionFrame, sticky="w")
#    tkgrid(bestFrame, sticky="w")

    tkgrid(directionFrame, labelRcmdr(optionsFrame, text="   "), bestFrame, sticky="w")
	tkgrid(optionsFrame, sticky="nw")

    tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","Supply weights if false positive and false negative predictions are not equivalent"), fg="blue"), sticky="w")  
    tkgrid(labelRcmdr(costFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Cost of of false negative classification")), costField, sticky="w")
    tkgrid(costFrame, sticky="w")
    tkgrid(labelRcmdr(prevalenceFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Prevalence")), prevalenceField, sticky="w")
    tkgrid(prevalenceFrame, sticky="w")
	StatMedSubsetBox(model=TRUE)
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
    }


StatMedROCtest <- function(){
    Library("pROC")
defaults <- list(response=NULL, predictor1=NULL, predictor2=NULL, subset = "")
dialog.values <- getDialog("StatMedROCtest", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Compare two ROC curves"))
    variablesFrame <- tkframe(top)
    responseBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Response (encoded as 0 or 1) (pick one)"), listHeight=12, initialSelection=varPosn(dialog.values$response, "all"))
    variables2Frame <- tkframe(top)
    predictor1Box <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Predictor1 (pick one)"), listHeight=12, initialSelection=varPosn(dialog.values$predictor1, "all"))
    predictor2Box <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Predictor2 (pick one)"), listHeight=12, initialSelection=varPosn(dialog.values$predictor2, "all"))
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Compare two ROC curves"), "#####", sep=""))
        response <- getSelection(responseBox)
        predictor1 <- getSelection(predictor1Box)
        predictor2 <- getSelection(predictor2Box)
			subset <- tclvalue(subsetVariable)
			if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
				subset1 <- ""
				subset2 <- ""
				subset <- ""
			} else {
				subset1 <- "subset("
				subset2 <- paste(", ", subset, ")", sep="")
				subset <- paste(", subset=", subset, sep="")
			}
putDialog("StatMedROCtest", list(response=response, predictor1=predictor1, predictor2=predictor2, subset = tclvalue(subsetVariable)))
        if (length(response) == 0 || length(predictor1) == 0 || length(predictor2) == 0){
            errorCondition(recall=StatMedROCtest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select three variables."))
            return()
            }
        closeDialog()
		command <- paste("ROC1 <- roc(", response, "~", predictor1, ", data=", subset1, ActiveDataSet(), subset2, 
            	", ci=TRUE)", sep="")
        doItAndPrint(command)
		command <- paste("ROC2 <- roc(", response, "~", predictor2, ", data=", subset1, ActiveDataSet(), subset2, 
            	", ci=TRUE)", sep="")
        doItAndPrint(command)
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		doItAndPrint("plot(ROC1, lty=1)")
		doItAndPrint("plot(ROC2, lty=2, add=TRUE)")
		doItAndPrint(paste('legend("bottomright", c("', predictor1, '", "', predictor2, '"), lty=1:2, box.lty=0)', sep=""))
		command <- paste("(res <- roc.test(", response, "~", predictor1, "+", predictor2, ", data=", subset1, ActiveDataSet(), subset2, 
            	"))", sep="")
        doItAndPrint(command)
		doItAndPrint("roc.table <-  signif(cbind(res$estimate, res$p.value), digits=3)")
		doItAndPrint(paste('rownames(roc.table) <- c("', predictor1, '", "', predictor2, '")', sep=""))
		doItAndPrint('colnames(roc.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Area under the curve", "p.value"))')
		doItAndPrint('roc.table[2,2] <- ""')
		doItAndPrint("data.frame(roc.table)")
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="roc.test", apply="StatMedROCtest", reset="StatMedROCtest")
    tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text="    "), sticky="nw")
    tkgrid(getFrame(predictor1Box), labelRcmdr(variables2Frame, text="    "), getFrame(predictor2Box), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(variables2Frame, sticky="nw")
	StatMedSubsetBox(model=TRUE)
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
    }
		
	
StatMedTest <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Accuracy of qualitative test"))
	textFrame <- tkframe(top)
	variableFrame <- tkframe(top)
	pospos <- tclVar("")
	posposEntry <- ttkentry(variableFrame, width="10", textvariable=pospos)
	posneg <- tclVar("")
	posnegEntry <- ttkentry(variableFrame, width="10", textvariable=posneg)
	variable2Frame <- tkframe(top)
	negpos <- tclVar("")
	negposEntry <- ttkentry(variable2Frame, width="10", textvariable=negpos)
	negneg <- tclVar("")
	negnegEntry <- ttkentry(variable2Frame, width="10", textvariable=negneg)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Accuracy of qualitative test"), "#####", sep=""))
		pospos <- tclvalue(pospos)
		posneg <- tclvalue(posneg)
		negpos <- tclvalue(negpos)
		negneg <- tclvalue(negneg)
		closeDialog()
	if (length(pospos) == 0 || length(posneg) == 0 || length(negpos) == 0 || length(negneg) == 0){
			errorCondition(recall=StatMedTest, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
#		library(epiR, quietly=TRUE)
#		command <- paste("res <- epi.tests(", pospos, ", ", posneg, ", ", negpos, ", ", negneg, ", conf.level = 0.95)", sep="")
#		result <- doItAndPrint(command)
#		doItAndPrint("summary.test <- round(rbind(res$se, res$sp, res$ppv, res$npv, res$da, res$lr.pos, res$lr.neg), 3)")
#		doItAndPrint('rownames(summary.test) <- c("Sensitivity", "Specificity", "Positive predictive value", "Negative predictive value", "Diagnstic accuracy", "Likelihood ratio of a positive test", "Likelihood ratio of a negative test")')
#		doItAndPrint('colnames(summary.test) <- c("Estimation", "Lower 95%CI", "Upper 95%CI")')
		doItAndPrint(paste(".Table <- matrix(c(", pospos, ", ", posneg, ", ", negpos, ", ", negneg, "), 2, 2, byrow=TRUE)", sep=""))
#		doItAndPrint('colnames(.Table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Disease positive", "Disease negative"))')
#		doItAndPrint('rownames(.Table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Test positive", "Test negative"))')
#		doItAndPrint(".Table")
		command <- "epi.tests(.Table, conf.level = 0.95)"
		doItAndPrint(command)
#		doItAndPrint("summary.test")
#		doItAndPrint("remove(summary.test)")
#		doItAndPrint("remove(res)")
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="epi.tests")
	tkgrid(labelRcmdr(textFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number      Disease (+)         (-)")), sticky="w")  
	tkgrid(textFrame, sticky="w")
	tkgrid(tklabel(variableFrame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Test (+)"), "         ", sep="")), posposEntry, posnegEntry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Test  (-)"), "         ", sep="")), negposEntry, negnegEntry, sticky="w")
	tkgrid.configure(posposEntry, sticky="w")
	tkgrid.configure(posnegEntry, sticky="w")
	tkgrid.configure(negposEntry, sticky="w")
	tkgrid.configure(negnegEntry, sticky="w")
	tkgrid(variableFrame, sticky="nw")
	tkgrid(variable2Frame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedPredictiveValue <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Compute positive and negative predictive values"))
	preprob <- tclVar("")
	preprobEntry <- ttkentry(top, width="20", textvariable=preprob)
	sens <- tclVar("")
	sensEntry <- ttkentry(top, width="20", textvariable=sens)
	spec <- tclVar("")
	specEntry <- ttkentry(top, width="20", textvariable=spec)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Compute positive and negative predictive values"), "#####", sep=""))
		preprob <- as.numeric(tclvalue(preprob))
		sens <- as.numeric(tclvalue(sens))
		spec <- as.numeric(tclvalue(spec))
		closeDialog()
		if (length(preprob) == 0 || length(sens) == 0 || length(spec) == 0){
			errorCondition(recall=StatMedPredictiveValue, message=gettext(domain="R-RcmdrPlugin.EZR","You
must select a variable."))
			return()
		}
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		doItAndPrint("x <- seq(0, 1, 0.01)")
		doItAndPrint(paste("plot(x, x*", sens, "/(x*", sens, "+(1-x)*(1-", spec, ')), ylim=c(0,1), type="l", ylab="Predictive value", xlab="Pretest probability")', sep=""))
		doItAndPrint("par(new=T)")
		doItAndPrint(paste("plot(x, (1-x)*", spec, "/(x*(1-", sens, ")+(1-x)*", spec, '), ylim=c(0,1), type="l", lty=2, ylab="", xlab="")', sep=""))
		doItAndPrint('legend("bottom", c("Positive predictive value", "Negative predictive value"), lty=1:2, box.lty=0)')		
		doItAndPrint(paste("PPT <- ", preprob, "*", sens, "/(", preprob, "*", sens, "+(1-", preprob, ")*(1-", spec, "))", sep=""))
		doItAndPrint(paste("NPT <- (1-", preprob, ")*", spec, "/(", preprob, "*(1-", sens, ")+(1-", preprob, ")*", spec, ")", sep=""))
		doItAndPrint(paste("predictive.value <- data.frame(c(", preprob, ", ", sens, ", ", spec, ', " ", gettext(domain="R-RcmdrPlugin.EZR","Estimated"), round(PPT, 3), round(NPT,3)))', sep=""))
		doItAndPrint('colnames(predictive.value) <- gettext(domain="R-RcmdrPlugin.EZR","Assumptions")')
		doItAndPrint('rownames(predictive.value) <- gettext(domain="R-RcmdrPlugin.EZR",c("Pretest probability", "Sensitivity", "Specificity", " ", "  ", "Positive predictive value", "Negative predictive value"))')
		doItAndPrint("predictive.value")
		doItAndPrint("remove(predictive.value)")
		tkfocus(CommanderWindow())
		}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Pretest probability")), preprobEntry,
sticky="w")
	tkgrid.configure(preprobEntry, sticky="w")
	tkgrid(tklabel(top, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Sensitivity"), "(0-1)", sep="")), sensEntry, sticky="w")
	tkgrid.configure(sensEntry, sticky="w")
	tkgrid(tklabel(top, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Specificity"), "(0-1)", sep="")), specEntry, sticky="w")
	tkgrid.configure(specEntry, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedKappa <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Kappa statistics for agreement of two tests"))
	textFrame <- tkframe(top)
	variableFrame <- tkframe(top)
	variableFrame <- tkframe(top)
	pospos <- tclVar("")
	posposEntry <- ttkentry(variableFrame, width="10", textvariable=pospos)
	posneg <- tclVar("")
	posnegEntry <- ttkentry(variableFrame, width="10", textvariable=posneg)
	variable2Frame <- tkframe(top)
	negpos <- tclVar("")
	negposEntry <- ttkentry(variable2Frame, width="10", textvariable=negpos)
	negneg <- tclVar("")
	negnegEntry <- ttkentry(variable2Frame, width="10", textvariable=negneg)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Kappa statistics for agreement of two tests"), "#####", sep=""))
		pospos <- tclvalue(pospos)
		posneg <- tclvalue(posneg)
		negpos <- tclvalue(negpos)
		negneg <- tclvalue(negneg)
		closeDialog()
	if (length(pospos) == 0 || length(posneg) == 0 || length(negpos) == 0 || length(negneg) == 0){
			errorCondition(recall=StatMedKappa, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}		
#		library(epiR, quietly=TRUE)
		doItAndPrint(paste(".Table <- matrix(c(", pospos, ", ", posneg, ", ", negpos, ", ", negneg, "), 2, 2, byrow=TRUE)", sep=""))
		doItAndPrint('colnames(.Table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Test2 (+)", "Test2 (-)"))')
		doItAndPrint('rownames(.Table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Test1 (+)", "Test1 (-)"))')
		doItAndPrint(".Table")
		command <- "res <- epi.kappa(.Table, conf.level = 0.95)"
		doItAndPrint(command)
		doItAndPrint("remove(.Table)")
		doItAndPrint('colnames(res$kappa) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(res$kappa))')
#		doItAndPrint('colnames(res$mcnemar) <- gettext(domain="R-RcmdrPlugin.EZR", colnames(res$mcnemar))')
		doItAndPrint("res[1]")
		doItAndPrint("remove(res)")
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="epi.kappa")
	tkgrid(labelRcmdr(textFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number      Test2   (+)         (-)")), sticky="w")  
	tkgrid(textFrame, sticky="w")
	tkgrid(tklabel(variableFrame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Test1 (+)"), "         ", sep="")), posposEntry, posnegEntry, sticky="w")
	tkgrid(tklabel(variable2Frame, text=paste(gettext(domain="R-RcmdrPlugin.EZR","Test1  (-)"), "         ", sep="")), negposEntry, negnegEntry, sticky="w")
	tkgrid.configure(posposEntry, sticky="w")
	tkgrid.configure(posnegEntry, sticky="w")
	tkgrid.configure(negposEntry, sticky="w")
	tkgrid.configure(negnegEntry, sticky="w")
	tkgrid(variableFrame, sticky="nw")
	tkgrid(variable2Frame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

		    
StatMedReliability <- function(){
defaults <- list(x=NULL)
dialog.values <- getDialog("StatMedReliability", defaults)
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Cronbach's alpha coefficient for reliability"))
    xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables (pick three or more)"), initialSelection=varPosn(dialog.values$x, "numeric"))
    onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Cronbach's alpha coefficient for reliability"), "#####", sep=""))
        x <- getSelection(xBox)
putDialog("StatMedReliability", list(x=x))
        closeDialog()
        if (3 > length(x)) {
            errorCondition(recall=StatMedReliability, message=gettext(domain="R-RcmdrPlugin.EZR","Fewer than 3 variables selected."))
            return()
            }
        x <- paste('"', x, '"', sep="")
        doItAndPrint(paste("res <- reliability(cov(", ActiveDataSet(), "[,c(", paste(x, collapse=","),
            ')], use="complete.obs"))', sep=""))
		doItAndPrint("res$rel.matrix <- signif(res$rel.matrix, digits=4)")
		doItAndPrint('colnames(res$rel.matrix) <- gettext(domain="R-RcmdrPlugin.EZR", c("Alpha reliability", "Standardized alpha", "r(item, total)"))')
		doItAndPrint("res$rel.matrix <- cbind(rownames(res$rel.matrix), res$rel.matrix)")
		doItAndPrint("rownames(res$rel.matrix) <- NULL")
		doItAndPrint('colnames(res$rel.matrix)[1] <- gettext(domain="R-RcmdrPlugin.EZR","Deleted item")')
		doItAndPrint('cat("\n", gettext(domain="R-RcmdrPlugin.EZR","Alpha reliability"), "=", signif(res$alpha, digits=4), ", ", gettext(domain="R-RcmdrPlugin.EZR","Standardized alpha"), "=", signif(res$st.alpha, digits=4), "\n\n", gettext(domain="R-RcmdrPlugin.EZR","Reliability deleting each item in turn:"), "\n\n"); data.frame(res$rel.matrix)')
#		doItAndPrint("res$rel.matrix")
		doItAndPrint("remove(res)")
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="reliability", apply="StatMedReliability", reset="StatMedReliability")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
    }


StatMedOptMatch <- function(){
defaults <- list(group=NULL, strata=NULL, matchnumber="1", unmatch="FALSE", newDataSetName="Add _MP at the end of original name")
dialog.values <- getDialog("StatMedOptMatch", defaults)

	dataSet <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Extract matched controls"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (control=0, case=1) (pick one)"), listHeight=15, initialSelection=varPosn(dialog.values$group, "all"))
    strataBox <- variableListBox(variablesFrame, Variables(),selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Matching variables (pick at least one)"), listHeight=15, initialSelection=varPosn(dialog.values$strata, "all"))
#	newDataSetName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR","Add _MP at the end of original name"))

	dataSetNameFrame <- tkframe(top)
	dataSetName <- tclVar(gettext(domain="R-RcmdrPlugin.EZR",dialog.values$newDataSetName))
	dataSetNameField <- ttkentry(dataSetNameFrame, width="25", textvariable=dataSetName)
	
    optionsFrame <- tkframe(top)	
    matchnumberFrame <- tkframe(optionsFrame)
    matchnumberLevel <- tclVar(dialog.values$matchnumber)
    matchnumberField <- ttkentry(matchnumberFrame, width="6", textvariable=matchnumberLevel)
	
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Extract matched controls"), "#####", sep=""))
	    group <- getSelection(groupBox)
        strata <- getSelection(strataBox)
	    matchnumber <- tclvalue(matchnumberLevel)
        unmatch <- as.character(tclvalue(unmatchVariable))
		newName <- trim.blanks(tclvalue(dataSetName))
		if (newName == gettext(domain="R-RcmdrPlugin.EZR","Add _MP at the end of original name")) newName <- paste(ActiveDataSet(), "_MP", sep="")
putDialog("StatMedOptMatch", list(group=group, strata=strata, matchnumber=matchnumber, unmatch=unmatch, newDataSetName=newName))
        if (length(group) == 0) {
            errorCondition(recall=StatMedOptMatch, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups variable."))
            return()
            }
        if (length(strata) == 0) {
            errorCondition(recall=StatMedOptMatch, message=gettext(domain="R-RcmdrPlugin.EZR","Pick at least one matching variable"))
            return()
            }
		if (!is.valid.name(newName)){
			errorCondition(recall=StatMedOptMatch,
				message=paste('"', newName, '" ', gettext(domain="R-RcmdrPlugin.EZR","is not a valid name."), sep=""))
			return()
		}
		if (is.element(newName, listDataSets())) {
			if ("no" == tclvalue(checkReplace(newName, type=gettext(domain="R-RcmdrPlugin.EZR","Data set")))){
				closeDialog()
				StatMedOptMatch()
				return()
			}
		}
		closeDialog()
		.activeDataSet <- ActiveDataSet()
		doItAndPrint("library(optmatch, quietly=TRUE)")
		nacheck.command <- paste("TempDataSet <- ", .activeDataSet, "[complete.cases(", .activeDataSet, "$", group, ", ", .activeDataSet, "$", strata[1], sep="")
		mdist.command <- paste("match.distance <- mdist(", group, " ~ ", strata[1], sep="")
		if (length(strata) >1 ){
			for (i in 2:length(strata)){
				nacheck.command <- paste(nacheck.command, ", ", .activeDataSet, "$", strata[i], sep="")
				mdist.command <- paste(mdist.command, " + ", strata[i], sep="")
			}
		}
		nacheck.command <- paste(nacheck.command, "),]", sep="")
		mdist.command <- paste(mdist.command, ", data=TempDataSet)", sep="")
		doItAndPrint(nacheck.command)
		doItAndPrint(mdist.command)
		match.command <- paste("pairmatch.results <- pairmatch(match.distance, control=", matchnumber, ", remove.unmatchables=", unmatch, ", data=TempDataSet)", sep="")
		logger(match.command)
		result <- justDoIt(match.command)
		if (class(result)[1] ==  "try-error"){
				errorCondition(recall=StatMedOptMatch, message=gettext(domain="R-RcmdrPlugin.EZR","Matching failed"))
				return()
		}
		doItAndPrint("TempDataSet$pairmatch <- as.integer(substring(pairmatch.results, 3))")
		command <- paste(newName, " <- TempDataSet[!is.na(TempDataSet$pairmatch),]", sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(newName)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="pairmatch", apply="StatMedOptMatch", reset="StatMedOptMatch")
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
	tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text="    "), getFrame(strataBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(matchnumberFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Number of controls matched to one case"), fg="blue"), sticky="w")
    tkgrid(matchnumberField, sticky="w")
    radioButtons(optionsFrame, name="unmatch", buttons=c("No", "Yes"), values=c("FALSE", "TRUE"),
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Not remove", "Remove")), title=gettext(domain="R-RcmdrPlugin.EZR","Remove unmatched cases?"))
 	tkgrid(matchnumberFrame, labelRcmdr(optionsFrame, text="    "), unmatchFrame, sticky="nw")		
    tkgrid(optionsFrame, sticky="nw")
	tkgrid(labelRcmdr(dataSetNameFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Name for new data set")), sticky="w")
	tkgrid(dataSetNameField, sticky="w")
	tkgrid(dataSetNameFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}


StatMedMH <- function(){
defaults <- list(group=NULL, var=NULL, strata=NULL, continuity="TRUE")
dialog.values <- getDialog("StatMedMH", defaults)
	dataSet <- activeDataSet()
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Mantel-Haenzel test for matched proportions"))
    variablesFrame <- tkframe(top)
    groupBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable (control=0, case=1) (pick one)"), listHeight=8, initialSelection=varPosn(dialog.values$group, "all"))
    varBox <- variableListBox(variablesFrame, Variables(),selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Binary response variable (pick at least one)"), listHeight=8, initialSelection=varPosn(dialog.values$var, "all"))
    variables2Frame <- tkframe(top)
    strataBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Stratifying variable for matching (pairmatch)"), listHeight=8, initialSelection=varPosn(dialog.values$strata, "all"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Mantel-Haenzel test for matched proportions"), "#####", sep=""))
	    group <- getSelection(groupBox)
        var <- getSelection(varBox)
        strata <- getSelection(strataBox)
		continuity <- tclvalue(continuityVariable)
putDialog("StatMedMH", list(group=group, var=var, strata=strata, continuity=continuity))
        if (length(group) == 0) {
            errorCondition(recall=StatMedMH, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a groups variable."))
            return()
            }
        if (length(var) == 0) {
            errorCondition(recall=StatMedMH, message=gettext(domain="R-RcmdrPlugin.EZR","Pick at least one binary response variable."))
            return()
            }
        if (length(strata) == 0) {
            errorCondition(recall=StatMedMH, message=gettext(domain="R-RcmdrPlugin.EZR",""))
            return()
            }
		closeDialog()
		.activeDataSet <- ActiveDataSet()
		nvar = length(var)
		doItAndPrint("MH.summary.table <- NULL")
		for (i in 1:nvar) {
        	if (var[i] == group) {
            	errorCondition(recall=StatMedMH, message=gettext(domain="R-RcmdrPlugin.EZR","Row and column variables are the same."))
 	           return()
            }
			command <- paste("xtabs(~", var[i], "+", group, ", data=", .activeDataSet, ")", sep="")
#       	 	logger(paste(".Table <- ", command, sep=""))
#        	assign(".Table", justDoIt(command), envir=.GlobalEnv)
       	 	doItAndPrint(paste(".Table <- ", command, sep=""))
        	doItAndPrint(".Table")
			command <- paste("(res <- mantelhaen.test(", .activeDataSet, "$", group, ", ", .activeDataSet, "$", var[i], ", ", .activeDataSet, "$", strata, ", correct=", continuity, "))", sep="")			
			doItAndPrint(command)
			doItAndPrint("MH.summary.table <- rbind(MH.summary.table, summary.table.MH(table=.Table, res=res))")		
			doItAndPrint("remove(res)")	
		}        
		doItAndPrint("MH.summary.table")
#		doItAndPrint("remove(MH.summary.table)")				
		logger("remove(.Table)")
        remove(.Table, envir=.GlobalEnv)
        tkfocus(CommanderWindow())
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="mantelhaen.test", apply="StatMedMH", reset="StatMedMH")
	tkgrid(labelRcmdr(variablesFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
	tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text="    "), getFrame(varBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(getFrame(strataBox), labelRcmdr(variables2Frame, text="    "), sticky="nw")
    tkgrid(variables2Frame, sticky="nw")	
	analysisFrame <- tkframe(top)
    radioButtons(analysisFrame, name="continuity",
        buttons=c("yes", "no"),
        values=c("TRUE", "FALSE"), initialValue=dialog.values$continuity,
        labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes", "No")), title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction of chi-square test"))
    tkgrid(continuityFrame, sticky="w")
	tkgrid(analysisFrame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}


StatMedCLogistic <- function(){
defaults <- list(lhs = "", rhs = "", actmodelVariable = 0, strata = NULL)
dialog.values <- getDialog("StatMedCLogistic", defaults)
currentFields$lhs <- dialog.values$lhs	
currentFields$rhs <- dialog.values$rhs
currentFields$subset <- dialog.values$subset	

#    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Generalized Linear Model"))
    initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Conditional logistic regression for matched-pair analysis"))
    .activeModel <- ActiveModel()
    currentModel <- if (!is.null(.activeModel))
        class(get(.activeModel, envir=.GlobalEnv))[1] == "glm"
#        eval(parse(text=paste("class(", .activeModel, ")[1] == 'glm'", sep="")),
#            envir=.GlobalEnv)
        else FALSE
    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv), glm=TRUE)
#        currentFields <- formulaFields(eval(parse(text=.activeModel),
#            envir=.GlobalEnv), glm=TRUE)
        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
        }
		
	currentModel <- TRUE
	
    StatMedModelFormula()
    variables2Frame <- tkframe(top)
    strataBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Stratifying variable for matching (pairmatch)"), listHeight=8, initialSelection=varPosn(dialog.values$strata, "all"))
    UpdateModelNumber()
    modelName <- tclVar(paste("GLM.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
	optionsFrame <- tkframe(top)
	
	checkBoxes(frame="checkboxFrame", boxes=c("actmodel"), initialValues=c(dialog.values$actmodelVariable),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Keep results as active model for further analyses")))	
	
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
#	stepwise1Variable <- tclVar("0")
#	stepwise1CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise1Variable)
#	stepwise2Variable <- tclVar("0")
#	stepwise2CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise2Variable)
    onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Conditional logistic regression for matched-pair analysis"), "#####", sep=""))
        strata <- getSelection(strataBox)
        modelValue <- trim.blanks(tclvalue(modelName))
        formula <- paste(tclvalue(lhsVariable), " ~ ", tclvalue(rhsVariable), " + strata(", strata, ")", sep="")
		actmodel <- tclvalue(actmodelVariable)
#		stepwise1 <- tclvalue(stepwise1Variable)
#       stepwise2 <- tclvalue(stepwise2Variable)
putDialog("StatMedCLogistic", list(lhs = tclvalue(lhsVariable), rhs = tclvalue(rhsVariable), actmodelVariable = actmodel, strata = strata))
        closeDialog()		
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedCLogistic, model=TRUE, message=gettext(domain="R-RcmdrPlugin.EZR","Left-hand side of model empty."))
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedCLogistic, model=TRUE, message=gettext(domain="R-RcmdrPlugin.EZR","Right-hand side of model empty."))
            return()
            }
        if (length(strata) == 0) {
            errorCondition(recall=StatMedCLogistic, message=gettext(domain="R-RcmdrPlugin.EZR","Pick one stratifying variable for matching."))
            return()
            }
        if (!is.valid.name(modelValue)){
            errorCondition(recall=StatMedCLogistic, model=TRUE, message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue))
            return()
            }
        if (is.element(modelValue, listGeneralizedLinearModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
                UpdateModelNumber(-1)
                closeDialog()
                StatMedCLogistic()
                return()
                }
           }		
		Library("survival")
        command <- paste("clogit(", formula, ", data=", ActiveDataSet(), ")", sep="")
#        logger(paste(modelValue, " <- ", command, sep=""))
#        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(paste(modelValue, " <- ", command, sep=""))
        doItAndPrint(paste("(res <- summary(", modelValue, "))", sep=""))
		
	if(eval(parse(text="length(res$coefficients[,1])"))==1){
			doItAndPrint(paste("odds <- signif(c(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)", sep=""))
			doItAndPrint("odds <- t(odds)")
			doItAndPrint("rownames(odds) <- rownames(res$coefficients)")
	} else {
			doItAndPrint(paste("odds <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)", sep=""))
	}
		
		doItAndPrint("odds <- data.frame(odds)")
		doItAndPrint("odds <- signif(odds, digits=3)")
		doItAndPrint('names(odds) <- gettext(domain="R-RcmdrPlugin.EZR",c("odds ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
		doItAndPrint("odds")
#		if (stepwise1 == 1 | stepwise2 == 1){
#			x <- strsplit(tclvalue(rhsVariable), split="\\+")
#			command <- paste("TempDF <- with(", ActiveDataSet(), ", ", ActiveDataSet(), "[complete.cases(", paste(x[[1]], collapse=","), "),])", sep="")
#			doItAndPrint(command)
#			command <- paste("clogit(", formula, ", data=TempDF)", sep="")
#			doItAndPrint(paste(modelValue, " <- ", command, sep=""))
#			}
#		if (stepwise1 == 1){
#			doItAndPrint("odds <- data.frame(exp( summary(res)$coef[,1:2] %*% rbind(c(1,1,1), 1.96*c(0,-1,1))))")
#		doItAndPrint(paste("odds <- cbind(odds, summary(res)$coefficients[,4])", sep=""))
#			doItAndPrint("odds <- signif(odds, digits=3)")
#			doItAndPrint('names(odds) <- c("odds ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
#			doItAndPrint("summary(res)")
#			doItAndPrint("odds")
#			doItAndPrint("remove(res)")			
#		}
#		if (stepwise2 == 1){
#			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="BIC")', sep=""))
#			doItAndPrint("odds <- data.frame(exp( summary(res)$coef[,1:2] %*% rbind(c(1,1,1), 1.96*c(0,-1,1))))")
#			doItAndPrint(paste("odds <- cbind(odds, summary(res)$coefficients[,4])", sep=""))
#			doItAndPrint("odds <- signif(odds, digits=3)")
#			doItAndPrint('names(odds) <- c("odds ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
#			doItAndPrint("summary(res)")
#			doItAndPrint("odds")
#			doItAndPrint("remove(res)")
#			}		
#		doItAndPrint("remove(odds)")
		doItAndPrint("remove(res)")
		if (actmodel==1) activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="clogit", apply="StatMedCLogistic", reset="StatMedCLogistic")
    helpButton <- buttonRcmdr(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(labelRcmdr(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(getFrame(xBox), sticky="w")
    tkgrid(outerOperatorsFrame, sticky="w")
    tkgrid(formulaFrame, sticky="w")
    tkgrid(getFrame(strataBox), labelRcmdr(variables2Frame, text="    "), sticky="nw")
    tkgrid(variables2Frame, sticky="nw")	
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on AIC")), stepwise1CheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on BIC")), stepwise2CheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on AIC/BIC not performed when missing data included.")), sticky="w")  
	tkgrid(optionsFrame, sticky="w", columnspan=2)
	
 tkgrid(checkboxFrame, sticky="w")
	
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
}

		
StatMedStCox  <- function(){
    xx <- getRcmdr("modelClasses")
    bolCoxphExists = FALSE
    for(ii in 1:length(xx)){if (xx[ii] == "coxph") bolCoxphExists = TRUE}
    if (bolCoxphExists == FALSE) putRcmdr("modelClasses", c(getRcmdr("modelClasses"), "coxph"))

defaults <- list(SurvivalTimeVariable = "", StatusVariable = "", rhs = "", waldVariable = 0,  prophazVariable = 0, basecurveVariable = 0, actmodelVariable = 0, stepwise1Variable = 0, stepwise2Variable = 0, stepwise3Variable = 0, strata = NULL)
dialog.values <- getDialog("StatMedStCox", defaults)
currentFields$SurvivalTimeVariable <- dialog.values$SurvivalTimeVariable	
currentFields$StatusVariable <- dialog.values$StatusVariable
currentFields$rhs <- dialog.values$rhs

  initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Stratified Cox proportional hazard regression for matched-pair analysis"))
  .activeModel <- ActiveModel()
  currentModel <- if (!is.null(.activeModel)) 
        class(get(.activeModel, envir=.GlobalEnv))[1] == "coxph"
#    eval(parse(text=paste("class(", .activeModel, ")[1] == 'coxph'", sep="")), 
#         envir=.GlobalEnv) 
    else FALSE
	
	currentModel <- TRUE
	#  if(currentModel){
#    currentFields <- formulaFields(eval(parse(text=.activeModel), 
#     envir=.GlobalEnv))
#    if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
#  }
    variables2Frame <- tkframe(top)
    strataBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Stratifying variable for matching (pairmatch)"), listHeight=8, initialSelection=varPosn(dialog.values$strata, "all"))
  UpdateModelNumber()
  modelName <- tclVar(paste("CoxModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="30", textvariable=modelName)
  	optionsFrame <- tkframe(top)
	
	checkBoxes(frame="checkboxFrame", boxes=c("wald", "prophaz", "actmodel", "stepwise1", "stepwise2"), initialValues=c(dialog.values$waldVariable, dialog.values$prophazVariable, dialog.values$actmodelVariable, dialog.values$stepwise1Variabl, dialog.values$stepwise2Variabl),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Wald test for overall p-value for factors with >2 levels", "Test proportional hazards assumption", "Keep results as active model for further analyses", "Stepwise selection based on AIC", "Stepwise selection based on BIC")))	
	
#	waldVariable <- tclVar("0")
#	waldCheckBox <- tkcheckbutton(optionsFrame, variable=waldVariable)
#	prophazVariable <- tclVar("0")
#	prophazCheckBox <- tkcheckbutton(optionsFrame, variable=prophazVariable)
#	actmodelVariable <- tclVar("0")
#	actmodelCheckBox <- tkcheckbutton(optionsFrame, variable=actmodelVariable)
#	stepwise1Variable <- tclVar("0")
#	stepwise1CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise1Variable)
#	stepwise2Variable <- tclVar("0")
#	stepwise2CheckBox <- tkcheckbutton(optionsFrame, variable=stepwise2Variable)
  onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Stratified Cox proportional hazard regression for matched-pair analysis"), "#####", sep=""))
#    XXX <- getSelection(timeBox)
    modelValue <- trim.blanks(tclvalue(modelName))
    strata <- getSelection(strataBox)
		prophaz <- tclvalue(prophazVariable)
		wald <- tclvalue(waldVariable)
		actmodel <- tclvalue(actmodelVariable)
		stepwise1 <- tclvalue(stepwise1Variable)
		stepwise2 <- tclvalue(stepwise2Variable)
#    library(survival, quietly=TRUE)
#    formula <- paste("Surv(", XXX, ", ", tclvalue(lhsVariable), ") ~ ", tclvalue(rhsVariable), sep="")
putDialog("StatMedStCox", list(SurvivalTimeVariable = tclvalue(SurvivalTimeVariable), StatusVariable = tclvalue(StatusVariable), rhs = tclvalue(rhsVariable), waldVariable = wald,  prophazVariable = prophaz, actmodelVariable = actmodel, stepwise1Variable = stepwise1, stepwise2Variable = stepwise2, strata = strata))
    closeDialog()
#    check.empty <- gsub(" ", "", tclvalue(lhsVariable))
#    if ("" == check.empty) {
#      errorCondition(recall=StatMedCoxRegression,
#        message=gettext(domain="R-RcmdrPlugin.EZR","Left-hand side of model empty."), model=TRUE) 
#      return()
#    }
        check.empty <- gsub(" ", "", tclvalue(SurvivalTimeVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedStCox, message=gettext(domain="R-RcmdrPlugin.EZR","Survival time variable of model empty."), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(StatusVariable))
        if ("" == check.empty) {
            errorCondition(recall=StatMedStCox, message=gettext(domain="R-RcmdrPlugin.EZR","Status variable of model empty."), model=TRUE)
            return()
            }

    check.empty <- gsub(" ", "", tclvalue(rhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=StatMedStCox,
        message=gettext(domain="R-RcmdrPlugin.EZR","Right-hand side of model empty."), model=TRUE)
      return()
    }
    if (is.element(modelValue, listCoxModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettext(domain="R-RcmdrPlugin.EZR","Model")))){
        UpdateModelNumber(-1)
        StatMedCoxRegression()
        return()
      }
    }
    if (!is.valid.name(modelValue)){
      errorCondition(recall=StatMedStCox, 
        message=sprintf(gettext(domain="R-RcmdrPlugin.EZR",'"%s" is not a valid name.'), modelValue), model=TRUE)
      return()
    }
        if (length(strata) == 0) {
            errorCondition(recall=StatMedStCox, message=gettext(domain="R-RcmdrPlugin.EZR","Pick one stratifying variable for matching."))
            return()
            }

    Library("survival")
		Library("aod")
     formula <- paste("Surv(", tclvalue(SurvivalTimeVariable), ", ", tclvalue(StatusVariable), "==1)~ ", tclvalue(rhsVariable), " + strata(", strata, ")", sep="")

    command <- paste("coxph(", formula,
      ", data=", ActiveDataSet(), ', method="breslow")', sep="")
#    logger(paste(modelValue, " <- ", command, sep=""))
#    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
    doItAndPrint(paste("summary(", modelValue, ")", sep=""))
	doItAndPrint(paste("res <- ", command, sep=""))
	doItAndPrint("res <- summary(res)")
#	if(eval(parse(text="length(res$coefficients[,1])"))==1){
#		doItAndPrint("cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)")
#		doItAndPrint(paste('rownames(cox.table) <- "', tclvalue(rhsVariable), '"', sep=""))
#		doItAndPrint('colnames(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')	
#	} else {
#		doItAndPrint("cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)")
#		doItAndPrint("cox.table <- data.frame(cox.table)")
#		doItAndPrint('names(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
#	}

	if(eval(parse(text="length(res$coefficients[,1])"))==1){
		doItAndPrint("cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), p.value=res$coefficients[,5]), digits=4)")
		doItAndPrint("rownames(cox.table) <- rownames(res$coefficients)")
		doItAndPrint('colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	} else {
		doItAndPrint("cox.table <- signif(cbind(res$conf.int[,c(1,3,4)], res$coefficients[,5]), digits=4)")
		doItAndPrint("cox.table <- data.frame(cox.table)")
		doItAndPrint('colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value"))')
	}

#	doItAndPrint("cox.table <- signif(cox.table, digits=3)")
	doItAndPrint("cox.table")
	if (wald==1) doItAndPrint(paste("waldtest(", modelValue, ")", sep=""))
	if (prophaz == 1){
			doItAndPrint(paste("print(cox.zph(", modelValue, "))", sep=""))
	}

	if (stepwise1 == 1 | stepwise2 == 1){
		x <- strsplit(tclvalue(rhsVariable), split="\\+")
		command <- paste("TempDF <- with(", ActiveDataSet(), ", ", ActiveDataSet(), "[complete.cases(", paste(x[[1]], collapse=","), "),])", sep="")
		doItAndPrint(command)
		command <- paste("coxph(", formula, ', data=TempDF, method="breslow")', sep="")
		doItAndPrint(paste(modelValue, " <- ", command, sep=""))
		}
	if (stepwise1 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="AIC")', sep=""))
			doItAndPrint("summary(res)")
			doItAndPrint("res2 <- summary(res)")
			if(eval(parse(text="length(res2$coefficients[,1])"))==1){
				doItAndPrint("cox.table <- signif(cbind(t(res2$conf.int[,c(1,3,4)]), p.value=res2$coefficients[,5]), digits=4)")
				doItAndPrint("rownames(cox.table) <- rownames(res2$coefficients)")
				doItAndPrint('colnames(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')	
				doItAndPrint("cox.table")
			} else if(eval(parse(text="length(res2$coefficients[,1])"))>1){
				doItAndPrint("cox.table <- signif(cbind(res2$conf.int[,c(1,3,4)], res2$coefficients[,5]), digits=4)")
				doItAndPrint("cox.table <- data.frame(cox.table)")
				doItAndPrint('names(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
				doItAndPrint("cox.table")
			}
			if (wald==1) doItAndPrint("waldtest(res)")
	}
	if (stepwise2 == 1){
			doItAndPrint(paste("res <- stepwise(", modelValue, ', direction="backward/forward", criterion="BIC")', sep=""))
			doItAndPrint("summary(res)")
			doItAndPrint("res2 <- summary(res)")
			if(eval(parse(text="length(res2$coefficients[,1])"))==1){
				doItAndPrint("cox.table <- signif(cbind(t(res2$conf.int[,c(1,3,4)]), p.value=res2$coefficients[,5]), digits=4)")
				doItAndPrint("rownames(cox.table) <- rownames(res2$coefficients)")
				doItAndPrint('colnames(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')	
				doItAndPrint("cox.table")
			} else if(eval(parse(text="length(res2$coefficients[,1])"))>1){
				doItAndPrint("cox.table <- signif(cbind(res2$conf.int[,c(1,3,4)], res2$coefficients[,5]), digits=4)")
				doItAndPrint("cox.table <- data.frame(cox.table)")
				doItAndPrint('names(cox.table) <- c("Hazard ratio", "Lower 95%CI", "Upper 95%CI", "p.value")')
				doItAndPrint("cox.table")
			}
			if (wald==1) doItAndPrint("waldtest(res)")
	}

	doItAndPrint("remove(res)")
#	doItAndPrint("remove(cox.table)")
	if (actmodel==1) activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="coxph", model=TRUE, apply="StatMedStCox", reset="StatMedStCox")
  tkgrid(tklabel(modelFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
#  StatMedModelFormula()
  modelFormulaCox()
 tkgrid(getFrame(xBox), sticky="w")
#  tkgrid(getFrame(xBox), getFrame(timeBox), sticky="w")
  tkgrid(outerOperatorsFrame, sticky="w")
  tkgrid(formulaFrame, sticky="w")
    tkgrid(getFrame(strataBox), labelRcmdr(variables2Frame, text="    "), sticky="nw")
    tkgrid(variables2Frame, sticky="nw")	

 tkgrid(checkboxFrame, sticky="w")

#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Wald test for overall p-value for factors with >2 levels")), waldCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Test proportional hazards assumption")), prophazCheckBox, sticky="w")
#	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Keep results as active model for further analyses")), actmodelCheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on AIC")), stepwise1CheckBox, sticky="w")
#  	tkgrid(labelRcmdr(optionsFrame, text=gettext(domain="R-RcmdrPlugin.EZR","Stepwise selection based on BIC")), stepwise2CheckBox, sticky="w")
	tkgrid(optionsFrame, sticky="w", columnspan=2)
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
}


StatMedSampleProportionsSingle <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison with specified proportion"))
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1),
labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	radioButtons(name="continuity", buttons=c("Yes", "No"), values=c(1, 0), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes (or exact test)", "No correction")),title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction of chi-square test"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison with specified proportion"), "#####", sep=""))
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		method <- tclvalue(methodVariable)
		continuity <- tclvalue(continuityVariable)
		closeDialog()
		if (length(group1) == 0 || length(group2) == 0){
			errorCondition(recall=StatMedSampleProportionsSingle,
message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0 ){
			errorCondition(recall=StatMedSampleProportionsSingle,
message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleProportionSingleArm(", group1, ", ", group2,
", ", alpha, ", ", power, ", ", method, ", ", continuity, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion (control)")),
group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion (test)")),
group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry,
sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")),
powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(continuityFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

	
StatMedPowerProportionsSingle <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison with specified proportion"))
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	sample <- tclVar("")
	sampleEntry <- ttkentry(top, width="20", textvariable=sample)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	radioButtons(name="continuity", buttons=c("Yes", "No"), values=c(1, 0), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes (or exact test)", "No correction")),title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction of chi-square test"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison with specified proportion"), "#####", sep=""))
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		alpha <- tclvalue(alpha)
		sample <- tclvalue(sample)
		method <- tclvalue(methodVariable)
		continuity <- tclvalue(continuityVariable)
		closeDialog()
		if (length(group1) == 0 || length(group2) == 0){
			errorCondition(recall=StatMedPowerProportionsSingle, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(sample) == 0 ){
			errorCondition(recall=StatMedPowerProportionsSingle, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("PowerProportionSingleArm(", group1, ", ", group2, ", ", alpha, ", ", sample, ", ", method, ", ", continuity, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion (control)")), group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion (test)")), group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size")), sampleEntry, sticky="w")
	tkgrid.configure(sampleEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(continuityFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

	
StatMedSampleProportionsCI <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size from proportion and confidence interval"))
	p1 <- tclVar("")
	p1Entry <- ttkentry(top, width="20", textvariable=p1)
	delta <- tclVar("")
	deltaEntry <- ttkentry(top, width="20", textvariable=delta)
	ci <- tclVar("95")
	ciEntry <- ttkentry(top, width="20", textvariable=ci)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size from proportion and confidence interval"), "#####", sep=""))
		p1 <- tclvalue(p1)
		delta <- tclvalue(delta)
		ci <- tclvalue(ci)
		closeDialog()
		if (length(p1) == 0 || length(delta) == 0 || length(ci) == 0){
			errorCondition(recall=StatMedSampleProportionsCI, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleProportionCI(", p1, ", ", delta, ", ", ci, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion")), p1Entry, sticky="w")
	tkgrid.configure(p1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval width")), deltaEntry, sticky="w")
	tkgrid.configure(deltaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence level")), ciEntry, sticky="w")
	tkgrid.configure(ciEntry, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

	
StatMedSampleMeansCI <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size from standard deviation and confidence interval"))
	sd <- tclVar("")
	sdEntry <- ttkentry(top, width="20", textvariable=sd)
	delta <- tclVar("")
	deltaEntry <- ttkentry(top, width="20", textvariable=delta)
	ci <- tclVar("95")
	ciEntry <- ttkentry(top, width="20", textvariable=ci)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size from standard deviation and confidence interval"), "#####", sep=""))
		sd <- tclvalue(sd)
		delta <- tclvalue(delta)
		ci <- tclvalue(ci)
		closeDialog()
		if (length(sd) == 0 || length(delta) == 0 || length(ci) == 0){
			errorCondition(recall=StatMedSampleMeansCI, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleMeanCI(", sd, ", ", delta, ", ", ci, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation (expected)")), sdEntry, sticky="w")
	tkgrid.configure(sdEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence interval width")), deltaEntry, sticky="w")
	tkgrid.configure(deltaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Confidence level")), ciEntry, sticky="w")
	tkgrid.configure(ciEntry, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}

	
StatMedSamplePhaseII <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size from control and desired response rates"))
	p1 <- tclVar("")
	p1Entry <- ttkentry(top, width="20", textvariable=p1)
	p2 <- tclVar("")
	p2Entry <- ttkentry(top, width="20", textvariable=p2)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	checkBoxes(frame="twostage", boxes=c("twostage"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Calculate two-stage model")))
#	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size from control and desired response rates"), "#####", sep=""))
		p1 <- tclvalue(p1)
		p2 <- tclvalue(p2)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		twostage <- tclvalue(twostageVariable)
		closeDialog()
		if (length(p1) == 0 || length(p2) == 0){
			errorCondition(recall=StatMedSamplePhaseII, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (as.numeric(p1) >= as.numeric(p2)){
			errorCondition(recall=StatMedSamplePhaseII, message=gettext(domain="R-RcmdrPlugin.EZR","Desirable response rate must be higher than unacceptable response rate."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0 ){
			errorCondition(recall=StatMedSamplePhaseII, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		doItAndPrint("library(clinfun, quietly=TRUE)")
		command <- paste("ph2single(", p1, ", ", p2, ", ", alpha, ", (1-", power, "), nsoln=1)", sep="")		
		doItAndPrint(command)
		if (twostage==1){
			command <- paste("ph2simon(", p1, ", ", p2, ", ", alpha, ", (1-", power, "), nmax=200)", sep="") #Two-stage		
			doItAndPrint(command)			
		}
		logger(gettext(domain="R-RcmdrPlugin.EZR","# r: if the number of response is equal to or fewer than r, the treatment is rejected."))
		if (twostage==1){
		logger(gettext(domain="R-RcmdrPlugin.EZR","# r1, n1: numbers in the first stage, r, n: total numbers in the study."))
		}
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Unacceptable response rate")), p1Entry, sticky="w")
	tkgrid.configure(p1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Desirable response rate")), p2Entry, sticky="w")
	tkgrid.configure(p2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")), powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
    tkgrid(twostage, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedSampleMeans <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two means"))
	difference <- tclVar("")
	differenceEntry <- ttkentry(top, width="20", textvariable=difference)
	stddevi <- tclVar("")
	stddeviEntry <- ttkentry(top, width="20", textvariable=stddevi)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	ratio <- tclVar("1")
	ratioEntry <- ttkentry(top, width="20", textvariable=ratio)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1),
labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two means"), "#####", sep=""))
		difference <- tclvalue(difference)
		stddevi <- tclvalue(stddevi)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		ratio <- tclvalue(ratio)
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(difference) == 0 || length(stddevi) == 0){
			errorCondition(recall=StatMedSampleMeans, message=gettext(domain="R-RcmdrPlugin.EZR","You
must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0 || length(ratio) == 0){
			errorCondition(recall=StatMedSampleMeans, message=gettext(domain="R-RcmdrPlugin.EZR","You
must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleMean(", difference, ", ", stddevi, ", ",
alpha, ", ", power, ", ", method, ", ", ratio, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Difference in means")), differenceEntry,
sticky="w")
	tkgrid.configure(differenceEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation in each group")),
stddeviEntry, sticky="w")
	tkgrid.configure(stddeviEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry,
sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")),
powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size ratio (1:X)")),
ratioEntry, sticky="w")
	tkgrid.configure(ratioEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedPowerMeans <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two means"))
	difference <- tclVar("")
	differenceEntry <- ttkentry(top, width="20", textvariable=difference)
	stddevi <- tclVar("")
	stddeviEntry <- ttkentry(top, width="20", textvariable=stddevi)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	sample1 <- tclVar("")
	sample1Entry <- ttkentry(top, width="20", textvariable=sample1)
	sample2 <- tclVar("")
	sample2Entry <- ttkentry(top, width="20", textvariable=sample2)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two means"), "#####", sep=""))
		difference <- tclvalue(difference)
		stddevi <- tclvalue(stddevi)
		alpha <- tclvalue(alpha)
		sample1 <- as.numeric(tclvalue(sample1))
		sample2 <- as.numeric(tclvalue(sample2))
		if (sample1 >= sample2){
			sample <- sample2
			ratio <- sample1/sample2
		} else {
			sample <- sample1
			ratio <- sample2/sample1
		}
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(difference) == 0 || length(stddevi) == 0){
			errorCondition(recall=StatMedPowerMeans, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(sample) == 0 || length(ratio) == 0){
			errorCondition(recall=StatMedPowerMeans, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("PowerMean(", difference, ", ", stddevi, ", ", alpha, ", ", sample, ", ", method, ", ", ratio, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Difference in means")), differenceEntry, sticky="w")
	tkgrid.configure(differenceEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation in each group")), stddeviEntry, sticky="w")
	tkgrid.configure(stddeviEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size of group 1")), sample1Entry, sticky="w")
	tkgrid.configure(sample1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size of group 2")), sample2Entry, sticky="w")
	tkgrid.configure(sample2Entry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedSampleMeansNonInf <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for non-inferiority trial of two means"))
	difference <- tclVar("")
	differenceEntry <- ttkentry(top, width="20", textvariable=difference)
	delta <- tclVar("")
	deltaEntry <- ttkentry(top, width="20", textvariable=delta)
	stddevi <- tclVar("")
	stddeviEntry <- ttkentry(top, width="20", textvariable=stddevi)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
#	ratio <- tclVar("1")
#	ratioEntry <- ttkentry(top, width="20", textvariable=ratio)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), initialValue=1, 
labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
		logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for non-inferiority trial of two means"), "#####", sep=""))
		difference <- tclvalue(difference)
		delta <- tclvalue(delta)
		stddevi <- tclvalue(stddevi)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
#		ratio <- tclvalue(ratio)
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(difference) == 0 || length(delta) == 0 || length(stddevi) == 0){
			errorCondition(recall=StatMedSampleMeans, message=gettext(domain="R-RcmdrPlugin.EZR","You
must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0){
			errorCondition(recall=StatMedSampleMeans, message=gettext(domain="R-RcmdrPlugin.EZR","You
must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleMeanNonInf(", difference, ", ", delta, ", ", stddevi, ", ",
alpha, ", ", power, ", ", method, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Difference in means (test - control)")), differenceEntry, sticky="w")
	tkgrid.configure(differenceEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Meaningful difference in mean")), deltaEntry, sticky="w")
	tkgrid.configure(deltaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation in each group")), stddeviEntry, sticky="w")
	tkgrid.configure(stddeviEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")), powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
#	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size ratio (1:X)")), ratioEntry, sticky="w")
#	tkgrid.configure(ratioEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedSampleMeansPaired <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two paired means"))
	difference <- tclVar("")
	differenceEntry <- ttkentry(top, width="20", textvariable=difference)
	stddevi <- tclVar("")
	stddeviEntry <- ttkentry(top, width="20", textvariable=stddevi)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1),
labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two paired means"), "#####", sep=""))
		difference <- tclvalue(difference)
		stddevi <- tclvalue(stddevi)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(difference) == 0 || length(stddevi) == 0){
			errorCondition(recall=StatMedSampleMeansPaired, message=gettext(domain="R-RcmdrPlugin.EZR","You
must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0){
			errorCondition(recall=StatMedSampleMeansPaired, message=gettext(domain="R-RcmdrPlugin.EZR","You
must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleMeanPaired(", difference, ", ", stddevi, ", ",
alpha, ", ", power, ", ", method, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Difference in means")), differenceEntry,
sticky="w")
	tkgrid.configure(differenceEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation in each group")),
stddeviEntry, sticky="w")
	tkgrid.configure(stddeviEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry,
sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")),
powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedPowerMeansPaired <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two paired means"))
	difference <- tclVar("")
	differenceEntry <- ttkentry(top, width="20", textvariable=difference)
	stddevi <- tclVar("")
	stddeviEntry <- ttkentry(top, width="20", textvariable=stddevi)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	sample <- tclVar("")
	sampleEntry <- ttkentry(top, width="20", textvariable=sample)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two paired means"), "#####", sep=""))
		difference <- tclvalue(difference)
		stddevi <- tclvalue(stddevi)
		alpha <- tclvalue(alpha)
		sample <- as.numeric(tclvalue(sample))
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(difference) == 0 || length(stddevi) == 0){
			errorCondition(recall=StatMedPowerMeansPaired, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(sample) == 0){
			errorCondition(recall=StatMedPowerMeansPaired, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("PowerMeanPaired(", difference, ", ", stddevi, ", ", alpha, ", ", sample, ", ", method, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Difference in means")), differenceEntry, sticky="w")
	tkgrid.configure(differenceEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation in each group")), stddeviEntry, sticky="w")
	tkgrid.configure(stddeviEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size")), sampleEntry, sticky="w")
	tkgrid.configure(sampleEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedSampleProportions <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two proportions"))#Chi-square test with continuity correnction
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	ratio <- tclVar("1")
	ratioEntry <- ttkentry(top, width="20", textvariable=ratio)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	radioButtons(name="continuity", buttons=c("Yes", "No"), values=c(1, 0), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes (or Fisher's exact test)", "No correction")),title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction of chi-square test"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two proportions"), "#####", sep=""))	
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		ratio <- tclvalue(ratio)
		method <- tclvalue(methodVariable)
		continuity <- tclvalue(continuityVariable)
		closeDialog()
		if (length(group1) == 0 || length(group2) == 0){
			errorCondition(recall=StatMedSampleProportions, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0 || length(ratio) == 0){
			errorCondition(recall=StatMedSampleProportions, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleProportion(", group1, ", ", group2, ", ", alpha, ", ", power, ", ", method, ", ", ratio, ", ", continuity, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion in group 1")), group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion in group 2")), group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")), powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size ratio (1:X)")), ratioEntry, sticky="w")
	tkgrid.configure(ratioEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(continuityFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedPowerProportions <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two proportions"))
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	sample1 <- tclVar("")
	sample1Entry <- ttkentry(top, width="20", textvariable=sample1)
	sample2 <- tclVar("")
	sample2Entry <- ttkentry(top, width="20", textvariable=sample2)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	radioButtons(name="continuity", buttons=c("Yes", "No"), values=c(1, 0), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Yes (or Fisher's exact test)", "No correction")),title=gettext(domain="R-RcmdrPlugin.EZR","Continuity correction of chi-square test"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two proportions"), "#####", sep=""))
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		alpha <- tclvalue(alpha)
		sample1 <- as.numeric(tclvalue(sample1))
		sample2 <- as.numeric(tclvalue(sample2))
		if (sample1 >= sample2){
			sample <- sample2
			ratio <- sample1/sample2
		} else {
			sample <- sample1
			ratio <- sample2/sample1
		}
		method <- tclvalue(methodVariable)
		continuity <- tclvalue(continuityVariable)
		closeDialog()
		if (length(group1) == 0 || length(group2) == 0){
			errorCondition(recall=StatMedPowerProportions, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(sample1) == 0 || length(sample2) == 0){
			errorCondition(recall=StatMedPowerProportions, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
#		library(statmod)
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("PowerProportion(", group1, ", ", group2, ", ", alpha, ", ", sample, ", ", method, ", ",  ratio, ", ", continuity, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion in group 1")), group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion in group 2")), group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size of group 1")), sample1Entry, sticky="w")
	tkgrid.configure(sample1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size of group 2")), sample2Entry, sticky="w")
	tkgrid.configure(sample2Entry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(continuityFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedSampleProportionsNonInf <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for non-inferiority trial of two proportions"))
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	delta <- tclVar("")
	deltaEntry <- ttkentry(top, width="20", textvariable=delta)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), initialValue=1, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for non-inferiority trial of two proportions"), "#####", sep=""))
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		delta <- tclvalue(delta)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(group1) == 0 || length(group2) == 0 || length(delta) == 0){
			errorCondition(recall=StatMedSampleProportionsNonInf, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0){
			errorCondition(recall=StatMedSampleProportionsNonInf, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleProportionNonInf(", group1, ", ", group2, ", ", delta, ", ", alpha, ", ", power, ", ", method, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion in control group")), group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Proportion in test group")), group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Meaningful difference in proportion")), deltaEntry, sticky="w")
	tkgrid.configure(deltaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")), powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedSampleHazard <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two survival curves"))
	enrol <- tclVar("")
	enrolEntry <- ttkentry(top, width="20", textvariable=enrol)
	studyperiod <- tclVar("")
	studyperiodEntry <- ttkentry(top, width="20", textvariable=studyperiod)
	followup <- tclVar("")
	followupEntry <- ttkentry(top, width="20", textvariable=followup)
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	ratio <- tclVar("1")
	ratioEntry <- ttkentry(top, width="20", textvariable=ratio)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for comparison between two survival curves"), "#####", sep=""))
		enrol <- tclvalue(enrol)
		studyperiod <- tclvalue(studyperiod)
		followup <- tclvalue(followup)
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		ratio <- tclvalue(ratio)
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(enrol) == 0 || length(studyperiod) == 0 || length(followup) == 0){
			errorCondition(recall=StatMedSampleHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(group1) == 0 || length(group2) == 0){
			errorCondition(recall=StatMedSampleHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0 || length(ratio) == 0){
			errorCondition(recall=StatMedSampleHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleHazard(", enrol, ", ", studyperiod, ", ", followup, ", ", group1, ", ", group2, ", ", alpha, ", ", power, ", ", method, ", ", ratio, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Accrual duration")), enrolEntry, sticky="w")
	tkgrid.configure(enrolEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Total (accrual + follow-up) duration")), studyperiodEntry, sticky="w")
	tkgrid.configure(studyperiodEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival ratio at n year in each group")), followupEntry, sticky="w")
	tkgrid.configure(followupEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival rate in group 1")), group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival rate in group 2")), group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")), powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size ratio (1:X)")), ratioEntry, sticky="w")
	tkgrid.configure(ratioEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedPowerHazard <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two survival curves"))
	enrol <- tclVar("")
	enrolEntry <- ttkentry(top, width="20", textvariable=enrol)
	studyperiod <- tclVar("")
	studyperiodEntry <- ttkentry(top, width="20", textvariable=studyperiod)
	followup <- tclVar("")
	followupEntry <- ttkentry(top, width="20", textvariable=followup)
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	sample1 <- tclVar("")
	sample1Entry <- ttkentry(top, width="20", textvariable=sample1)
	sample2 <- tclVar("")
	sample2Entry <- ttkentry(top, width="20", textvariable=sample2)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate power for comparison between two survival curves"), "#####", sep=""))
		enrol <- tclvalue(enrol)
		studyperiod <- tclvalue(studyperiod)
		followup <- tclvalue(followup)
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		alpha <- tclvalue(alpha)
		sample1 <- as.numeric(tclvalue(sample1))
		sample2 <- as.numeric(tclvalue(sample2))
		if (sample1 >= sample2){
			sample <- sample2
			ratio <- sample1/sample2
		} else {
			sample <- sample1
			ratio <- sample2/sample1
		}
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(enrol) == 0 || length(studyperiod) == 0 || length(followup) == 0){
			errorCondition(recall=StatMedPowerHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(group1) == 0 || length(group2) == 0){
			errorCondition(recall=StatMedPowerHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(sample1) == 0 || length(sample2) == 0){
			errorCondition(recall=StatMedPowerHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("PowerHazard(", enrol, ", ", studyperiod, ", ", followup, ", ", group1, ", ", group2, ", ", alpha, ", ", sample, ", ", method, ", ", ratio, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Accrual duration")), enrolEntry, sticky="w")
	tkgrid.configure(enrolEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Total (accrual + follow-up) duration")), studyperiodEntry, sticky="w")
	tkgrid.configure(studyperiodEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival ratio at n year in each group")), followupEntry, sticky="w")
	tkgrid.configure(followupEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival rate in group 1")), group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival rate in group 2")), group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size of group 1")), sample1Entry, sticky="w")
	tkgrid.configure(sample1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size of group 2")), sample2Entry, sticky="w")
	tkgrid.configure(sample2Entry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedSampleHazardNonInf <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for non-inferiority trial of two survival curves"))
	enrol <- tclVar("")
	enrolEntry <- ttkentry(top, width="20", textvariable=enrol)
	studyperiod <- tclVar("")
	studyperiodEntry <- ttkentry(top, width="20", textvariable=studyperiod)
	followup <- tclVar("")
	followupEntry <- ttkentry(top, width="20", textvariable=followup)
	group1 <- tclVar("")
	group1Entry <- ttkentry(top, width="20", textvariable=group1)
	group2 <- tclVar("")
	group2Entry <- ttkentry(top, width="20", textvariable=group2)
	lowerlimit <- tclVar("")
	lowerlimitEntry <- ttkentry(top, width="20", textvariable=lowerlimit)
	alpha <- tclVar("0.05")
	alphaEntry <- ttkentry(top, width="20", textvariable=alpha)
	power <- tclVar("0.80")
	powerEntry <- ttkentry(top, width="20", textvariable=power)
	ratio <- tclVar("1")
	ratioEntry <- ttkentry(top, width="20", textvariable=ratio)
	radioButtons(name="method", buttons=c("Two.sided", "One.sided"), values=c(2, 1), initialValue=1, labels=gettext(domain="R-RcmdrPlugin.EZR",c("Two-sided", "One-sided")),title=gettext(domain="R-RcmdrPlugin.EZR","Method"))
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Calculate sample size for non-inferiority trial of two survival curves"), "#####", sep=""))
		enrol <- tclvalue(enrol)
		studyperiod <- tclvalue(studyperiod)
		followup <- tclvalue(followup)
		group1 <- tclvalue(group1)
		group2 <- tclvalue(group2)
		lowerlimit <- tclvalue(lowerlimit)
		alpha <- tclvalue(alpha)
		power <- tclvalue(power)
		ratio <- tclvalue(ratio)
		method <- tclvalue(methodVariable)
		closeDialog()
		if (length(enrol) == 0 || length(studyperiod) == 0 || length(followup) == 0){
			errorCondition(recall=StatMedSampleHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(group1) == 0 || length(group2) == 0 || length(lowerlimit) == 0){
			errorCondition(recall=StatMedSampleHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
		if (length(alpha) == 0 || length(power) == 0 || length(ratio) == 0){
			errorCondition(recall=StatMedSampleHazard, message=gettext(domain="R-RcmdrPlugin.EZR","You must select a variable."))
			return()
		}
        if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		command <- paste("SampleHazardNonInf(", enrol, ", ", studyperiod, ", ", followup, ", ", group1, ", ", group2, ", ", lowerlimit, ", ", alpha, ", ", power, ", ", method, ", ", ratio, ")", sep="")
		result <- doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Accrual duration")), enrolEntry, sticky="w")
	tkgrid.configure(enrolEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Total (accrual + follow-up) duration")), studyperiodEntry, sticky="w")
	tkgrid.configure(studyperiodEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival ratio at n year in each group")), followupEntry, sticky="w")
	tkgrid.configure(followupEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival rate in control group")), group1Entry, sticky="w")
	tkgrid.configure(group1Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Survival rate in test group")), group2Entry, sticky="w")
	tkgrid.configure(group2Entry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Non-inferiority lower limit")), lowerlimitEntry, sticky="w")
	tkgrid.configure(lowerlimitEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Alpha error")), alphaEntry, sticky="w")
	tkgrid.configure(alphaEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Power (1 - beta error)")), powerEntry, sticky="w")
	tkgrid.configure(powerEntry, sticky="w")
	tkgrid(tklabel(top, text=gettext(domain="R-RcmdrPlugin.EZR","Sample size ratio (1:X)")), ratioEntry, sticky="w")
	tkgrid.configure(ratioEntry, sticky="w")
	tkgrid(methodFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)
}


StatMedMeta  <- function(){
defaults <- list(studyname=NULL, testpositive=NULL, testnumber=NULL, controlpositive=NULL, controlnumber=NULL, group=NULL, reg=NULL, endpoint="OR", dsl=1, detail=1, funnel=0, subset = "")
dialog.values <- getDialog("StatMedMeta", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Metaanalysis and metaregression for proportions"))
    studynameBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable to identify studies (pick 0 or 1)"), initialSelection=varPosn(dialog.values$studyname, "all"))
    variablesFrame <- tkframe(top)
    testpositiveBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Number of events in test group"), initialSelection=varPosn(dialog.values$testpositive, "all"))
    testnumberBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in test group"), initialSelection=varPosn(dialog.values$testnumber, "all"))
    variables2Frame <- tkframe(top)
    controlpositiveBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Number of events in control group"), initialSelection=varPosn(dialog.values$controlpositive, "all"))
    controlnumberBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in control group"), initialSelection=varPosn(dialog.values$controlnumber, "all"))
    variables3Frame <- tkframe(top)
    groupBox <- variableListBox(variables3Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), initialSelection=varPosn(dialog.values$group, "all"))
    regBox <- variableListBox(variables3Frame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables for meta-regression"), initialSelection=varPosn(dialog.values$reg, "all"))
    radioButtons(name="endpoint", buttons=c("OR", "RR", "RD"), initialValue=dialog.values$endpoint, 
    values=c("OR", "RR", "RD"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Odds ratio", "Relative risk", "Risk difference")),title=gettext(domain="R-RcmdrPlugin.EZR","Summary measure"))
	
    optionsFrame <- tkframe(top)
	checkBoxes(frame="optionsFrame", boxes=c("dsl", "detail", "funnel"), initialValues=c(dialog.values$dsl, dialog.values$detail, dialog.values$funnel),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Conduct random effects meta-analysis", "Show detailed data in forest plot", "Evaluate publication bias with funnel plot")))	
#    checkBoxes(frame="dsl", boxes=c("dsl"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Conduct random effects meta-analysis")))
#    checkBoxes(frame="detail", boxes=c("detail"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show detailed data in forest plot")))
#    checkBoxes(frame="funnel", boxes=c("funnel"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Evaluate publication bias with funnel plot")))
	StatMedSubsetBox(model=TRUE)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Metaanalysis and metaregression for proportions"), "#####", sep=""))
    studyname <- getSelection(studynameBox)
    testpositive <- getSelection(testpositiveBox)
    testnumber <- getSelection(testnumberBox)
    controlpositive <- getSelection(controlpositiveBox)
    controlnumber <- getSelection(controlnumberBox)
    group <- getSelection(groupBox)
    reg <- getSelection(regBox)
	dataSet <- ActiveDataSet()
	subset <- tclvalue(subsetVariable)
	if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
		subset <- ""
	} 
	endpoint <- tclvalue(endpointVariable)
    dsl <- tclvalue(dslVariable)
	detail <- tclvalue(detailVariable)
	funnel <- tclvalue(funnelVariable)	
putDialog("StatMedMeta", list(studyname=studyname, testpositive=testpositive, testnumber=testnumber, controlpositive=controlpositive, controlnumber=controlnumber, group=group, reg=reg, endpoint=endpoint, dsl=dsl, detail=detail, funnel=funnel, subset = tclvalue(subsetVariable)))
    closeDialog()
    if (length(testpositive) == 0 || length(testnumber) == 0 || length(controlpositive) == 0 || length(controlnumber) == 0) {
      errorCondition(recall=StatMedMeta, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick all required variables"))
      return()
    }
	if (length(studyname) == 0 ){
		studyname <- NULL
	}
	if (length(group) == 0 ){
		group1 <- NULL
		group2 <- NULL
	}
	else {
		group1 <- paste(", byvar=", group, ', bylab="', group, '"')
		group2 <- paste(', bylab="', group, '"')
	}
	if (subset==""){
		doItAndPrint(paste("TempDF <- ", dataSet, "[complete.cases(", dataSet, "$", testpositive, ", ", dataSet, "$", testnumber, ", ", dataSet, "$", controlpositive, ", ", dataSet, "$", controlnumber, "),]", sep=""))	
	}else{
		doItAndPrint(paste("TempDF <- subset(", dataSet, ", subset=", subset, ")[complete.cases(subset(", dataSet, ", subset=", subset, ")$", testpositive, ", subset(", dataSet, ", subset=", subset, ")$", testnumber, ", subset(", dataSet, ", subset=", subset, ")$", controlpositive, ", subset(", dataSet, ", subset=", subset, ")$", controlnumber, "),]", sep=""))
	}
#    library(meta, quietly=TRUE)
	Library("meta")
	if (dsl==0) {
		command <- paste("res <- metabin(", testpositive, ", ", testnumber, ", ", controlpositive, ", ", controlnumber, ', data=TempDF, sm="', endpoint, '", studlab=', studyname, group1, ", comb.fixed=TRUE, comb.random=FALSE)", sep="")
	} else {
		command <- paste("res <- metabin(", testpositive, ", ", testnumber, ", ", controlpositive, ", ", controlnumber, ', data=TempDF, sm="', endpoint, '", studlab=', studyname, group1, ", comb.fixed=TRUE, comb.random=TRUE)", sep="")
	}
	doItAndPrint(command)
	doItAndPrint("res")
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
	if (detail == 0){
		doItAndPrint("plot(res)")
	} else{
		doItAndPrint(paste("forest.meta(res", group2, ")", sep=""))
	}
	if (funnel == 1) {
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		doItAndPrint("funnel(res)")
		doItAndPrint("metabias(res)")
	}
	if (length(reg) > 0) {
		doItAndPrint("Var <- (res$seTE)^2")
		doItAndPrint("library(metatest, quietly=TRUE)")
		for (i in 1:length(reg)){
			doItAndPrint("y <- exp(res$TE)")
			doItAndPrint(paste("(metareg <- metatest(res$TE~TempDF$", reg[i], ", Var))", sep=""))
			doItAndPrint(paste("x <- TempDF$", reg[i], sep=""))
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}			
			doItAndPrint("y.L <- exp(res$TE-qnorm(0.975)*res$seTE)")
			doItAndPrint("y.H <- exp(res$TE+qnorm(0.975)*res$seTE)")
			doItAndPrint("max.weight <- sqrt(max(res$w.fixed))")
			doItAndPrint(paste('plot(y ~ x, ylab="Effect size", xlab="', reg[i], '", log="y", pch=15, cex=sqrt(res$w.fixed)*2.5/max.weight, ylim=c(min(y.L), max(y.H)))', sep=""))
			doItAndPrint("if(is.numeric(x)) arrows(x, y.L, x, y.H, code=3, angle=90, length=0.1)")
			doItAndPrint("metareg.table <- signif(cbind(metareg$coefficients, metareg$se, metareg$coef-qnorm(0.975)*metareg$se, metareg$coef+qnorm(0.975)*metareg$se, metareg$pZtest), digits=3)")
			doItAndPrint(paste('rownames(metareg.table) <- c("(Intercept)", "', reg[i], '")', sep=""))
			doItAndPrint('colnames(metareg.table) <- c("Coef", "SE", "Lower 95%CI", "Upper 95%CI", "p.value")')
			doItAndPrint("metareg.table<- data.frame(metareg.table)")
			doItAndPrint("metareg.table")
#			doItAndPrint("remove(metareg.table)")
		}
	}
	doItAndPrint("remove(res)")
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="metabin", apply="StatMedMeta", reset="StatMedMeta")
    tkgrid(getFrame(studynameBox), sticky="nw")
    tkgrid(getFrame(testpositiveBox), labelRcmdr(variablesFrame, text="    "), getFrame(testnumberBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(getFrame(controlpositiveBox), labelRcmdr(variables2Frame, text="    "), getFrame(controlnumberBox), sticky="nw")
	tkgrid(variables2Frame, sticky="nw")
	tkgrid(labelRcmdr(variables3Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
    tkgrid(getFrame(groupBox), labelRcmdr(variables3Frame, text="    "), getFrame(regBox), sticky="nw")
	tkgrid(variables3Frame, sticky="nw")
#    tkgrid(endpointFrame, sticky="w")
#    tkgrid(optionsFrame, sticky="w")

tkgrid(optionsFrame, endpointFrame, sticky="w")

#    tkgrid(dsl, sticky="w")
#    tkgrid(detail, sticky="w")
#    tkgrid(funnel, sticky="w")
	tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedMetaHazard  <- function(){
defaults <- list(input="CI", studyname=NULL, hazard=NULL, ci=NULL, group=NULL, reg=NULL, dsl=1, detail=1, funnel=0, subset = "")
dialog.values <- getDialog("StatMedMetaHazard", defaults)
currentFields$subset <- dialog.values$subset
currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Metaanalysis and metaregression for hazard ratios"))
    radioButtons(name="input", buttons=c("CI", "SE"), initialValue=dialog.values$input, 
    values=c("CI", "SE"), labels=gettext(domain="R-RcmdrPlugin.EZR",c("Combine hazard ratio and 95% confidence interval", "Combine log hazard ratio and standard error")),title=gettext(domain="R-RcmdrPlugin.EZR","Choose data to combine"))
    studynameBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable to identify studies (pick 0 or 1)"), initialSelection=varPosn(dialog.values$studyname, "all"))
    variablesFrame <- tkframe(top)
    hazardBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Hazard ratio or log hazard ratio (pick one)"), initialSelection=varPosn(dialog.values$hazard, "all"))
    ciBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Upper limit of 95% confidence interval or standard error (pick one)"), initialSelection=varPosn(dialog.values$ci, "all"))
    variables2Frame <- tkframe(top)
    groupBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), initialSelection=varPosn(dialog.values$group, "all"))
	regBox <- variableListBox(variables2Frame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables for meta-regression"), initialSelection=varPosn(dialog.values$reg, "all"))
    optionsFrame <- tkframe(top)
	checkBoxes(frame="optionsFrame", boxes=c("dsl", "detail", "funnel"), initialValues=c(dialog.values$dsl, dialog.values$detail, dialog.values$funnel),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Conduct random effects meta-analysis", "Show detailed data in forest plot", "Evaluate publication bias with funnel plot")))		
#    checkBoxes(frame="dsl", boxes=c("dsl"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Conduct random effects meta-analysis")))
#    checkBoxes(frame="detail", boxes=c("detail"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show detailed data in forest plot")))
#    checkBoxes(frame="funnel", boxes=c("funnel"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Evaluate publication bias with funnel plot.")))
	StatMedSubsetBox(model=TRUE)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Metaanalysis and metaregression for hazard ratios"), "#####", sep=""))
    studyname <- getSelection(studynameBox)
    hazard <- getSelection(hazardBox)
    upperci <- getSelection(ciBox)
    group <- getSelection(groupBox)
	reg <- getSelection(regBox)
	dataSet <- ActiveDataSet()
	subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")){
		subset <- ""
	}
	input <- tclvalue(inputVariable)
    dsl <- tclvalue(dslVariable)
	detail <- tclvalue(detailVariable)
	funnel <- tclvalue(funnelVariable)	
putDialog("StatMedMetaHazard",list(input=input, studyname=studyname, hazard=hazard, ci=upperci, group=group, reg=reg, dsl=dsl, detail=detail, funnel=funnel, subset = tclvalue(subsetVariable)))
	closeDialog()
    if (length(hazard) == 0 || length(upperci) == 0) {
      errorCondition(recall=StatMedMetaHazard, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick all required variables"))
      return()
    }
	if (length(studyname) == 0 ){
		studyname <- NULL
	}
	if (length(group) == 0 ){
		group1 <- NULL
		group2 <- NULL
	}
	else {
		group1 <- paste(", byvar=", group, ', bylab="', group, '"')
		group2 <- paste(', bylab="', group, '"')
	}

	if (subset==""){
		doItAndPrint(paste("TempDF <- ", dataSet, "[complete.cases(", dataSet, "$", hazard, ", ", dataSet, "$", upperci, "),]", sep=""))	
	} else {
		doItAndPrint(paste("TempDF <- subset(", dataSet, ", subset=", subset, ")[complete.cases(subset(", dataSet, ", subset=", subset, ")$", hazard, ", subset(", dataSet, ", subset=", subset, ")$", upperci, "),]", sep=""))
	}
#    library(meta, quietly=TRUE)
	Library("meta")
	if (input == "CI"){
		doItAndPrint(paste("logHR <- log(TempDF$", hazard, ")", sep=""))
		doItAndPrint(paste("logSE <- (log(TempDF$", upperci, ")-log(TempDF$", hazard, ")) / qnorm(0.975)", sep=""))				
	} else {
		doItAndPrint(paste("logHR <- TempDF$", hazard, sep=""))
		doItAndPrint(paste("logSE <- TempDF$", upperci, sep=""))		
	}
	if (dsl==0) {
		command <- paste('res <- metagen(logHR, logSE, data=TempDF, sm="HR", studlab=', studyname, group1, ", comb.fixed=TRUE, comb.random=FALSE)", sep="")
	} else {
		command <- paste('res <- metagen(logHR, logSE, data=TempDF, sm="HR", studlab=', studyname, group1, ", comb.fixed=TRUE, comb.random=TRUE)", sep="")
	}
	doItAndPrint(command)
	doItAndPrint("res")
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
	if (detail == 0){
		doItAndPrint("plot(res)")
	}
	else{
		doItAndPrint(paste("forest.meta(res", group2, ")", sep=""))
	}
	if (funnel == 1) {
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		doItAndPrint("funnel(res)")
		doItAndPrint("metabias(res)")
	}
	if (length(reg) > 0) {
		doItAndPrint("Var <- (res$seTE)^2")
		doItAndPrint("library(metatest, quietly=TRUE)")
		for (i in 1:length(reg)){
			doItAndPrint("y <- exp(res$TE)")
			doItAndPrint(paste("(metareg <- metatest(res$TE~TempDF$", reg[i], ", Var))", sep=""))
			doItAndPrint(paste("x <- TempDF$", reg[i], sep=""))
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint("y.L <- exp(res$TE-qnorm(0.975)*res$seTE)")
			doItAndPrint("y.H <- exp(res$TE+qnorm(0.975)*res$seTE)")
			doItAndPrint("max.weight <- sqrt(max(res$w.fixed))")
			doItAndPrint(paste('plot(y ~ x, ylab="Effect size", xlab="', reg[i], '", log="y", pch=15, cex=sqrt(res$w.fixed)*2.5/max.weight, ylim=c(min(y.L), max(y.H)))', sep=""))
			doItAndPrint("if(is.numeric(x)) arrows(x, y.L, x, y.H, code=3, angle=90, length=0.1)")
			doItAndPrint("metareg.table <- signif(cbind(metareg$coefficients, metareg$se, metareg$coef-qnorm(0.975)*metareg$se, metareg$coef+qnorm(0.975)*metareg$se, metareg$pZtest), digits=3)")
			doItAndPrint(paste('rownames(metareg.table) <- c("(Intercept)", "', reg[i], '")', sep=""))
			doItAndPrint('colnames(metareg.table) <- c("Coef", "SE", "Lower 95%CI", "Upper 95%CI", "p.value")')
			doItAndPrint("metareg.table<- data.frame(metareg.table)")
			doItAndPrint("metareg.table")
#			doItAndPrint("remove(metareg.table)")
		}
	}
	doItAndPrint("remove(res)")
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="metagen", apply="StatMedMetaHazard", reset="StatMedMetaHazard")
    tkgrid(inputFrame, sticky="w")
    tkgrid(getFrame(studynameBox), sticky="nw")
    tkgrid(getFrame(hazardBox), labelRcmdr(variablesFrame, text="    "), getFrame(ciBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
	tkgrid(labelRcmdr(variables2Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
    tkgrid(getFrame(groupBox), labelRcmdr(variables2Frame, text="    "), getFrame(regBox), sticky="nw")
	tkgrid(variables2Frame, sticky="nw")
	tkgrid(optionsFrame, sticky="nw")
 #   tkgrid(dsl, sticky="w")
 #   tkgrid(funnel, sticky="w")
	tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


StatMedMetaCont  <- function(){
defaults <- list(studyname=NULL, testmean=NULL, testnumber=NULL, testsd=NULL, controlmean=NULL, controlnumber=NULL, controlsd=NULL, group=NULL, reg=NULL, smd=0, dsl=1, detail=1, funnel=0, smd=0, subset = "")
dialog.values <- getDialog("StatMedMetaCont", defaults)
currentFields$subset <- dialog.values$subset	
currentModel <- TRUE
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","Metaanalysis and metaregression for means"))
    studynameBox <- variableListBox(top, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Variable to identify studies (pick 0 or 1)"), initialSelection=varPosn(dialog.values$studyname, "all"))
    variablesFrame <- tkframe(top)
    testmeanBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Mean in test group"), initialSelection=varPosn(dialog.values$testmean, "all"))
    testnumberBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in test group"), initialSelection=varPosn(dialog.values$testnumber, "all"))
    testsdBox <- variableListBox(variablesFrame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation in test group"), initialSelection=varPosn(dialog.values$testsd, "all"))
    variables2Frame <- tkframe(top)
    controlmeanBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Mean in control group"), initialSelection=varPosn(dialog.values$controlmean, "all"))
    controlnumberBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Number of samples in control group"), initialSelection=varPosn(dialog.values$controlnumber, "all"))
    controlsdBox <- variableListBox(variables2Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Standard deviation in control group"), initialSelection=varPosn(dialog.values$controlsd, "all"))
    variables3Frame <- tkframe(top)
    groupBox <- variableListBox(variables3Frame, Variables(), title=gettext(domain="R-RcmdrPlugin.EZR","Grouping variable(pick 0 or 1)"), initialSelection=varPosn(dialog.values$group, "all"))
	regBox <- variableListBox(variables3Frame, Variables(), selectmode="multiple", title=gettext(domain="R-RcmdrPlugin.EZR","Variables for meta-regression"), initialSelection=varPosn(dialog.values$reg, "all"))
	
    optionsFrame <- tkframe(top)
	checkBoxes(frame="optionsFrame", boxes=c("smd", "dsl", "detail", "funnel"), initialValues=c(dialog.values$smd, dialog.values$dsl, dialog.values$detail, dialog.values$funnel),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Pool standard mead difference", "Conduct random effects meta-analysis", "Show detailed data in forest plot", "Evaluate publication bias with funnel plot")))		
#    checkBoxes(frame="dsl", boxes=c("dsl"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Conduct random effects meta-analysis")))
#    checkBoxes(frame="detail", boxes=c("detail"),initialValues=c(1),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Show detailed data in forest plot")))
#    checkBoxes(frame="funnel", boxes=c("funnel"),initialValues=c(0),labels=gettext(domain="R-RcmdrPlugin.EZR",c("Evaluate publication bias with funnel plot")))
	StatMedSubsetBox(model=TRUE)
	onOK <- function(){
	logger(paste("#####", gettext(domain="R-RcmdrPlugin.EZR","Metaanalysis and metaregression for means"), "#####", sep=""))
    studyname <- getSelection(studynameBox)
    testmean <- getSelection(testmeanBox)
    testnumber <- getSelection(testnumberBox)
    testsd <- getSelection(testsdBox)
    controlmean <- getSelection(controlmeanBox)
    controlnumber <- getSelection(controlnumberBox)
    controlsd <- getSelection(controlsdBox)
	group <- getSelection(groupBox)
	reg <- getSelection(regBox)
	dataSet <- ActiveDataSet()
	subset <- tclvalue(subsetVariable)
	if (trim.blanks(subset) == gettext(domain="R-RcmdrPlugin.EZR","<all valid cases>")) {
		subset <- ""
	} 
	if (length(studyname) == 0 ){
		studyname <- NULL
	}
	if (length(group) == 0 ){
		group1 <- NULL
		group2 <- NULL
	}
	else {
		group1 <- paste(", byvar=", group, ', bylab="', group, '"')
		group2 <- paste(', bylab="', group, '"')
	}
    smd <- tclvalue(smdVariable)
    dsl <- tclvalue(dslVariable)
	detail <- tclvalue(detailVariable)
	funnel <- tclvalue(funnelVariable)
	if (subset==""){
		doItAndPrint(paste("TempDF <- ", dataSet, "[complete.cases(", dataSet, "$", testnumber, ", ", dataSet, "$", testmean, ", ", dataSet, "$", testsd, ", ", dataSet, "$", controlnumber, ", ", dataSet, "$", controlmean, ", ", dataSet, "$", controlsd, "),]", sep=""))	
	}else{
		doItAndPrint(paste("TempDF <- subset(", dataSet, ", subset=", subset, ")[complete.cases(subset(", dataSet, ", subset=", subset, ")$", testnumber, ", subset(", dataSet, ", subset=", subset, ")$", testmean, ", subset(", dataSet, ", subset=", subset, ")$", testsd, ", subset(", dataSet, ", subset=", subset, ")$", controlnumber, ", subset(", dataSet, ", subset=", subset, ")$", controlmean, ", subset(", dataSet, ", subset=", subset, ")$", controlsd, "),]", sep=""))
	}
putDialog("StatMedMetaCont", list(studyname=studyname, testmean=testmean, testnumber=testnumber, testsd=testsd, controlmean=controlmean, controlnumber=controlnumber, controlsd=controlsd, group=group, reg=reg, dsl=dsl, detail=detail, funnel=funnel, smd=smd, subset = tclvalue(subsetVariable)))
    closeDialog()
    if (length(testmean) == 0 || length(testnumber) == 0 || length(testsd) == 0 || length(controlmean) == 0 || length(controlnumber) == 0 || length(controlsd) == 0) {
      errorCondition(recall=StatMedMetaCont, 
        message=gettext(domain="R-RcmdrPlugin.EZR","Pick all required variables"))
      return()
    }	
#    library(meta, quietly=TRUE)
	Library("meta")
	smd <- ifelse(smd==0, ', sm="MD"', ', sm="SMD"')
	if (dsl==0) {
		command <- paste("res <- metacont(", testnumber, ", ", testmean, ", ", testsd, ", ", controlnumber, ", ", controlmean, ", ", controlsd, ", data=TempDF, studlab=", studyname, group1, ", comb.fixed=TRUE, comb.random=FALSE", smd, ")", sep="")
	} else {
		command <- paste("res <- metacont(", testnumber, ", ", testmean, ", ", testsd, ", ", controlnumber, ", ", controlmean, ", ", controlsd, ", data=TempDF, studlab=", studyname, group1, ", comb.fixed=TRUE, comb.random=TRUE", smd, ")", sep="")
	}
	doItAndPrint(command)
	doItAndPrint("res")
	if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
	if (detail == 0){
		doItAndPrint("plot(res)")
	}
	else{
		doItAndPrint(paste("forest.meta(res", group2, ")", sep=""))
	}
	if (funnel == 1) {
		if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
		doItAndPrint("funnel(res)")
		doItAndPrint("metabias(res)")
	}
		if (length(reg) > 0) {
		doItAndPrint("Var <- (res$seTE)^2")
		doItAndPrint("library(metatest, quietly=TRUE)")
		for (i in 1:length(reg)){
			doItAndPrint("y <- res$TE")
			doItAndPrint(paste("(metareg <- metatest(res$TE~TempDF$", reg[i], ", Var))", sep=""))
			doItAndPrint(paste("x <- TempDF$", reg[i], sep=""))
			if (.Platform$OS.type == 'windows'){doItAndPrint(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else if (MacOSXP()==TRUE) {doItAndPrint(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))} else {doItAndPrint(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))}
			doItAndPrint("y.L <- res$TE-qnorm(0.975)*res$seTE")
			doItAndPrint("y.H <- res$TE+qnorm(0.975)*res$seTE")
			doItAndPrint("max.weight <- sqrt(max(res$w.fixed))")
			doItAndPrint(paste('plot(y ~ x, ylab="Effect size", xlab="', reg[i], '", pch=15, cex=sqrt(res$w.fixed)*2.5/max.weight, ylim=c(min(y.L), max(y.H)))', sep=""))
			doItAndPrint("if(is.numeric(x)) arrows(x, y.L, x, y.H, code=3, angle=90, length=0.1)")
			doItAndPrint("metareg.table <- signif(cbind(metareg$coefficients, metareg$se, metareg$coef-qnorm(0.975)*metareg$se, metareg$coef+qnorm(0.975)*metareg$se, metareg$pZtest), digits=3)")
			doItAndPrint(paste('rownames(metareg.table) <- c("(Intercept)", "', reg[i], '")', sep=""))
			doItAndPrint('colnames(metareg.table) <- c("Coef", "SE", "Lower 95%CI", "Upper 95%CI", "p.value")')
			doItAndPrint("metareg.table<- data.frame(metareg.table)")
			doItAndPrint("metareg.table")
#			doItAndPrint("remove(metareg.table)")
		}
	}
	doItAndPrint("remove(res)")
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="metacont", apply="StatMedMetaCont", reset="StatMedMetaCont")
    tkgrid(getFrame(studynameBox), sticky="nw")
    tkgrid(getFrame(testmeanBox), labelRcmdr(variablesFrame, text="    "), getFrame(testnumberBox), getFrame(testsdBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(getFrame(controlmeanBox), labelRcmdr(variables2Frame, text="    "), getFrame(controlnumberBox), getFrame(controlsdBox), sticky="nw")
	tkgrid(variables2Frame, sticky="nw")
	tkgrid(labelRcmdr(variables3Frame, text=gettext(domain="R-RcmdrPlugin.EZR","Click pressing Ctrl key to select multiple variables."), fg="blue"), sticky="w")
    tkgrid(getFrame(groupBox), labelRcmdr(variables3Frame, text="    "), getFrame(regBox), sticky="nw")
	tkgrid(variables3Frame, sticky="nw")
	tkgrid(optionsFrame, sticky="nw")
#    tkgrid(dsl, sticky="w")
#    tkgrid(funnel, sticky="w")
	tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}


EZRVersion <- function(){
	initializeDialog(title=gettext(domain="R-RcmdrPlugin.EZR","EZR version"))
	onOK <- function(){
		closeDialog()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="Rcmdr")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR","  EZR on R commander (programmed by Y.Kanda) "), fg="blue"), sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR"," "), fg="blue"), sticky="w")
	tkgrid(labelRcmdr(top, text=paste("      ", gettext(domain="R-RcmdrPlugin.EZR","Current version:"), " 1.32", sep="")), sticky="w")
	tkgrid(labelRcmdr(top, text=paste("        ", gettext(domain="R-RcmdrPlugin.EZR","February 1, 2016"), sep="")), sticky="w")
	tkgrid(labelRcmdr(top, text=gettext(domain="R-RcmdrPlugin.EZR"," "), fg="blue"), sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}


StatMedCloseCommander <- function() StatMedcloseCommander(ask=getRcmdr("ask.to.exit"), ask.save=getRcmdr("ask.on.exit"))


StatMedcloseCommanderAndR <- function(){
	response <- StatMedCloseCommander()
	if (response == "cancel") return()
	cat("\n")
	quit(save="no")
}


#StatMedcloseCommander <- function(ask=TRUE, ask.save=ask){
#	if (ask){
#		response <- tclvalue(RcmdrTkmessageBox(message=gettext(domain="R-RcmdrPlugin.EZR","Exit?"),
#						icon="question", type="okcancel", default="cancel"))
#		if (response == "cancel") return(invisible(response))
#	}
#	else {
#		ask.save=FALSE
#		response <- "ok"
#	}
#	sink(type="message")

###add save data function
#	if (ask.save && !is.null(ActiveDataSet())){	
#
#		logger("Active dataset")
#		response1 <- RcmdrTkmessageBox(message=gettext(domain="R-RcmdrPlugin.EZR","Save active dataset?"),
#				icon="question", type="yesno", default="yes")
#		if ("yes" == tclvalue(response1)){
#			file <- tclvalue(tkgetSaveFile(filetypes=
#				gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"R Data Files" {".rda" ".Rda" ".RDA" ".RData"}}'),
#				defaultextension="rda", initialfile=paste(activeDataSet(), "rda", sep=".")))
#			if (file == "") return()
#			command <- paste('save("', activeDataSet(), '", file="', file, '")', sep="")
#			justDoIt(command)
#			logger(command)
#		}
#	}
	
#	if (!is.null(ActiveDataSet()) && getRcmdr("attach.data.set"))
#		justDoIt(logger(paste("detach(", ActiveDataSet(), ")", sep="")))
#	putRcmdr(".activeDataSet", NULL)
#	putRcmdr(".activeModel", NULL)
#	if (ask.save && getRcmdr("log.commands") && tclvalue(tkget(LogWindow(), "1.0", "end")) != "\n"){
#		response2 <- RcmdrTkmessageBox(message=gettext(domain="R-RcmdrPlugin.EZR","Save script file?"),
#				icon="question", type="yesno", default="yes")
#		if ("yes" == tclvalue(response2)) saveLog()
#
#	if (ask.save && !getRcmdr("console.output") && tclvalue(tkget(OutputWindow(), "1.0", "end")) != "\n"){
#		response3 <- RcmdrTkmessageBox(message=gettext(domain="R-RcmdrPlugin.EZR","Save output file?"),
#				icon="question", type="yesno", default="yes")
#		if ("yes" == tclvalue(response3)) saveOutput()
#	}
#	if (.Platform$OS.type != "windows") options(getRcmdr("oldPager"))
#	if (getRcmdr("suppress.X11.warnings")) {
#		sink(type = "message")
#		close(getRcmdr("messages.connection"))
#	}
#	options(getRcmdr("saveOptions"))
#	tkdestroy(CommanderWindow())
#	putRcmdr("commanderWindow", NULL)
#	putRcmdr("logWindow", NULL)
#	putRcmdr("messagesWindow", NULL)
#	putRcmdr("outputWindow", NULL)
#	options(getRcmdr("quotes"))
#	tkwait <- options("Rcmdr")[[1]]$tkwait  # to address problem in Debian Linux
#	if ((!is.null(tkwait)) && tkwait) putRcmdr(".commander.done", tclVar("1"))
#	return(invisible(response))
#}


StatMedcloseCommander <- function(ask=TRUE, ask.save=ask){

###add save data function
	if (!is.null(ActiveDataSet())){	
		logger("Active dataset")
		response1 <- RcmdrTkmessageBox(message=gettext(domain="R-RcmdrPlugin.EZR","Save active dataset?"),
				icon="question", type="yesno", default="yes")
		if ("yes" == tclvalue(response1)){
			file <- tclvalue(tkgetSaveFile(filetypes=
				gettext(domain="R-RcmdrPlugin.EZR",'{"All Files" {"*"}} {"R Data Files" {".rda" ".Rda" ".RDA" ".RData"}}'),
				defaultextension="rda", initialfile=paste(activeDataSet(), "rda", sep=".")))
			if (file == "") return()
			command <- paste('save("', activeDataSet(), '", file="', file, '")', sep="")
			justDoIt(command)
			logger(command)
		}
	}
		
	closeCommander(ask=TRUE, ask.save=ask)
#	closeCommander()
}


EZRhelp <- function(){
	flag <- 0
	for(i in search()) if(i=="package:RcmdrPlugin.EZR")flag <- 1
	if(flag==0){
		doItAndPrint('browseURL(paste(file.path(path.package(package="Rcmdr"), "doc"), "/", "EZR.htm", sep=""))')
	}else{
		doItAndPrint("help(EZR)")
	}
}


EZR <- function(){
	cat(gettext(domain="R-RcmdrPlugin.EZR","EZR on R commander (programmed by Y.Kanda) Version 1.32", "\n"))
}

if (getRversion() >= '2.15.1') globalVariables(c('top', 'buttonsFrame',
'TempTD', 'actmodelVariable', 'subsetVariable', 
'subsetFrame', 'oneWayAnova', 'graphVariable', 'pairwiseVariable',
'dunnettVariable', 'bonferroniVariable', 'holmVariable', 'graphFrame',
'lineVariable', 'placeVariable', 'censorVariable', 'atriskVariable',
'xscaleVariable', 'lineFrame', 'placeFrame', 'xscaleFrame', 'censor',
'atrisk', 'colorVariable', 'besideVariable', 'percentVariable',
'errorBarsVariable', 'errorBarsFrame', 'levelsVariable', 'binVariable',
'methodVariable', 'subdialog', 'subButtonsFrame', 'entry1', 'onCancel',
'levelNames', 'levelsFrame', 'methodFrame', 'logyVariable',
'whiskerVariable', 'logy', 'whiskerFrame', 'lhsVariable',
'rhsVariable', 'onHelp', 'xBox', 'outerOperatorsFrame', 'formulaFrame',
'checkboxFrame', 'lhsEntry', 'paletteVariable', 'paletteFrame',
'alternativeVariable', 'alternativeFrame', 'waldVariable',
'prophazVariable', 'basecurveVariable', 'stepwise1Variable',
'stepwise2Variable', 'stepwise3Variable', 'SurvivalTimeVariable',
'StatusVariable', 'posthocVariable', 'posthocFrame', 'ymdVariable',
'ymdFrame', 'percentsVariable', 'chisqVariable',
'chisqComponentsVariable', 'expFreqVariable', 'fisherVariable',
'.Test', '.Table', 'percentsFrame', 'testsFrame', 'optionsFrame',
'delimiterFrame', 'delimiterVariable', 'colnamesVariable',
'rownamesVariable', 'quotesVariable', 'numericToFactor', 'filterNA',
'subwin', '.Probs', '.Responses', 'window.sizeVariable',
'window.typeVariable', 'lwdVariable', 'lasVariable', 'familyVariable',
'cexVariable', 'window.sizeFrame', 'window.typeFrame', 'lwdFrame',
'lasFrame', 'familyFrame', 'cexFrame', 'scaleVariable', 'color',
'scaleFrame', 'importMinitab', 'importRODBCtable', 'importSPSS',
'importSTATA', 'ciVariable', 'separatestrataVariable', 'testVariable',
'testFrame', 'steeldwassVariable', 'steelVariable', 'logVariable',
'multiVariable', 'y', 'linearRegressionModel', 'helpButton',
'baseVariable', 'baseFrame', 'continuityVariable', 'continuityFrame',
'endpointVariable', 'dslVariable', 'detailVariable', 'funnelVariable',
'endpointFrame', 'inputVariable', 'inputFrame', 'interactionVariable',
'numbersButton', 'namesButton', 'meanVariable', 'sdVariable',
'.groups', 'checkBoxFrame', 'groupsFrame', 'unmatchVariable',
'unmatchFrame', 'typeVariable', 'trendVariable', 'typeFrame',
'trendFrame', 'chisqTestVariable', 'exactTestVariable',
'directionVariable', 'bestVariable', 'thresholdVariable',
'directionFrame', 'bestFrame', 'locationVariable', 'decimalVariable',
'locationFrame', 'decimalFrame', 'renameVariables', 'chrtofacVariable',
'chrtofac', 'reorderFactor', 'removeVariable', 'removeFrame',
'twostageVariable', 'twostage', 'jitterXVariable', 'jitterYVariable',
'logXVariable', 'logYVariable', 'identifyVariable', 'boxplotsVariable',
'lsLineVariable', 'smoothLineVariable', 'spreadVariable',
'scatterPlot', 'diagonalVariable', 'diagonalFrame',
'StaMedSetContrasts', 'contrastsVariable', 'contrastsFrame', '..',
'hex.1', 'hex.2', 'hex.3', 'hex.4', 'hex.5', 'hex.6', 'hex.7', 'hex.8',
'decreasingVariable', 'decreasingFrame', 'Stack', 'partsVariable',
'styleVariable', 'trimOutliersVariable', 'showDepthsVariable',
'reverseNegativeVariable', 'partsFrame', 'styleFrame',
'variancesVariable', 'variancesFrame', 'fisherTestVariable', 'saveLog',
'saveOutput', '.commander.done', 'ci.summary.table', 'cox.table',
'km.summary.table', 'summary.ttest', 'Fisher.summary.table', 
'StatMedcloseCommander', 'hist2', 'separatestrata', 'diagnosisVariable',
'martinVariable', 'res', 'HistEZR', 'QQPlot', '.Workbook', 'par.lwd', 'par.cex',
'getSheets', 'analysisVariable', 'outputVariable', 'languageVariable',
'analysisFrame', 'outputFrame', 'languageFrame', 'exactVariable', 
'rangeVariable', 'explainVariable', 'exactFrame', 'rangeFrame',
'explainFrame', 'multireg.table', 'smdVariable', 'survfit', 'survdiff',
'odbcCloseAll', 'odbcConnectExcel', 'odbcConnectExcel2007', 'odbcConnectAccess',
'odbcConnectAccess2007', 'odbcConnectDbase', 'sqlTables', '.Tcl.args',
'cuminc', 'Anova', 'pmvt', 'wald.test', 'timepoints', 'ci', 'sqlQuery',
'groupingVariable', 'groupingFrame', 'othervarVariable', 'rocVariable',
'columnmergeVariable', 'column.name1', 'column.name2', 'columnmergeFrame',
'deleteVariable', 'RecodeDialog'))