## ----setup, include=FALSE------------------------------------------------
library(matconv)
library(knitr)

knitr::opts_chunk$set(fig.pos='center', echo=FALSE, comment='>', results='asis')

showExampleConv <- function(matIn){
	lout <- mat2r(matIn, verbose = 0)
	diff <- length(lout[[2]]) - length(lout[[1]])
	if (diff > 0){
		lout[[1]] <- c(lout[[1]], rep("", diff))
	}
	output <- c("<table><tr>",
    paste0("<th>",names(lout)[1],"</th>"),
		paste0("<th>",names(lout)[2],"</th>"),
		"<tr> <td valign = top>",
		paste('|', lout[[1]], "\n"),
		"</td>", "<td valign = top>",
		paste('|', lout[[2]], "\n"),
		"</td> </table>")
	
	cat(paste(output, collapse = "\n"))
}


## ----basSyntax-----------------------------------------------------------
example <- c( "thing = 5 * 3;", "thing2 = (thing ~= 14);")
showExampleConv(example)


## ----flow----------------------------------------------------------------
example <- c("if argLen == 1", "  doThing = 9999;","else", "  doThing = 1;","end")
showExampleConv(example)


## ----userFunctions-------------------------------------------------------
example <- c("function [ dat ] = xlsReadPretty(varargin))", "  didThing = 1*3;", "  dat = didThing / 3;", "end")
showExampleConv(example)


