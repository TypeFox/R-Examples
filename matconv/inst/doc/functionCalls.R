## ----setup, include=FALSE------------------------------------------------
library(matconv)
library(knitr)

knitr::opts_chunk$set(fig.pos='center', echo=FALSE, comment='>', results='asis')

showExampleConv <- function(matIn, dict = ""){
	maps <- makeFuncMaps(dict)
	lout <- mat2r(matIn, funcConverters = maps, verbose = 0)

	lout[[2]] <- gsub("\\$", "\\\\$", lout[[2]])
	diff <- length(lout[[2]]) - length(lout[[1]])
	if (diff > 0){
		lout[[1]] <- c(lout[[1]], rep("", diff))
	}
	output <- c("<table><tr>",
    paste0("<th> Dictionary Used </th>",
		"<th>",names(lout)[1],"</th>"),
		paste0("<th>",names(lout)[2],"</th>"),
		"<tr> <td valign = top>",
		paste('|', dict, "\n"),
		"</td> <td valign = top>",
		paste('|', lout[[1]], "\n"),
		"</td>", "<td valign = top>",
		paste('|', lout[[2]], "\n"),
		"</td> </table>")

	cat(paste(output, collapse = "\n"))
}


## ----first_example-------------------------------------------------------

showExampleConv(
	"thing = linspace(first, second)",
	"linspace:seq, 2, 1")


## ----Function_literals---------------------------------------------------
showExampleConv(
	c("thing = linspace(1, 2, 3)",
	"hjkl = binornd(2.3, 1.5)",
	"asdf = erf(2)"),
	c("linspace:seq, 1, 2, len = %3",
	"binornd:rbinom, 1L, 1, 2",
	"erf: , 2 * pnorm(%1 * sqrt(2)) - 1"))


## ----Func_switcher-------------------------------------------------------
showExampleConv(
	c("thing = rand(1, 5)",
	 "thing = rand(5, 1)",
	 "thing = rand(5, 5)"),
	c("rand--if 1 == 1L:runif, 2",
    "rand--if 2 == 1L:runif, 1",
    "rand--if finally:matrix, runif(%1 * %2), %1)"))


## ----mult_out------------------------------------------------------------
showExampleConv(
  c("[myL myU myP] = lu(badMatrix)"),
  c("lu: , expand(lu(Matrix::Matrix(%1))) --out L U P")
  )


