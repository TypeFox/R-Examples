## ----set-options, echo=FALSE, cache=FALSE-----------------------------------------------------------------------------
options(width=120)
#opts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")
#read_chunk("some/script/I/want/to/load.R")

## ---- echo = FALSE, eval = FALSE--------------------------------------------------------------------------------------
#  library(DiagrammeR)
#  
#  grViz('
#  digraph AHP {
#          rankdir=LR;
#  	node [shape = box, style = filled, color = cornflowerblue, fontname="helvetica"];
#          fun1 [label="function(a1, a2)\nin ahp file"];
#          fun2 [label="function(a) in\nahp file", shape = none];
#          fun3 [label="scoreFun param\nof Calculate()", shape = none];
#          fun4 [label="pairwiseFun param\nof Calculate()", shape = none];
#          node [shape = none, style = empty];
#  
#    pairwiseFunction -> fun1 [arrowhead=none]
#    fun1 -> pairwise
#    pairwise -> fun4 [arrowhead=none]
#    fun4 -> priority
#    priority -> weightContribution
#    scoreFunction -> fun2 [arrowhead=none]
#    fun2 -> score
#    score -> fun3 [arrowhead=none]
#    fun3 -> priority
#  }
#  
#  ')
#  

## ---- comment = NA----------------------------------------------------------------------------------------------------
ahpFile <- system.file("extdata", "vacation.ahp", package="ahp")
cat(readChar(ahpFile, file.info(ahpFile)$size))

