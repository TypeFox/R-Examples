MST.plot.latex <-
function(tree, file="tree-code.tex", digits=5){
  n.node <- nrow(tree)
  sink(file=file)
  # DEFINE SYMBOLS
  cat("\\begin{figure} \\centering \n")
  cat("\\newcommand{\\Tgcircle}{\\Tcircle[linecolor=black, fillcolor=gray]} \n")
  cat("\\newcommand{\\Tgoval}{\\Toval[linecolor=black, fillcolor=gray]} \n")
  cat("\\newcommand{\\Tgframe}[1]{\\Tr{\\psframebox[linecolor=black, fillstyle=solid, fillcolor=orange]{#1}}} \n")
  # OPTION
  cat("\\psset{nodesep=0.7pt, linecolor=black, treesep=1.2cm, levelsep=1.8cm} \n")
  I0 <- i0 <- NULL
  for (i in 1:n.node) {
    node.i <- as.character(tree[i, 1])
    de.i <- de(node.i, tree)
    blanks <- paste(rep(" ", (nchar(node.i)-1)*8), sep="")  # 8 SPACES IN ONE TAB
    n.i <- tree$size[i]
    mean.i <- ifelse(!is.null(tree$mean.y), tree$mean.y[i], " "); 	##### THIS MEASURE MAY BE DIFFERENT FOR OTHER TYPES OF TREES
    if (!is.na(de.i[1])) {	# INTERNAL NODE
      if (nchar(node.i)==1 ||  substr(node.i, nchar(node.i), nchar(node.i))=="2")
        cat(blanks, "\\pstree{\\Tgcircle{~~}} \n", sep = "")
      else cat(blanks, "\\pstree{\\Tgcircle{~~} \\tlput{\\color{blue}", rule.i, "\\hspace{-.6in}}} \n", sep = "")
      cat(blanks, "{ \n", sep = "")
      I0 <- c(I0, i)
      i0 <- c(i0, i + length(de.i))
      # UPDATE THE SPLITTING RULE
      vname.i <- tree$vname[i]; 
      op <- as.character(tree$operator[i])
      if (op == ">") {operator.i <- ">"}
      else if (op=="<=")  {operator.i <- "\\leq "}
      else if (op=="in") {operator.i <- "\\in "}
      else if (op == "not in") {operator.i <- " \\not \\in "}
      else {stop("Wrong operate!")}
      cut.i <- strtrim(as.character(tree$cut[i]), width=digits);
      if (op=="in"|| op=="not in")  rule.i <- paste("\\texttt{", vname.i, "} ", "$", operator.i, "\\{", cut.i, "\\}", "$", sep="")
      else rule.i <- paste("\\texttt{", vname.i, "} ", "$", operator.i, cut.i, "$", sep="")
    } else if (substr(node.i, nchar(node.i), nchar(node.i))=="1") { # TERMINAL NODE
      cat(blanks, "\\Tgframe{",  mean.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{", n.i, "}}",
          "\\tlput{\\color{blue} ", rule.i, " \\hspace{-.3in}} \n", sep = "")
    } else cat(blanks, "\\Tgframe{",  mean.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{",  n.i, "}} \n", sep = "")
    if (is.element(i, i0)) {
      rep0 <- rep("}", sum(i0==i))
      node.i0 <- as.character(tree[I0[i0==i][1] , 1])
      blanks0 <- paste(rep(" ", (nchar(node.i0)-1)*8), sep="")
      cat(blanks0, rep0, "\n", sep = "")
    }
  }
  cat("\\end{figure} \n")
  sink()
}
