latex.cover <-
function(y, label="", caption="", scale=c(1, 1), ...) {
  
  vertices <- -vertices(y)
  lab <- rownames(y)
  
  head <-  "\\documentclass{article}\n\\usepackage{tikz}\n\\tikzstyle{node}=[circle, draw, fill=white, inner sep=1pt, minimum width=1pt]\n\\begin{document}\n"
  
  begin <- paste("\\begin{figure}[!h]\n\\label{", label,"}\n\\centering\n\\begin{tikzpicture}[scale=1.2]\n", sep="")
  
  n <- nrow(y)
  nodes <- rep("", n)
  for(i in 1:n) {
    nodes[i] <- paste("\\node(", i, ") at (", vertices$x[i]*scale[1], ",", vertices$y[i]*scale[2], ")[node]{", lab[i],"};\n", sep="")
  }
  
  k <- sum(y)/2
  lines <- rep("", k)
  count <- 1
  for(a in 1:(n-1)) for(b in (a+1):n) if(y[a, b]) {
    lines[count] <- paste("\\draw[-, very thin] (", a,") to (", b,");\n", sep="")
    count <- count + 1
  }
  
  end <- paste("\\end{tikzpicture}\n\\caption{", caption, "}\n\\end{figure}\n\\end{document}", sep="")
  return(cat(paste(head, begin, paste(nodes, collapse="\n"), paste(lines, collapse="\n"), end, sep="\n")))
}
