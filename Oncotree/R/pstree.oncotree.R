"pstree.oncotree" <-
function(x, edge.weights=c("none","observed","estimated"), edge.digits=2,
         shape=c("none","oval", "circle", "triangle", "diamond"),
         pstree.options=list(arrows="->", treefit="loose", arrowscale="1.5 0.8", nodesep="3pt")) {
  edge.weights <- match.arg(edge.weights)
  if ((edge.weights=="estimated") & is.null(x$parent$est.weight))
  	stop("No estimated weights are available in the tree.")
 	
 	shape <- match.arg(shape)
  
	put.vertex <- function(v.num, PA, labs){
	   if (is.null(labs)){  #no edge weights
	      res <- paste("\\pstree{%\n\\lab{", PA$child[v.num], "}}{", sep="")
	      }
     else {
       res <- paste("\\pstree{%\n\\lab{", PA$child[v.num], "}\\ncput*{",
			              labs[v.num],"}}{", sep="")
       }
	   children <- which(PA$parent.num == v.num)
	   for (ch in children){
     	  res <- paste(res, put.vertex(ch, PA, labs), sep="")
 	   }
		 res <- paste(res, "}", sep="")
		 res
  }
 
  labs <- switch(edge.weights,
	         none = NULL,
	         observed = formatC(x$parent$obs.weight, edge.digits, flag="#"),
	         estimated = formatC(x$parent$est.weight, edge.digits, flag="#"))
  node <- switch(shape,
           none = "TR",
					 oval = "Toval",
					 circle = "Tcircle",
					 triangle = "Ttri",
					 diamond = "Tdia")
           
  
  pstr <- put.vertex(1, x$parent, labs)
  cat("\\providecommand{\\lab}[1]{\\", node, "[name=#1]{$#1$}}\n", sep="")
	#cat("\\renewcommand{\\psedge}{\\ncline{",arrows,"}}\n", sep="")
	cat("\\psset{")
	 cat(paste(names(pstree.options), pstree.options, sep="="), sep=", ")
	 cat("}\n")
  cat(pstr, "\n")
  
 }
 
    