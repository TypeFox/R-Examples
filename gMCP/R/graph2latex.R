#' Graph2LaTeX
#' 
#' Creates LaTeX code that represents the given graph.
#' 
#' For details see the given references.
#' 
#' @param graph A graph of class \code{\link{graphMCP}}.
#' @param package A character string specifying the LaTeX package that should
#' be used.  Up to now only \code{TikZ} is available.
#' @param scale A numeric scalar specifying a possible scaling of the graph.
#' It is only used if \code{tikzEnv==TRUE}.
#' Note that this does only effect the fontsize of the graph if \code{scaleText==FALSE}.
#' (Coordinates are interpreted in big points: 72 bp = 1 inch).
#' @param showAlpha Logical whether local alpha levels or weights should be shown.
#' @param alpha An optional numeric argument to specify the type I error rate.
#' @param pvalues If the optional numeric argument pvalues is given, nodes that
#' can be rejected, will be marked.
#' @param fontsize An optional character vector specifying the fontsize for the
#' graph, must be one of \code{"tiny"}, \code{"scriptsize"},
#' \code{"footnotesize"}, \code{"small"}, \code{"normalsize"}, \code{"large"},
#' \code{"Large"}, \code{"LARGE"}, \code{"huge"} or \code{"Huge"}.
#' @param nodeTikZ A character string with additional arguments for the TikZ
#' \code{node} command like for example \code{nodeTikZ="minimum size=2cm"}.
#' @param labelTikZ A character string with arguments for the TikZ \code{node}
#' command within an edge.
#' @param tikzEnv Logical whether the LaTeX code should be wrapped in a TikZ
#' environment.
#' @param offset A numeric of length 2 specifying the x and y offset in the
#' TikZ environment.
#' @param fill A list containing 2 elements \code{reject} and \code{retain}
#' specifying node fill colour of rejected and retained (or not yet rejected)
#' nodes.
#' @param fig Logical whether a figure environment should be created.
#' @param fig.label Label for figure environment (if \code{fig==TRUE}).
#' @param fig.caption Caption for figure environment (if \code{fig==TRUE}).
#' @param fig.caption.short Optional short version of fig.caption for list of figures (if \code{fig==TRUE}).
#' @param nodeR Radius of nodes (pixel in Java, bp in LaTeX).
#' @param scaleText Only used if scale is unequal 1 and \code{tikzEnv==TRUE}. 
#' If \code{scaleText} is \code{TRUE} (the default) a scalebox environment is used.
#' If it is \code{FALSE} the optional parameter \code{scale} from the
#' tikzpicture environment is used and font size will not change. 
#' Note that while you easily can change the scale in the scalebox environment,
#' it is more problematic to adjust the scale in the tikzpicture environment
#' afterwards in the LaTeX document, since for curved edges the parameters
#' are calculated for a certain relative node size which changes if the graph
#' is scaled but the text size stays the same.
#' @return A character string that contains LaTeX code representing the given
#' graph.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{graphMCP}}, \code{\link{gMCPReport}}
#' @references The TikZ and PGF Packages Manual for version 2.00, Till Tantau,
#' \url{https://www.ctan.org/pkg/pgf/}
#' @keywords print graphs
#' @examples
#' 
#' 
#' g <- BonferroniHolm(5)
#' 
#' graph2latex(g)
#' 
#' 
#' @export graph2latex
graph2latex <- function(graph, package="TikZ", scale=1, showAlpha=FALSE, alpha=0.05, pvalues,
		fontsize,	nodeTikZ, labelTikZ="near start,above,fill=blue!20",
		tikzEnv=TRUE, offset=c(0,0),fill=list(reject="red!80",retain="green!80"),
		fig=FALSE, fig.label=NULL, fig.caption=NULL, fig.caption.short=NULL, nodeR=25, scaleText=TRUE) {
	graph <- placeNodes(graph)
	colors <- c("yellow","black","blue","red","green")
	if (tikzEnv) {        
		tikz <- paste("\\begin{tikzpicture}",
                  ifelse(scaleText, 
                         "",
                         paste("[scale=",scale,"]", sep=""))
                  ,"\n", sep="")    
    if (!scaleText) nodeR <- nodeR / scale
	} else {
		tikz <- ""
	}
	nodes2 <- getUsableNames(getNodes(graph))
	names(nodes2) <- getNodes(graph)
	#tikz <- paste(tikz, "\\tikzset{help lines/.style=very thin}", paste="\n")	
	for (node in getNodes(graph)) {
		nodeColor <- ifelse(getRejected(graph, node),fill$reject, fill$retain)
		x <- getXCoordinates(graph, node) + nodeR # Actually I'm not sure whether this "+ nodeR" in this function is really necessary.
		y <- getYCoordinates(graph, node) + nodeR # It is more important to use nodeR to reduce the arc that should be drawn.
		#alpha <- format(getWeights(graph,node), digits=3, drop0trailing=TRUE)
		weight <- paste(getLaTeXFraction(getWeights(graph,node)), collapse=" ")
		if (showAlpha) {
		  if (weight == 1) {
		    weight <- "\\alpha"
		  } else if (weight != "0") {
		    weight <- paste(weight, "\\alpha", sep="")
		  }
		}
		double <- ""
		if (!missing(pvalues)) {
			if (is.null(names(pvalues))) {
				names(pvalues) <- getNodes(graph)
			}
			if (canBeRejected(graph, node, alpha, pvalues)) { double <- "double," }
		}		
		nodeLine <- paste("\\node (",nodes2[node],")",
				" at (",x+offset[1],"bp,",-y-offset[2],"bp)",
				"[draw,circle split,",ifelse(missing(nodeTikZ),"",paste(nodeTikZ,", ",sep="")),double,"fill=",nodeColor,"]",
				" {$",node,"$ \\nodepart{lower} $",weight,"$};",sep="")
		tikz <- paste(tikz, nodeLine,sep="\n")			
	}
	# A second loop for the edges is necessary:
	if ("entangledMCP" %in% class(graph)) {
		for(k in 1:length(graph@subgraphs)) {
			subgraph <- graph@subgraphs[[k]]
			for (i in getNodes(subgraph)) {
				for (j in getNodes(subgraph)) {			
					if (subgraph@m[i,j]!=0) {
						# The following to lines test whether the edge in opposite direction exists:				
						to <- ifelse(subgraph@m[j,i]==0, "auto", "bend left=15")
						#weight <- ifelse(edgeL[i]==0, "\\epsilon", getLaTeXFraction(edgeL[i])) # format(edgeL[i], digits=3, drop0trailing=TRUE))
						weight <- getWeightStr(subgraph, i, j, LaTeX=TRUE) 
						edgeLine <- paste("\\draw [draw=",colors[k%%length(colors)+1],",->,line width=1pt] (",nodes2[i],") to[",to,"] node[",labelTikZ,"] {$",weight,"$} (",nodes2[j],");",sep="")
						tikz <- paste(tikz, edgeLine,sep="\n")
					}
				}
			}
		}
	} else {
		for (i in getNodes(graph)) {
			for (j in getNodes(graph)) {			
				if (graph@m[i,j]!=0) {
				  
				  weight <- getWeightStr(graph, i, j, LaTeX=TRUE) 
				  edgeNode <- paste("node[",labelTikZ,"] {$",weight,"$}",sep="")
					# The following line tests whether the edge in opposite direction exists:				          
					to <- paste(") to[",ifelse(graph@m[j,i]==0, "auto", "bend left=15"),"] ", edgeNode, sep="")
          
          edgeNode <- paste("node[","fill=blue!20","] {$",weight,"$}",sep="") # TODO labelTikZ is ignored in this case
          # New arc function:
					x <- try(unlist(edgeAttr(graph, i, j, "labelX")), silent = TRUE)          
					y <- try(unlist(edgeAttr(graph, i, j, "labelY")), silent = TRUE)
					if (!("try-error" %in% class(x)) && !is.null(x) && !is.na(x) && class(y)!="try-error" && !is.null(y) && !is.na(y) && x>-10 && y>-10) {
					  b <- c(x+offset[1],y-offset[2]) + nodeR
					  x <- getXCoordinates(graph, c(i,j)) + nodeR + offset[1]
					  y <- getYCoordinates(graph, c(i,j)) + nodeR - offset[2]
					  to2 <- try(getArc(c(x[1],y[1]),b,c(x[2],y[2]), edgeNode, nodeR=nodeR))
					  if (!("try-error" %in% class(to2))) to <- to2
					}          
					#weight <- ifelse(edgeL[i]==0, "\\epsilon", getLaTeXFraction(edgeL[i])) # format(edgeL[i], digits=3, drop0trailing=TRUE))					
					edgeLine <- paste("\\draw [->,line width=1pt] (",nodes2[i], to," (",nodes2[j],");",sep="")
					tikz <- paste(tikz, edgeLine,sep="\n")
				}
			}
		}
	}
	if (tikzEnv) tikz <- paste(tikz, "\\end{tikzpicture}\n",sep="\n")
	if (!missing(fontsize)) {
		tikz <- paste(paste("{\\", fontsize, sep=""), tikz, "}",sep="\n")
	}
  if ( !isTRUE(all.equal(scale,1, check.attributes=FALSE, check.names=FALSE)) && scaleText && tikzEnv) {
    tikz <- paste(paste("\\scalebox{",scale,"}{",sep=""),tikz,"}",sep="\n")
  }
  if (fig) {
    label <- " "
    if (!is.null(fig.label)) {
      label <- paste("\\label{", fig.label, "} ", sep="")
    }
    caption <- ""
    if (!is.null(fig.caption)) {
      short.caption <- ""
      if (!is.null(fig.caption.short)) {
        short.caption <- paste("[",fig.caption.short,"]", sep="")
      }
      caption <- paste("\\caption",short.caption,"{",label,fig.caption,"}", sep="")
    }
    tikz <- paste("\\begin{figure}[ht]\n\\begin{center}", tikz,
                  "\\end{center}", caption, "\\end{figure}", sep="\n")
    
  }
	return(tikz)
}

# Arc from a to b and from b to c.
getArc <- function(a, b, c, edgeNode, col="black", nodeR=25) {
  #a <- invertY(a)
  #b <- invertY(b)
  #c <- invertY(c)
  m <- try(getCenter(a,b,c,0.001), silent=TRUE)
  if ("try-error" %in% class(m)) {
    return(getLine(a, b, c, edgeNode, col="black"))    
  }
  r <- sqrt(sum((m-a)^2))
  if (r>500) return(getLine(a, b, c, edgeNode, col="black"))
  phi <- getAngle(a,b,c,m, nodeR)  
  #cat("a: ",a,", b: ",b,", c:", c,"m: ",m,"r: ",r,", phi: ",phi,"\n")
  return(paste(".",round(phi[1]+ifelse(phi[1]>phi[2],-90,90)),") arc(",round(phi[1]),":",round(phi[3]),":",round(r),"bp) ",edgeNode," arc(",round(phi[3]),":",round(phi[2]),":",round(r),"bp) to",sep=""))
}

# Line from a to b and from b to c.
getLine <- function(a, b, c, edgeNode, col="black") {
  return(paste(") to (",round(b[1]),"bp, ",-round(b[2]),"bp) ",edgeNode," to",sep=""))
}

invertY <- function(x) {
  return(c(x[1],-x[2]))
}

getAngle <- function(a,b,c,m, nodeR=25) {
  # phi correction factor:
  r <- sqrt((m[1]-a[1])*(m[1]-a[1])+(m[2]-a[2])*(m[2]-a[2]))
  #phiCF <- (nodeR*360)/(2*pi*r)
  phiCF <- 2*asin(nodeR/(2*r))/(2*pi)*360
  
  if ((a[1]-m[1])==0) {
    phi1 <- 90 + ifelse((m[2]-a[2]>0),0,180)
  } else {
    phi1 <- atan((-a[2]+m[2])/(a[1]-m[1]))*360/(2*pi)+ifelse((a[1]-m[1]<0),180,0)
  }
  if ((c[1]-m[1])==0) {
    phi2 <- 90 + ifelse((m[2]-c[2]>0),0,180)
  } else {
    phi2 <- atan((-c[2]+m[2])/(c[1]-m[1]))*360/(2*pi)+ifelse((c[1]-m[1]<0),180,0)
  }
  if ((b[1]-m[1])==0) {
    phi3 <- 90 + ifelse((m[2]-b[2]>0),0,180)
  } else {
    phi3 <- atan((-b[2]+m[2])/(b[1]-m[1]))*360/(2*pi)+ifelse((b[1]-m[1]<0),180,0)
  }		
  phi1 <- (phi1 + 360) %% 360 # phi for a
  phi2 <- (phi2 + 360) %% 360 # phi for c
  phi3 <- (phi3 + 360) %% 360 # phi for b
  #return(c(phi1, phi2, phi3))
  if ((phi1 > phi2 && phi1 > phi3 && phi3 > phi2) || (phi2 > phi1 && (phi3>phi2 || phi3<phi1))) {  
    # Clockwise direction: phi2 < phi1
    phi1 <- phi1 - phiCF
    phi2 <- phi2 + phiCF + 2
    if (phi3>phi1) phi3 <- phi3 -360
    if (phi2<phi1){
      return(c(phi1, phi2, phi3))
    } else {
      return(c(phi1, phi2-360, phi3))
    }
  } else {
    # Counter clockwise: phi1 < phi2
    phi1 <- phi1 + phiCF
    phi2 <- phi2 - phiCF - 2
    if (phi3<phi1) phi3 <- phi3 + 360
    if (phi1<phi2) {      
      return(c(phi1, phi2, phi3))
    } else {
      return(c(phi1, phi2+360, phi3))
    }
  }
}

getCenter <- function(a,b,c, eps=0.05) {  
  if((b[2]-c[2])==0) {
    x <- c(0,1)
  } else {
    x <- c(1,-(b[1]-c[1])/(b[2]-c[2]))
  }
  if ((a[2]-b[2])==0) {
    z <- c(0,1)
  } else {
    z <- c(1,-(a[1]-b[1])/(a[2]-b[2]))
  }
  if (abs((b[1]-a[1])/(b[2]-a[2])-(c[1]-b[1])/(c[2]-b[2]))<eps && sign(b[1]-a[1])==sign(c[1]-b[1])) {
    stop("Slopes are to similar")
  }
  if (z[1]!=0 && x[1]==0) {			
    c <- (c[1]-a[1])/(2*z[1])
    return(c((a[1]+b[1])/2+c*z[1], (a[2]+b[2])/2+c*z[2]))
  } else if (x[1]!=0 && z[1]==0) {
    d <- (a[1]-c[1])/(2*x[1])
    return(c((b[1]+c[1])/2+d*x[1], (b[2]+c[2])/2+d*x[2]))
  } else if ((x[1]==0 && z[1]==0)||(x[2]==0 && z[2]==0)) {
    stop("Slopes are too similar.")
  } else if (z[2]!=0 && x[2]==0) {			
    c <- (c[2]-a[2])/(2*z[2])
    return(c((a[1]+b[1])/2+c*z[1], (a[2]+b[2])/2+c*z[2]))
  } else if (x[2]!=0 && z[2]==0) {			
    d <- (a[2]-c[2])/(2*x[2])
    return(c((b[1]+c[1])/2+d*x[1], (b[2]+c[2])/2+d*x[2]))
  } else {
    if ((x[2]-x[1]*z[2]/z[1])==0) {
      if ((z[2]-z[1]*x[2]/x[1])==0) stop("Can this happen?")
      c <- ((c[2]-a[2])/2+((a[1]-c[1])/2*x[1])*x[2])/(z[2]-z[1]*x[2]/x[1])
      return(c((a[1]+b[1])/2+c*z[1], (a[2]+b[2])/2+c*z[2]))			
    }
    d <- ((a[2]-c[2])/2+((c[1]-a[1])/2*z[1])*z[2])/(x[2]-x[1]*z[2]/z[1])		
  }  
  m <- c((b[1]+c[1])/2+d*x[1], (b[2]+c[2])/2+d*x[2])
  return(m)
}

# x <- c("H+1","H-1","H/1","H1","H2","H1+")
# getUsableNames(x)
# [1] "H1"  "H12" "H13" "H14" "H2"  "H15"
getUsableNames <- function(x) {
	x <- removeSymbols(x)
	for (i in which(duplicated(x))) {
		name <- x[i]
		j <- 2
		while (any(x==paste(x[i],j,sep=""))) {
			j <- j + 1
		}
		x[i] <- paste(x[i],j,sep="")
	}
	return(x)
}

getLaTeXFraction <- function(x) {
	result <- c()
	for (nom in strsplit(as.character(getFractionString(x)),split="/")) {		
		if (length(nom)==1) {
			result <- c(result, nom)
		} else {
			result <- c(result, paste("\\frac{",nom[1],"}{",nom[2],"}", sep=""))
		}
	}
	return(result)
}

#' Automatic Generation of gMCP Reports
#' 
#' Creates a LaTeX file with a gMCP Report.
#' 
#' This function uses \code{cat} and \code{graph2latex}.
#' 
#' @param object A graph of class \code{\link{graphMCP}} or an object of class
#' \code{\link{gMCPResult}}.
#' @param file A connection, or a character string naming the file to print to.
#' If \code{""} (the default), the report is printed to the standard output
#' connection (the console unless redirected by \code{sink}).  If it is
#' \code{"|cmd"}, the output is piped to the command given by \code{cmd}, by
#' opening a pipe connection [taken from the manual page of \code{cat}, which
#' is called in this function].
#' @param ...  Arguments to be passed to method \code{\link{graph2latex}} like
#' \code{package} and \code{scale}.
#' @return None (invisible \code{NULL}).
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{cat}} \code{\link{graph2latex}}
#' @references The TikZ and PGF Packages Manual for version 2.00, Till Tantau,
#' \url{https://www.ctan.org/pkg/pgf/}
#' @keywords print IO file graphs
#' @examples
#' 
#' g <- BretzEtAl2011()
#' 
#' result <- gMCP(g, pvalues=c(0.1, 0.008, 0.005, 0.15, 0.04, 0.006))
#' 
#' gMCPReport(result)
#' 
#' @export gMCPReport
#' 
gMCPReport <- function(object, file="", ...) {
	report <- LaTeXHeader()
	if (class(object)=="gMCPResult") {
		report <- paste(report, "\\subsection*{Initial graph}", sep="\n")
		report <- paste(report, graph2latex(object@graphs[[1]], ..., pvalues=object@pvalues), sep="\n")
		report <- paste(report, "\\subsection*{P-Values}", sep="\n")
		report <- paste(report, createTable(object@pvalues), sep="\n")	
		if (length(object@adjPValues)>0) {
			report <- paste(report, "\\subsection*{Adjusted p-values}", sep="\n")
			report <- paste(report, createTable(object@adjPValues), sep="\n")	
		}
		if (length(object@rejected)>0) {
			report <- paste(report, paste("\\subsection*{Rejected Hypotheses with $\\alpha=",object@alpha,"$}", sep=""), sep="\n")
			report <- paste(report, createTable(object@rejected), sep="\n")
		}
		if (length(object@graphs)>1) {
			for(i in 2:length(object@graphs)) {
				report <- paste(report, paste("\\subsection*{Graph in Step ",i,"}", sep=""), sep="\n")
				report <- paste(report, graph2latex(object@graphs[[i]], ..., pvalues=object@pvalues), sep="\n")
			}		
		}
	} else if (class(object)=="graphMCP") {		
		report <- paste(report, "\\subsection*{Graph for SRMTP}", sep="\n")
		report <- paste(report, graph2latex(object, ...), sep="\n")
	} else {
		stop("object has to be of class gMCPResult or graphMCP.")
	} 
	report <- paste(report, "\\end{document}", sep="\n")
	cat(report, file=file)
}

createTable <- function(vector) {
	table <- paste("\\begin{table}[ht]",
	"\\begin{center}", sep="\n")
	table <- paste(table, 
		"\n\\begin{tabular}{",paste(rep("r",length(vector)),collapse=""),"}\n",
		"\\hline\n", sep="")
    values <- paste(vector, collapse="&")
	if (is.numeric(vector)) values <- paste(sprintf("%.5f", vector), collapse="&")
	table <- paste(table, "\n", paste(names(vector), collapse="&"), " \\\\\n\\hline\n ", values, "\\\\\n\\hline\n ", sep="") 
	table <- paste(table, 
		"\\end{tabular}",
		"\\end{center}",
		"\\end{table}", sep="\n");
	return(table)
}

LaTeXHeader <- function() {
	report <- "\\documentclass[11pt]{article}"
	report <- paste(report, "\\usepackage{tikz}", sep="\n")
	report <- paste(report, "\\usetikzlibrary{decorations,arrows,shapes}", sep="\n")
	report <- paste(report, "\\begin{document}", sep="\n")
	return(report)
}
