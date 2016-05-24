test.LaTeX <- function() {
  #TODO Create a function getAllGraphs(variables=TRUE) for all these occcurences.
	graphs <- list(BonferroniHolm(5),
			parallelGatekeeping(),
			improvedParallelGatekeeping(),
			BretzEtAl2011(),
			HungEtWang2010(),
			HuqueAloshEtBhore2011(),
			HommelEtAl2007(),
			HommelEtAl2007Simple(),
			MaurerEtAl1995(),
			improvedFallbackI(weights=rep(1/3, 3)),
			improvedFallbackII(weights=rep(1/3, 3)),
			cycleGraph(nodes=paste("H",1:4,sep=""), weights=rep(1/4, 4)),
			fixedSequence(5),
			fallback(weights=rep(1/4, 4)),
			generalSuccessive(weights = c(1/2, 1/2)),
			simpleSuccessiveI(),
			simpleSuccessiveII(),
			truncatedHolm(),
			BauerEtAl2001(),
			BretzEtAl2009a(),
			BretzEtAl2009b(),
			BretzEtAl2009c(),
			Ferber2011())
	
	report <- gMCP:::LaTeXHeader()
	for (graph in graphs) {
		report <- paste(report, graph2latex(graph), sep="\n")
	}
	report <- paste(report, "\\end{document}", sep="\n")
	if (interactive() && gMCP:::tests("interactive")) {
		file <- paste(Sys.getenv("GMCP_UNIT_TEST_OPATH"), "report_test.tex", sep="/")
		cat(report, file=file)
    wd <- getwd()
    setwd(Sys.getenv("GMCP_UNIT_TEST_OPATH"))
		exitValue <- system(paste("pdflatex", file))
    setwd(wd)
		if (exitValue != 0) {
			stop("Error calling pdflatex.")
		}
		answer <- readline(paste("Does report_test.pdf look fine (Y/n)? "))
		if (substr(answer, 1, 1) %in% c("n","N")) {
			stop("User reported error for report_test.pdf.")
		}
	} else {
		cat("Skipping interactive control of LaTeX output.\n")
	}
}

test.fractionStrings <- function() {
	checkEquals(gMCP:::getLaTeXFraction(1/9), "\\frac{1}{9}")
	checkTrue(gMCP:::getFractionString(1/9+0.000001) !="1/9")
}