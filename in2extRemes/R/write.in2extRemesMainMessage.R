write.in2extRemesMainMessage <- function(txt) {

    tkconfigure( txt, state="normal")
    nl <- paste("\n", "\n", sep="")

    msg01 <- paste("Into the extRemes Package (in2extRemes): ",
			"Weather and Climate Applications of Extreme-Value Statistics", " ", " ", sep="\n")

    msg02 <- paste("Type \'help(extRemes)\' for more information, and",
			"to get started, please see the tutorial at:", " ", sep="\n")

    msg03 <- paste("http://www.assessment.ucar.edu/toolkit/", sep="")

    msg04 <- paste("To see how to cite this GUI package use:")

    msg05 <- paste("To see how to cite extRemes use:")

    tkinsert( txt, "end", nl)
    tkinsert( txt, "end", msg01)
    tkinsert( txt, "end", nl)
    tkinsert( txt, "end", msg02)
    tkinsert( txt, "end", nl)
    tkinsert( txt, "end", msg03)
    tkinsert( txt, "end", nl)
    tkinsert(txt, "end", msg04)
    tkinsert( txt, "end", nl)
    tkinsert(txt, "end", "citation(\"in2extRemes\")")
    tkinsert(txt, "end", nl)
    # tkinsert(txt, "end", print(citation("in2extRemes"), style="text"))
    tkinsert( txt, "end", nl)
    tkinsert(txt, "end", msg05)
    tkinsert( txt, "end", nl)
    tkinsert(txt, "end", "citation(\"extRemes\")\n\n")
    tkinsert( txt, "end", nl)
    tkconfigure(txt, state = "disabled")

invisible()
}
