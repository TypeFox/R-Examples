s.plot.convex.hull.titles <-
function# outside titles Convex Hull matrix plot
### internal function for sisus
(n.sources
### internal variable
, n.isotopes
### internal variable
, names.isotopes
### internal variable
, analysis.name
### internal variable
, title.mixture
### internal variable
, tit
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # box("outer", col="black")  # draw a box around the entire figure
  mtext(paste("Isotopic Mixing Convex Hull"), side=1, line=1, outer=TRUE);                                                # bottom
  mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=3, outer=TRUE);                              # bottom (line specifies)
  mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                                                        # left
  mtext(paste(paste(title.mixture, ": ", tit)), cex = 1.2, side=3, line=1.5, outer=TRUE);                                   # top
  mtext(paste(analysis.name), cex = 0.8, side=3, line=0.25, outer=TRUE);                                                   # top
  mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=1, outer=TRUE);  # right
  # mtext(expression(paste("Isotopes: ", delta^{13}, "C")), side=4, line=-1, outer=TRUE);  # right
      # having trouble getting the n.isotopes in the front of the line, also trouble feeding the delta^13 as a variable.

  ### internal variable
}
