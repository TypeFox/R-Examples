#
# Genetic code table as in Text Books
#

tablecode <- function(numcode = 1, urn.rna = s2c("TCAG"), dia = FALSE,
latexfile = NULL, label = latexfile, size = "normalsize", caption = NULL,
preaa = rep("", 64), postaa = rep("", 64), 
precodon = preaa, postcodon = postaa)
{
  aa1 <- a()
  aa3 <- aaa()
  codename <- SEQINR.UTIL$CODES.NCBI[numcode, "ORGANISMES"]
  urn <- s2c("tcag") # internal
#
# Make default caption for LaTeX table:
#
  if( is.null(caption) ){
    caption <- paste("\\caption{Genetic code number ", 
                     numcode, ": ", codename, ".}", sep = "")
  }
#
# As a LaTex table:
#  
  if( ! is.null(latexfile) ) {
    Tfile <- file(latexfile, open = "w")
    cat("\\begin{table}", file = Tfile, sep = "
")
    cat("\\begin{center}", file = Tfile, sep = "
")
#
# Character size:
#
    cat(paste("{\\", size, sep = ""), file = Tfile, sep = "
")
    
    cat("\\begin{tabular}{*{13}{l}}", file = Tfile, sep = "
")
    cat("\\hline", file = Tfile, sep = "
")
    cat("\\\\", file = Tfile, sep = "
")
    
    ncodon <- 1 # codon rank as in
    #paste(paste(rep(s2c("tcag"), each = 16), s2c("tcag"), sep = ""), rep(s2c("tcag"), each = 4), sep = "")
    for( i in 0:3 )
    {
      for( j in 0:3 )
      {
        for( k in 0:3 )
        {
          codon <- c(urn[i+1], urn[k+1], urn[j+1])
          codon.urn <- paste(urn.rna[i+1], urn.rna[k+1], urn.rna[j+1], sep = "", collapse = "")
          codon.urn <- paste(precodon[ncodon], codon.urn, postcodon[ncodon], sep = "")

          aminoacid <- aa3[which(aa1 == translate(codon, numcode = numcode))]
          aminoacid <- paste(preaa[ncodon], aminoacid, postaa[ncodon], sep = "")

          cat(paste(codon.urn, aminoacid, " &", sep = " & "), file = Tfile, sep = "
")
          ncodon <- ncodon + 1
        }
        cat("\\\\", file = Tfile, sep = "
")
      }
      cat("\\\\", file = Tfile, sep = "
")
    }
    cat("\\hline", file = Tfile, sep = "
")
    cat("\\end{tabular}", file = Tfile, sep = "
")
#
# Caption:
#
    cat(caption, file = Tfile, sep = "
")
#
# LaTeX label:
#
    cat(paste("\\label{", label, "}", sep = ""), file = Tfile, sep = "
")
#
# End character size:
#
    cat("}", file = Tfile, sep = "
")

    cat("\\end{center}", file = Tfile, sep = "
")
    cat("\\end{table}", file = Tfile, sep = "
")
    close(Tfile)
    return(invisible(NULL))
  }
#
# END LATEX
#
  if( dia )
  {  
    op <- par(no.readonly = TRUE)
    par(bg = "blue")
    par(fg = "yellow")
    par(col = "yellow")
    par(col.axis = "yellow")
    par(col.lab = "yellow")
    par(col.main = "yellow")
    par(col.sub = "yellow")
  }

  plot.new()
  plot.window(xlim=c(0,100),ylim=c(0,100))

  segments( 0, 102, 100, 102, lwd = 2)
  segments( 0, 0, 100, 0, lwd = 2)
  segments( 0, 97, 100, 97)


  text(x=0, y = 98.5, font = 2, adj = c(0, 0),
    lab = paste("Genetic code", numcode,":",codename))

  urn <- c("t","c","a","g") # internal
  for( i in 0:3 )
  {
    for( j in 0:3 )
    {
      for( k in 0:3 )
      {
        codon <- c(urn[i+1], urn[j+1], urn[k+1])

        text( x = 100*j/4, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab = urn.rna[i+1] )

        text( x = 100*j/4 + 3, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab = urn.rna[j+1] )

        text( x = 100*j/4 + 6, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab = urn.rna[k+1] )

        aminoacid <- aa3[which(aa1 == translate(codon, numcode = numcode))]
        text( x = 100*j/4 + 12, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab =  aminoacid )
      }
    }
  }
  if(dia)
    par(op)
}
