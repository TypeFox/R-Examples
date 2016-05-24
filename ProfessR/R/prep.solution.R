`prep.solution` <-
function(ofile)
  {


 if(file.exists(ofile)) { file.remove(ofile) }
   ##   system(paste("rm ", ofile))
    
    cat(file=ofile, "\\documentclass[10pt]{article}", sep="\n", append = TRUE)
    cat(file=ofile, "\\usepackage{amsmath} %Never write a paper without using amsmath for its many new commands", sep="\n", append = TRUE)
    cat(file=ofile, "\\usepackage{amssymb} %Some extra symbols", sep="\n", append = TRUE)
    cat(file=ofile, "\\usepackage{makeidx} %If you want to generate an index, automatically", sep="\n", append = TRUE)
    cat(file=ofile, "\\usepackage{graphicx} %If you want to include postscript graphics", sep="\n", append = TRUE)
    cat(file=ofile, "%%%  \\usepackage{mystyle} %Create your own file, mystyle.sty where you put all your own \\newcommand statements", sep="\n", append = TRUE)
    cat(file=ofile, "\\usepackage{amscd}", sep="\n", append = TRUE)
    ##    cat(file=ofile, "\\usepackage{/home/lees/Class/TESTBANK/UOFTEXAM}", sep="\n", append = TRUE)
    cat(file=ofile, "\\newcommand{\\mydegree}{ \\ensuremath{^\\circ} }", sep="\n", append = TRUE)
    cat(file=ofile, "\\newcommand{\\Rivpt}{\\rule{.1pt}{1pt}}", sep="\n", append = TRUE)
    cat(file=ofile, "\\clubpenalty10000", sep="\n", append = TRUE)
    cat(file=ofile, "\\widowpenalty10000", sep="\n", append = TRUE)
    cat(file=ofile, "\\raggedbottom", sep="\n", append = TRUE)
    cat(file=ofile, "\\addtolength{\\topskip}{0pt plus 10pt}", sep="\n", append = TRUE)

    cat(file=ofile, "\\let\\oldsubsubsection=\\subsubsection", sep="\n", append = TRUE)
    cat(file=ofile, "\\renewcommand{\\subsubsection}{%", sep="\n", append = TRUE)
    cat(file=ofile, "\\filbreak", sep="\n", append = TRUE)
    cat(file=ofile, "\\oldsubsubsection", sep="\n", append = TRUE)
    cat(file=ofile, "}", sep="\n", append = TRUE)

    cat(file=ofile, "\\setlength\\topmargin{0in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\headheight{0in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\headsep{0in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\textheight{9in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\textwidth{7in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\oddsidemargin{-.3in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\evensidemargin{-.3in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\parindent{0.25in}", sep="\n", append = TRUE)
    cat(file=ofile, "\\setlength\\parskip{0.25in}", sep="\n", append = TRUE)

    cat(file=ofile, "\\begin{document}", sep="\n", append = TRUE)
    cat(file=ofile, "\\begin{enumerate}", sep="\n", append = TRUE)



  }

