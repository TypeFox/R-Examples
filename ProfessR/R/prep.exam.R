
prep.exam<-function(OF, incfile, instructor="", examdate=" ", course="",  examname="", instructions="", ncol=2)
  {

    if(missing(instructor)) instructor=""
    if(missing(examdate)) examdate=" "
    if(missing(course)) course=""
    if(missing(examname)) examname=""
    if(missing(instructions)) instructions=""
    if(missing(ncol))  ncol =  2


   
    
    docclass = "\\documentclass[12pt]{article}"
    
    if(ncol == 2 )
      {
        docclass = "\\documentclass[11pt,twocolumn]{article}"
      }

    ##################   head of latex document
    cat(file=OF, docclass, sep="\n", append=FALSE)

    cat(file=OF, "\\usepackage{amsmath} %use  amsmath ", sep="\n", append=TRUE)
    cat(file=OF, "\\usepackage{amssymb} %Some extra symbols", sep="\n", append=TRUE)
    cat(file=OF, "\\usepackage{makeidx} % to generate an index, automatically", sep="\n", append=TRUE)
    cat(file=OF, "\\usepackage{graphicx} % for postscript graphics", sep="\n", append=TRUE)
    cat(file=OF, "%%%  \\usepackage{mystyle} %Create your own file, mystyle.sty ", sep="\n", append=TRUE)
    cat(file=OF, "\\usepackage{amscd}", sep="\n", append=TRUE)

    cat(file=OF, "\\newcommand{\\mydegree}{ \\ensuremath{^\\circ} }", sep="\n", append=TRUE)
    cat(file=OF, "\\newcommand{\\Rivpt}{\\rule{.1pt}{1pt}}", sep="\n", append=TRUE)
    cat(file=OF, "\\clubpenalty10000", sep="\n", append=TRUE)
    cat(file=OF, "\\widowpenalty10000", sep="\n", append=TRUE)
    cat(file=OF, "\\raggedbottom", sep="\n", append=TRUE)
    cat(file=OF, "\\addtolength{\\topskip}{0pt plus 10pt}", sep="\n", append=TRUE)
    cat(file=OF, "\\let\\oldsubsubsection=\\subsubsection", sep="\n", append=TRUE)
    cat(file=OF, "\\renewcommand{\\subsubsection}{%", sep="\n", append=TRUE)
    cat(file=OF, "\\filbreak", sep="\n", append=TRUE)
    cat(file=OF, "\\oldsubsubsection", sep="\n", append=TRUE)
    cat(file=OF, "}", sep="\n", append=TRUE)



    cat(file=OF, "\\def\\instructions{\\begin{center}", sep="\n", append=TRUE)
    cat(file=OF, "{\\bf INTRUCTIONS\\vspace{-.5em}\\vspace{0pt}}\\end{center}}", sep="\n", append=TRUE)
   
    
    cat(file=OF, "\\def\\endinstructions{\\par}", sep="\n", append=TRUE)
    
    cat(file=OF, "\\def\\studentinfo{\\noindent\\hspace*{0pt}{\\bfseries Given name}:\\rule{46mm}{0.6pt}", sep="\n", append=TRUE)

    cat(file=OF, " \\par\\vspace{5mm}", sep="\n", append=TRUE)
    
    cat(file=OF, "\\hspace*{0mm}{\\bfseries Family name}:\\rule{43mm}{0.6pt}", sep="\n", append=TRUE)
    
    cat(file=OF, "\\hspace*{0pt}\\par\\vspace{5mm}\\noindent\\hspace*{0pt}{\\bfseries Student number}:\\rule{45mm}{0.6pt}", sep="\n", append=TRUE)
     cat(file=OF, " \\par\\vspace{5mm}", sep="\n", append=TRUE)
    
    cat(file=OF, "\\hspace*{0mm}{\\bfseries Signature}:\\rule{60mm}{0.6pt}\\hspace*{0pt}}", sep="\n", append=TRUE)

    cat(file=OF, "\\newcommand{\\WriteOnTest}{", sep="\n", append=TRUE)
    cat(file=OF, "\\studentinfo", sep="\n", append=TRUE)
    cat(file=OF, "\\vspace*{4mm}}", sep="\n", append=TRUE)

    
 ###   cat(file=OF, "\\textwidth 6.5in", sep="\n", append=TRUE)
  ###  cat(file=OF, "\\textheight 9in", sep="\n", append=TRUE)

############  % Move the % down one line, as above, if necessary
  ###  cat(file=OF, "\\topmargin 0in", sep="\n", append=TRUE)
####%\topmargin .5in



cat(file=OF, "\\setlength\\topmargin{0in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\headheight{0in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\headsep{0in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\textheight{9in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\textwidth{7in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\oddsidemargin{-.3in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\evensidemargin{-.3in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\parindent{0.25in}", sep="\n", append=TRUE)
cat(file=OF, "\\setlength\\parskip{0.25in}", sep="\n", append=TRUE)


    


    cat(file=OF, "\\begin{document}", sep="\n", append=TRUE)

    cat(file=OF, "\\author{", sep="", append=TRUE)
    cat(file=OF, instructor, sep="", append=TRUE)
    if(course!=""){
      cat(file=OF, "\\\\", sep="\n", append=TRUE)
      cat(file=OF, course, sep="", append=TRUE)
    }
    cat(file=OF, "}", sep="\n", append=TRUE)

    tit = paste(sep="", "\\title{", examname, "}")

    cat(file=OF, tit, sep="\n", append=TRUE)

    edate =  paste(sep="", "\\date{", examdate, "}")

    cat(file=OF, edate, sep="\n", append=TRUE)

    cat(file=OF, "\\maketitle", sep="\n", append=TRUE)

    cat(file=OF, "\\WriteOnTest", sep="\n", append=TRUE)

    
    if(length(instructions)>0)
      {

         cat(file=OF, "{\\bfseries", sep="\n", append=TRUE)

        if(instructions[1]!=""){
          cat(file=OF, "\\begin{instructions}", sep="\n", append=TRUE)
          for(i in 1:length(instructions))
            {
              cat(file=OF,instructions[i], sep="\n", append=TRUE)
              
            }
          cat(file=OF, "\\end{instructions}", sep="\n", append=TRUE)
        }

             cat(file=OF, "}", sep="\n", append=TRUE)
      }
    
    cat(file=OF, "\\begin{enumerate}", sep="\n", append=TRUE)
    cat(file=OF, "\\setlength{\\parskip}{1cm}%", sep="\n", append=TRUE)

    zout = paste(sep="", "\\include{", incfile, "}")

    cat(file=OF, zout, sep="\n", append=TRUE)
    cat(file=OF, "\\end{enumerate}", sep="\n", append=TRUE)
    cat(file=OF, "\\end{document}", sep="\n", append=TRUE)


  }
