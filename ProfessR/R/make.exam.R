`make.exam` <-
function(Qbank, ofile="examq.tex", ncol=2)
{
  if(missing(ofile)) { ofile = "examq.tex" }
     if(missing(ncol))  ncol =  2


  ansfile = paste(sep=".", ofile, "ANS")

  if(file.exists(ofile)) { file.remove(ofile) }
  if(file.exists(ansfile)) { file.remove(ansfile) }

  ##  system(paste(sep=" ", "rm", ofile))
  
  ## system(paste(sep=" ", "rm", ansfile))

  for(i in 1:length(Qbank))
    {
      z = Qbank[[i]]
      A1 = substring(z$A, 3, 10000)
      A2 = paste("\\item {", A1, "}")

      THEQ = paste(z$Q)

      cat(file=ofile, "\\item {", sep="\n", append = TRUE)
      ## cat(file=ofile, "\\begin{minipage}{1\\linewidth}", sep="\n", append = TRUE)

       cat(file=ofile, "\\setlength{\\itemsep}{0cm}", sep="\n", append = TRUE)
     cat(file=ofile, "\\setlength{\\parskip}{.2cm}", sep="\n", append = TRUE)
   
      cat(file=ofile, "\\begin{samepage}", sep="\n", append = TRUE)

      ###########   questions are in Bold Font
      cat(file=ofile,  "\\textbf{", sep="\n", append = TRUE)
      cat(file=ofile,  z$Q , sep="\n", append = TRUE)
      cat(file=ofile,  "}", sep="\n", append = TRUE)
      cat(file=ofile, "\\begin{enumerate}", sep="\n", append = TRUE)
        cat(file=ofile, A2 , sep="\n", append = TRUE)

      cat(file=ofile, "\\end{enumerate}", sep="\n", append = TRUE)


      if(!is.null(z$FIG))
         {
           cat(file=ofile, "\\begin{figure}[htp]", sep="\n", append = TRUE)
           if(ncol==2)
             {
               colwidth = 0.5

             }
           else
             {

               colwidth = 0.85
             }

           
           cat(file=ofile, "\\centering", sep="\n", append = TRUE)
           ofig = paste(sep="", "\\includegraphics[width=",colwidth,"\\textwidth]{", z$FIG$fn, "}")
           cat(file=ofile, ofig, sep="\n", append = TRUE)
           labfig = paste(sep="", "\\caption{}\\label{",z$FIG$tag ,"}")
           cat(file=ofile,labfig, sep="\n", append = TRUE)
           cat(file=ofile, "\\end{figure}", sep="\n", append = TRUE)
         }
      ##   cat(file=ofile, "}", sep="\n", append = TRUE)
      ##  cat(file=ofile, "\\end{minipage}", sep="\n", append = TRUE)
      cat(file=ofile, "\\end{samepage}", sep="\n", append = TRUE)

      cat(file=ofile, "}", sep="\n", append = TRUE)
      cat(file=ansfile, paste(i, z$a) , sep="\n", append = TRUE)
    }
  
}

