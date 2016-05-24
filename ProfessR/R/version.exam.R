`version.exam` <-
function(Qbank, V, exnumber="Exam 1", seqnum="2", examdate='', instructor="", course="", instructions="", SAMP=TRUE , ncol=2)
{
  if(missing(exnumber)) {  exnumber="Exam 1" }
  if(missing(seqnum)) { seqnum="1" }
  if(missing(examdate)) {   examdate=''}
    if(missing(instructor)) instructor=""
   
    if(missing(course)) course=""
    
    if(missing(instructions)) instructions=""
    if(missing(SAMP)) SAMP=TRUE
    if(missing(ncol))  ncol =  2

  


if(SAMP==TRUE)   {  QTEMP = ran.exam(Qbank)  } else { QTEMP = Qbank ; attr(QTEMP, "QBord")<-c(1:length(Qbank))   }


  QBord  = attr(QTEMP, "QBord")

    
  outtex = paste(sep=".",V, "tex" )
  outMAST  = paste(sep="", V, "MAST" )
  outSOLTN = paste(sep="", V, "SOLTN" )


  examname=paste(sep=" ", exnumber, "Seq", seqnum)



  
  MASTtex  = paste(sep=".", outMAST , "tex" )


  
  MASTdvi  = paste(sep=".", outMAST , "dvi" )
  MASTps  = paste(sep=".", outMAST , "ps" )
  MASTpdf  = paste(sep=".", outMAST , "pdf" )


SOLTNtex  = paste(sep=".", outSOLTN , "tex" )

  
  SOLTNdvi  = paste(sep=".", outSOLTN , "dvi" )
  SOLTNps  = paste(sep=".", outSOLTN , "ps" )
  SOLTNpdf  = paste(sep=".", outSOLTN , "pdf" )


  
  make.exam(QTEMP, ofile=outtex, ncol =ncol )

  
  make.solution(QTEMP, ofile=outSOLTN)


  prep.exam(MASTtex, V , instructor=instructor, examdate=examdate, course=course,
            examname=examname, instructions=instructions, ncol =ncol)


 
#### OLD:    prep.exam(MASTtex, V, exnumber=exnumber, seqnum=seqnum , examdate=examdate)

####  system(paste(sep=" ", "latex", outMAST ))

####  system(paste(sep=" ", "latex", outMAST ))

  
 
 #### system(paste(sep=" ", "dvips -Ppdf",  MASTdvi , " >" , MASTps ))
  
####  system(paste(sep=" ", "ps2pdf", MASTps,"  >", MASTpdf))

  cat("####  To get the final output, change directory to current directory.", sep="\n")
  cat("####  Execute the following system commands outside of R:", sep="\n")
  cat(paste(sep=" ", "latex", outMAST ), sep="\n")
  cat(paste(sep=" ", "latex", outMAST ), sep="\n")
  
 cat(paste(sep=" ", "dvips -Ppdf",  MASTdvi , " >" , MASTps ), sep="\n")
 cat(paste(sep=" ", "ps2pdf", MASTps,"  >", MASTpdf), sep="\n")



 #### cat(file="latex.run", "####  To get the final output, change directory to current directory.", sep="\n")
 #### cat(file="latex.run","####  Execute the following system commands outside of R:", sep="\n", append = TRUE)

  cat(file="latex.run", sep="\n", append = TRUE)
  cat(file="latex.run",paste(sep=" ", "latex", outMAST ), sep="\n", append = TRUE)
  cat(file="latex.run",paste(sep=" ", "latex", outMAST ), sep="\n", append = TRUE)
  
 cat(file="latex.run",paste(sep=" ", "dvips -Ppdf",  MASTdvi , " >" , MASTps ), sep="\n", append = TRUE)
 cat(file="latex.run",paste(sep=" ", "ps2pdf", MASTps,"  >", MASTpdf), sep="\n", append = TRUE)

cat(file="latex.run", sep="\n", append = TRUE)
  cat(file="latex.run",paste(sep=" ", "latex", outSOLTN ), sep="\n", append = TRUE)
  cat(file="latex.run",paste(sep=" ", "latex", outSOLTN ), sep="\n", append = TRUE)
  
 cat(file="latex.run",paste(sep=" ", "dvips -Ppdf",  SOLTNdvi , " >" , SOLTNps ), sep="\n", append = TRUE)
 cat(file="latex.run",paste(sep=" ", "ps2pdf", SOLTNps,"  >", SOLTNpdf), sep="\n", append = TRUE)



 cat(file="clean.run","/bin/rm *.aux *.dvi *.pdf *.ps *.eps *.log *.key *.ANS *SOLTN *.tex latex.run" , sep="\n", append = FALSE)








  
###########   make a key:
IDfier = paste(sep=" ",  V, exnumber, seqnum, examdate, instructor, course)

  
  keyfile = paste(sep=".", outMAST, "key") 
  cat(file=keyfile , sep="\n", append = FALSE)
  cat(file=keyfile ,IDfier,  sep="\n", append=TRUE)
  for(i in 1:length(QTEMP)) {
     txtout = QTEMP[[i]]$a
     txtout = substr(txtout, 1, min(nchar(txtout), 50)   )
    putout  =  paste(sep="\t", i, "ORIG:",  QBord[i], "BUBBLE:",  QTEMP[[i]]$numANS, txtout)
    cat(file=keyfile, putout  , sep="\n",append = TRUE)

  }

  
 ####
 ####  latex examBsolutions.tex


  

  
}

