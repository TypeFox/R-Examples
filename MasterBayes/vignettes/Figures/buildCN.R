#source("~/Work/MasterBayes/MasterBayes_2.50/inst/doc/Figures/buildCN.R")
set.seed(33)
forCRAN=TRUE
LINUX=TRUE
options(width=80)

MCpath="~/Work/MasterBayes/MasterBayes_2.50/inst/doc/"

  setwd(MCpath)	
  system("rm Tutorial.pdf")
  system("rm Tutorial.tex")
  system("rm Tutorial.log")
  system("rm Tutorial.aux")
  system("rm Tutorial.blg")
  system("rm Tutorial.bbl")
  system("rm Tutorial.toc")

  system(paste("cp ", MCpath, "Figures/Tutorial.Rnw ", MCpath, sep=""))
  		
  Sweave("Tutorial.Rnw")
  
  system(paste("pdflatex", "Tutorial.tex"))
  system(paste("bibtex", "Tutorial"))
  system(paste("pdflatex", "Tutorial.tex"))
  system(paste("pdflatex", "Tutorial.tex"))
	
  system(paste("evince", "Tutorial.pdf&"))

  system("cp Tutorial.tex Tutorial.Rnw")
  system("rm Tutorial.tex")

  system(paste("rm *\\.aux", sep=""))
  system(paste("rm *\\.bbl", sep="")) 
  system(paste("rm *\\.blg", sep=""))
  system(paste("rm *\\.eps", sep=""))
  system(paste("rm *\\.log", sep=""))
  system(paste("rm *\\.toc", sep=""))

  system(paste("rm Tutorial*-[A-Z]*\\.pdf", sep=""))





