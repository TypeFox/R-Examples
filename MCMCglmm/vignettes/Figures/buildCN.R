#source("~/Work/AManal/MCMCglmm_2.19/vignettes/Figures/buildCN.R")
#set.seed(32)
set.seed(34)
alone=FALSE
lectures<-c(7,9)
forCRAN=TRUE
JSS<-TRUE
LINUX=TRUE


options(SweaveHooks=list(fig=function() par(bty="l")))
options(useFancyQuotes=FALSE)
options(width=60)



# Things to check add up

# Need to make sure the data are more likely under N(0,0.5) than N(0,1) in Chapter 1 
# Need to make sure the binary data are completely separated in Chapter 2 

#

UPpath="~/Work/AManal/UP_course/Tex/"
MCpath="~/Work/AManal/MCMCglmm_2.19/vignettes/"
MCpath.old="~/Work/AManal/MCMCglmm_2.19/inst/doc/"

library(MCMCglmm)
library(MCMCpack)
library(R2WinBUGS)
library(kinship2)
library(tools)

data(BTdata)

if(length(lectures)>0){

  if(alone==FALSE){
	setwd(MCpath)	
	for(i in lectures){
	   system(paste("rm -f Lecture", i, "*", sep=""))
	}
	system("rm -f CourseNotes.pdf")
	system("rm -f CourseNotes.tex")
	system("rm -f CourseNotes.log")
	system("rm -f CourseNotes.aux")
	system("rm -f CourseNotes.blg")
  }   # remove files from master

  for(i in lectures){
	
   UPpath_tmp<-paste(UPpath,  "Lecture", i, "/", sep="")

   setwd(UPpath_tmp)
  
   system("rm -f *.pdf")
   system("rm -f *.eps")
   system("rm -f *.aux")
   system("rm -f *.log")  # remove old files

   Rnw<-paste(UPpath_tmp, "Lecture", i, ".Rnw", sep="")
   Tex<-paste(UPpath_tmp, "Lecture", i, ".tex", sep="")
   Bib<-paste(UPpath_tmp, "Lecture", i, sep="")
   Pdf<-paste(UPpath_tmp, "Lecture", i, ".pdf", sep="")

   if(alone){
  	 system(paste("ex", Rnw, "-s +/alonefalse +:s/alonefalse/alonetrue +:x"))
   }else{ 
     system(paste("ex", Rnw, "-s +/alonetrue +:s/alonetrue/alonefalse +:x"))
   } 
  		
   Sweave(Rnw)
  
#	system(paste("ex",Tex, "+:s/alonefalse/alonetrue +:x"))

   if(alone){
      if(LINUX){
	system(paste("cp ~/Work/Tex/library/JarLib.bib ", UPpath_tmp, sep=""))
        if(lectures[1]==3){
          system(paste("cp ", MCpath, "/rgl1.pdf ", UPpath_tmp, sep=""))
          system(paste("cp ", MCpath, "/rgl2.pdf ", UPpath_tmp, sep=""))
        }

        system(paste("pdflatex", Tex))
	system(paste("bibtex", Bib))
	system(paste("pdflatex", Tex))
	system(paste("pdflatex", Tex))
	
	if(length(lectures)==1){system(paste("evince", Pdf, "&"))}

      }else{

	system(paste("cp ~/Work/Tex/library/JarLib.bib ", UPpath_tmp, sep=""))

	system(paste("~/Library/TeXShop/bin/pdflatexc", Tex))
	system(paste("~/Library/TeXShop/bin/bibtexc", Bib))
	system(paste("~/Library/TeXShop/bin/pdflatexc", Tex))
	system(paste("~/Library/TeXShop/bin/pdflatexc", Tex))
	
	if(length(lectures)==1){system(paste("open -a Preview", Pdf))}
      }
    }else{
      
      system(paste("cp *.tex", MCpath))
      system(paste("cp *.pdf", MCpath))  # copy tex and pdf files over to master	
    }	
  
  }
			   
  if(alone==FALSE){	
    setwd(MCpath)
    system(paste("cp ", MCpath, "Figures/CourseNotes.Rnw ", MCpath, sep="")) # copy master file back out of figures
    system(paste("cp ~/Work/Tex/library/JarLib.bib ", MCpath, "CourseNotes.bib", sep="")) 
    system(paste("cp ~/Work/Tex/library/JarLib.bib ", MCpath, "Overview.bib", sep="")) 
    system(paste("cp ~/Work/AManal/Tex/jss.cls " , MCpath, sep="")) 
    system(paste("cp ~/Work/AManal/Tex/jss.bst " , MCpath, sep="")) 

    Sweave("CourseNotes.Rnw") 

    if(LINUX){
      system("pdflatex CourseNotes.tex")
      system("bibtex CourseNotes")
      system("pdflatex CourseNotes.tex")
      system("pdflatex CourseNotes.tex")
      system("xdg-open CourseNotes.pdf&")
    }else{
      system("~/Library/TeXShop/bin/pdflatexc CourseNotes.tex")
      system("~/Library/TeXShop/bin/bibtexc CourseNotes")
      system("~/Library/TeXShop/bin/pdflatexc CourseNotes.tex")
      system("~/Library/TeXShop/bin/pdflatexc CourseNotes.tex")
      system("open -a Preview CourseNotes.pdf")
    }
  }
}


if(JSS){
      if(LINUX){
        setwd(MCpath)
        system(paste("cp ", MCpath, "Figures/Overview.Rnw ", MCpath, sep="")) # copy master file back out of figures
        Sweave("Overview.Rnw")
	system("pdflatex Overview.tex")
	system("bibtex Overview")
	system("pdflatex Overview.tex")
	system("pdflatex Overview.tex")
        system("xdg-open Overview.pdf&")
      }else{
	setwd("~/Desktop/Work/Tex")
	Sweave("article.Rnw")
	system("~/Library/TexShop/bin/pdflatexc article.tex")
	system("~/Library/TexShop/bin/bibtexc article")
	system("~/Library/TexShop/bin/pdflatexc article.tex")
	system("~/Library/TexShop/bin/pdflatexc article.tex")
        system("open -a Preview article.pdf")
      }
    }

if(forCRAN){
	system(paste("bibtool -x ", MCpath, "CourseNotes.aux ", "-o ", MCpath, "CourseNotes.bib", sep=""))
	system(paste("bibtool -x ", MCpath, "Overview.aux ", "-o ", MCpath, "Overview.bib", sep=""))
	system(paste("cp ", MCpath, "CourseNotes.Rnw ", MCpath, "Figures", sep=""))
	system(paste("cp ", MCpath, "CourseNotes.tex ", MCpath, "CourseNotes.Rnw", sep=""))
	system(paste("cp ", MCpath, "Overview.Rnw ", MCpath, "Figures", sep=""))
	system(paste("cp ", MCpath, "Overview.tex ", MCpath, "Overview.Rnw", sep=""))	
	system(paste("rm -f ", MCpath, "*\\.aux", sep=""))
	system(paste("rm -f ", MCpath, "*\\.bbl", sep=""))
	system(paste("rm -f ", MCpath, "*\\.blg", sep=""))
	system(paste("rm -f ", MCpath, "*\\.idx", sep=""))
	system(paste("rm -f ", MCpath, "*\\.toc", sep=""))
	system(paste("rm -f ", MCpath, "*\\.out", sep=""))
	system(paste("rm -f ", MCpath, "*\\.log", sep=""))
	system(paste("rm -f ", MCpath, "CourseNotes.tex", sep=""))
	system(paste("rm -f ", MCpath, "Overview.tex", sep=""))
	system(paste("rm -f ", MCpath, "jss.cls", sep=""))
 	system(paste("rm -f ", MCpath, "jss.bst", sep=""))
        tools::compactPDF(MCpath)
        tools::compactPDF(MCpath.old)
	system(paste("cp ", MCpath, "CourseNotes.Rnw ", MCpath, "Figures", sep=""))
	system(paste("cp ", MCpath, "Overview.Rnw ", MCpath, "Figures", sep=""))
	system(paste("cp ", MCpath, "CourseNotes.pdf ", MCpath.old, sep=""))
	system(paste("cp ", MCpath, "Overview.pdf ", MCpath.old, sep=""))
 	system(paste("rm -f ", MCpath, "Overview.ps", sep=""))
}




