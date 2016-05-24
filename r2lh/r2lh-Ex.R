pkgname <- "r2lh"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('r2lh')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("examCheating")
### * examCheating

flush(stderr()); flush(stdout())

### Name: examCheating
### Title: ~ Data: Exam Cheating at French University ~
### Aliases: examCheating
### Keywords: datasets print interface utilities univar

### ** Examples

data(examCheating)
str(examCheating)



cleanEx()
nameEx("rthMainFile")
### * rthMainFile

flush(stderr()); flush(stdout())

### Name: r2hMainFile
### Title: ~ Generation of HTML main document ~
### Aliases: r2hMainFile
### Keywords: print interface utilities

### ** Examples

 # # # # # # # # # # # # # # # # # # #
#   R to HTML, Main file generation   #
 #             Examples              #
  #           r2hMainFile           #
   # # # # # # # # # # # # # # # # #


### Creates a file named "main.html" that includes "univ.html"
if(!file.exists("tmp")){dir.create("tmp")};setwd("tmp")
r2hMainFile()
setwd("../..")



cleanEx()
nameEx("rthb")
### * rthb

flush(stderr()); flush(stdout())

### Name: r2hb
### Title: ~ Function: R to HTML, Bivariate analysis ~
### Aliases: r2hb r2htmlBiv
### Keywords: print classes programming interface utilities univar

### ** Examples

 # # # # # # # # # # # # # # # # # # #
#    R to HTML, Bivariate Analyses    #
 #        Artificial examples        #
  #         Single variable         #
   # # # # # # # # # # # # # # # # #

### Create some data
V1 <- factor(LETTERS[floor(runif(50,1,4))])
V2 <- rnorm(50,1,1)<0
V3 <- ordered(LETTERS[floor(runif(50,1,4))])

### Create a directory for the output
if(!file.exists("tmp/r2hbExample",recursive=TRUE)){dir.create("tmp/r2hbExample",recursive=TRUE)}else{}
setwd("tmp/r2hbExample")

### Execute r2hb
r2hb(V1~V2,fileOut="first.html",textBefore="<H2>Variables V1, V2, V3</H2>",graphName="Gr1",type="png")
r2hb(V2~V1,fileOut="second.html",graphName="Gr2",type="png")
r2hb(V3~V1,fileOut="third.html",textBefore="This is V3 vs. V1",graphDir="P",graphName="Gr3",type="png",displayStyle=2)
r2hMainFile(text="
<LU>
<LI><A HREF='first.html'>First example</A></LI>
<LI><A HREF='second.html'>Second example</A></LI>
<LI><A HREF='third.html'>Third example</A></LI>
</LU>
")
setwd("..")




cleanEx()
nameEx("rthu")
### * rthu

flush(stderr()); flush(stdout())

### Name: r2hu
### Title: ~ Function: R to HTML, Univariate analysis ~
### Aliases: r2hu r2htmlUniv
### Keywords: print classes programming interface utilities univar

### ** Examples

 # # # # # # # # # # # # # # # # # # #
#   R to HTML,  Univariate Analyses   #
 #             Examples              #
  #       r2hu single variable      #
   # # # # # # # # # # # # # # # # #

########################

### Create some data
V1 <- factor(LETTERS[floor(runif(50,1,4))])
V2 <- rnorm(50,1,1)<0
V3 <- ordered(LETTERS[floor(runif(50,1,4))])

### Create a directory for the output
if(!file.exists("tmp")){dir.create("tmp")};setwd("tmp")
if(!file.exists("r2huExample")){dir.create("r2huExample")};setwd("r2huExample")


### Execute r2hu
r2hu(V1,fileOut="first.html",textBefore="<h2>Variable 1 to 3</h2>",graphName="V1")
r2hu(V2,fileOut="second.html",graphName="V2")
r2hu(V3,fileOut="third.html",textBefore="This is variable 3",graphDir="P")
r2hMainFile(text="
<LU>
<LI><A HREF='first.html'>First example</A></LI>
<LI><A HREF='second.html'>Second example</A></LI>
<LI><A HREF='third.html'>Third example</A></LI>
</LU>
")


 # # # # # # # # # # # # # # # # # # #
#   R to HTML,  Univariate Analyses   #
 #          Real examples            #
  #        r2hu data.frame          #
   # # # # # # # # # # # # # # # # #

########################
###### Step 1: Create the data

data(examCheating)
str(examCheating)

########################
###### Step 2: ordering variable

examCheating$YearOfStudy <- ordered(examCheating$YearOfStudy,levels=c("L1","L2","L3","M1","M2"))
examCheating$Bac <- ordered(examCheating$Bac,levels=c("Remedial exam","Pass","Fairly good","Good","Very good","Summa cum laude"))
for(iColumn in 8:17){
    examCheating[,iColumn] <- ordered(examCheating[,iColumn],levels=c("Never","Rarely","Sometimes","Often","Always"))
}
str(examCheating)


########################
###### Step 3: running r2hu

### Preparation of textBefore, for transition between variable

textBefore <- paste("<h3>",names(examCheating)[c(2:5,18:20)],"</h3>",sep="")

text <- "
<h2>Survey</h2>
  <ul>
    <li> What is your age?</li>
    <li> What is your gender?<li>
    <li> What is your level?<li>
    <li> What is your field?<li>
    <li> Did you cheat at Bac?<li>
    <li> Did you cheat high scool?<li>
    <li> Cheating score<li>
  </ul>

<h2>Univariate analysis</h2>
  <OBJECT data = 'ExamCheat-univ.html' type = 'text/html'></OBJECT>


<h2>More information?</h2>
For a detailled analysis, see
http://christophe.genolini.free.fr/EPO/2007 Fraude/EPO2007-Fraude-Rapport.pdf"

### We can run r2hu
r2hu(examCheating[,c(2:5,18:20)],fileOut="ExamCheat-univ.html",textBefore=textBefore)
r2hMainFile("ExamCheat-main.html",text=text)
setwd("../..")

### Then open ExamCheat-main.html in your browser. It is ready !



cleanEx()
nameEx("rtlMainFile")
### * rtlMainFile

flush(stderr()); flush(stdout())

### Name: r2lMainFile
### Title: ~ Generation of LaTeX main document ~
### Aliases: r2lMainFile
### Keywords: print interface utilities

### ** Examples

 # # # # # # # # # # # # # # # # # # #
#   R to LaTeX, Main file generation  #
 #             Examples              #
  #           r2lMainFile           #
   # # # # # # # # # # # # # # # # #

### Creates a Sweave file
text <- "
\\maketitle
\\tableofcontents

<<>>=
data(examCheating)
@

\\section{Univariate analysis}

<<>>=
r2lu(examCheating$CheatScore,fileOut='ExamCheat-univ1.tex')
@
\\input{ExamCheat-univ1.tex}


\\section{Bivariate analysis}

<<>>=
#r2lb(examCheating$CheatScore~examCheating$Sexe,fileOut='ExamCheat-biv1.tex')
@
\\input{ExamCheat-biv1.tex}
"

if(!file.exists("tmp/sweave",recursive=TRUE)){dir.create("tmp/sweave",recursive=TRUE)}else{}
setwd("tmp/sweave")

r2lMainFile(fileOut="main.Rnw",text=text,sweave=TRUE)
Sweave("main.Rnw")
setwd("../..")



cleanEx()
nameEx("rtlb")
### * rtlb

flush(stderr()); flush(stdout())

### Name: r2lb
### Title: ~ Function: R to LaTeX, Bivariate analysis ~
### Aliases: r2lb r2latexBiv
### Keywords: print classes programming interface utilities univar

### ** Examples

 # # # # # # # # # # # # # # # # # # #
#    R to LaTeX, Bivariate Analyses   #
 #        Artificial examples        #
  #         Single variable         #
   # # # # # # # # # # # # # # # # #

### Create some data
V1 <- factor(LETTERS[floor(runif(50,1,4))])
V2 <- rnorm(50,1,1)<0
V3 <- ordered(LETTERS[floor(runif(50,1,4))])

### Create a directory for the output
if(!file.exists("tmp/r2lbExample",recursive=TRUE)){dir.create("tmp/r2lbExample",recursive=TRUE)}else{}
setwd("tmp/r2lbExample")

### Execute r2lb
r2lb(V1~V2,fileOut="first.tex",textBefore="\\section{Variables V1, V2, V3}",graphName="Gr1",type="postscript")
r2lb(V2~V1,fileOut="second.tex",graphName="Gr2",type="postscript")
r2lb(V3~V1,fileOut="third.tex",textBefore="This is V3 vs. V1",graphDir="P",graphName="Gr3",type="postscript",displayStyle=2)
r2lMainFile(text="\\input{first.tex}\n\\input{second.tex}\n\\input{third.tex}")
setwd("../..")




cleanEx()
nameEx("rtlu")
### * rtlu

flush(stderr()); flush(stdout())

### Name: r2lu
### Title: ~ Function: R to LaTeX, Univariate analysis ~ ~ Function: R to
###   HTML, Univariate analysis ~
### Aliases: r2lu r2latexUniv
### Keywords: print classes programming interface utilities univar

### ** Examples

 # # # # # # # # # # # # # # # # # # #
#   R to LaTeX, Univariate Analyses   #
 #        Artificial examples        #
  #         Single variable         #
   # # # # # # # # # # # # # # # # #

### Create some data
V1 <- factor(LETTERS[floor(runif(50,1,4))])
V2 <- rnorm(50,1,1)<0
V3 <- ordered(LETTERS[floor(runif(50,1,4))])

### Create a directory for the output
if(!file.exists("tmp/r2luExample",recursive=TRUE)){dir.create("tmp/r2luExample",recursive=TRUE)}else{}
setwd("tmp/r2luExample")

### Execute r2lu
r2lu(V1,fileOut="first.tex",textBefore="\\section{Variable 1 to 3}",graphName="V1")
r2lu(V2,fileOut="second.tex",graphName="V2")
r2lu(V3,fileOut="third.tex",textBefore="This is variable 3",graphDir="P")
r2lMainFile(text="\\input{first.tex}\n\\input{second.tex}\n\\input{third.tex}")



 # # # # # # # # # # # # # # # # # # #
#   R to LaTeX, Univariate Analyses   #
 #          Real examples            #
  #        r2lu data.frame          #
   # # # # # # # # # # # # # # # # #

########################
###### Step 1: Create the data

data(examCheating)
str(examCheating)

########################
###### Step 2: ordering variable

examCheating$YearOfStudy <- ordered(examCheating$YearOfStudy,levels=c("L1","L2","L3","M1","M2"))
examCheating$Bac <- ordered(examCheating$Bac,levels=c("Remedial exam","Pass","Fairly good","Good","Very good","Summa cum laude"))
for(iColumn in 8:17){
    examCheating[,iColumn] <- ordered(examCheating[,iColumn],levels=c("Never","Rarely","Sometimes","Often","Always"))
}
str(examCheating)


########################
###### Step 3: running r2lu

### Preparation of textBefore, for transition between variables

textBefore <- paste("\\subsection{",names(examCheating)[c(2:5,18:20)],"}",sep="")

text <- "\\maketitle
\\tableofcontents
\\section{Survey}
  \\begin{enumerate}
    \\item What is your age?
    \\item What is your gender?
    \\item What is your level?
    \\item What is your field?
    \\item Did you cheat at Bac?
    \\item Did you cheat high scool?
    \\item Cheating score
  \\end{enumerate}
\\section{Univariate analysis}
  \\input{ExamCheat-univ.tex}

\\section{More information?}
For a detailled analysis, see
http://christophe.genolini.free.fr/EPO/2007 Fraude/EPO2007-Fraude-Rapport.pdf"

### We can run r2lu
r2lu(examCheating[,c(2:5,18:20)],fileOut="ExamCheat-univ.tex",textBefore=textBefore)
r2lMainFile("ExamCheat-main.tex",text=text)
setwd("../..")

### Then compile ExamCheat-main.tex twice. It is ready !



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
