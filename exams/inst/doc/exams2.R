### R code from vignette source 'exams2.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 76, prompt = "R> ", continue = "+  ")
library("exams")
library("base64enc")
library("tth")


###################################################
### code chunk number 2: exams2.Rnw:228-231
###################################################
invisible(file.copy(system.file("exercises", "tstat.Rnw", package = "exams"), "tstat.Rnw"))
Rnw <- readLines("tstat.Rnw")
cat(c("\\begin{verbatim}", Rnw, "\\end{verbatim}"), sep = "\n")


###################################################
### code chunk number 3: exams2.Rnw:243-247
###################################################
set.seed(1090)
Sweave("tstat.Rnw")
tex <- readLines("tstat.tex")
file.remove(c("tstat.Rnw", "tstat.tex"))


###################################################
### code chunk number 4: exams2.Rnw:249-250
###################################################
cat(c("\\begin{verbatim}", tex, "\\end{verbatim}"), sep = "\n")


###################################################
### code chunk number 5: exams2.Rnw:263-264
###################################################
cat(tex, sep = "\n")


###################################################
### code chunk number 6: exams1 (eval = FALSE)
###################################################
## library("exams")
## set.seed(1090)
## exams("tstat.Rnw")


###################################################
### code chunk number 7: myexam1
###################################################
myexam <- list(
  "boxplots",
  c("confint", "ttest", "tstat"),
  c("anova", "regression"),
  "scatterplot",
  "relfreq")
odir <- tempfile()
set.seed(1090)
exams(myexam, n = 3, dir = odir, template = c("exam", "solution"))


###################################################
### code chunk number 8: myexam1-odir
###################################################
dir(odir)


###################################################
### code chunk number 9: exams2pdf (eval = FALSE)
###################################################
## set.seed(1090)
## exams2pdf("tstat.Rnw")


###################################################
### code chunk number 10: exams2html (eval = FALSE)
###################################################
## set.seed(1090)
## exams2html("tstat.Rnw")


###################################################
### code chunk number 11: exams2moodle
###################################################
set.seed(1090)
exams2moodle(myexam, n = 3, dir = odir)


###################################################
### code chunk number 12: myexam2-odir
###################################################
dir(odir)


###################################################
### code chunk number 13: exams2qti12
###################################################
set.seed(1090)
exams2qti12(myexam, n = 3, dir = odir)


###################################################
### code chunk number 14: myexam3-odir
###################################################
dir(odir)


###################################################
### code chunk number 15: exams_skeleton
###################################################
mydir <- file.path(tempdir(), "myexam")
exams_skeleton(dir = mydir, absolute = TRUE,
  writer = c("exams2html", "exams2pdf", "exams2moodle"))
dir(mydir)


###################################################
### code chunk number 16: answerlist-question
###################################################
qu <- c("Zurich is the capital of Switzerland.",
        "Italian is an official language in Switzerland.",
        "Switzerland is part of the European Union (EU).")
answerlist(qu)


###################################################
### code chunk number 17: answerlist-solution
###################################################
sol <- c(FALSE, TRUE, FALSE)
ex <- c("The capital of Switzerland is Bern.",
        "The official languages are: German, French, Italian, Romansh.",
        "Switzerland is part of the Schengen Area but not the EU.")
answerlist(ifelse(sol, "True", "False"), ex)


###################################################
### code chunk number 18: myexam4-xexams1
###################################################
set.seed(1090)
x <- xexams(myexam, n = 3)


###################################################
### code chunk number 19: myexam4-xexams2
###################################################
writeLines(x[[1]][[1]]$question)
x[[1]][[1]]$supplements


###################################################
### code chunk number 20: myexam4-xexams3
###################################################
x[[1]][[1]]$questionlist
x[[1]][[1]]$metainfo[c("file", "type", "solution")]


###################################################
### code chunk number 21: myexam4-xexams4
###################################################
trafo <- make_exercise_transform_html(converter = "ttm", base64 = FALSE)
writeLines(trafo(x[[1]][[1]])$question)


###################################################
### code chunk number 22: ttm-tth
###################################################
tex <- c("This is \\textbf{bold} and this \\textit{italic}.",
  "Points on the unit circle: $x^2 + y^2 = 1$.")
ttm(tex)
tth(tex)


###################################################
### code chunk number 23: tex2image (eval = FALSE)
###################################################
## (tex2image(tex, dir = odir, show = FALSE))


###################################################
### code chunk number 24: tex2image-out
###################################################
t2i <- file.path(odir, "tex2image_1.png")
foo <- file.create(t2i)
print(t2i)


###################################################
### code chunk number 25: exams2qti12-boxhist (eval = FALSE)
###################################################
## set.seed(1090)
## exams2qti12("boxhist", n = 1, name = "boxhist")


