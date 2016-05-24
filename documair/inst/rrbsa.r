#
# 14_09_06 14_09_07 14_09_12 14_09_15
#
## loading the libraries
#  (in fact /rbsa/ is not available in /documair/
#   you can get it from the R forge
#   https://r-forge.r-project.org/projects/riskassessment/)
library(rbsa); library(documair);
#
# /rbsa/ code file path
rbsapath <- "/home/jbdenis/attente.liana/inra/paquets/rbsa/pro/perso";
rbsapath <- "/home/jbdenis/bananier/inra/p/r/paquets/rbsa/pro/perso";
#
# list of the necessary objects
nefu <- c("places4text6tags","explore8list","get8comp7list","list2text",
          "vma2text",
          "rbsa0","text3places8word","object9","form3crop","file2text",
          "fidi9","now","bc","bf","belong9","text3replace","text3brackets",
          "text3preparation","dipa","erreur","text3places8brackets","void9",
          "form3title","parse8text","file2list","text3translate","form3titre",
          "form3display","text2list","texts4text","text3ij2n","pause",
          "filter8text","form3repeat","form3justify","list4text",
          "interv7belonging","bd","text2vma","text3n2ij","text3stext",
          "text3acceptance","text2file","text3interval","form3parag",
          "display8k","rbsa7list9");
#
# files to receive the functions
form3title("Preparing The Intermediate File",4);
newint <- "tutu.inte.r";
newfi2 <- "rrbsa.code.r";
if (fidi9(newint) == "f") { unlink(newint);}
if (fidi9(newfi2) == "f") { unlink(newfi2);}
#
# getting the objects from /rbsa/ directory
rcopy <- copy8objects(Sys.glob(paste0(rbsapath,"/*.code.r")),nefu,newint);
if (length(rcopy)>0) {
  print(rcopy);
  form3title("Found difficulties",6);
  stop("Have a look");
}
#
# translating the objects names
form3title("Translating Object Names",4);
newnames <- cbind(nefu,paste0("rrr",nefu));
faire <- change8names(c(newint,newfi2),newnames);
if (length(faire)>0) {
  print(faire);
  form3title("The translations of the function names went wrong!",5);
  stop("Better to see what happened");
}
#
# checking by building the subpackage
form3title("Checking the sub-package by R itself",4);
rrbbii <- "rrbsa.essairrbsa.r";
if (fidi9(rrbbii)=="f") { unlink(rrbbii);}
if (!file.copy(newfi2,rrbbii)) {
  erreur(list(newfi2,rrbbii),"Not possible to copy 'rrbbii'?");
}
construire <- build8pkg("rrbsa",".","rrbsa",what="pzl");
if (length(construire)>0) {
  print(construire);
  form3title("The Building of the package was not good!",4);
  stop("Look also at the R messages");
}
#
## the end
if (fidi9(newint) == "f") { unlink(newint);}
form3title("rrbsa.r* finished its job!",8);
