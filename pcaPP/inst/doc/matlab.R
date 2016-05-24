### R code from vignette source 'matlab.rnw'

###################################################
### code chunk number 1: matlab.rnw:9-28
###################################################

source ("load.package.name.R")

library (package.name, character.only = TRUE)

vt <- eval (parse (text = paste (package.name, ":::", ".getVtext", sep = "")))

cat (sep = "",
#	"%\\VignetteIndexEntry{Compiling ", vt (1), " for Matlab}\n",	##	these lines cannot be created automatically - unfortunately..
#	"%\\VignetteDepends{", vt (1), "}\n",
#	"%\\VignetteKeywords{Matlab}\n",
#	"%\\VignettePackage{", vt (1), "}\n",
#	"\n",
	"\n",
	"\\newcommand{\\dapck}{", vt (1), "}\n",
	"\\newcommand{\\daver}{", vt (2), "}\n",
	"\n",
	"\n"
)


###################################################
### code chunk number 2: matlab.rnw:94-96
###################################################
cat (sep = "",
	">> cd ('C:/work/", vt(1), "/matlab')")


###################################################
### code chunk number 3: matlab.rnw:104-109
###################################################
cat (sep = "",
	"Compiling the ", vt(1), " package ... ok\n",
	"Copying the '", vt(1), ".mex*' file(s) to '../matlab' ... ok\n",
	"Changing the current directory back to '../matlab' ... ok\n\n",
	"  Successfully compiled the ", vt(1), " package for Matlab!")


###################################################
### code chunk number 4: matlab.rnw:120-121
###################################################
cat (paste ("\\code{", vt(3), "}", sep = "", collapse = ", "))


###################################################
### code chunk number 5: matlab.rnw:127-128
###################################################
	cat (vt (4))


