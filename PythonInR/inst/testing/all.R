q("no")
R

library(testthat)
require(PythonInR)

pyConnect()
#pyConnect("C:\\Python34_64bit\\python.exe")

test_this <- c("Basics.R", "PyAttach.R", "PyCall.R", "PyExec.R", "PyFunction.R", 
               "PyGetSet.R", "PyObject.R", "PyOptions.R", ##"PySource.R", 
               "Utf8.R")#, "PyImport.R")

dir()
fpath <- "."
fpath <- file.path(path.package("PythonInR"), "testing")

pyOptions()

for (i in 1:2){    
    for (i in test_this){
        cat(i, "\n")
        fname <- file.path(fpath, i)
        test_file(fname)
        cat("\n\n")
    }    
}

print("Finished Testing!")

pyExec("print('Hallo R')")
pyExecp("'Hallo R'")

traceback()

