## Ian Kopacka
## 2010-08-14
##
## Function: createStyleFile
## Package: FFD
## 
## The function checks if the css-file exists and creates 
## one if it doesn't.
##
## This internal function is called by the method 'HTML'
##
## Input parameters:
##     css.file...file name + path of css-file
##
## Return value: Creates a style file.

createStyleFile <- function(css.file){
    if (!file.exists(css.file)){
        file.create(css.file, showWarnings = TRUE)
        cat("BODY {\n  margin-top: 1em;\n  margin-left: 1em;\  nmargin-right: 1em;\n  ",
            "margin-bottom: 1em;\n  color: rgb(68, 68, 68);\n  background-color: white;\n  ",
            "font-family: arial, helvetica, sans-serif;\n  padding: 5px;\nborder: 0px;\n}\n\n",
            "p {\n  font-size: small;\n}\n\n",
            "img {\n  border: 0;\n}\n\n",
            "h1,h2,h3,h4,h5,h6 {\n  font-family: verdana, arial, helvetica, sans-serif;\n}\n\n",
            "h1 {\n  font-size: x-large;\n}\n\n",
            "h2 {\n  font-size: large;\n  color:#6e7e96;\n}\n\n",
            "h3 {\n  font-size: medium;\n}\n\n",
            "h4 {\n  font-size: small;\n}\n\n",
            "h5 {\n  font-size: x-small;\n}\n\n",
            "h6 {\n  font-size: xx-small;\n}\n\n",
            "table.dataframe td {\n  padding: 5px 20px 5px 20px;\n}", file = css.file, append = TRUE)        
    } 
}
