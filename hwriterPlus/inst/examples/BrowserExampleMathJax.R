### Browser Example
### This example is intended to show all the capability of hwriterPlus
### to produce an HTML document for display in a browser
###
### This version produces the header information for use with MathJax
###
### DJS 15/12/2012
date()
options(width=80)

### Set some state variables
opSys <- Sys.info()["sysname"]
if (opSys == "Windows"){
  windows <- TRUE
} else {
  windows <- FALSE
}


reportName <- "BrowserExampleMathJax.html"

### Packages required
require(hwriterPlus)
require(Cairo)
require(xtable)
require(MASS)
require(lattice)

### Break function
br <- function(page){
  hwrite("", page, br = TRUE)
}
### Open file for writing
pg <- newPage(reportName,
              title = "Example of a Document for Display in a Browser",
              link.css = c("../css/BrowserExample.css"),
              link.javascript =
              c("http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"))
br(pg)

### Check if values have been set
hwriterEquation
hwriterEquationList

hwrite("<span class = 'title'>
       Example of a Document for Display in a Browser</span>",
       pg, center = TRUE, br = TRUE)
hwrite("<span class = 'subtitle'>
       Dr David J Scott</span>",
       pg, center = TRUE, br = TRUE)
hwrite(hwrite("Entering Text", name = "intro"), pg,
       heading = 1, center = FALSE, br = TRUE)
hwrite("Ordinary paragraph text can be entering using
<font face = 'monospace' >hwrite</font>.
Aspects of the text such as the font family, face, size and colour can
be altered in line by using tags or styles.", pg,
       center = FALSE, br = TRUE)
br(pg)
hwrite("For example the font can be changed to sans-serif
<font face = 'Arial, sans-serif'> like this, </font>
italic <i> like this, </i> bold <b> like this,</b> coloured
<font color = 'blue'> like this,</font> or any combination of these
<b><i><font color = 'blue'>like this.</font></i></b>", pg, br = TRUE)
br(pg)

hwrite("Some text variations may not be rendered appropriately as can be
seen in the examples above. Classes can be defined in the
<span class = 'pkg'>.css</span> file
to deal with this problem.", pg, br = TRUE)
br(pg)

hwrite("Incorporating Special Symbols",
       pg, heading = 1, center = FALSE, br = TRUE)
hwrite("Using appropriate tags and a list of symbols it is possible to
output mathematical and other symbols. Superscripts and subscripts are
also possible.", pg, br = TRUE)
br(pg)
hwrite("Here are some Greek letters: &alpha;
&beta; and &Gamma;,
followed by some examples of superscripts and subscripts:
&alpha;<sup>2</sup> and &beta;<sub>n</sub>. As examples
of other symbols, here are some arrows:
&larr;  and &rArr;, some set symbols:
&exist;  and &cap;, and some operators:
&int;  and &prod;.", pg, br = TRUE)
br(pg)
hwrite("More than one approach is
available to produce mathematical symbols.
Greek letters can be produced by using symbol font:
<font face='symbol, fantasy'>a, b, G</font>, although this may not
work in some browsers without the presence of additional fonts.
Codes can be used instead of \\(\\LaTeX\\)-like names:
&#8592;, &#8747; and &#8719;",
       pg, br = TRUE)
br(pg)
hwrite("Quotes can be incorporated by escaping with \"\\\",
so for example we can obtain \" and \' in the middle of some text.",
       pg, br = TRUE)

### Mathematics
hwrite("Rendering Mathematics", pg, id = "mathematics",
       heading = 1, center = FALSE, br = TRUE)
hwrite("Entering more complex mathematical expressions including
displayed mathematical expressions is difficult.
One approach is to insert
an image produced by other means, for example by \\(\\LaTeX\\) or Word's
equation editor. The following expression was produced by Word's
equation editor then saved as a jpeg.",
       pg, br = TRUE)
br(pg)
hwriteImage("DisplayedEquation.jpg", pg, center = TRUE, br = TRUE)

hwrite("To enable the rendering of complex mathematical expressions
using \\(\\LaTeX\\) syntax, JavaScript is required, the JavaScript being
run when the file is loaded into the browser. The package
<span class = 'pkg'> R2HTML</span> uses MathPlayer to do the
rendering. An alternative is MathJax which is used by
<span class = 'pkg'> org-babel.</span> MathJax is supported by a number of
scientific societies and is sponsored by the American Mathematical Society
and the Society for Industrial and Applied Mathematics:
see <a href = 'http://www.mathjax.org/sponsors/'></a>.
MathPlayer does not work reliably in all browsers and may require an
add-on to function. MathJax is preferred for these reasons in extending
the capabilities of <span class = 'pkg'> hwriter</span> in the
<span class = 'pkg'> hwriterPlus</span> package",
       pg, br = TRUE)
br(pg)
hwrite("Here is an inline expression:\\(\\int_{-\\infty}^{1}f(x)dx\\),
followed by two displayed expressions, one of which is numbered and
has also been assigned a label.
The first example also has a box around it, by assigning the value
<font face = 'monospace'> \"border = '1'\"</font> to the argument
<font face = 'monospace'> table.attributes</font>.",
       pg, br = TRUE)
br(pg)


hwriteLatex(as.latex("\\int_{-\\infty}^{1}f(x)dx",
                     inline = FALSE, count = FALSE),
            page = pg,
            table.attributes = "border = '1'",
            tr.attributes = "bgcolor = 'white'")

hwriteLatex(as.latex("\\{ 26.119 < \\sum_{i=1}^n(X_i-\\bar{X})^2\\}
\\bigcup\\ \\{ 5.629 > \\sum_{i=1}^n (X_i-\\bar{X})^2 \\}.",
                     inline = FALSE, label = "equation1"),
            page = pg,
            tr.attributes = "bgcolor = 'white'",
            td.attributes = c("width = '50'", "align = 'center'",
                              "align = 'right' width = '50'"))
br(pg)

### Check again on equation numbering
hwriterEquation
hwriterEquationList

### Images
hwrite("Incorporating Images", pg,
       heading = 1, center = FALSE, br = TRUE)
hwrite("Incorporating images is complicated by the inability of
different browsers to display images in different formats.
Firefox will not display <span class ='pkg'>.wmf</span> images for example.
The only scalable image format which can be displayed by all the
common browsers (Internet Explorer, Firefox, Safari and Chrome)
appears to be <span class ='pkg'>.svg</span>,
scalable vector graphics.",
        pg, center = FALSE, br = TRUE)
br(pg)
hwrite(paste("Here is a windows metafile, <span class ='pkg'>.wmf</span> image.
Unfortunately this will only display in Internet Explorer (or Microsoft Word),
and there will be no output produced at all if the operating system is not
Microsoft Windows.
The example is the cats data used by Leisch as an Sweave example,
taken from Venables and Ripley (1987).
The data frame contains measurements of heart and body weight
of ",
             nrow(cats), " cats (",
             sum(cats$Sex=='F')," female, ",
             sum(cats$Sex=='F')," male).", sep = ""),
       pg, center = FALSE, br = TRUE)
br(pg)
hwrite("A linear regression model of heart weight by sex and gender was
fitted to this data. The graph is a scatter plot of the data including
the regression lines.",
       pg, center = FALSE, br = TRUE)
br(pg)
### Fit linear model
lm1 <- lm(Hwt~Bwt*Sex, data=cats)
lm1
lattice.options(theme = "col.whitebg")

### Draw graph
if (windows) {
  win.metafile("cats.wmf", height = 4)
  print(xyplot(Hwt~Bwt|Sex, data=cats, type=c("p", "r")))
  dev.off()
  hwriteImage("cats.wmf", pg, center = FALSE, br = TRUE)
}




hwrite("Here is the cats data plot in svg format.
This uses an extension to hwriter to produce HTML code which enables
display in up to date versions of all common browsers.",
       pg, center = FALSE, br = TRUE)
CairoSVG("cats.svg", width = 4, height = 4)
lattice.options(theme = "col.whitebg")
print(xyplot(Hwt~Bwt|Sex, data=cats, type=c("p", "r")))
dev.off()
hwriteSVG("cats.svg", pg, height = 600, width = 600, id = "catsSVG",
          center = FALSE, br = TRUE)

hwrite("A further format is png (portable network graphics).
This should display in up to date versions of all common browsers.
It is a bitmap format however so not scalable.",
       pg, center = FALSE, br = TRUE)
CairoPNG("cats.png", width = 4*80, height = 4*80)

lattice.options(theme = "col.whitebg")
print(xyplot(Hwt~Bwt|Sex, data=cats, type=c("p", "r")))
dev.off()
hwriteImage("cats.png", pg, center = FALSE, br = TRUE)

hwrite("One problem with images is that the size of the image displayed
can vary widely from browser to browser. Obtaining the right sized image
for a particular browser may require a lot of trial and error.",
       pg, center = FALSE, br = TRUE)

### Vectors, matrices and dataframes
hwrite("Vectors, Matrices and Dataframes", pg,
       heading = 1, center = FALSE, br = TRUE)
hwrite("Here is some code producing a character vector", pg,
       center = FALSE, br = TRUE)
form <- y ~ a + b + c
example <- as.character(form)
hwriteOutput(
"form <- y ~ a + b + c
example <- as.character(form)", pg, br = TRUE)
hwrite("This is the result of printing the vector", pg,
       center = FALSE, br = TRUE)
hwrite(example, pg,
       center = TRUE, br = TRUE)

hwrite("Here is some code producing a numeric vector", pg,
       center = FALSE, br = TRUE)
y <- 3*(1:5)
hwriteOutput("y <- 3*(1:5)", pg, br = TRUE)
hwrite("This is the result of printing the vector", pg,
       center = FALSE, br = TRUE)
hwrite(y, pg,
       center = TRUE, br = TRUE)
hwrite("Here is some code producing a matrix", pg,
       center = FALSE, br = TRUE)
mdat <- matrix(c(1,2,3, 11,12,13),
               nrow = 2, ncol = 3, byrow = TRUE,
               dimnames = list(c('row1', 'row2'),
                               c('C.1', 'C.2', 'C.3')))
hwriteOutput(
"mdat <- matrix(c(1,2,3, 11,12,13),
               nrow = 2, ncol = 3, byrow = TRUE,
               dimnames = list(c('row1', 'row2'),
                               c('C.1', 'C.2', 'C.3')))",
             pg, br = TRUE)
hwrite("This is the result of printing the matrix", pg,
       center = FALSE, br = TRUE)
hwrite(mdat, pg,
       center = TRUE, br = TRUE)
hwrite("Here is some code producing a dataframe", pg,
       center = FALSE, br = TRUE)
L3 <- LETTERS[1:3]
d <- data.frame(cbind(x = 1, y = 1:10),
                fac = sample(L3, 10, replace = TRUE))
hwriteOutput(
"L3 <- LETTERS[1:3]
d <- data.frame(cbind(x = 1, y = 1:10),
                fac = sample(L3, 10, replace = TRUE))",
             pg, center = FALSE, br = TRUE)
hwrite("This is the result of printing the dataframe", pg,
       center = FALSE, br = TRUE)
hwrite(d, pg,
       center = TRUE, br = TRUE)


### Tables
hwrite("Tables", pg, heading = 2, center = FALSE, br = TRUE)
hwrite("Here is an example taken from the documentation of
<font face='monospace'>xtable</font>. It shows the output of an
ANOVA table.",
       pg, br = TRUE)
br(pg)
require(xtable)
data(tli)
fm1 <- aov(tlimth ~ sex + ethnicty + grade + disadvg, data=tli)
fm1.table <- print(xtable(fm1), type="html")

hwrite(fm1.table, pg, center = TRUE, br = TRUE)
hwrite("Here is an example taken from the documentation of
<font face='monospace'>table</font>. It shows the output of a
two-dimensional table.",
       pg, br = TRUE)
br(pg)
require(stats)
tbl <- print(xtable(table(state.division, state.region)), type = "html")
hwrite(tbl, pg, center = TRUE, br = TRUE)

### Output
hwrite("R Output", pg, heading = 2, center = FALSE, br = FALSE)
hwrite("To include output from R in a file, the
command <code>capture.output</code> is used to record the output.
Then the output is included in the HTML file by using the command
<code>hwriteOutput</code>, which is not in the package
 <span class = 'pkg'>hwriter</span>.", pg, center = FALSE, br = TRUE)
hwrite("Here is an example from the <code>glm</code> help.",
       pg, center = FALSE, br = TRUE)
aggOut <-
    capture.output(data(iris),
                   str(iris),
                   aggregate(Sepal.Length~Species, iris, mean)
                   )
aggOut
hwriteOutput(
"aggOut <-
    capture.output(data(iris),
                   str(iris),
                   aggregate(Sepal.Length~Species, iris, mean)
                   )",
             pg, center = FALSE, br = TRUE)
hwrite("This produces the following result.", pg, center = FALSE, br = TRUE)
hwriteOutput(aggOut, pg, center = FALSE, br = TRUE)
hwrite("Note the commas separating the parts of the output to be captured.
See the help and examples for <code>capture.output</code> for details.",
       pg, br = TRUE)

### R Session
hwrite("An R Session", pg, heading = 2, center = FALSE, br = FALSE)
hwrite("To capture an R session, or part of one, including both commands
and output, the command <code>txtStart</code> from the package
<span class = 'pkg'>TeachingDemos</span> can be used then
<code>hwriteOutput</code>.
This requires writing to a file, and reading the results back
from the file. A temporary file can used for holding the output.
See <code>?tempfile</code>.",
       pg, center = FALSE, br = TRUE)
hwrite("Here is an example.",
       pg, center = FALSE, br = TRUE)
tmpFile <- tempfile("Session")
#txtStart("Temp1.txt")
txtStart(tmpFile)
clotting <-
    data.frame(
               u = c(5,10,15,20,30,40,60,80,100),
               lot1 = c(118,58,42,35,27,25,21,19,18),
               lot2 = c(69,35,26,21,18,16,13,12,12)
               )
clotting
coef(glm(lot1 ~ log(u), data=clotting, family=Gamma))
txtStop()
#sessionOut <- readLines("Temp1.txt")
sessionOut <- readLines(tmpFile)
sessionOut
hwriteOutput(sessionOut, pg, br = TRUE)
hwrite("An alternative approach which should not mangle line breaks is provided by the <code>script</code> function written by Ross Ihaka",
       pg, center = FALSE, br = TRUE)
tmpFile <- tempfile("Session")
script(tmpFile)
clotting <-
    data.frame(
               u = c(5,10,15,20,30,40,60,80,100),
               lot1 = c(118,58,42,35,27,25,21,19,18),
               lot2 = c(69,35,26,21,18,16,13,12,12)
               )
clotting
coef(glm(lot1 ~ log(u), data=clotting, family=Gamma))
q()
sessionOut <- readLines(tmpFile)
sessionOut
hwriteOutput(sessionOut, pg, br = TRUE)
hwrite("The new function <code>hwriteScript</code> can be used to drop off lines at the beginning and end of the file created by <code>script</code>.",
       pg, center = FALSE, br = TRUE)
hwriteScript(tmpFile, pg, br = TRUE)


### Links
hwrite("Creating Links", pg, heading = 2, center = FALSE, br = FALSE)
hwrite("Since HTML is being produced, it is easy to create links to
other websites. Here is an example of a link to the Statistics Department
website: ",
       pg, br = FALSE)
hwrite("The Department of Statistics.", pg, br = TRUE,
       link = 'http://www.stat.auckland.ac.nz/uoa/')
br(pg)
hwrite("Links may be created within pages using anchors.
Destination anchors in HTML documents may be specified either by the
<font face = 'monospace'>a</font> element (naming it with the
<font face = 'monospace'>name </font> attribute), or by any other
element (naming with the <font face = 'monospace'>id</font> attribute).
Here is a link to the first section created using the
<font face = 'monospace'>a</font> element: ",
       pg, br = FALSE)
hwrite("Entering Text.", pg, br = TRUE,
       link = "#intro")
br(pg)
hwrite("Here is a link to the third section created by naming with
the <font face = 'monospace'>id</font> attribute: ",
       pg, br = FALSE)
hwrite("Rendering Mathematics.", pg, br = TRUE,
       link = "#mathematics")
br(pg)
hwrite("This use of named or identified elements of a documents is
how cross-referencing can be implemented.",
       pg, br = TRUE)
br(pg)
hwrite("The numbered equation entered previously can be cross-referenced
using the link argument to <font face = 'monospace'>hwrite</font>.
That equation is equation 1 at present. Here is the link to the equation: ",
       pg, br = FALSE)
hwrite("Numbered Equation.", pg, br = TRUE,
       link = "#eq:equation1")
### Test paste
paste("We can also retrieve the number of a labeled equation.
For example the numbered equation was the first equation and had the label
\"equation1\" which in full is \"eq:equation1\". We can retrieve the number
using <code> which(hwriterEquationList == \"eq:equation1\")
</code>",which(hwriterEquationList == "eq:equation1"), "or more simply using
the convenience function <code>eqRef</code>, by <code>eqRef(\"equation1\") ")

hwrite(paste("We can also retrieve the number of a labeled equation.
For example the numbered equation was the first equation and had the label
\"equation1\" which in full is \"eq:equation1\". We can retrieve the number
using <code> which(hwriterEquationList == \"eq:equation1\")</code>
which produces ",which(hwriterEquationList == "eq:equation1"),
"or more easily using the convenience function <code>eqRef</code>,
using simply <code>eqRef(\"equation1\")</code>."),
       pg, br = TRUE)
br(pg)
hwrite("The graph in SVG format included previously can also be
cross-referenced since it was assigned an identifier using the
<font face = 'monospace'>id</font> attribute.
Here is a link to that graph: ",
       pg, br = FALSE)
hwrite("Graph in SVG format.", pg, br = TRUE,
       link = "#catsSVG")
br(pg)
### Close file
closePage(pg)

directory <- getwd()
reportName <- paste("file://", directory, "/", reportName, sep = "")
reportName

browseURL(reportName)

q(save = "no")



