## ----setup,echo=FALSE,results='hide'-------------------------------------
  library(tikzDevice)
  if( !file.exists('figs') ){dir.create( 'figs' )}

  options(tikzMetricsDictionary='.tikzMetrics')

  knitr::opts_chunk$set(
    echo=FALSE,
    fig.path='figs/fig',
    message=FALSE,
    width = 100,
    comment = NA
  )
  knitr::knit_hooks$set(
    source = function(x, options) {
      paste("\\vspace{-1ex}",
            "\\begin{tikzCodeBlock}[listing style=sweavechunk]",
            paste(x, collapse = '\n'),
            "\\end{tikzCodeBlock}\n\n", sep = '\n')
    },
    output = function(x, options) {
      paste("\\vspace{-2ex}",
            "\\begin{Verbatim}[frame=single]",
            sub("\n+$", "", paste(x, collapse = '\n')),
            "\\end{Verbatim}",
            "",
            sep = '\n')
    },
    warning = function(x, options) {
      paste("\\vspace{-2ex}",
            "\\begin{Verbatim}[frame=single,formatcom=\\color{warningcolor}]",
            sub("\n+$", "", paste(strwrap(x, width = options$width),
                                  collapse = '\n')),
            "\\end{Verbatim}",
            "",
            sep = '\n')
    },
    chunk = function(x, options) x
  )

## ----tikzTitlePlot,results='hide'----------------------------------------
  tikz('figs/titlePlot.tex',width=4,height=4)

  x <- seq(-4.5,4.5,length.out=100)
  y <- dnorm(x)

  xi <- seq(-2,2,length.out=30)
  yi <- dnorm(xi)

  plot(x,y,type='l',col='blue',ylab='$p(x)$',xlab='$x$')
  lines(xi,yi,type='s')
  lines(range(xi),c(0,0))
  lines(xi,yi,type='h')
  title(main="$p(x)=\\frac{1}{\\sqrt{2\\pi}}e^{-\\frac{x^2}{2}}$")
  int <- integrate(dnorm,min(xi),max(xi),subdivisions=length(xi))
  text(2.8,0.3,paste(
    "\\small$\\displaystyle\\int_{",min(xi),"}^{",max(xi),"}p(x)dx",
    "\\approx",round(int[['value']],3),'$',sep=''))

  dev.off()

## ----pdf-example,echo=FALSE,results='hide'-------------------------------
  pdf('figs/pdf-example.pdf', width = 3.25, height = 3.25)
  plot(1, 1, main = 'Hello!', ps = 10)
  dev.off()

## ----tikz-example,echo=FALSE,results='hide'------------------------------
  tikz('figs/tikz-example.tex', width = 3.25, height = 3.25)
  plot(1, 1, main = 'Hello \\TeX !')
  dev.off()

## ----tikzArgs, code=formatR::usage(tikz), eval=FALSE---------------------
#  tikz(file = ifelse(onefile, "./Rplots.tex", "./Rplot%03d.tex"), width = 7, height = 7, 
#      onefile = TRUE, bg = "transparent", fg = "black", pointsize = 10, lwdUnit = getOption("tikzLwdUnit"), 
#      standAlone = FALSE, bareBones = FALSE, console = FALSE, sanitize = FALSE, 
#      engine = getOption("tikzDefaultEngine"), documentDeclaration = getOption("tikzDocumentDeclaration"), 
#      packages, footer = getOption("tikzFooter"), symbolicColors = getOption("tikzSymbolicColors"), 
#      colorFileName = "%s_colors.tex", maxSymbolicColors = getOption("tikzMaxSymbolicColors"), 
#      timestamp = TRUE, verbose = interactive())

## ----simpleEx,echo=TRUE,results='hide'-----------------------------------
library(tikzDevice)
tikz('figs/simpleEx.tex',width=3.5,height=3.5)
plot(1,main='Hello World!')
dev.off()

## ----latexEx,echo=TRUE,results='hide',tidy=FALSE-------------------------
library(tikzDevice)
tikz('figs/latexEx.tex',
  width=3.5,height=3.5)

x <- rnorm(10)
y <- x + rnorm(5,sd=0.25)

model <- lm(y ~ x)
rsq <- summary( model )$r.squared
rsq <- signif(rsq,4)

plot(x, y, main='Hello \\LaTeX!')
abline(model, col='red')

mtext(paste("Linear model: $R^{2}=",
  rsq, "$"), line=0.5)

legend('bottomright', legend =
  paste("$y = ",
    round(coef(model)[2],3), 'x +',
    round(coef(model)[1],3), '$',
    sep = ''), bty = 'n')

dev.off()

## ----bareBonesExample,echo=TRUE,results='hide',tidy=FALSE----------------
library(tikzDevice)
library(maps)

tikz('figs/westCoast.tex', bareBones=TRUE)

map('state', regions=c('california', 'oregon', 'washington'),
    lwd=4, col='grey40')

# Insert some named coordinates into the picture that will
# be available once the picture is included into the
# TeX document.
tikzCoord(-124.161, 40.786, 'humBay')
tikzCoord(-122.962, 46.148, 'longView')
tikzCoord(-124.237, 43.378, 'coosBay')
tikzCoord(-122.419, 37.775, 'sfBay')

dev.off()

## ----standAloneExample,echo=TRUE,results='hide',tidy=FALSE---------------
library(tikzDevice)
tikz('standAloneExample.tex',standAlone=TRUE)
plot(sin,-pi,2*pi,main="A Stand Alone TikZ Plot")
dev.off()

## ----standAloneCompileExample, results='hide', eval=FALSE----------------
#  
#    library(tools)
#  
#    catch <- system(paste(Sys.which('pdflatex'),
#      '-interaction=batchmode -output-directory figs/ figs/standAloneExample.tex'),
#      ignore.stderr=T)
#  
#    # If compiling the example failed, we don't want to include a broken link.
#    if( catch == 0 ){
#      pdfLink <- "The file \\\\code{standAloneExample.tex} may then be compiled to produce
#        \\\\href{./figs/standAloneExample.pdf}{standAloneExample.pdf}. "
#    }else{
#      pdfLink <- ""
#    }
#      #%\Sexpr{print(pdfLink)}

## ----xelatexFontVariantExample,tidy=FALSE,echo=TRUE,eval=FALSE,results='hide'----
#  # Set options for using XeLaTeX font variants.
#  options(tikzXelatexPackages = c(
#    getOption('tikzXelatexPackages'),
#    "\\usepackage[colorlinks, breaklinks]{hyperref}",
#    "\\usepackage{color}",
#    "\\definecolor{Gray}{rgb}{.7,.7,.7}",
#    "\\definecolor{lightblue}{rgb}{.2,.5,1}",
#    "\\definecolor{myred}{rgb}{1,0,0}",
#    "\\newcommand{\\red}[1]{\\color{myred} #1}",
#    "\\newcommand{\\reda}[1]{\\color{myred}\\fontspec[Variant=2]{Zapfino}#1}",
#    "\\newcommand{\\redb}[1]{\\color{myred}\\fontspec[Variant=3]{Zapfino}#1}",
#    "\\newcommand{\\redc}[1]{\\color{myred}\\fontspec[Variant=4]{Zapfino}#1}",
#    "\\newcommand{\\redd}[1]{\\color{myred}\\fontspec[Variant=5]{Zapfino}#1}",
#    "\\newcommand{\\rede}[1]{\\color{myred}\\fontspec[Variant=6]{Zapfino}#1}",
#    "\\newcommand{\\redf}[1]{\\color{myred}\\fontspec[Variant=7]{Zapfino}#1}",
#    "\\newcommand{\\redg}[1]{\\color{myred}\\fontspec[Variant=8]{Zapfino}#1}",
#    "\\newcommand{\\lbl}[1]{\\color{lightblue} #1}",
#    "\\newcommand{\\lbla}[1]{\\color{lightblue}\\fontspec[Variant=2]{Zapfino}#1}",
#    "\\newcommand{\\lblb}[1]{\\color{lightblue}\\fontspec[Variant=3]{Zapfino}#1}",
#    "\\newcommand{\\lblc}[1]{\\color{lightblue}\\fontspec[Variant=4]{Zapfino}#1}",
#    "\\newcommand{\\lbld}[1]{\\color{lightblue}\\fontspec[Variant=5]{Zapfino}#1}",
#    "\\newcommand{\\lble}[1]{\\color{lightblue}\\fontspec[Variant=6]{Zapfino}#1}",
#    "\\newcommand{\\lblf}[1]{\\color{lightblue}\\fontspec[Variant=7]{Zapfino}#1}",
#    "\\newcommand{\\lblg}[1]{\\color{lightblue}\\fontspec[Variant=8]{Zapfino}#1}",
#    "\\newcommand{\\old}[1]{",
#    "\\fontspec[Ligatures={Common, Rare},Variant=1,Swashes={LineInitial, LineFinal}]{Zapfino}",
#    "\\fontsize{25pt}{30pt}\\selectfont #1}%",
#    "\\newcommand{\\smallprint}[1]{\\fontspec{Hoefler Text}
#      \\fontsize{10pt}{13pt}\\color{Gray}\\selectfont #1}"
#  ))
#  
#  # Set the content using custom defined commands
#  label <- c(
#    "\\noindent{\\red d}roo{\\lbl g}",
#    "\\noindent{\\reda d}roo{\\lbla g}",
#    "\\noindent{\\redb d}roo{\\lblb g}",
#    "\\noindent{\\redf d}roo{\\lblf g}\\\\[.3cm]",
#    "\\noindent{\\redc d}roo{\\lblc g}",
#    "\\noindent{\\redd d}roo{\\lbld g}",
#    "\\noindent{\\rede d}roo{\\lble g}",
#    "\\noindent{\\redg d}roo{\\lblg g}\\\\[.2cm]"
#  )
#  
#  # Set the titles using custom defined commands, and hyperlinks
#  title <- c(
#  paste(
#    "\\smallprint{D. Taraborelli (2008),",
#    "\\href{http://nitens.org/taraborelli/latex}",
#    "{The Beauty of \\LaTeX}}"
#  ), paste(
#    "\\smallprint{\\\\\\emph{Some rights reserved}.",
#    "\\href{http://creativecommons.org/licenses/by-sa/3.0/}",
#    "{\\textsc{cc-by-sa}}}"
#  ))
#  
#  # Draw the graphic
#  tikz('xelatexEx.tex',
#    standAlone=TRUE,width=5,height=5,
#    engine = 'xetex')
#  lim <- 0:(length(label)+1)
#  plot(lim,lim,cex=0,pch='.',xlab = title[2],ylab='', main = title[1])
#  for(i in 1:length(label))
#    text(i,i,label[i])
#  dev.off()

## ----annotation,echo=TRUE,results='hide',tidy=FALSE----------------------
library(tikzDevice)

# Load some additional TikZ libraries
tikz("figs/annotation.tex",width=4,height=4,
  packages = c(getOption('tikzLatexPackages'),
    "\\usetikzlibrary{decorations.pathreplacing}",
    "\\usetikzlibrary{positioning}",
    "\\usetikzlibrary{shapes.arrows,shapes.symbols}")
)

p <- rgamma (300 ,1)
outliers <- which( p > quantile(p,.75)+1.5*IQR(p) )
boxplot(p)

# Add named coordinates that other TikZ commands can hook onto
tikzCoord(1, min(p[outliers]), 'min outlier')
tikzCoord(1, max(p[outliers]), 'max outlier')

# Use tikzAnnotate to insert arbitrary code, such as drawing a
# fancy path between min outlier and max outlier.
tikzAnnotate(c("\\draw[very thick,red,",
  # Turn the path into a brace.
  'decorate,decoration={brace,amplitude=12pt},',
  # Shift it 1em to the left of the coordinates
  'transform canvas={xshift=-1em}]',
  '(min outlier) --',
  # Add a node with some text in the middle of the path
  'node[single arrow,anchor=tip,fill=white,draw=green,',
  'left=14pt,text width=0.70in,align=center]',
  '{Holy Outliers Batman!}', '(max outlier);'))

# tikzNode can be used to place nodes with customized options and content
tikzNode(
  opts='starburst,fill=green,draw=blue,very thick,right=of max outlier',
  content='Wow!'
)

dev.off()

## ----strWidthDemo,echo=T-------------------------------------------------
getLatexStrWidth( "The symbol: alpha" )
getLatexStrWidth( "The symbol: $\\alpha$" )

## ----charMetricDemo,echo=T,tidy=FALSE------------------------------------
# Get metrics for 'y'
getLatexCharMetrics(121)

# Get metrics for 'x' - the second value is the descent
# and should be zero or very close to zero.
getLatexCharMetrics(120)

## ----charMetricErrors,echo=T,tidy=FALSE----------------------------------
getLatexCharMetrics('y')
getLatexCharMetrics(20)

# Will return metrics for 'y'
getLatexCharMetrics(121.99)

