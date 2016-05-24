#########################################################################/**
# @RdocPackage R.rsp
#
# \description{
#   @eval "gsub('(\\\\|%)', '\\\\\\1', getDescription(R.rsp))"
# }
#
# \section{Installation}{
#   To install this package, call \code{install.packages("R.rsp")}.
# }
#
# \section{To get started}{
#   We recommend that you start by reading one of the '\href{../doc/index.html}{vignettes}';
#   \enumerate{
#     \item A 5 minute slideshow covering the basics of RSP.
#     \item Detailed description of the RSP markup language.
#     \item A one-page RSP reference card.
#     \item How to use RSP for package vignettes.
#     \item How to use plain LaTeX for package vignettes.
#     \item How to use static PDF or HTML package vignettes.
#   }
#
#   Then, when you're ready to try it yourself, these are commands you can start with:
#   \enumerate{
#     \item Play with @see "rcat", which works like @see "base::cat" but also processed RSP expressions, e.g. \code{rcat("A random number: <\%=sample(100, size=1)\%>\n")}.
#     \item To source a RSP document as you do with R scripts, use @see "rsource", e.g. \code{rsource("report.md.rsp")} which will run the RSP and display the output as it appears.
#     \item To compile a RSP document to a final document, use @see "rfile", e.g. \code{rfile("report.md.rsp")} outputs Markdown file 'report.md' which is automatically compiled into a final 'report.html'.
#   }
# }
#
# \section{Acknowledgments}{
#   Several of the post-processing features of this package utilize
#   packages such as \pkg{base64enc}, \pkg{ascii}, \pkg{knitr}, and
#   \pkg{markdown}.
#   Not enough credit can be given to the authors and contributors
#   of those packages.  Thank you for your great work.
# }
#
# \section{License}{
#   The releases of this package is licensed under
#   LGPL version 2.1 or newer.
#
#   The development code of the packages is under a private licence
#   (where applicable) and patches sent to the author fall under the
#   latter license, but will be, if incorporated, released under the
#   "release" license above.
# }
#
# \section{How to cite this package}{
#  @eval "paste(capture.output(print(citation('R.rsp'), style='latex')), collapse='\n')"
# }
#
# @author
#*/#########################################################################
