
lns <-
c("1. Instructions to install packages provided on p.8 of our book", 
  "   were modified as follows (Mar. 2013): \n",
  'pckgs <- c("plyr", "reshape", "RLRsim", "WWGbook", "ellipse")',
 "install.packages(pckgs) \n",
 "2. To install nlmeU package go to book website at:",
 " http://www-personal.umich.edu/~agalecki/Web/book2.html \n",
 "3. To install lme4.0 package (revision 1780  dated 2012-06-26)",
 'install.packages("lme4.0", repos = "http://R-Forge.R-project.org")'
)

writeLines(lns)  


