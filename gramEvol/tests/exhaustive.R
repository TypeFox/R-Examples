# This program performs regression testing for package "gramEvol"

# Author: Farzad Noorian <farzad.noorian@sydney.edu.au>

# This program is free software, distributed under the terms of
# the GNU General Public License version 2.
# Refer to <http://www.gnu.org/licenses/> for full terms.
################################################################################
library("gramEvol")

ruleDef <- list(expr = gsrule("<var><op><var>"),
                op   = gsrule("+", "-", "*"),
                var  = gsrule("A", "B"))

# Create a grammar object
grammarDef <- CreateGrammar(ruleDef)  		   


# use exhaustive search to find the 
evalFunc <- function(expr) {
  if (as.character(expr) == "B - A") {
    return(0)
  } else {
    return(1)
  }
}

res = GrammaticalExhaustiveSearch(grammarDef, evalFunc)

stopifnot(all(res$bestSequence == c(0,1,1,0)))
