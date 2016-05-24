# This program performs regression testing for package "gramEvol"

# Author: Farzad Noorian <farzad.noorian@sydney.edu.au>

# This program is free software, distributed under the terms of
# the GNU General Public License version 2.
# Refer to <http://www.gnu.org/licenses/> for full terms.
################################################################################

library("gramEvol")

set.seed(5)

ruleDef <- list(list("<expr>",     list("<der-expr><op><der-expr>")),
                list("<der-expr>", list("<func>(<var>)", "<var>")),
                list("<func>",     list("log", "exp", "sin", "cos")),
                list("<op>",       list("+", "-", "*")),
                list("<var>",      list("A", "B", "<n>")),
                list("<n>",        list("1", "2", "3", "4")))

grammarDef <- CreateGrammar(ruleDef, startSymb = "<expr>")

costFunction <- function(expr) {
  A <- 1:6
  B <- c(2, 5, 8, 3, 4, 1)
  
  result <- eval(expr)
  
  X <- log(A) * B
  err <- sum((result - X)^2)
  
  return(err)
}

ge <- GrammaticalEvolution(grammarDef, costFunction, terminationCost = 1e-3)

stopifnot(as.character(ge$best$expressions) == "log(A) * B")


