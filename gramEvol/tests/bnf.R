# This program performs regression testing for package "gramEvol"

# Author: Farzad Noorian <farzad.noorian@sydney.edu.au>

# This program is free software, distributed under the terms of
# the GNU General Public License version 2.
# Refer to <http://www.gnu.org/licenses/> for full terms.
################################################################################
library("gramEvol")

bnfGrammarDef <- CreateGrammar("test.bnf", startSymb = '<expr>')
stopifnot(GetGrammarNumOfExpressions(bnfGrammarDef) == 18500)

genome = c(2, 1, 0, 0, 3, 3, 3, 1)
expr = GrammarMap(genome, bnfGrammarDef)
stopifnot(as.character(expr) == "c1 * v1/c2 * v2")
stopifnot(expr$type == "T")

