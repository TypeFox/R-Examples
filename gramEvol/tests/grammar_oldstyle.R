# This program performs regression testing for package "gramEvol"

# Author: Farzad Noorian <farzad.noorian@sydney.edu.au>

# This program is free software, distributed under the terms of
# the GNU General Public License version 2.
# Refer to <http://www.gnu.org/licenses/> for full terms.
################################################################################

library("gramEvol")

ruleDef <- list(list('<expr>',     list('<expr><op><expr>', '<coef>*<var>')),
                list('<op>',       list('+', '-', '*', '/')),
                list('<coef>',     list('c1', 'c2')),
                list('<var>',      list('v1', 'v2')))

grammarDef <- CreateGrammar(ruleDef, startSymb = '<expr>')			   

grammarDef
################
stopifnot(GetGrammarDepth(grammarDef) == 4)
stopifnot(GetGrammarDepth(grammarDef, max.depth=10) == 10)
################

ruleDef2 <- list(list('<expr>',    list('<subexpr><op><subexpr>')),
                 list('<subexpr>', list('<coef>*<var>', 'sqrt(<var>)')),
                 list('<op>',      list('+', '-', '*', '/')),
                 list('<coef>',    list('c1', 'c2')),
                 list('<var>',     list('v1', 'v2')))

grammarDef2 <- CreateGrammar(ruleDef2, startSymb = '<expr>')			   

stopifnot(GetGrammarDepth(grammarDef2) == 3)
################
stopifnot(GetGrammarDepth(grammarDef2, startSymb = "<subexpr>") == 2)
stopifnot(GetGrammarDepth(grammarDef2, startSymb = "<coef>") == 1)
################
stopifnot(GetGrammarMaxRuleSize(grammarDef) == 4)
################
stopifnot(GetGrammarNumOfExpressions(grammarDef) == 18500)
stopifnot(GetGrammarNumOfExpressions(grammarDef, max.depth = 2) == 4)
stopifnot(GetGrammarNumOfExpressions(grammarDef, max.depth = 3) == 68)
stopifnot(GetGrammarNumOfExpressions(grammarDef, startSymb = "<coef>") == 2)
################
stopifnot(GetGrammarMaxSequenceLen(grammarDef) == 18)
stopifnot(GetGrammarMaxSequenceLen(grammarDef2, startSymb = "<subexpr>") == 3)
################
stopifnot(GetGrammarMaxSequenceLen(grammarDef, max.depth = 3) == 8)
################
genome = c(2, 1, 0, 0, 3, 3, 3, 1)
expr = GrammarGenotypeToPhenotype(genome, grammarDef)
stopifnot(expr$expr == "c1*v1/c2*v2")
stopifnot(expr$type == "T")
################
inputString = c(0) # a recursive example
expr2 = GrammarGenotypeToPhenotype(inputString, grammarDef, wrapping=2)
stopifnot(expr2$expr == "<expr><op><expr><op><expr>")
stopifnot(expr2$type == "NT")
################
string = expr$expr

c1 = 1
c2 = 2
v1 = 3
v2 = 4
r1 = eval(parse(text = string))
stopifnot(r1 == 6)
################
df = data.frame(c1 = c(1,2),
                c2 = c(2,3),
                v1 = c(3,4),
                v2 = c(4,5))

expressions = list(expr, "c1 * v2")
r2 = EvalExpressions(expressions, envir = df)
stopifnot(dim(r2) == c(2,2))
stopifnot(colnames(r2) == c("expr1","expr2"))
stopifnot(r2$expr1[1] == 6)
stopifnot(r2$expr1[2] * 3 == 40)
stopifnot(r2$expr2 == c(4,10))

