# This program performs regression testing for package "gramEvol"

# Author: Farzad Noorian <farzad.noorian@sydney.edu.au>

# This program is free software, distributed under the terms of
# the GNU General Public License version 2.
# Refer to <http://www.gnu.org/licenses/> for full terms.
################################################################################

library("gramEvol")

ruleDef <- list(
  expr = gsrule("<expr><op><expr>", "<coef> * <var>"),
  op   = gsrule("+", "-", "*", "/"),
  coef = grule(c1, c2),
  var  = grule(v1, v2)
)

grammarDef <- CreateGrammar(ruleDef)

stopifnot(grammarDef$def[[1]][[2]][[1]] == "<expr><op><expr>")
stopifnot(grammarDef$def[[2]][[2]][[1]] == "+")

grammarDef
################
rule.cat = c(ruleDef$expr, ruleDef$op, ruleDef$coef)
stopifnot("GERule" %in% class(rule.cat) )
stopifnot(length(rule.cat) == 8)
################
stopifnot(GrammarGetDepth(grammarDef) == 4)
stopifnot(GrammarGetDepth(grammarDef, max.depth=10) == 10)
################
ruleDef2 <- list(
  expr =    grule(op(subexpr, subexpr)),
  subexpr = grule(coef * var, sqrt(var)),
  op   =    grule(`+`, `-`, `*`, `/`),
  coef =    grule(c1, c2),
  var  =    grule(v1, v2)
)

grammarDef2 <- CreateGrammar(ruleDef2, startSymb = '<expr>')  		   

stopifnot(GrammarGetDepth(grammarDef2) == 3)
################
ruleDef3 <- list(
  start = grule(rbind(cbind(var, var), cbind(var, var)), cbind(var, var)),
  op =    grule(`+`, 'a', `>`),
  var =   grule(x, y, z)
)  

grammarDef3 = CreateGrammar(ruleDef3)
stopifnot(grammarDef3$def[[1]][[2]][[1]] == "rbind(cbind(<var>, <var>), cbind(<var>, <var>))")
stopifnot(grammarDef3$def[[2]][[2]][[1]] == "`+`")
stopifnot(grammarDef3$def[[2]][[2]][[2]] == "\"a\"")
stopifnot(grammarDef3$startSymb == "<start>")
################
stopifnot(GrammarGetDepth(grammarDef2, startSymb = "<subexpr>") == 2)
stopifnot(GrammarGetDepth(grammarDef2, startSymb = "<coef>") == 1)
################
stopifnot(GrammarMaxRuleSize(grammarDef) == 4)
################
stopifnot(GrammarNumOfExpressions(grammarDef) == 18500)
stopifnot(GrammarNumOfExpressions(grammarDef, max.depth = 2) == 4)
stopifnot(GrammarNumOfExpressions(grammarDef, max.depth = 3) == 68)
stopifnot(GrammarNumOfExpressions(grammarDef, startSymb = "<coef>") == 2)
################
stopifnot(GrammarMaxSequenceLen(grammarDef) == 18)
stopifnot(GrammarMaxSequenceLen(grammarDef2, startSymb = "<subexpr>") == 3)
################
stopifnot(GrammarMaxSequenceLen(grammarDef, max.depth = 3) == 8)
################
genome = c(2, 1, 0, 0, 3, 3, 3, 1)
expr = GrammarMap(genome, grammarDef)
stopifnot(expr$type == "T")
stopifnot(GrammarIsTerminal(expr))
stopifnot(as.character(expr) == "c1 * v1/c2 * v2")
################
inputString = c(0) # a recursive example
expr2 = GrammarMap(inputString, grammarDef, wrapping=2)
stopifnot(expr2$expr == "<expr><op><expr><op><expr>")
stopifnot(expr2$type == "NT")
stopifnot(GrammarIsTerminal(expr2) == FALSE)
################
stopifnot(GrammarIsRecursive(grammarDef, "<expr>"))
stopifnot(!GrammarIsRecursive(grammarDef, "<op>"))
stopifnot(!GrammarIsRecursive(grammarDef2))
################
string = as.expression(expr)

c1 = 1
c2 = 2
v1 = 3
v2 = 4
r1 = eval(string)
stopifnot(r1 == 6)
################
df = data.frame(c1 = c(1,2),
                c2 = c(2,3),
                v1 = c(3,4),
                v2 = c(4,5))

expressions = list(as.expression(expr), quote(c1 * v2))
r2 = lapply(expressions, eval, envir = df)
stopifnot(dim(r2) == c(2,2))
stopifnot(colnames(r2) == c("expr1","expr2"))
stopifnot(r2$expr1[1] == 6)
stopifnot(r2$expr1[2] * 3 == 40)
stopifnot(r2$expr2 == c(4,10))
################
vrule = gvrule(1:10)
stopifnot(length(vrule) == 10)
stopifnot(all(as.vector(vrule) == 1:10))


