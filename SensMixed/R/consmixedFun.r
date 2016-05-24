##############################################################################
# performs  Conjoint Analysis for lmer object
##############################################################################
consmixedFun <- function(response, Prod_effects, Cons_effects=NULL, Cons, data, structure = 3, alpha.random = 0.1, alpha.fixed = 0.05 )
{
#structure=1  (default structure) : Analysis of main effects, Random consumer effect AND interaction 
#         between consumer and the main effects. 
#        (Automized reduction in random part, NO reduction in fixed part).
#structure=2 : Main effects AND all 2-factor interactions. 
#         Random consumer effect AND interaction between consumer and 
#         all fixed effects (both main and interaction ones). 
#         (Automized reduction in random part, NO reduction in fixed part). 
#structure=3 : Full factorial model with ALL possible fixed and random effects. 
#         (Automized reduction in random part, AND automized reduction in fixed part).
  
#attach(data)

model <- createLMERmodel(structure, data, response, 
                         fixed = list(Product=Prod_effects, 
                                      Consumer=Cons_effects), 
                         random = Cons, FALSE)
  
#check if reduction of the fixed part is required
if(structure==1 || structure==2)
 isFixReduce <- FALSE
else
 isFixReduce <- TRUE
isRandReduce <- TRUE
isLsmeans <- TRUE
  

#check if there are correlations between intercepts and slopes
checkCorr <- function(model)
{
   corr.intsl <- FALSE
   lnST <- length(getME(model, "ST"))
   for(i in 1:lnST)
   {    
      if(nrow(getME(model, "ST")[[i]])>1)
         corr.intsl <- TRUE
   } 
   return(corr.intsl) 
}

if(checkCorr(model))
  isRandReduce <- FALSE

t <- step(model, reduce.fixed = isFixReduce, reduce.random = isRandReduce, 
          alpha.random = alpha.random, alpha.fixed = alpha.fixed)

#detach(data)

return(t)
}

