# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
#   
#   essay-scoring.R
#   fridolin.wild@wu-wien.ac.at, June 5th 2006
#   
#   Written for a tutorial at the 
#   ProLearn Summer School 2006, Bled, Slowenia
#   


# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# PREPARE TRAINING DATA

# files have been generated with 
# textmatrix(stemming=FALSE, minWordLength=3, minDocFreq=1)

data(corpus_training) 

weighted_training = corpus_training * gw_entropy(corpus_training)
space = lsa( weighted_training, dims=dimcalc_share(share=0.5) )


# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# FOLD IN ESSAYS

# files have been prepared with
# textmatrix( stemming=FALSE, minWordLength=3, 
# vocabulary=rownames(training) )

data(corpus_essays)

weighted_essays = corpus_essays * gw_entropy(corpus_training)
lsaEssays = fold_in( weighted_essays, space )


# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# TEST THEM, BENCHMARK

essay2essay = cor(lsaEssays, method="spearman")
goldstandard = c( "data6_golden_01.txt", "data6_golden_02.txt", "data6_golden_03.txt" )
machinescores = colSums( essay2essay[goldstandard, ] ) / 3

data(corpus_scores)
humanscores = corpus_scores

cor.test(humanscores[names(machinescores),], machinescores, exact=FALSE, method="spearman", alternative="two.sided")


# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# COMPARE TO PURE VECTOR SPACE MODEL

essay2essay = cor(corpus_essays, method="spearman")
machinescores = colSums( essay2essay[goldstandard, ] ) / 3
cor.test(humanscores[names(machinescores),], machinescores, exact=FALSE, method="spearman", alternative="two.sided")

# => impressingly good!
# => in other experiments at the Vienna University
#    of Economics and Business Administration, the 
#    interrater correlation was in the best case .88, 
#    but going down to -0.17 with unfamiliar topics/raters
