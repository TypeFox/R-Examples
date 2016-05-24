## ----Setup,hidden=TRUE,echo=FALSE---------------------------------------------
library(knitr)
options(width = 80)
opts_chunk$set(fig.align="center", fig.width=5, fig.height=5, size="scriptsize")

## ----install.package, echo=TRUE, eval=FALSE-----------------------------------
#  install.packages("textreg_0.1.tar.gz", repos = NULL, type="source")

## ----LoadPackage, echo=TRUE, message=FALSE------------------------------------
library( textreg )
library( tm )
data( bathtub )
bathtub

## ----GetLabeling, echo=TRUE---------------------------------------------------
mth.lab = meta(bathtub)$meth.chl
table( mth.lab )

## ----Getbanwords,echo=TRUE----------------------------------------------------
banwords = c( "methylene", "chloride")

## ----DoRegression, echo=TRUE, out.width="0.5\\textwidth"----------------------
corpus = build.corpus(bathtub, mth.lab, banwords, C=4, gap=1, min.support = 1, 
                      verbosity=0, convergence.threshold=0.00001, maxIter=100 )
rs = textreg(corpus, C=4, gap=1, min.support = 1, 
            verbosity=0, convergence.threshold=0.00001, maxIter=100 )
rs

## ----See.results, echo=TRUE---------------------------------------------------
print( reformat.textreg.model( rs ), row.names=FALSE )

## ----plot_results, echo=TRUE--------------------------------------------------
plot( rs )

## ----Play.with.parameters, echo=TRUE------------------------------------------
corpus = build.corpus( bathtub, mth.lab, banwords, C = 5, gap=1, min.support = 10, 
                      verbosity=0, convergence.threshold=0.00001, maxIter=100 )
rs5 = textreg( corpus, C = 5, gap=1, min.support = 1, 
            verbosity=0, convergence.threshold=0.00001, maxIter=100 )
rsLq5 = textreg( corpus, C = 3, Lq=5, gap=1, min.support = 1, 
               verbosity=0, convergence.threshold=0.00001, maxIter=100 )
rsMinSup10 = textreg( corpus, C = 3, Lq=5, gap=1, min.support = 10,
                    verbosity=0, positive.only=TRUE, convergence.threshold=0.00001, maxIter=100 )
rsMinPat2 = textreg( corpus, C = 3, Lq=5, gap=1, min.support = 1, 
                   min.pattern=2, verbosity=0, convergence.threshold=0.00001, maxIter=100 )

## ----show.different.models, results='asis', echo=TRUE-------------------------
library(xtable)
lst = list( rs5, rsLq5, rsMinSup10, rsMinPat2 )
names(lst) = c("C=5", "Lq=5","sup=10", "pat=2")
tbl = make.list.table( lst, topic="Misc Models" )
print( xtable( tbl, caption="Table from the make.list.table call" ), 
       latex.environments="tiny" )

## ----plot_different_models, echo=TRUE-----------------------------------------
list.table.chart( tbl )

## ----FindC, echo=TRUE---------------------------------------------------------
corpus = build.corpus( bathtub, mth.lab, banwords, gap=1, min.support = 5, 
                       verbosity=0, convergence.threshold=0.00001  )
Cs = find.threshold.C( corpus, R = 100, gap=1, min.support = 5, 
                       verbosity=0, convergence.threshold=0.00001 )

Cs[1]
summary( Cs[-1] )

C = quantile( Cs, 0.95 )
C

## ----dropDocs, echo=TRUE------------------------------------------------------
mth.lab.lit = mth.lab
mth.lab.lit[20:length(mth.lab)] = 0

corpus = build.corpus( bathtub, mth.lab.lit, banwords, C = 4, gap=1, min.support = 1, verbosity=0 )
rs.lit = textreg( corpus, C = 4, gap=1, min.support = 1, verbosity=0 )
rs.lit
rs.lit$labeling

## ----See.results.loc, echo=TRUE-----------------------------------------------
hits = phrase.matrix( rs )
dim( hits )
t( hits[ 1:10, ] )
hits.lit = phrase.matrix( rs.lit )
dim(hits.lit)

## ----See.results.loc2, echo=TRUE----------------------------------------------
apply( hits[ mth.lab == 1, ], 1, sum )
apply( hits[ mth.lab == 1, ], 2, sum )

## ----phrase.count.demo, echo=TRUE---------------------------------------------
tt2 = phrase.count( "tub * a", bathtub )
head( tt2 )
table( tt2, dnn="Counts for tub * a" )

## ----appearance.pat, echo=TRUE------------------------------------------------
tab = make.phrase.matrix( c( "bathtub", "tub * a" ), bathtub )
head( tab )
table( tab[,2] )

## ----make.count.table.demo, echo=TRUE-----------------------------------------
ct = make.count.table( c( "bathtub", "tub * a", "bath" ), mth.lab, bathtub )
ct

## ----grab.frag.demo, echo=TRUE------------------------------------------------
tmp = grab.fragments( "bathtub", bathtub, char.before=30, char.after=30, clean=TRUE )
tmp[1:3]

## ----sample.frag.demo, echo=TRUE----------------------------------------------
frags = sample.fragments( "tub * a", mth.lab, bathtub, 20, char.before=30, char.after=30 )
frags

## ----ClusterPhraes, echo=TRUE, out.width="0.5\\textwidth"---------------------
cluster.phrases( rs, num.groups=3 )

## ----Make_phrase_cor_chart, echo=TRUE, out.width="0.5\\textwidth"-------------
make.phrase.correlation.chart( rs, count=TRUE, num.groups=3 )

## ----CalcLoss, echo=TRUE------------------------------------------------------
calc.loss( rs )

## ----Prediction, echo=TRUE, out.width="0.5\\textwidth"------------------------
pds = predict( rs )
labs = rs$labeling
table( labs )
boxplot( pds ~ labs, ylim=c(-1,1) ) 
abline( h=c(-1,1), col="red" )

## ----Outofsample, echo=TRUE, out.width="0.5\\textwidth"-----------------------
  smp = sample( length(bathtub), length(bathtub)*0.5 )
  corpus = build.corpus( bathtub[smp], mth.lab[smp], C = 3, gap=1, min.support = 5, 
              verbosity=0, convergence.threshold=0.00001, maxIter=100 )
	rs = textreg( corpus, C = 3, gap=1, min.support = 5, 
              verbosity=0, convergence.threshold=0.00001, maxIter=100 )
	rs
	train.pred = predict( rs )
	test.pred = predict( rs, bathtub[-smp] )

	train.loss = calc.loss( rs )
	train.loss
	test.loss = calc.loss( rs, bathtub[-smp], mth.lab[-smp] )
	test.loss

## ----Cross Validation, echo=TRUE----------------------------------------------
  tbl = find.CV.C( bathtub, mth.lab, c("methylene","chloride"), 4, 8, verbosity=0 )
  print( round( tbl, digits=3 ) )

## ----CrossValidationPlot, echo=TRUE, out.width="0.5\\textwidth"---------------
  rs = make.CV.chart( tbl )
  rs

## ----CleanAndStem, echo=TRUE--------------------------------------------------
data( dirtyBathtub )
strwrap( dirtyBathtub$text[[1]] )
bc = Corpus( VectorSource( dirtyBathtub$text ) )

bc.clean = clean.text( bc )
strwrap( bc.clean[[1]] )
  
bc.stem = stem.corpus(bc.clean, verbose=FALSE)
strwrap( bc.stem[[1]] )

## ----CleanAndStem2, echo=TRUE-------------------------------------------------
  corpus = build.corpus( bc.stem, mth.lab, c("chlorid+", "methylen+"), C=4, verbosity=0 )
  res.stm = textreg( corpus, C=4, verbosity=0 )
  res.stm

  sample.fragments( "that contain+", res.stm$labeling, bc.stem, 5, char.before=10 )
  sample.fragments( "that contain+", res.stm$labeling, bc.clean, 5, char.before=10 )

