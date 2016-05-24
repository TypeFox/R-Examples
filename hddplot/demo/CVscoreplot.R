## High-dimensional data, classification, and plots 

##                    What groups are of interest?
data(Golub)
data(golubInfo)      # 7129 rows by 72 columns
with(golubInfo, table(cancer, tissue.mf)) 
attach(golubInfo) 
## Identify allB samples for that are BM:f or BM:m or PB:m 
subsetB <- cancer=="allB" & tissue.mf%in%c("BM:f","BM:m","PB:m") 
## Form vector that identifies these as BM:f or BM:m or PB:m 
tissue.mfB <- tissue.mf[subsetB, drop=TRUE] 
## Separate off the relevant columns of the matrix Golub 
GolubB <- Golub[, subsetB] 
detach(golubInfo) 

##    Cross-validation to determine the optimum number of features
tissue.mfB.cv <- cvdisc(GolubB, cl=tissue.mfB, nfeatures=1:20,
                         nfold=c(10,4))  # 10-fold CV, repeat 4 times
 # Accuracy measures is: tissue.mfB.cv$acc.cv

## Calculations for random normal data: 
set.seed(43) 
rGolubB <- matrix(rnorm(prod(dim(GolubB))), nrow=dim(GolubB)[1])
rtissue.mfB.cv <- cvdisc(rGolubB, cl=tissue.mfB, nfeatures=1:20,
                          nfold=c(10,4))

##      Cross-validation: bone marrow (BM) samples only
attach(golubInfo) 
Golub.BM <- Golub[, BM.PB=="BM"] 
cancer.BM <- cancer[BM.PB=="BM"]
detach(golubInfo)
BMonly.cv <- cvdisc(Golub.BM, cl=cancer.BM, nfeatures=1:20,
                    nfold=c(10,4))

## Which genes appear most frequently in the first 3 positions?
genelist <- matrix(tissue.mfB.cv$genelist[1:3, ,], nrow=3)
tab <- table(genelist, row(genelist))
ord <- order(apply(tab,1,sum), decreasing=TRUE)
tab[ord,]
                          
## Panel A: Uses tissue.mfB.cv from above 
tissue.mfB.scores <- 
  cvscores(cvlist = tissue.mfB.cv, nfeatures = 3, cl.other = NULL)

scoreplot(scorelist = tissue.mfB.scores, cl.circle=NULL, 
       prefix="B-cell subset -")

## Panel B: Uses BMonly.cv from above
BMonly.scores <- cvscores(cvlist=BMonly.cv, nfeatures=13, 
                          cl.other=NULL)
scoreplot(scorelist=BMonly.scores, cl.circle=tissue.mfB, 
       circle=tissue.mfB%in%c("BM:f","BM:m"), 
       params=list(circle=list(col=c("cyan","gray"))), 
       prefix="B: BM samples -") 

