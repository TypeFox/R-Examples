library(prabclus)

# example(prabclust)
data(kykladspecreg)
data(nb)
set.seed(1234)
x <- prabinit(prabmatrix=kykladspecreg, neighborhood=nb)
# If you want to use your own ASCII data files, use
# x <- prabinit(file="path/prabmatrixfile",
# neighborhood="path/neighborhoodfile")
print(prabclust(x))

# Here is an example for species delimitation with codominant markers;
# only 50 individuals were used in order to have a fast example. 
data(tetragonula)
ta <- alleleconvert(strmatrix=tetragonula[1:50,])
tai <- alleleinit(allelematrix=ta)
print(prabclust(tai))

# Here is an example for species delimitation with dominant markers;
# only 50 individuals were used in order to have a fast example.
# You may want to use stressvals to choose mdsdim.
data(veronica)
vei <- prabinit(prabmatrix=veronica[1:50,],distance="jaccard")
print(prabclust(vei,mdsmethod="kruskal",mdsdim=3))

# example(crmatrix)
  options(digits=3)
  data(kykladspecreg)
  data(nb)
  set.seed(1234)
  x <- prabinit(prabmatrix=kykladspecreg, neighborhood=nb)
  xc <- prabclust(x)

  crmatrix(x,xc)
  crmatrix(x,xc, percentages=TRUE)


# example(lociplots)
  options(digits=4)
  data(veronica)
  vei <- prabinit(prabmatrix=veronica[1:50,],distance="jaccard")
  ppv <- prabclust(vei)
  veloci <- prabinit(prabmatrix=veronica[1:50,],rows.are.species=FALSE)
  velociclust <- prabclust(veloci,nnk=0)
  lociplots(ppv,velociclust$clustering,veloci,lcluster=3)


# Results:

# R version 3.1.2 (2014-10-31) -- "Pumpkin Helmet"
# Copyright (C) 2014 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
# 
# > library(prabclus)
# Loading required package: MASS
# Loading required package: mclust
# Package 'mclust' version 4.4
# Type 'citation("mclust")' for citing this R package in publications.
# > 
# > # example(prabclust)
# > data(kykladspecreg)
# > data(nb)
# > set.seed(1234)
# > x <- prabinit(prabmatrix=kykladspecreg, neighborhood=nb)
# > # If you want to use your own ASCII data files, use
# > # x <- prabinit(file="path/prabmatrixfile",
# > # neighborhood="path/neighborhoodfile")
# > print(prabclust(x))
# * Clustered presence-absence matrix * 
# 
# Clustered:  4 -dim. MDS result from method  classical 
# 
# Noise-detector NNclean has been used with k= 2 
# NNclean is explained in S. Byers and A. E. Raftery, JASA 95 (1998), 781-794
# A Normal mixture model with noise component (mclust) has been used.
# Mixture component memberships:
#  [1] 0 1 0 2 2 8 6 0 7 0 2 0 0 4 1 6 6 8 4 0 0 0 4 1 4 0 6 5 3 1 3 5 0 6 1 0 0 1
# [39] 0 8 1 2 3 3 5 0 1 3 2 1 7 0 0 4 5 3 7 4 0 0 4 1 5 7 0 3 2 0 2 3 0 1 7 4 0 0
# [77] 2 5 0 6
# 
# Clustering (N denotes noise or one-point components):
#  [1] "N" "1" "N" "2" "2" "8" "6" "N" "7" "N" "2" "N" "N" "4" "1" "6" "6" "8" "4"
# [20] "N" "N" "N" "4" "1" "4" "N" "6" "5" "3" "1" "3" "5" "N" "6" "1" "N" "N" "1"
# [39] "N" "8" "1" "2" "3" "3" "5" "N" "1" "3" "2" "1" "7" "N" "N" "4" "5" "3" "7"
# [58] "4" "N" "N" "4" "1" "5" "7" "N" "3" "2" "N" "2" "3" "N" "1" "7" "4" "N" "N"
# [77] "2" "5" "N" "6"
# > 
# > # Here is an example for species delimitation with codominant markers;
# > # only 50 individuals were used in order to have a fast example. 
# > data(tetragonula)
# > ta <- alleleconvert(strmatrix=tetragonula[1:50,])
# > tai <- alleleinit(allelematrix=ta)
# > print(prabclust(tai))
# * Clustered presence-absence matrix * 
# 
# Clustered:  4 -dim. MDS result from method  classical 
# 
# Noise-detector NNclean has been used with k= 2 
# NNclean is explained in S. Byers and A. E. Raftery, JASA 95 (1998), 781-794
# A Normal mixture model with noise component (mclust) has been used.
# Mixture component memberships:
#  [1] 2 2 1 1 1 1 1 2 2 1 1 1 2 2 2 1 2 1 2 1 2 2 1 2 1 2 2 2 1 1 1 1 2 1 2 3 0 3
# [39] 0 0 0 3 3 3 3 0 3 3 3 3
# 
# Clustering (N denotes noise or one-point components):
#  [1] "2" "2" "1" "1" "1" "1" "1" "2" "2" "1" "1" "1" "2" "2" "2" "1" "2" "1" "2"
# [20] "1" "2" "2" "1" "2" "1" "2" "2" "2" "1" "1" "1" "1" "2" "1" "2" "3" "N" "3"
# [39] "N" "N" "N" "3" "3" "3" "3" "N" "3" "3" "3" "3"
# > 
# > # Here is an example for species delimitation with dominant markers;
# > # only 50 individuals were used in order to have a fast example.
# > # You may want to use stressvals to choose mdsdim.
# > data(veronica)
# > vei <- prabinit(prabmatrix=veronica[1:50,],distance="jaccard")
# > print(prabclust(vei,mdsmethod="kruskal",mdsdim=3))
# initial  value 28.163173 
# iter   5 value 20.897590
# iter  10 value 19.154545
# iter  15 value 18.814679
# iter  20 value 18.493361
# iter  20 value 18.475223
# final  value 18.228921 
# converged
# * Clustered presence-absence matrix * 
# 
# Clustered:  3 -dim. MDS result from method  kruskal 
# 
# Noise-detector NNclean has been used with k= 2 
# NNclean is explained in S. Byers and A. E. Raftery, JASA 95 (1998), 781-794
# A Normal mixture model with noise component (mclust) has been used.
# Mixture component memberships:
#  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1
# [39] 1 0 0 1 1 1 1 1 1 0 1 0
# 
# Clustering (N denotes noise or one-point components):
#  [1] "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1"
# [20] "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "N" "1" "1" "1" "1" "1" "1" "1" "1"
# [39] "1" "N" "N" "1" "1" "1" "1" "1" "1" "N" "1" "N"
# > 
# > # example(crmatrix)
# >   options(digits=3)
# >   data(kykladspecreg)
# >   data(nb)
# >   set.seed(1234)
# >   x <- prabinit(prabmatrix=kykladspecreg, neighborhood=nb)
# >   xc <- prabclust(x)
# > 
# >   crmatrix(x,xc)
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#  [1,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#  [2,]    0    0    0    0    0    0    0    1    0     1     0     0     0
#  [3,]    0    0    0    0    0    0    0    0    0     1     2     2     2
#  [4,]    1    0    0    0    2    7    3    0    1     2     0     0     1
#  [5,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#  [6,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#  [7,]    5    4    1    0    0    0    0    0    2     0     0     0     0
#  [8,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#  [9,]    9   10    3    3    4    4    4    9    7     9     3     1     4
#       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
#  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
#  [2,]     2     3     2     2     0     1     6     1     3     6     2     2
#  [3,]     1     1     0     1     0     0     0     0     0     0     2     0
#  [4,]     0     0     0     0     0     1     0     0     0     1     0     0
#  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
#  [6,]     0     0     0     0     0     0     0     0     0     0     0     0
#  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
#  [8,]     0     0     0     0     0     0     0     0     0     0     0     0
#  [9,]     8     7     6     5     6     8    10     3     4     7     6     8
#       [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34]
#  [1,]     0     0     7     8     3     0     0     0     0
#  [2,]     2     0     0     0     0     0     0     0     0
#  [3,]     2     0     0     0     0     0     0     0     0
#  [4,]     0     0     0     0     0     0     0     0     0
#  [5,]     0     6     0     0     0     0     0     0     0
#  [6,]     0     0     0     0     3     4     2     6     0
#  [7,]     0     0     0     0     0     0     0     0     0
#  [8,]     0     3     3     3     0     0     0     0     0
#  [9,]     5    10     6    10     5     6     2     7     4
# >   crmatrix(x,xc, percentages=TRUE)
#        [,1] [,2] [,3]   [,4] [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#  [1,] 0.000  0.0 0.00 0.0909 0.00 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [2,] 0.000  0.0 0.00 0.0000 0.00 0.000 0.000 0.125 0.000 0.125  0.00  0.00
#  [3,] 0.000  0.0 0.00 0.0000 0.00 0.000 0.000 0.000 0.000 0.125  0.25  0.25
#  [4,] 0.125  0.0 0.00 0.0000 0.25 0.875 0.375 0.000 0.125 0.250  0.00  0.00
#  [5,] 0.000  0.0 0.00 0.0000 0.00 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [6,] 0.000  0.0 0.00 0.0000 0.00 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [7,] 1.000  0.8 0.20 0.0000 0.00 0.000 0.000 0.000 0.400 0.000  0.00  0.00
#  [8,] 0.000  0.0 0.00 0.0000 0.00 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [9,] 0.360  0.4 0.12 0.1200 0.16 0.160 0.160 0.360 0.280 0.360  0.12  0.04
#       [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#  [1,] 0.000 0.000 0.000  0.00 0.000  0.00 0.000  0.00 0.000 0.000 0.000  0.00
#  [2,] 0.000 0.250 0.375  0.25 0.250  0.00 0.125  0.75 0.125 0.375 0.750  0.25
#  [3,] 0.250 0.125 0.125  0.00 0.125  0.00 0.000  0.00 0.000 0.000 0.000  0.25
#  [4,] 0.125 0.000 0.000  0.00 0.000  0.00 0.125  0.00 0.000 0.000 0.125  0.00
#  [5,] 0.000 0.000 0.000  0.00 0.000  0.00 0.000  0.00 0.000 0.000 0.000  0.00
#  [6,] 0.000 0.000 0.000  0.00 0.000  0.00 0.000  0.00 0.000 0.000 0.000  0.00
#  [7,] 0.000 0.000 0.000  0.00 0.000  0.00 0.000  0.00 0.000 0.000 0.000  0.00
#  [8,] 0.000 0.000 0.000  0.00 0.000  0.00 0.000  0.00 0.000 0.000 0.000  0.00
#  [9,] 0.160 0.320 0.280  0.24 0.200  0.24 0.320  0.40 0.120 0.160 0.280  0.24
#       [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34]
#  [1,]  0.00  0.00   0.0 0.636 0.727 0.273 0.000 0.000  0.00  0.00
#  [2,]  0.25  0.25   0.0 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [3,]  0.00  0.25   0.0 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [4,]  0.00  0.00   0.0 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [5,]  0.00  0.00   1.0 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [6,]  0.00  0.00   0.0 0.000 0.000 0.500 0.667 0.333  1.00  0.00
#  [7,]  0.00  0.00   0.0 0.000 0.000 0.000 0.000 0.000  0.00  0.00
#  [8,]  0.00  0.00   1.0 1.000 1.000 0.000 0.000 0.000  0.00  0.00
#  [9,]  0.32  0.20   0.4 0.240 0.400 0.200 0.240 0.080  0.28  0.16
# > 
# > 
# > # example(lociplots)
# >   options(digits=4)
# >   data(veronica)
# >   vei <- prabinit(prabmatrix=veronica[1:50,],distance="jaccard")
# >   ppv <- prabclust(vei)
# >   veloci <- prabinit(prabmatrix=veronica[1:50,],rows.are.species=FALSE)
# >   velociclust <- prabclust(veloci,nnk=0)
# >   lociplots(ppv,velociclust$clustering,veloci,lcluster=3)
# $locfreq
#  [1] 0.4737 0.3684 0.4737 0.5263 0.3684 0.4211 0.3684 0.5263 0.5263 0.5263
# [11] 0.6842 0.4211 0.4211 0.5263 0.4211 0.2632 0.4211 0.4737 0.5263 0.4211
# [21] 0.4737 0.4211 0.5263 0.4211 0.6316 0.5263 0.4211 0.4737 0.4211 0.3684
# [31] 0.5263 0.4737 0.4737 0.4211 0.4737 0.5789 0.5263 0.5263 0.4211 0.3684
# [41] 0.2632 0.4211 0.3684 0.4737 0.5263 0.4211 0.4211 0.2632 0.5263 0.5263
# 
# $locfreqmin
# [1] 0.2632 0.3684 0.2632 0.4211
# 
# $locfreqmax
# [1] 0.5263 0.5263 0.5789 0.6842
# 
# $locfreqmean
# [1] 0.3947 0.4520 0.4575 0.5044
# 
# > 
# > proc.time()
#    user  system elapsed 
#   1.708   0.020   1.729 
