p.value.correcture <-  function(Dv.pairwise){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#          empirical.value,bootstrapped.values <- pair.pops.Gst(), pair.pops.Dest.Chao(), pair.pops.Dest();

# Output:
#          Dv.pairwise.adjusted -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------  

          # Function, that adjusts the p-values obtained for the values of genetic differentiation
          # when the populations were compared pairwise.

          # The argument of the function is a list that equals the result from
          # the functions pair.pops.Dest(), pair.pops.Dest(Chao) and pair.pops.Gst.
          # Example from pair.pops.Dest():
          
          #       $Dest.loci.pairwise.comparison
          #         Dest.for.locus      locus population1 population2     populationpair 0.90.confidence.level 0.95.confidence.level 0.99.confidence.level p.values
          #      1       0.02380959        L12      Borkum    Langeoog     BorkumLangeoog     0.397058823529413     0.581632653061224     0.811968774427113    0.377
          #      2       0.02380959 LgleichL12      Borkum    Langeoog     BorkumLangeoog     0.435039370078737     0.568645833333334     0.882908352849339   0.3679
          #      3      -0.37980780     LoPi89      Borkum    Langeoog     BorkumLangeoog     0.445945945945951     0.570634420813546     0.795169266055046     0.76
          #      4      -1.62345483    ZuPop56      Borkum    Langeoog     BorkumLangeoog     0.457676703930741      0.58603896103893     0.811968833045222  0.98999
          #      5      -0.38055556        L12      Borkum  Wangerooge   BorkumWangerooge     0.366245107632092     0.494946666036428     0.720114393549218     0.86
          #      6      -0.38055556 LgleichL12      Borkum  Wangerooge   BorkumWangerooge     0.377626994315056     0.489318402011915     0.697002091010118     0.87
          #      7       0.36363636     LoPi89      Borkum  Wangerooge   BorkumWangerooge     0.375735294117647                   0.5     0.691348039215686    0.108
          #      8       0.05092607    ZuPop56      Borkum  Wangerooge   BorkumWangerooge                   0.5      0.62518115942029     0.860641307155773    0.365
          #      9      -0.40833324        L12    Langeoog  Wangerooge LangeoogWangerooge     0.353817629784486     0.511339285714285     0.753258482168656    0.877
          #      10     -0.40833324 LgleichL12    Langeoog  Wangerooge LangeoogWangerooge      0.38461300309598     0.513626104647465     0.734965986394557    0.872
          #      11     -0.42721505     LoPi89    Langeoog  Wangerooge LangeoogWangerooge     0.352728781412997     0.511904761904764      0.83483064516129   0.9009
          #      12     -0.09848485    ZuPop56    Langeoog  Wangerooge LangeoogWangerooge     0.312668918918919     0.458360777014127     0.732727077773542    0.588
          #
          #      $Dest.mean.pairwise.comparison
          #                 Mean.Dest population1 population2     populationpair 0.90.confidence.level 0.95.confidence.level 0.99.confidence.level p.values
          #      1   -0.48891086317704      Borkum    Langeoog     BorkumLangeoog     0.173210036584964     0.243220947634115     0.417173218728602   0.9622
          #      2 -0.0866371684787738      Borkum  Wangerooge   BorkumWangerooge     0.167947625574689     0.224381882908422     0.338647995015044    0.565
          #      3  -0.335591596175755    Langeoog  Wangerooge LangeoogWangerooge     0.136355453244588     0.211930415868027     0.339991808971277   0.9544

          # A Bonferroni, Holm, Hommel and (Benjamini and Hochberg) correction of the p-values per locus
          # is calculated because of the multiple comparison from one data set.
          # Only the p-values referring to the same locus are adjusted at once (to
          # one another).
          # The p-values of the mean values of genetic differentiation measurements
          # (over all loci) for the several population pairs are adjusted in the same way.

Dv.locis.correction <- split(Dv.pairwise[[1]],Dv.pairwise[[1]]$locus)

          # The table is splitted according to the locis that have been
          # examined.
          
Dv.locis.corrected1 <- lapply(Dv.locis.correction,Dv.locis.corrected.calc)
Dv.locis.corrected <- do.call(rbind,Dv.locis.corrected1)

Dv.locis.corrected <- as.data.frame(Dv.locis.corrected)
                          
                                    # Now the correction is carried out for the p-values that were
                                    # calculated for the measures of genetic distance over all loci (the mean values).
                                    
p=as.numeric(as.vector(Dv.pairwise[[2]]$p.values))
      
p.bonferroni <- p.adjust(p,method="bonferroni")
p.holm <- p.adjust(p,method="holm")
p.hommel <- p.adjust(p,method="hommel")
pBH <- p.adjust(p,method="BH")
         
Dv.means.corrected <- cbind(as.data.frame(as.matrix(Dv.pairwise[[2]])),p.bonferroni,p.holm,p.hommel,pBH)
Dv.means.corrected <- as.data.frame(Dv.means.corrected)
      
Dv.pairwise.adjusted=list(Dv.locis.corrected,Dv.means.corrected)
      
assign("Dv.pairwise.adjusted",Dv.pairwise.adjusted,pos = DEMEtics.env)

}
