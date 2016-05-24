twoStageGwasPower <-
function(pD, pG, grr, inheritance="multiplicative", pi.samples, 
                              pi.markers, alpha.marker, n.cases, n.controls) {
 # Purpose: Calculate critical values and powers for a two-stage genome-wide association study
 #          Results for one-stage study, replication study, and joint analysis are given
 # Reference: Skol AD, Scott, LJ, Abecasis GR, Boehnke M (2006) Nature Genetics
 #              doi:10.1038/ng1706
 # R implementation of CaTS computer program:   http://www.sph.umich.edu/csg/abecasis/CaTS/
 # Variables
 # pD           probability of disease in population (prevalence)
 # pG           allele frequency of disease gene in population
 # grr          genotypic relative risk
 # c1           critical value for stage I
 # c2           critical value for stage II only
 # c.joint      critical value for combined (stages I and II) statistic
 # c.onestage   critial value for one stage analyis
 # pi.samples   proportion of samples in Stage I
 # pi.markers   probabilty of selecting a gene in stage I
 # M            number of markers
 # pi.markers*M  number of markers selected for Stage II
 # alpha.genome overall Type I error over all M markers
 # alpha.marker alpha.genome/M, level of test for each marker
 # z1           normalized statistic from Stage I data only
 # z2           normalized statistics from Stage II data only
 # z.joint      z1*sqrt(pi.samples) + z2*sqrt(1 - pi.samples)
 # savings      genotyping savings for joint analysis design as compared to a single-stage design
 # 
  p0p1.out <- freqConvert(pD=pD, pG=pG, grr=grr, inheritance=inheritance)
                            
  constants.out <- twoStageNull(pi.samples=pi.samples, pi.markers=pi.markers, alpha.marker=alpha.marker)
  power.outA <- twoStagePower(pi.samples=pi.samples, pi.markers=pi.markers, p0=p0p1.out$p0, p1=p0p1.out$p1, 
                           c1=constants.out$c1, c2=constants.out$c2, 
                           c.joint=constants.out$c.joint, c.singleStage=constants.out$c.singleStage, 
                           n.cases=n.cases, n.controls=n.controls)
  savings <- genotypingSavings(pi.samples, pi.markers)
  result <- list(power.singleStage=power.outA$power.singleStage, power.joint=power.outA$power.joint, 
                 power.rep=power.outA$power.rep, 
                 c1=constants.out$c1, c2=constants.out$c2, c.joint=constants.out$c.joint, c.singleStage=constants.out$c.singleStage,
                 penetrance.GG=p0p1.out$pD.given.GG, penetrance.Gg=p0p1.out$pD.given.Gg, penetrance.gg=p0p1.out$pD.given.gg,
                 p0=p0p1.out$p0, p1=p0p1.out$p1, p.stageOne=power.outA$p.stageOne,
                 savings=savings) 
  class(result) <- "twoStageGwasPower" 
  result                                                  
  }
