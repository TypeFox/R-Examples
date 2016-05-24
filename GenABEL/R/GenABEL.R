#' GenABEL: an R package for Genome Wide Association Analysis
#' 
#' Genome-wide association (GWA) analysis is a tool of choice 
#' for identification of genes for complex traits. Effective 
#' storage, handling and analysis of GWA data represent a 
#' challenge to modern computational genetics. GWA studies 
#' generate large amount of data: hundreds of thousands of 
#' single nucleotide polymorphisms (SNPs) are genotyped in 
#' hundreds or thousands of patients and controls. Data on 
#' each SNP undergoes several types of analysis: 
#' characterization of frequency distribution, testing of 
#' Hardy-Weinberg equilibrium, analysis of association between 
#' single SNPs and haplotypes and different traits, and so on. 
#' Because SNP genotypes in dense marker sets are correlated, 
#' significance testing in GWA analysis is preferably performed 
#' using computationally intensive permutation test procedures, 
#' further increasing the computational burden.
#' 
#' To make GWA analysis possible on standard desktop computers 
#' we developed GenABEL library which addresses the following
#' objectives: 
#' 
#' (1) Minimization of the amount of rapid access memory (RAM) used 
#' and the time required for data transactions. For this, we developed 
#' an effective data storage and manipulation model.
#' 
#' (2) Maximization of the throughput of GWA analysis. For this, 
#' we designed optimal fast procedures for specific genetic tests. 
#' 
#' Embedding GenABEL into R environment allows for easy data 
#' characterization, exploration and presentation of the results 
#' and gives access to a wide range of standard and special 
#' statistical analysis functions available in base R and specific 
#' R packages, such as "haplo.stats", "genetics", etc.
#' 
#' To see (more or less complete) functionality of GenABEL, try running
#' 
#' demo(ge03d2).
#' 
#' Other demo of interest could be run with demo(srdta). 
#' Depending on your user priveleges in Windows, it may well not run. 
#' In this case, try demo(srdtawin).
#' 
#' The most important functions and classes are:
#' 
#' For converting data from other formats, see
#' 
#' \code{\link{convert.snp.illumina}} (Illumina/Affymetrix-like format). This is 
#' our preferred converting function, very extensively tested. Other conversion 
#' functions include: 
#' \code{\link{convert.snp.text}} (conversion from human-readable GenABEL format),
#' \code{\link{convert.snp.ped}} (Linkage, Merlin, Mach, and similar files),
#' \code{\link{convert.snp.mach}} (Mach-format),
#' \code{\link{convert.snp.tped}} (from PLINK TPED format),
#' \code{\link{convert.snp.affymetrix}} (BRML-style files).
#' 
#' For converting of GenABEL's data to other formats, see
#' \code{\link{export.merlin}} (MERLIN and MACH formats), 
#' \code{\link{export.impute}} (IMPUTE, SNPTEST and CHIAMO formats),
#' \code{\link{export.plink}} (PLINK format, also exports phenotypic data). 
#' 
#' To load the data, see \code{\link{load.gwaa.data}}.
#' 
#' For conversion to DatABEL format (used by ProbABEL and some other 
#' GenABEL suite packages), see 
#' \code{\link{impute2databel}}, 
#' \code{\link{impute2mach}}, 
#' \code{\link{mach2databel}}. 
#' 
#' For data managment and manipulations see
#' \code{\link{merge.gwaa.data}},
#' \code{\link{merge.snp.data}},
#' \code{\link{gwaa.data-class}},
#' \code{\link{snp.data-class}},
#' \code{\link{snp.names}},
#' \code{\link{snp.subset}}.
#' 
#' For merging extra data to the phenotypic part of \code{\link{gwaa.data-class}} object, 
#' see \code{\link{add.phdata}}.
#' 
#' For traits manipulations see 
#' \code{\link{ztransform}} (transformation to standard Normal),
#' \code{\link{rntransform}} (rank-transformation to normality),
#' \code{\link{npsubtreated}} (non-parametric routine to "impute" trait's values in these medicated).
#' 
#' 
#' For quality control, see
#' \code{\link{check.trait}},
#' \code{\link{check.marker}},
#' \code{\link{HWE.show}},
#' \code{\link{summary.snp.data}},
#' \code{\link{perid.summary}},
#' \code{\link{ibs}},
#' \code{\link{hom}}.
#' 
#' For fast analysis function, see
#' \code{\link{scan.gwaa-class}},
#' \code{\link{ccfast}},
#' \code{\link{qtscore}},
#' \code{\link{mmscore}},
#' \code{\link{egscore}},
#' \code{\link{ibs}},
#' \code{\link{r2fast}} (estimate linkage disequilibrium using R2),
#' \code{\link{dprfast}} (estimate linkage disequilibrium using D'),
#' \code{\link{rhofast}}  (estimate linkage disequilibrium using 'rho')
#' 
#' For specific tools facilitating analysis of the data with stratification
#' (population stratification or (possibly unknown) pedigree structure), see
#' \code{\link{qtscore}} (implements basic Genomic Control),
#' \code{\link{ibs}} (computations of IBS / genomic IBD),
#' \code{\link{egscore}} (stratification adjustment following Price et al.),
#' \code{\link{polygenic}} (heritability analysis),
#' \code{\link{polygenic_hglm}} (another function for heritability analysis),
#' \code{\link{mmscore}} (score test of Chen and Abecasis),
#' \code{\link{grammar}} (grammar, grammar-gc, and garmmar-gamma tests of 
#' Aulchenko et al., Amin et al., and Svishcheva et al.).
#' 
#' For functions facilitating construction of tables for your manuscript, see
#' \code{\link{descriptives.marker}},
#' \code{\link{descriptives.trait}},
#' \code{\link{descriptives.scan}}.
#' 
#' For functions recunstructing relationships from genomic data, 
#' see 
#' \code{\link{findRelatives}}, \code{\link{reconstructNPs}}. 
#' 
#' For meta-analysis and related, see help on
#' \code{\link{formetascore}}.
#' 
#' For link to WEB databases, see
#' \code{\link{show.ncbi}}.
#' 
#' For interfaces to other packages and standard R functions, 
#' also for 2D scans, see
#' \code{\link{scan.glm}},
#' \code{\link{scan.glm.2D}},
#' \code{\link{scan.haplo}},
#' \code{\link{scan.haplo.2D}},
#' \code{\link{scan.gwaa-class}},
#' \code{\link{scan.gwaa.2D-class}}.
#' 
#' For graphical facilities, see
#' \code{\link{plot.scan.gwaa}},
#' \code{\link{plot.check.marker}}.
#' 
#' @author Yurii Aulchenko et al. 
#' (see help pages for specific functions)
#' 
#' @references 
#' If you use GenABEL package in your analysis, please cite the following work:
#' 
#' Aulchenko Y.S., Ripke S., Isaacs A., van Duijn C.M. GenABEL: an R package 
#' for genome-wide association analysis. Bioinformatics. 2007 23(10):1294-6.
#' 
#' If you used \code{\link{polygenic}}, please cite
#' 
#' Thompson EA, Shaw RG (1990) Pedigree analysis for quantitative 
#' traits: variance components without matrix inversion. Biometrics 
#' 46, 399-413.
#' 
#' If you used environmental residuals from \code{\link{polygenic}} for 
#' \code{\link{qtscore}}, or used \code{\link{grammar}}, please cite
#' 
#' for original GRAMMAR
#' 
#' Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid association using mixed model 
#' and regression: a fast and simple method for genome-wide pedigree-based quantitative 
#' trait loci association analysis. Genetics. 2007 177(1):577-85.
#' 
#' for GRAMMAR-GC
#' 
#' Amin N, van Duijn CM, Aulchenko YS. A genomic background based method for 
#' association analysis in related individuals. PLoS ONE. 2007 Dec 5;2(12):e1274.
#' 
#' for GRAMMAR-Gamma
#' 
#' Svischeva G, Axenovich TI, Belonogova NM, van Duijn CM, Aulchenko YS. Rapid 
#' variance components-based method for whole-genome association analysis. 
#' Nature Genetics. 2012 44:1166-1170. doi:10.1038/ng.2410 
#' 
#' for GRAMMAR+ transformation
#' 
#' Belonogova NM, Svishcheva GR, van Duijn CM, Aulchenko YS, Axenovich TI (2013) 
#' Region-Based Association Analysis of Human Quantitative Traits in Related Individuals. 
#' PLoS ONE 8(6): e65395. doi:10.1371/journal.pone.0065395 
#' 
#' If you used \code{\link{mmscore}}, please cite
#' 
#' Chen WM, Abecasis GR. Family-based association tests for genome-wide association 
#' scans. Am J Hum Genet. 2007 Nov;81(5):913-26. 
#' 
#' For exact HWE (used in \code{\link{summary.snp.data}}), please cite:
#' 
#' Wigginton G.E., Cutler D.J., Abecasis G.R. A note on exact tests of 
#' Hardy-Weinberg equilibrium. Am J Hum Genet. 2005 76: 887-893.
#' 
#' For haplo.stats (\code{\link{scan.haplo}}, \code{\link{scan.haplo.2D}}), please cite:
#' 
#' Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland GA. Score tests for 
#' association between traits and haplotypes when linkage phase is ambiguous. 
#' Am J Hum Genet. 2002 70:425-434.
#' 
#' For fast LD computations (function \code{\link{dprfast}}, \code{\link{r2fast}}), please cite:
#' 
#' Hao K, Di X, Cawley S. LdCompare: rapid computation of single- and 
#' multiple-marker r2 and genetic coverage. Bioinformatics. 2006 23:252-254.
#' 
#' If you used \code{\link{npsubtreated}}, please cite
#' 
#' Levy D, DeStefano AL, Larson MG, O'Donnell CJ, Lifton RP, Gavras H, Cupples LA, 
#' Myers RH. Evidence for a gene influencing blood pressure on chromosome 17. 
#' Genome scan linkage results for longitudinal blood pressure phenotypes in 
#' subjects from the framingham heart study. Hypertension. 2000 Oct;36(4):477-83.
#' 
#' @keywords package
#' 
#' @seealso \code{DatABEL}, \code{genetics}, \code{haplo.stats}, \code{qvalue}
#'
#' @examples
#' \dontrun{
#' demo(ge03d2)
#' demo(srdta)
#' demo(srdtawin)
#' }
#'
#' @name GenABEL
#' @docType package
#' @title GWAS in R
#' @aliases GenABEL genabel
#' @keywords package
#'
NULL