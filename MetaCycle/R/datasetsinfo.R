##' @name cycMouseLiverRNA
##' @title cycMouseLiverRNA
##' @description This data set lists expression profiles of 20 circadian
##'   transcripts with 1h-resolution covering two days.
##' @docType data
##' @usage cycMouseLiverRNA
##' @format A dataframe containing 49 columns(column 1 = transcript name,
##'   column 2 to 49 = time points from CT18 to CT65).
##' @source Hughes M. E., et al. (2009). Harmonics of circadian gene
##'   transcription in mammals. \emph{PLoS Genet}, \bold{5(4)}, e1000442.
NULL

##' @name cycMouseLiverProtein
##' @title cycMouseLiverProtein
##' @description This data set lists expression profiles of 5 circadian
##'   proteins with 3h-resolution covering two days.
##' @docType data
##' @usage cycMouseLiverProtein
##' @format A dataframe containing 49 columns(column 1 = protein name,
##'   column 2 to 49 = time points from CT0 to CT45 with three replicates
##'   at each time point).
##' @source Robles M. S., Cox J., Mann M. (2014). In-vivo quantitative
##'   proteomics reveals a key contribution of post-transcriptional
##'   mechanisms to the circadian regulation of liver metabolism.
##'   \emph{PLoS Genet}, \bold{10(1)}, e1004047.
NULL

##' @name cycSimu4h2d
##' @title cycSimu4h2d
##' @description This data set lists 20 simulated profiles(periodic and
##'   non-periodic) with 4h-resolution covering two periods.
##' @docType data
##' @usage cycSimu4h2d
##' @format A dataframe containing 13 columns(column 1 = curve ID, column 2 to
##'   13 = time points from 0 to 44).
##' @source Wu G., Zhu J., Yu J., Zhou L., Huang J. Z. and  Zhang Z. (2014).
##'   Evaluation of five methods for genome-wide circadian gene identification.
##'   \emph{Journal of Biological Rhythms}, \bold{29(4)}, 231--242.
NULL

##' @name cycYeastCycle
##' @title cycYeastCycle
##' @description This data set lists expression profiles of 10 cycling
##'   transcripts with 16-minutes resolution covering about two yeast cell
##'   cycles.
##' @docType data
##' @usage cycYeastCycle
##' @format A dataframe containing 12 columns(column 1 = transcript name,
##'   column 2 to 12 = time points from 2 minutes to 162 minutes after
##'   recovery phase).
##' @source Orlando D. A., et al. (2008). Global control of cell-cycle
##'   transcription by coupled CDK and network oscillators. \emph{Nature},
##'   \bold{453(7197)}, 944--947.
NULL

##' @name cycVignettesAMP
##' @title cycVignettesAMP
##' @description This data set lists meta2d's analysis results of three
##'   circadian transcripts selected from the same source dataset used by
##'   cycMouseLiverRNA.
##' @docType data
##' @usage cycVignettesAMP
##' @format A dataframe containing 71 columns described as below:
##'   \tabular{rlll}{
##'         [,1] \tab CycID                 \tab character  \tab transcript name\cr
##'         [,2] \tab ARS_pvalue            \tab numeric    \tab pvalue from ARS\cr
##'         [,3] \tab ARS_BH.Q              \tab numeric    \tab FDR from ARS\cr
##'         [,4] \tab ARS_period            \tab numeric    \tab period from ARS\cr
##'         [,5] \tab ARS_adjphase          \tab numeric    \tab adjusted phase from ARS\cr
##'         [,6] \tab ARS_amplitude         \tab numeric    \tab amplitude from ARS\cr
##'         [,7] \tab JTK_pvalue            \tab numeric    \tab pvalue from JTK\cr
##'         [,8] \tab JTK_BH.Q              \tab numeric    \tab FDR from JTK\cr
##'         [,9] \tab JTK_period            \tab numeric    \tab period from JTK\cr
##'        [,10] \tab JTK_adjphase          \tab numeric    \tab adjusted phase from JTK\cr
##'        [,11] \tab JTK_amplitude         \tab numeric    \tab amplitude from JTK\cr
##'        [,12] \tab LS_pvalue             \tab numeric    \tab pvalue from LS\cr
##'        [,13] \tab LS_BH.Q               \tab numeric    \tab FDR from JTK\cr
##'        [,14] \tab LS_period             \tab numeric    \tab period from LS\cr
##'        [,15] \tab LS_adjphase           \tab numeric    \tab adjusted phase from LS\cr
##'        [,16] \tab LS_amplitude          \tab numeric    \tab amplitude from LS\cr
##'        [,17] \tab meta2d_pvalue         \tab numeric    \tab integrated pvalue\cr
##'        [,18] \tab meta2d_BH.Q           \tab numeric    \tab FDR based on integrated pvalue\cr
##'        [,19] \tab meta2d_period         \tab numeric    \tab averaged period of three methods\cr
##'        [,20] \tab meta2d_phase          \tab numeric    \tab integrated phase\cr
##'        [,21] \tab meta2d_Base           \tab numeric    \tab baseline value given by meta2d\cr
##'        [,22] \tab meta2d_AMP            \tab numeric    \tab amplitude given by meta2d\cr
##'        [,23] \tab meta2d_rAMP           \tab numeric    \tab relative amplitude\cr
##'     [,24:71] \tab CT18 to CT65          \tab numeric    \tab sampling time point
##'   }
##' @source Hughes M. E., et al. (2009). Harmonics of circadian gene
##'   transcription in mammals. \emph{PLoS Genet}, \bold{5(4)}, e1000442.
NULL

##' @name cycHumanBloodData
##' @title cycHumanBloodData
##' @description This data set lists time-series profiles of 10 transcripts
##'   sampled from  multiple individuals under different sleep conditions.
##' @docType data
##' @usage cycHumanBloodData
##' @format A dataframe containing 439 columns (column 1 = transcript name,
##'   column 2 to 439 = samples from individuals at different time points and
##'   sleep conditions).
##' @source Moller-Levet C. S., et al. (2013). Effects of insufficient sleep on
##'   circadian rhythmicity and expression amplitude of the human blood
##'   transcriptome. \emph{Proc Natl Acad Sci U S A}, \bold{110(12)},
##'   E1132--1141.
NULL

##' @name cycHumanBloodDesign
##' @title cycHumanBloodDesign
##' @description This data set describes individual information, sleep
##'   condition and sampling time corresponding to each sample in
##'   'cycHumanBloodData'.
##' @docType data
##' @usage cycHumanBloodDesign
##' @format A dataframe containing 4 columns described as below:
##'   \tabular{rlll}{
##'     [,1] \tab sample_library  \tab character  \tab sample ID\cr
##'     [,2] \tab subject         \tab character  \tab individual ID\cr
##'     [,3] \tab group           \tab character  \tab sleep condition\cr
##'     [,4] \tab time_hoursawake \tab numeric    \tab hours after awake
##'   }
##' @source Moller-Levet C. S., et al. (2013). Effects of insufficient sleep on
##'   circadian rhythmicity and expression amplitude of the human blood
##'   transcriptome. \emph{Proc Natl Acad Sci U S A}, \bold{110(12)},
##'   E1132--1141.
NULL
