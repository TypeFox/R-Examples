#' @details
#' Bolded tasks, followed by their respective models, are itemized below.  
#' 
#' \describe{
#'  \item{\strong{Bandit}}{2-Armed Bandit (Rescorla-Wagner (delta)) --- \link{bandit2arm}}
#'  \item{\strong{Delay Discounting}}{Constant Sensitivity --- \link{dd_cs} \cr
#'                                    Exponential          --- \link{dd_exp} \cr
#'                                    Hyperbolic           --- \link{dd_hyperbolic}}
#'  \item{\strong{Orthogonalized Go/Nogo}}{RW + Noise                                   --- \link{gng_m1} \cr
#'                                         RW + Noise + Bias                            --- \link{gng_m2} \cr
#'                                         RW + Noise + Bias + Pavlovian Bias           --- \link{gng_m3} \cr
#'                                         RW(modified) + Noise + Bias + Pavlovian Bias --- \link{gng_m4}}
#'  \item{\strong{Iowa Gambling}}{Prospect Valence Learning-DecayRI --- \link{igt_pvl_decay} \cr
#'                                Prospect Valence Learning-Delta   --- \link{igt_pvl_delta} \cr
#'                                Value-Plus_Perseverance           --- \link{igt_vpp}}
#'  \item{\strong{Probabilistic Reversal Learning}}{Fictitious Update --- \link{prl_fictitious}}
#'  \item{\strong{Risk Aversion}}{Prospect Theory --- \link{ra_prospect}} \cr
#' }
#' 
#' 
#' For tutorials and further readings, visit : \url{http://u.osu.edu/ccsl/codedata/hbayesdm/}.
#' 
#' Please cite as: 
#' Ahn, W.-Y., Haines, N., & Zhang, L. (in preparation) Hierarchical Bayesian Modeling of Decision-Making Tasks with hBayesDM
#'
#' @author
#' Woo-Young Ahn \email{ahn.280@@osu.edu}
#'
#' Nathaniel Haines \email{haines.175@@osu.edu}
#'
#' Lei Zhang \email{bnuzhanglei2008@@gmail.com}
#'
#'
"_PACKAGE"
