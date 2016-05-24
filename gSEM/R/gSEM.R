##' Semi-supervised Generalized Structure Equation Modeling
##'
##' Conducts a semi-gSEM statistical analysis (semi-supervised generalized structural equation modeling) on a data frame of coincident observations of multiple predictive or intermediate variables and a final continuous, outcome variable, via two functions sgSEMp1() and sgSEMp2(), representing fittings based on two statistical principles. 
##' Principle 1 determines all sensible univariate relationships in the spirit of the Markovian process. The relationship between each pair of variables, including predictors and the final outcome variable, is determined with the Markovian property that the value of the current predictor is sufficient in relating to the next level variable, i.e., the relationship is independent of the specific value of the preceding-level variables to the current predictor, given the current value. 
##' Principle 2 resembles the multiple regression principle in the way multiple predictors are considered simultaneously. Specifically, the relationship of the first-level predictors (such as Time and irradiance etc) to the outcome variable (such as, module degradation or yellowing)  is fit by a supervised additive model. Then each significant intermediate variable is taken as the new outcome variable and the other variables (except the final outcome variable) as the predictors in investigating the next-level multivariate relationship by a supervised additive model. This fitting process is continued until all sensible models are investigated.
##' 
##' @import DiagrammeR
##' @import knitr
##' @import MASS
##' @import htmlwidgets
##' @references 1. Bruckman, Laura S., Nicholas R. Wheeler, Junheng Ma, Ethan Wang, Carl K. Wang, Ivan Chou, Jiayang Sun, and Roger H. French. "Statistical and Domain Analytics Applied to PV Module Lifetime and Degradation Science." IEEE Access 1 (2013): 384-403. doi:10.1109/ACCESS.2013.2267611
##'           
##' 2. Bruckman, Laura S., Nicholas R. Wheeler, Ian V. Kidd, Jiayang Sun, and Roger H. French. "Photovoltaic Lifetime and Degradation Science Statistical Pathway Development: Acrylic Degradation." In SPIE Solar Energy+ Technology, 8825:88250D-8. International Society for Optics and Photonics, 2013. doi:10.1117/12.2024717
##' @name gSEM
##' @docType package
##' @seealso sgSEMp1() for implementing principle 1 and sgSEMp2() for implementing principle 2.

NULL
