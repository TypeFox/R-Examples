#' Haberman survival data set
#' 
#' @description Training data for the Haberman dataset.
#' 
#' @details This data set contains cases from a study that was conducted between 1958 and 1970 
#'     at the University of Chicago's Billings Hospital on the survival of patients who had undergone 
#'     surgery for breast cancer. The task is to determine if the patient survived 5 years or longer 
#'     (positive) or if the patient died within 5 year (negative)
#'     
#' @format A keel class with 306 instances, 3 variables (without the target variable) and 2 values for the target variable.
#'     Three fuzzy labels for each numerical variable are defined.
#'     
#' @source Haberman, S. J. (1976). Generalized Residuals for Log-Linear Models, 
#'     Proceedings of the 9th International Biometrics Conference, Boston, pp. 104-122. 
#' @source Landwehr, J. M., Pregibon, D., and Shoemaker, A. C. (1984), Graphical Models for 
#'     Assessing Logistic Regression Models (with discussion), Journal of the American Statistical 
#'     Association 79: 61-83. 
#' @source Lo, W.-D. (1993). Logistic Regression Trees, PhD thesis, Department of Statistics,
#'     University of Wisconsin, Madison, WI. 
#'     
#' @examples 
#'     habermanTra$data 
#'     habermanTra$fuzzySets
#'     
#' @docType data
#' @name habermanTra
NULL



#' Haberman survival data set
#' 
#' @description Test data for the Haberman dataset.
#' 
#' @details This data set contains cases from a study that was conducted between 1958 and 1970 
#'     at the University of Chicago's Billings Hospital on the survival of patients who had undergone 
#'     surgery for breast cancer. The task is to determine if the patient survived 5 years or longer 
#'     (positive) or if the patient died within 5 year (negative)
#'     
#' @format A keel class with 62 instances, 3 variables (without the target variable) and 2 values for the target variable.
#'     Three fuzzy labels for each numerical variable are defined.
#'     
#' @source Haberman, S. J. (1976). Generalized Residuals for Log-Linear Models, 
#'     Proceedings of the 9th International Biometrics Conference, Boston, pp. 104-122. 
#' @source Landwehr, J. M., Pregibon, D., and Shoemaker, A. C. (1984), Graphical Models for 
#'     Assessing Logistic Regression Models (with discussion), Journal of the American Statistical 
#'     Association 79: 61-83. 
#' @source Lo, W.-D. (1993). Logistic Regression Trees, PhD thesis, Department of Statistics,
#'     University of Wisconsin, Madison, WI. 
#'     
#' @examples 
#'     habermanTra$data 
#'     habermanTra$fuzzySets
#'   @docType data
#'   @name habermanTst
NULL





#' Car evaluation dataset
#' 
#' @description Training data for the car dataset
#' 
#' @details Car Evaluation Database was derived from a simple hierarchical decision model. 
#'      The model evaluates cars according to six input attributes: buying, maint, doors, 
#'      persons, lug_boot, safety. 
#'      
#' @format A keel class with 1382 instances, 6 variables (without the target variable) and 4 values for the target Variable.
#'     Three labels for each variable are defined.
#'     
#' @source M. Bohanec and V. Rajkovic: Knowledge acquisition and explanation for multi-attribute 
#'     decision making. In 8th Intl Workshop on Expert Systems and their Applications, Avignon, 
#'     France. pages 59-78, 1988. 
#'     
#' @source B. Zupan, M. Bohanec, I. Bratko, J. Demsar: Machine learning by function decomposition. 
#'     ICML-97, Nashville, TN. 1997 (to appear).
#'  
#' @examples 
#'    carTra$data
#'    carTra$atributeNames
#'    
#' @docType data
#' @name carTra
NULL




#' Car evaluation dataset
#' 
#' @description Test data for the car dataset
#' 
#' @details Car Evaluation Database was derived from a simple hierarchical decision model. 
#'      The model evaluates cars according to six input attributes: buying, maint, doors, 
#'      persons, lug_boot, safety. 
#'      
#' @format A keel class with 346 instances, 6 variables (without the target variable) and 4 values for the target Variable.
#'     Three labels for each variable are defined.
#'     
#' @source M. Bohanec and V. Rajkovic: Knowledge acquisition and explanation for multi-attribute 
#'     decision making. In 8th Intl Workshop on Expert Systems and their Applications, Avignon, 
#'     France. pages 59-78, 1988. 
#'     
#' @source B. Zupan, M. Bohanec, I. Bratko, J. Demsar: Machine learning by function decomposition. 
#'     ICML-97, Nashville, TN. 1997 (to appear).
#'  
#' @examples 
#'    carTra$data
#'    carTra$atributeNames
#'    
#'    @docType data
#'    @name carTst
NULL






#' German Credit data set
#' 
#' @description Training data for the german dataset
#' 
#' @details A numerical version of the Statlog German Credit Data data set. 
#'     Here, the task is to clasify customers as good (1) or bad (2), 
#'     depending on 20 features about them and their bancary accounts.
#' 
#' @format A keel class with 800 instances, 20 variables (without the target variable)
#'     and 2 values for the target class.
#'     
#'  @source \url{http://sci2s.ugr.es/keel/dataset.php?cod=88}
#'  
#'  @examples 
#'      germanTra$data
#'      
#'      @docType data
#'      @name germanTra
NULL





#' German Credit data set
#' 
#' @description Test data for the german dataset
#' 
#' @details A numerical version of the Statlog German Credit Data data set. 
#'     Here, the task is to clasify customers as good (1) or bad (2), 
#'     depending on 20 features about them and their bancary accounts.
#' 
#' @format A keel class with 200 instances, 20 variables (without the target variable)
#'     and 2 values for the target class.
#'     
#'  @source \url{http://sci2s.ugr.es/keel/dataset.php?cod=88}
#'  
#'  @examples 
#'      germanTra$data
#'      
#'      @docType data
#'      @name germanTst
NULL