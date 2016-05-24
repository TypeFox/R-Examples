#' @title Basic PLS-PM algorithm
#' 
#' @description
#' Internal function. \code{get_pls_basic} is called by \code{pathmox}, 
#' \code{techmox}, \code{fix.pathmox}, \code{fix.techmox}, \code{treemox.pls}, 
#' and \code{treemox.boot}. 
#' 
#' @param DT data table
#' @param path_matrix inner design matrix
#' @param blocks blocks of manifest variables
#' @param specs list with pls algorithm specifications
#' @export
#' @keywords internal
get_pls_basic <- 
function(DT, path_matrix, blocks, specs)
{   
  ### Variable Names
  lvs.names = colnames(path_matrix)
  mvs.names = colnames(DT)
  
  # apply the selected scaling
  if (specs$scaled) {
    sd.X = sqrt((nrow(DT)-1)/nrow(DT)) * apply(DT, 2, sd)
    X = scale(DT, scale=sd.X)
  } else {
    X = scale(DT, scale = FALSE)
  }
  dimnames(X) = list(rownames(DT), mvs.names)
  
  # ==================== Stage 1: Iterative procedure ==================
  weights = get_weights(X, path_matrix, blocks, specs)
  ok_weights = test_null_weights(weights, specs)
  outer_weights = weights$w
  LV = get_scores(X, weights$W)
  xloads = cor(X, LV)
  loadings = rowSums(xloads * weights$ODM)
  
  # ============ Stage 2: Path coefficients and total effects ==========
  inner_results = get_paths(path_matrix, LV)
  inner_model = inner_results[[1]]
  Path = inner_results[[2]]
  R2 = inner_results[[3]]
  residuals <- inner_results[[4]]
  
  # ============================= Results ==============================
  model <- list(IDM = path_matrix, 
                blocks = blocks, 
                specs = specs)
  # output
  list(out.weights = outer_weights, 
       loadings = loadings, 
       scores = LV,   
       path.coefs = Path, 
       R2 = R2, 
       residuals = residuals, 
       model = model)
}
