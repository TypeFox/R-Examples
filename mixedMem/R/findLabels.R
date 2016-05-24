#' Mixed Membership Post-Processing
#'    
#' \code{findLabels} finds the optimal permutation of labels that minimizes 
#' the weighted squared difference between the arrays of subpopulation parameters from a fitted mixed membership
#' model, \eqn{\theta} and a given comparison model. 
#' 
#' 
#' Mixed Membership models are invariant to permutations of the sub-population labels; swapping the names of each sub-population yields an equivalent model. 
#' The ordering of the labels in a fitted model is dependent on the initialization points of the variational EM algorithim. The function \code{findLabels} selects a
#' permutation of the sub-population labels that best matches a given comparison model by minimizing the weighted squared difference between the 
#' \eqn{\theta} arrays. The weights are determined by the relative frequencies of each group.  
#' 
#' \eqn{Loss = \sum_j \sum_k \alpha_k/\alpha_0 [\sum_v (\hat\theta_{k,v} - \theta_{k,v})^2]}
#' where \eqn{\alpha_0 = \sum_k \alpha_k}
#' 
#' If K, number of sub-populations, is small, the method searches through all K! permutations of the sub-population labels and 
#' select the permutation which minimizes the loss. If K is large, a greedy algorithim can be used instead. This
#' algorithm selects the best match for each fitted sub-population starting with the group with the largest fitted 
#' relative frequency.
#'  
#' @param model the fitted \code{mixedMemModel} object.
#' @param comparison an array of the same dimensions as model$theta which contains the subpopulation parameters from another model.
#'  \code{findLabels} will return a permutation of the labels of \code{model} which match to \code{comparison} most closely.
#' @param exhaustive a boolean for whether an exhaustive search should be performed. If false, a greedy algorithim is used instead.
#' @return \code{findLabels} returns a list with two objects: \code{perm} and \code{loss}. \code{perm} is the optimal permutation of the labels with respect to the squared error loss.
#' \code{loss} is the calculated value of the weighted squared error loss (shown above) for the optimal permutation.
#' @seealso permuteLabels
#' @examples
#' 
#' \dontrun{
#' # See mixedMemModel documentation for how to generate data and instantiate a mixedMemModel object
#' # After the data as been generated, we initialize the array of sub-population parameters (theta) 
#' # according to a permutation of the true labeling
#' set.seed(123)
#' perm = sample.int(K, size = K, replace = FALSE)
#' theta1 = theta_truth[,perm,]
#' test_model <- mixedMemModel(Total = Total, J = J,Rj = Rj, Nijr= Nijr, K = K, Vj = Vj,dist = dist,
#'  obs = obs, alpha = alpha, theta = theta1)
#' out <- mmVarFit(test_model)
#' opt.perm <- findLabels(out, theta_truth)
#' opt.perm
#' 
#' # produce mixedMemModel object with sub-population labels permuted to best match
#' # the comparison model
#' out = permuteLabels(out, opt.perm$perm)
#' }
#' @export
findLabels = function(model, comparison,  exhaustive = TRUE)
{
    fitted.set = model$theta
    K = model$K
 
    optimal.perm = c(1:K)
    weight = model$alpha/sum(model$alpha)
    loss = 0
    if(exhaustive)
    {
    perms = gtools::permutations(K,K)
    loss = sum((fitted.set[,perms[1,],]-comparison[,c(1:K),])^2)
    
    for(i in 2:factorial(K))
      {
        if(dim(comparison)[1]==1) {
          diff = (fitted.set[,perms[i,],]-comparison[,c(1:K),])^2
        } else {
          diff = aperm((fitted.set[,perms[i,],]-comparison)^2, c(2,1,3))
        }
        loss.i = sum(diff*weight)
        if(loss.i<loss)
        {
          loss = loss.i
          optimal.perm = perms[i,] 
        }
      }
    } else {
      priority = order(weight, decreasing = T)
      selected = c()
      loss.j = rep(0,K)
      for(i in 1:K)
      {
        for(j in 1:K)
        {
          loss.j[j] = sum((fitted.set[,priority[i],] - comparison[,j,])^2)*weight[i] 
        }
        loss.j[selected] = max(loss.j)+1
        selected = c(selected, which.min(loss.j))
        optimal.perm[priority[i]] = which.min(loss.j)
        loss = loss + min(loss.j)
      }
    }
    return(list(perm = optimal.perm, loss = loss))
}

#' Mixed Membership Post-Processing
#' 
#' Mixed Membership models are invariant to permutations of the sub-population labels; swapping the names of each sub-population yields an equivalent model. 
#' The ordering of the labels in a fitted model is dependent on the initialization points of the variational EM algorithim.
#' The \code{permuteLabels} function returns a \code{mixedMemModel} object where the labels (for \eqn{\theta}, \eqn{\phi}, \eqn{\delta} and \eqn{\alpha}) have been permuted
#' according a given permutation of the integers 1 through K. The \code{findLabels} function can be used to find a permutation of the labels which
#' most closely matches another fitted model. 
#' 
#' @param model a fitted \code{mixedMemModel} object which will be relabeled.
#' @param perm a vector of length K with integers 1:K. This is the permutation by which to relabel the \code{mixedMemModel} object such that
#' group i in the returned mixedMemModel object corresponds to group \code{perm}[i] from the input mixedMemModel object.
#' @return \code{permuteLabels} returns a \code{mixedMemModel} object such that
#' group i in the returned \code{mixedMemModel} object corresponds to group perm[i] from the input \code{mixedMemModel} object
#' @seealso findLabels
#' @export
permuteLabels = function(model, perm)
{
  if(length(perm)!= model$K)
  {stop("Error: perm  must be of length model$K")}
  
  out = model
  
  out$alpha = out$alpha[perm]
  names(out$alpha) <- names(model$alpha)
  for(j in 1:out$J)
  {
    out$theta[j,,] = out$theta[j,perm,]
  }
  dimnames(out$theta) <- dimnames(model$theta)
  
  out$phi = out$phi[,perm]
  for(i in 1:out$Total)
  {
  for(j in 1:out$J)
  {
    for(r in 1:out$Rj[j])
    {
      for(n in 1:out$Nijr[i,j,r])
      {
        out$delta[i,j,r,n,] = out$delta[i,j,r,n,perm]
      }
    }
  }
  }
  
  dimnames(out$delta) <- dimnames(model$delta)
  dimnames(out$phi) <- dimnames(model$phi)
  
  return(out)
}