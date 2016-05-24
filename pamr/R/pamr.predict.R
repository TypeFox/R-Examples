pamr.predict <-  function(fit, newx, threshold, type = c("class", "posterior", "centroid", "nonzero"), 
                          prior = fit$prior,  threshold.scale = fit$
                          threshold.scale) {
  norm.cen <- fit$norm.cen
  if(!is.null(norm.cen)) {
    newx <- abs(t(scale(t(newx), center = norm.cen, scale = FALSE)))
  }
  type <- match.arg(type)
  sd <- fit$sd
  centroid.overall <- fit$centroid.overall
  centroids <- fit$centroids
  se.scale <- fit$se.scale
  delta <- scale((centroids - centroid.overall)/sd, FALSE, threshold.scale * 
                 se.scale)

  if(fit$sign.contrast=="positive"){delta <- delta*(delta>0)}
  if(fit$sign.contrast=="negative"){delta <- delta*(delta<0)}


  delta.shrunk <- scale(soft.shrink(delta, threshold), FALSE,
                        1/(  threshold.scale * se.scale))
  posid <- drop(abs(delta.shrunk) %*% rep(1, length(prior))) > 0
                
  if(!match(type, c("centroid", "nonzero"), FALSE))
    dd <- diag.disc((newx - centroid.overall)/sd, delta.shrunk, 
                    prior, posid)
  switch(type,
         class = softmax(dd),
         posterior = {
           dd <- safe.exp(dd)
           dd/drop(dd %*% rep(1, length(prior)))
         }
         ,
         centroid = centroid.overall + delta.shrunk * sd,
         nonzero = {
           nz <- drop(abs(delta.shrunk) %*% rep(1, ncol(centroids)
                                                )) > 0
           seq(nz)[nz]
         }
         )
}

safe.exp=function(x){
 xx=sign(x)*pmin(abs(x),500)
 return(exp(xx))
}

