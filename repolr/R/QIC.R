QIC <-
function(object, digits = 3){

 # independence model
 mod.indepen <- update(object, corr.mod = "independence")

 # quasilikelihood, QICu and QIC
 qlike <- object$y * object$linear.predictor + log(1 - object$fitted.values)
 QICu <- - 2 * sum(qlike) + 2 * length(object$coeff)
 QIC <- - 2 * sum(qlike) + 2 * sum(diag(solve(mod.indepen$naive.var) %*% object$robust.var))

 # output
 list(QLike = round(sum(qlike), digits = digits), QIC = round(QIC, digits = digits), 
                      QICu = round(QICu, digits = digits))
}
