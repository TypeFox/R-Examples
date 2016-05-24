##########################################################
#	Adaptative Gradient Descent (without momentum)
##########################################################
ADAPTgd.MLPnet <- function(net, P, T, n.epochs, n.threads=0L) {
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }
   net <- .Call("ADAPTgd_loop_MLPnet", net, t(P), t(T),as.integer(n.epochs), new.env(), PACKAGE="AMORE" )
   return(net)
}
##################################################


##########################################################
#	Adaptative Gradient Descent (with momentum)
##########################################################
ADAPTgdwm.MLPnet <- function(net,P,T, n.epochs, n.threads=0L) {
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }
   net <- .Call("ADAPTgdwm_loop_MLPnet", net, t(P), t(T),  as.integer(n.epochs), new.env(), PACKAGE="AMORE" )
   return(net)
}
##################################################


##############################################################
#	BATCHgd ( BATCH gradient descent without momentum )
##############################################################
BATCHgd.MLPnet <- function(net, P, T, n.epochs, n.threads=0L) { # Each pattern is a row of P, 
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }

#####  First Step: BATCHgd.Forward.MLPnet
   for (ind.MLPneuron in 1:length(net$neurons)) {
      net$neurons[[ind.MLPneuron]]$method.dep.variables$sum.delta.bias <- as.double(0)
      net$neurons[[ind.MLPneuron]]$method.dep.variables$sum.delta.x    <- as.double(numeric(length(net$neurons[[ind.MLPneuron]]$method.dep.variables$sum.delta.x)))
   }
   net <- .Call("BATCHgd_loop_MLPnet", net, t(P), t(T), as.integer(n.epochs), new.env(), as.integer(n.threads), PACKAGE="AMORE")
   return(net)
}
##############################################################
#	BATCHgdwm ( BATCH gradient descent with momentum )
##############################################################
BATCHgdwm.MLPnet <- function(net, P, T, n.epochs, n.threads=0L) { # Each pattern is a row of P, 
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }

##### First step: BATCHgdwm.Forward.MLPnet
   for (ind.MLPneuron in 1:length(net$neurons)) {
      net$neurons[[ind.MLPneuron]]$method.dep.variables$sum.delta.bias <- as.double(0)
      net$neurons[[ind.MLPneuron]]$method.dep.variables$sum.delta.x    <- as.double(numeric(length(net$neurons[[ind.MLPneuron]]$method.dep.variables$sum.delta.x)))
   }
   net <- .Call("BATCHgdwm_loop_MLPnet", net, t(P), t(T), as.integer(n.epochs), new.env(), as.integer(n.threads), PACKAGE="AMORE")
   return(net)
}
#######
