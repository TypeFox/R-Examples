#################################
# Creates a new MLPnet object
#############################
# OJO : FALTA completar la entrada de un deltae custom para aceptar como parametro la función custom

newff <- function (n.neurons, learning.rate.global, momentum.global=NA, error.criterium="LMS", Stao=NA, hidden.layer="tansig", output.layer="purelin", method="ADAPTgdwm") {

   net <- list( layers=list(), neurons=list(), input=as.double(numeric(n.neurons[1])), output=as.double(numeric(n.neurons[length(n.neurons)])), target=as.double(numeric(n.neurons[length(n.neurons)])), deltaE=list(fname=as.integer(0),f=function(){},Stao=as.double(NA)), other.elements=list() )

   if (length(n.neurons)<3) {
       stop("You should enter a vector containing the number of input neurons, the number of neurons of each hidden layer and the number of outputs.")
   }

   possible.activation.functions <- c("custom","tansig","sigmoid","purelin","hardlim")

   if ( is.na(hidden.activation.function.choice <- pmatch(hidden.layer,possible.activation.functions)) ) {
       stop("You should use a correct activation function for the hidden layers.")
   } 
   if ( is.na(output.activation.function.choice <- pmatch(output.layer,possible.activation.functions)) ) {
       stop("You should use a correct activation function for the output layer.")
   }

   possible.methods <- c("ADAPTgdwm","ADAPTgd","BATCHgd","BATCHgdwm")

   if ( is.na(method.choice <- pmatch(method,possible.methods)) ) {
       stop("You should use a correct training method: ADAPTgdwm, ADAPTgd, BATCHgdwm, BATCHgd. Read the help files.")
   } 

   layers.last.neuron  <- cumsum(n.neurons)
   layers.first.neuron <- c(1,1+layers.last.neuron[1:(length(layers.last.neuron)-1)])[-c(1)]-n.neurons[1]
   layers.last.neuron  <- layers.last.neuron[-c(1)]-n.neurons[1]

   net$layers[[1]] <- -c(1:n.neurons[1])
   for ( ind.layer in 1:length(layers.last.neuron) ) {
      net$layers[[ind.layer+1]] <- layers.first.neuron[ind.layer]:layers.last.neuron[ind.layer]
   }

   input.links  <- net$layers[1:(length(net$layers)-1)]
   output.links <- list()
   for ( ind.layer in 2:length(layers.last.neuron)) {
      output.links[[ind.layer-1]] <- layers.first.neuron[ind.layer]:layers.last.neuron[ind.layer]
   }
   output.links[[length(layers.last.neuron)]] <- NA

   for (ind.layer in 2:length(n.neurons)) {
      if (ind.layer == length(n.neurons)) {
         this.neuron.type="output"
         this.neuron.activation.function.choice <- output.activation.function.choice
      } else {
         this.neuron.type="hidden"
         this.neuron.activation.function.choice <- hidden.activation.function.choice
      }
      if (method == "ADAPTgd"  ) {
         method.dep.variables                      <- list()
         method.dep.variables$delta                <- as.double(0)
         method.dep.variables$learning.rate        <- as.double(learning.rate.global)
      } else if (method == "ADAPTgdwm") {
         method.dep.variables                      <- list()
         method.dep.variables$delta                <- as.double(0)
         method.dep.variables$learning.rate        <- as.double(learning.rate.global)
         method.dep.variables$momentum             <- as.double(momentum.global)
         method.dep.variables$former.weight.change <- as.double(numeric(n.neurons[ind.layer-1]))
         method.dep.variables$former.bias.change   <- as.double(0)
      } else if (method == "BATCHgd"  ) {
         method.dep.variables                      <- list()
         method.dep.variables$delta                <- as.double(0)
         method.dep.variables$learning.rate        <- as.double(learning.rate.global)
         method.dep.variables$sum.delta.x          <- as.double(numeric(n.neurons[ind.layer-1]))
         method.dep.variables$sum.delta.bias       <- as.double(0)
      } else if (method == "BATCHgdwm") {
         method.dep.variables                      <- list()
         method.dep.variables$delta                <- as.double(0)
         method.dep.variables$learning.rate        <- as.double(learning.rate.global)
         method.dep.variables$sum.delta.x          <- as.double(numeric(n.neurons[ind.layer-1]))
         method.dep.variables$sum.delta.bias       <- as.double(0)
         method.dep.variables$momentum             <- as.double(momentum.global)
         method.dep.variables$former.weight.change <- as.double(numeric(n.neurons[ind.layer-1]))
         method.dep.variables$former.bias.change   <- as.double(0)
      }

      for ( ind.MLPneuron.relative in 1:length(net$layers[[ind.layer]]) ) {
         ind.MLPneuron <- net$layers[[ind.layer]][[ind.MLPneuron.relative]]
         net$neurons[[ind.MLPneuron]] <- init.MLPneuron(id=ind.MLPneuron,type=this.neuron.type, activation.function=as.integer(this.neuron.activation.function.choice-1),output.links=output.links[[ind.layer-1]], output.aims=rep(ind.MLPneuron.relative,length(output.links[[ind.layer-1]])), input.links=input.links[[ind.layer-1]],weights=numeric(n.neurons[ind.layer-1]), bias=0, method, method.dep.variables )
      }
   }


   if (error.criterium == "LMS" ) {
      net$deltaE$fname <- as.integer(0)      # LMS_NAME  0
      net$deltaE$f <- deltaE.LMS
   } else if (error.criterium == "LMLS") {
      net$deltaE$fname <- as.integer(1)      # LMLS_NAME 1
      net$deltaE$f <- deltaE.LMLS
   } else if (error.criterium == "TAO") {   
      net$deltaE$fname <- as.integer(2)      # TAO_NAME  2
      net$deltaE$f <- deltaE.TAO 
      if (missing(Stao)){ 
         stop("You should enter the Stao value")
      } else {
         net$deltaE$Stao <-as.double(Stao)
      }
   } else {
      stop("You should enter either: \"LMS\", \"LMSL\" or \"TAO\". ")
   }

   class(net)              <- "MLPnet"
   net <- random.init.MLPnet(net)
return(net)
}
#################################
# Creates individual neurons
#########################
init.MLPneuron   <- function(id,type,activation.function,output.links, output.aims, input.links, weights, bias, method, method.dep.variables) {
aux <- select.activation.function(activation.function)
neuron                      <- list()
neuron$id                   <- as.integer(id)
neuron$type                 <- as.character(type)
neuron$activation.function  <- activation.function
neuron$output.links         <- as.integer(output.links)
neuron$output.aims          <- as.integer(output.aims)
neuron$input.links          <- as.integer(input.links)
neuron$weights              <- as.double(weights)
neuron$bias                 <- as.double(bias)
neuron$v0                   <- as.double(0)
neuron$v1                   <- as.double(0)
neuron$f0                   <- aux$f0
neuron$f1                   <- aux$f1
neuron$method               <- as.character(method)
neuron$method.dep.variables <- method.dep.variables

class(neuron) <- "neuron"
return(neuron)
}
#########################################
# Initialize the neuron bias and weights with random values according to the book:
# Neural Networks. A comprehensive foundation. 2nd Edition.
# Author: Simon Haykin.
# pages = 182, 183, 184.
#################################
random.init.MLPneuron <- function(net.number.weights, neuron) {
   extreme        <- sqrt(3/net.number.weights)
   n.weights      <- length(neuron$weights)
   neuron$weights <- runif(n.weights,min=-extreme,max=extreme)
   neuron$bias    <- runif(1,min=-extreme,max=extreme)
   return(neuron)
}
#################################################
# Runs random.init.MLPneuron upon each neuron.
###########################################
random.init.MLPnet <- function(net) {
   net.number.weights <- length(net$neurons)          #number of bias terms
   for (ind.MLPneuron in 1:length(net$neurons)) {
          net.number.weights <- net.number.weights + length(net$neurons[[ind.MLPneuron]]$weights)
       }

   for ( i in 1:length(net$neurons)) { 
      net$neurons[[i]] <- random.init.MLPneuron(net.number.weights,net$neurons[[i]] )
   }
return(net)
}

#########################################
# A simple function to bestow the neuron with the appropriate 
select.activation.function <- function(activation.function) {
   f0 <- NA
   f1 <- NA

# a.tansig  : 1/tanh(2/3)
# b.tansig  : 2/3
# a.sigmoid : 1.0

   if (activation.function == 1 ) { # TANSIG
     f0 <- function (v) {
                          a.tansig   <- 1.715904708575539
                          b.tansig   <- 0.6666666666666667
                          return ( a.tansig * tanh( v * b.tansig ) )
                        }
     f1 <- function (v) {         # realmente usaremos f1= b.tansig/a.tansig*(a.tansig-f0)*(a.tansig+f0)
                          a.tansig   <- 1.715904708575539
                          b.tansig   <- 0.6666666666666667
                          return( a.tansig * b.tansig * (1-tanh( v * b.tansig )^2)  )
                        } 
    } else if (activation.function == 2 ) { # SIGMOID
     f0 <- function (v) {
                          a.sigmoid  <- 1
                          return( 1/(1+exp(- a.sigmoid * v)) )
                        }
     f1 <- function (v) {           # realmente usaremos f1=a.sigmoid*f0*(1-f0)
                          a.sigmoid  <- 1
                          return ( a.sigmoid * exp(- a.sigmoid * v) / (1+exp(- a.sigmoid * v))^2 )
                        } 
   } else if (activation.function == 3 ) { # PURELIN
     f0 <- function (v) {
                          return( v )  
                        }
     f1 <- function (v) {
                          return( 1 ) 
                        }
   } else if (activation.function == 4 ) { # HARDLIM
     f0 <- function (v) {
                          if (v>=0) { return(1) } else { return(0) }
                        }
     f1 <- function (v) {
                          return ( NA )
                        }
   }

   return(list(f0=f0,f1=f1))
}

##############################################################
# Manually set the learning rate and momentum for each neuron
##############################################################
# deprecated. do not use. does not work
#set.learning.rate.and.momentum <- function(net, learning.rate, momentum) {
#   for (i in 1:length(net$neurons)) {
#      net$neurons[[i]]$learning.rate <- learning.rate
#      net$neurons[[i]]$momentum <- momentum
#   }
#   return(net)
#}
#
#

