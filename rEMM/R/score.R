#######################################################################
# rEMM - Extensible Markov Model (EMM) for Data Stream Clustering in R
# Copyrigth (C) 2011 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


### exponential
#.simil_weight <- function(d, th) {
#  w <- .5^(d/th - 1)
#  w[d<=th] <-1
#  w
#}

### use sigmoid curve instead
.simil_weight <- function(d, th) {
  w <- 1-plogis(d/th, 1.5, .2)
  #w <- 1- 1/(1+exp(-(d-1.5)/.2))
  w
}

## linear
#.simil_weight <- function(d, th) {
#  d <- d/th
#  w <- 2-d
#  w[d<=1] <- 1
#  w[d>2] <- 0
#  w
#}

## does newdata come from the EMM?
setMethod("score", signature(x = "EMM", newdata = "numeric"),
          function(x, newdata, method=NULL, 
                   match_cluster="exact", prior = TRUE, normalize=TRUE, 
                   initial_transition = FALSE, threshold = NA) 
            score(x, as.matrix(rbind(newdata)), method, 
                  match_cluster, prior, normalize, initial_transition, threshold)
)

setMethod("score", signature(x = "EMM", newdata = "data.frame"),
          function(x, newdata, method=NULL, 
                   match_cluster="exact", prior = TRUE, normalize = TRUE, 
                   initial_transition = FALSE, threshold = NA) 
            score(x, as.matrix(newdata), method, 
                  match_cluster, prior, normalize, initial_transition, threshold)
)

setMethod("score", signature(x = "EMM", newdata = "matrix"),
          function(x, newdata, method = c(
            "product", 
            "log_sum", 
            "sum", 
            "log_odds", 
            "supported_transitions",
            "supported_states",
            "sum_transitions",
            "log_loss",
            "likelihood",
            "log_likelihood",
            "AIC"
          ), 
                   match_cluster = "exact", 
                   prior = TRUE,
                   normalize = TRUE,
                   initial_transition = FALSE,
                   threshold = NA) {
            
            method <- match.arg(method)
            
            if(!is.numeric(match_cluster)) match_cluster <- match.arg(match_cluster, c("exact", "nn", "weighted"))

            ### deal with empty models
            if(nclusters(x) <1) return(if(method=="log_odds") -Inf else 0)

            if(method == "supported_transitions") {
              ###if(prior) warning("prior has no effect on supported transitions!")
              prior <- FALSE
            }

            ### FIXME: Don't know about log_odds!
            if(method == "log_odds") {
              log_odds <- transition_table(x, newdata, type="log_odds", 
                                           match_cluster, prior, 
                                           initial_transition)[,3]
              return(sum(log_odds, na.rm=TRUE))
            }            
            
            ### calc: tt with  count and (weighted) prob
            ###       state_weight ... (weighted) cluster membership
            ###       states ... matching states
            
            if(match_cluster != "weighted") {
              tt <- transition_table(x, newdata, type="probability", 
                                       match_cluster=match_cluster, prior=prior, 
                                       initial_transition=initial_transition)
              tt[["count"]] <- transition(x, tt, type="count", prior=prior)
              
              ### clusters
              states <- c(tt[,"from"], tail(tt[,"to"], n=1))
              state_weight <- as.integer(!is.na(states))
              
            }else{
              tt <- transition_table(x, newdata, type="probability",
                                         match_cluster="nn", prior=prior,
                                         initial_transition=initial_transition)
              tt[["count"]] <- transition(x, tt, type="count", prior=prior)
              
              states <- c(tt[,"from"], tail(tt[,"to"], n=1))
            
              n <- length(states)
              state_weight<- numeric(n)
              for(i in 1:n) {
                state_weight[i] <- as.numeric(.simil_weight(
                  dist(
                    newdata[i, , drop=FALSE],
                    cluster_centers(x)[states[i], , drop=FALSE],
                    measure=x@measure),
                  #x@threshold
                  x@tnn_d$var_thresholds[states[i]]
                ))
              }
       
              ### weight is the product of source and target weight
              weight <- state_weight[-n] * state_weight[-1] 
              tt[,"probability"] <- weight*tt[,"probability"]
              tt[,"count"] <- weight*tt[,"count"]
            }
            
            ### remove states and transitions with count < threshold
            if(!is.na(threshold)) {
              state_weight[cluster_counts(x)[states]<threshold] <- NA  
              rem <- tt[, "count"]<threshold
              tt[rem, "probability"] <- NA 
              tt[rem, "count"] <- NA 
              
            }
            
            if(method == "supported_states"){
              if(normalize) return(sum(state_weight, na.rm=TRUE)/length(state_weight))
              else return(sum(state_weight, na.rm=TRUE))
            }

            if(method == "supported_transitions") {
              ###if(prior) warning("prior has no effect on supported transitions!")
              if(match_cluster !="weighted") {
                if(normalize) return(sum(tt[["count"]]>0, na.rm=TRUE)/length(tt[["count"]]))
                else return(sum(tt[["count"]]>0, na.rm=TRUE))
              }else{
                if(normalize) return(sum((tt[["count"]]>0)*weight, na.rm=TRUE)/length(tt[["count"]]))
                else return(sum((tt[["count"]]>0)*weight, na.rm=TRUE))
              }
            }
            
            if(method == "sum_transitions") {
              if(normalize) return(sum(tt[["count"]], na.rm=TRUE)/
                                     sum(smc_countMatrix(x@tracds_d$mm))/nrow(tt))
              else return(sum(tt[["count"]], na.rm=TRUE)/
                            sum(smc_countMatrix(x@tracds_d$mm)))
            }
            
            if(method == "log_loss") {
              return(-sum(log2(tt[["probability"]]))/length(tt[["probability"]]))
            }
            
            if(method == "likelihood") { ### this is the unnormalized product
              return(prod(tt[["probability"]]))
            }
            
            if(method == "log_likelihood") { ### this is the unnormalized product
              return(sum(log(tt[["probability"]])))
            }
            
            if(method == "AIC") { 
              ### AICc = 2k - 2 ln(L) * 2k(k+1)/(n-k-1)  
              ### minimum AIC for model selection
              
              ### we use supported transitions a the Likelihood
              L <- sum(tt[["count"]]>0, na.rm=TRUE)/length(tt[["count"]])
              ### complexity is number of transitions
              k <- ntransitions(x, threshold=threshold)
              ### number of transitions
              n <- length(tt[["count"]])
              
              return(2*k-2*log(L) * 2*k*(k-1)/(n-k-1))
            }
            
            if(method == "product") {
              if(normalize) return(prod(tt[["probability"]])^(1/length(tt[["probability"]])))
              else return(prod(tt[["probability"]]))
            }
            
            if(method == "log_sum") {
              if(normalize) return(sum(log10(tt[["probability"]]))/length(tt[["probability"]]))
              else return(sum(log10(tt[["probability"]])))
            }
            
            if(method == "sum") {
              if(normalize) return(sum(tt[["probability"]])/length(tt[["probability"]]))
              else return(sum(tt[["probability"]]))
            }
       
            stop("Unknown method!")
            
          }
)


### score two models
setMethod("score", signature(x = "EMM", newdata = "EMM"),
          function(x, newdata, method=c("product", "log_sum", "sum",
                                        "supported_transitions"), 
                   match_cluster="exact", prior = TRUE, 
                   initial_transition = FALSE) {
            
            method <- match.arg(method)
            
            
            ### find transitions in newdata
            trans <- transitions(newdata)	
            
            ### match states in newdata to x
            cl <- find_clusters(x,cluster_centers(newdata), 
                                match_cluster=match_cluster)
            
            ### translate to states in x
            cl <- cbind(cl[as.integer(trans[,1])], cl[as.integer(trans[,2])])
            
            ### FIXME: add weighted versions. What weights should we use?
            
            if(method=="product") return(
              prod(transition(x, cl, type="probability", 
                              prior=prior)^(1/nrow(cl))))
            
            if(method=="sum") return(
              sum(transition(x, cl, type="probability", 
                             prior=prior)*(1/nrow(cl))))
            
            if(method=="log_sum") return(
              sum(log(transition(x, cl, type="probability", 
                                 prior=prior))*(1/nrow(cl))))
            
            if(method=="supported_transitions") return(
              (nrow(cl)-sum(transition(x, cl, type="count", 
                                       prior=FALSE)==0))/nrow(cl))
            
          }
)

