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


## predict next n states using P^n
setMethod("predict", signature(object = "TRACDS"),
	function(object, current_state=NULL, n=1,
		probabilities = FALSE, 
		randomized = FALSE,
		prior = FALSE) {

		## probabilistic max with random tie breaking
		.prob_max <- function(x) {
			m <- which(x==max(x))
			if(length(m)>1) m <- sample(m,1)
			m
		}

		## randomized
		.randomized <- function(x)
			sample((1:length(x))[x>0], 1, prob=x[x>0])


		
		if(is.null(current_state)) 
		    current_state <- current_state(object)
		else 
		    current_state <- as.character(current_state)

		current_state_i <- which(states(object) == current_state)
		
		## check is state exists!
		if(!is.element(current_state, states(object))) 
		    stop("State does not exist")


		P <- transition_matrix(object, prior=prior)
		## calculate P^n
		if(n>1) for(i in 1:(n-1)) P <- P%*%P

		prob <- P[current_state_i,]
		
		## create result
		if(probabilities) return(prob)
		if(randomized) return(states(object)[.randomized(prob)])
		
		return(states(object)[.prob_max(prob)])
	}
)
