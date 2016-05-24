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


setMethod("transition", signature(x = "TRACDS", 
		from = "matrix", to = "missing"),
	function(x, from, to, 
		type=c("probability", "counts", "log_odds"), prior = TRUE){

		to <- from[,2]
		from <- from[,1]

		transition(x, from, to, type, prior)
	}
)
	
setMethod("transition", signature(x = "TRACDS", 
          from = "data.frame", to = "missing"),
          function(x, from, to, 
                   type=c("probability", "counts", "log_odds"), prior = TRUE){
            
            to <- from[,2]
            from <- from[,1]
            
            transition(x, from, to, type, prior)
          }
)

setMethod("transition", signature(x = "TRACDS", from = "character", to =
                "character"), function(x, from, to, type=c("probability",
                        "counts", "log_odds"), prior = TRUE){ 
            type <- match.arg(type)

	    if(length(from) != length(to)) stop("vectors from and to are not of the same length!")
	    
	    ### deal with empty from/to
	    if(length(from) <1) return(numeric(0))	

            tm <- transition_matrix(x, type, prior)

            from <- match(from, states(x)) 
            to <- match(to, states(x)) 
            res <- sapply(1:length(from), FUN = function(i) tm[from[i], to[i]])
        
            ## handle missing states (NA) 
            res[is.na(res)] <- 0 
            res 
        })



setMethod("transition_matrix", signature(x = "TRACDS"),
	function(x,
		type=c("probability", "counts", "log_odds"), prior = TRUE){
		type <- match.arg(type)

		## get transition count matrix
		m <- smc_countMatrix(x@tracds_d$mm)
		
		if(prior) m <- m+1

		if(type=="counts") return(m)

		rs <- rowSums(m)
		prob <- m/rs

		## we have to handle absorbing states here (row sum is 0)
		absorbing <- which(rs==0)
		prob[absorbing,] <- 0
		for(i in absorbing) prob[i,i] <- 1

		switch(type,
			probability = prob,
			log_odds = log(prob*size(x))
		)
	}
)


setMethod("initial_transition", signature(x = "TRACDS"),
	function(x, 
		type=c("probability", "counts", "log_odds"), prior = TRUE){
		type <- match.arg(type)

		ic <- smc_initialCounts(x@tracds_d$mm)
		if(prior) ic <- ic+1

		switch(type,
			probability = ic / sum(ic),
			counts = ic,
			log_odds = log(ic / sum(ic)* size(x))
		)
	}
)


setMethod("transition_table", signature(x = "EMM", newdata = "numeric"),
	function(x, newdata, type= c("probability", "counts", "log_odds"), 
		match_cluster="exact", prior = TRUE, 
		initial_transition = FALSE) 
	transition_table(x, as.matrix(rbind(newdata)), type, 
		match_cluster, prior, initial_transition)
)

setMethod("transition_table", signature(x = "EMM", newdata = "data.frame"),
	function(x, newdata, type= c("probability", "counts", "log_odds"), 
		match_cluster="exact", prior = TRUE, 
		initial_transition = FALSE) 
	transition_table(x, as.matrix(newdata), type, 
		match_cluster, prior, initial_transition)
)

setMethod("transition_table", signature(x = "EMM", newdata = "matrix"),
        function(x, newdata, type= c("probability", "counts", "log_odds"), 
                match_cluster="exact", prior = TRUE, 
                initial_transition = FALSE) {

            type<- match.arg(type)

            ## make sure  newdata is a matrix (maybe a single row)
            if(!is.matrix(newdata)) newdata <- as.matrix(rbind(newdata))
            n <- nrow(newdata)

            ## empty EMM or single state?
            if(n<2) { 
                df <- data.frame(from=NA, to=NA, val=NA)
                names(df)[3] <- type
                return(df)
            }

            ## get sequence
            ssequence <- find_clusters(x, newdata, match_cluster=match_cluster, 
                    dist=FALSE)
            from <- ssequence[1:(n-1)]
            to <- ssequence[2:n]

            ## get values
            res <- transition(x, from, to, type=type, 
                    prior=prior)

            if(initial_transition) {
                from <- c(NA, from)
                to <- c(ssequence[1], to)
                res <- c(initial_transition(x, type=type, 
                                prior=prior)[ssequence[1]], 
                        res)
            }

            df <- data.frame(from=from, to=to, val=res, stringsAsFactors=FALSE)
            names(df)[3] <- type
            return(df)
        })
