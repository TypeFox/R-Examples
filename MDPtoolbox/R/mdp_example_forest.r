mdp_example_forest <- function(S, r1, r2, p) {

if ( nargs() >= 1 & ifelse(!missing(S), S <= 1, F) ) {
	print('----------------------------------------------------------')
	print('MDP Toolbox ERROR: Number of states S must be upper than 1')
	print('----------------------------------------------------------')
} else if ( nargs() >= 2 & ifelse(!missing(r1), r1 <= 0, F) ) {
	print('-----------------------------------------------------------')
	print('MDP Toolbox ERROR: The reward value r1 must be upper than 0')
	print('-----------------------------------------------------------')
} else if ( nargs() >= 3 & ifelse(!missing(r2), r1 <= 0, F) ) {
	print('-----------------------------------------------------------')
	print('MDP Toolbox ERROR: The reward value r2 must be upper than 0')
	print('-----------------------------------------------------------')
} else if ( nargs() >=4 & ifelse(!missing(p), p < 0 | p > 1, F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Probability p must be in [0; 1]')
	print('--------------------------------------------------------') 
} else {
	# initialization of optional arguments
	if (nargs() < 4) p <- 0.1
	if (nargs() < 3) r2 <- 2
	if (nargs() < 2) r1 <- 4
	if (nargs() < 1) S <- 3
	
	# Definition of Transition matrix P(:,:,1) associated to action Wait (action 1) and
	# P(:,:,2) associated to action Cut (action 2)
	#             | p 1-p 0.......0  |                  | 1 0..........0 |
	#             | .  0 1-p 0....0  |                  | . .          . |
	#  P(:,:,1) = | .  .  0  .       |  and P(:,:,2) =  | . .          . |
	#             | .  .        .    |                  | . .          . |
	#             | .  .         1-p |                  | . .          . |
	#             | p  0  0....0 1-p |                  | 1 0..........0 |

	P1 <- matrix(0,S,S)
	if (S > 2) diag(P1[-nrow(P1),-1]) <- (1-p) else P1[1,2] <- 1-p
	P1[,1] <- p
	P1[S,S] <- 1-p

	P2 <- matrix(0,S,S)
	P2[,1] <- 1

	P <- array(0, c(S,S,2))
	P[,,1] <- P1
	P[,,2] <- P2

	# Definition of Reward matrix R1 associated to action Wait and 
	# R2 associated to action Cut
	#           | 0  |                   | 0  |
	#           | .  |                   | 1  |
	#  R(:,1) = | .  |  and     R(:,2) = | .  |
	#           | .  |                   | .  |
	#           | 0  |                   | 1  |
	#           | r1 |                   | r2 |

	R1 <- numeric(S)
	R1[S] <- r1

	R2 <- rep(1,S)
	R2[1] <- 0
	R2[S] <- r2
	R <- cbind(R1,R2)

	return(list("P"=P, "R"=R))
}

}

