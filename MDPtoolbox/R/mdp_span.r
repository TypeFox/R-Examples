# function sp = mdp_span(W)

mdp_span <- function(W)

# mdp_span  Returns the span of W 
#           sp(W) = max W(s) - min W(s)

{
	return(max(W) - min(W))
}
