visitNodePreorder <- function(v, node_id = 1)
{
	left_child_id = node_id * 2
	right_child_id = node_id * 2 + 1
	left_subtree = c()
	right_subtree = c()
	if (any(left_child_id == v))
		left_subtree = visitNodePreorder(v, left_child_id)
	if (any(right_child_id == v))
		right_subtree = visitNodePreorder(v, right_child_id)
	return(c(which(node_id == v), left_subtree, right_subtree))
}