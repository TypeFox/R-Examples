EMSRb = function(Fare = Fare, Mean = Mean, Var = Var, p_up = numeric(length(Fare)), cap = cap) {
	return(RM2:::EMSRb_internal(Fare, Mean, Var, p_up, cap))
} # end function
