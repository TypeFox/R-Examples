create_cutoffs_for_K_fold_cv = function(pct_leave_out, n){
	begin_cutoffs_for_leave_outs = seq(1, n)[seq(1, n, n * pct_leave_out)] - 1
	begin_cutoffs_for_leave_outs[1] = 1
	end_cutoffs_for_leave_outs = seq(1, n)[seq(1, n, n * pct_leave_out)] - 2
	end_cutoffs_for_leave_outs = c(end_cutoffs_for_leave_outs[2 : length(end_cutoffs_for_leave_outs)], n)
	num_windows = length(begin_cutoffs_for_leave_outs)
	list(begin_cutoffs_for_leave_outs = begin_cutoffs_for_leave_outs, end_cutoffs_for_leave_outs = end_cutoffs_for_leave_outs, num_windows = num_windows)
}