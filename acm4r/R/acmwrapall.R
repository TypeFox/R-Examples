#' acmwrapall (acm warp all)
#'
#' This function takes only take two required parameters.
#' This function calls the function call_erra, call_acm, call_bablbs, call_gd1, call_gd2, and clusters with their default arguments.
#' 
#' @param replic  is the full name of the file containing the rflp result of replicate strains
#' @param patient is the full name of the file containing the rflp result of strains in patients
#' @param work_dir is where the datasets should be stored
#' @param dnum is the file number
#' @param delete logical value indicating if you want to delete any pre-existing files. Default is FALSE
#' @return none
#' @references 
#' Salamon et al. (1998) Accommodating Error Analysis in Comparison and Clustering of Molecular Fingerprints. Emerging Infectious Diseases Vol. 4, No. 2, April-June 1998
#' @author XiaoFei Zhao \email{xiaofei.zhao@@mail.mcgill.ca}
#' @export 
acmwrapall <- function(replic, patient, work_dir = ".", dnum = 1, delete = FALSE) {
	call_erra(replic, work_dir = work_dir, dnum = dnum, delete = delete)
	acm_table = call_acm(patient, work_dir = work_dir, dnum = dnum)
	bablbs_table = call_bablbs(acm_table)
	gd1_table = call_gd1(acm_table)
	print(gd1_table)
	gd2_table = call_gd2(acm_table)
	# clusters(acm_table)	
	print("acmwrapall finished")
}

