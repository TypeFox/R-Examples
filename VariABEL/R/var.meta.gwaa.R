#=====================================================================================
#
#       Filename:  var.meta.gwaa.R
#
#    Description:  Function for meta analysis ov variance. Read flat files in plink format and perform metanalysis.
#
#        Version:  1.0
#        Created:  26-Apr-2010
#       Revision:  none
#
#
#         Author:  Maksim V. Struchalin
#        Company:  ErasmusMC, Epidemiology, The Netherlands.
#          Email:  m.struchalin@erasmusmc.nl
#
#=====================================================================================


#
#
#"var.meta.gwaa" <-
#function(input_filenames, 
#				 testname,
#				 output_filename="output.variance.metaanalysis",
#				 exclude_snp_below_threshold=F,
#				 threshold=30,
#				 all_warnings=F,
#				 skip_first_lines_amount=0,
#				 delim=' ') {
#
#
#file_amount <- length(filenames_array)
#
#
#return_val <-  .C("var_meta_gwaa_C", input_filenames, as.integer(file_amount), output_filename, as.integer(skip_first_lines_amount), delim,
#									lambdas = double(file_amount+1),
#									lambdas_NA = integer(file_amount+1),
#									as.logical(exclude_snp_below_threshold),
#									as.integer(threshold),
#									as.logical(all_warnings),
#									as.character(testname))
#			
#
#return_val$lambdas[return_val$lambdas_NA==1] <- NA
#
#lambda_df <- data.frame(filename=filenames_array, lambda=return_val$lambdas[1:file_amount])
#
#lambda_df <- rbind(lambda_df, data.frame(filename=output_filename, lambda=return_val$lambdas[file_amount+1]))
#
#lambda_df
#}
