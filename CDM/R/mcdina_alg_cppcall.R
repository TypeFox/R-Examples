

###############################################################################
# calculation of individual likelihood
# ---------
# dat / dat_resp ... data and response indicator matrix
# group_ ... vector of group identifiers. Note that it must be integers
#       with zero as the lowest integer
# probs_ ... matrix with probabilities
#      original array in R: [items , categories , theta points , groups ]
#                           [ i , c , t , g ]
#      input in the Cpp routine is the 2-dimensional format [ i , c*t*g ]
# CC_ ... number of item response categories
# TP_ ... number of theta points (skill classes)
# ---
# INPUT:
# extern "C" {
# SEXP probs_pcm_groups_C( SEXP dat_, SEXP dat_resp_, SEXP group_, 
#     SEXP probs_, SEXP CC_, SEXP TP_) ;
#  }
probs_pcm_groups_Cpp <- 
function( dat_ , dat_resp_,  group_ , probs_,  CC_ ,  TP_ ){
	res <- .Call("probs_pcm_groups_C", 
					dat_, dat_resp_,  group_, probs_,  CC_,  TP_ ,
					PACKAGE = "CDM")
    return(res)			
			}
# ---			
# OUTPUT:			
#     return List::create(  
#     		_["fyiqk"] = fyiqk  
#     			) ;  			
##########################################################################

##########################################################################
# calculation of expected counts
# ---
# INPUT:
# extern "C" {
# SEXP calccounts_pcm_groups_C( SEXP dat_, SEXP dat_resp_, SEXP group_, 
#      SEXP fyiqk_, SEXP pik_, SEXP CC_, SEXP weights_) ;
# }
# fyiqk_  ... individual likelihood
# pik_    ... matrix with group distributions
# CC_     ... maximum number of response categories per item
# weight_ ... vector of sampling or frequency weights
calccounts_pcm_groups_Cpp <- 
function( dat_,  dat_resp_,  group_, fyiqk_,  pik_,  CC_,  weights_ ){
	res <- .Call("calccounts_pcm_groups_C", 
					dat_,  dat_resp_,  group_, fyiqk_,  pik_,  CC_,  weights_ ,
					PACKAGE = "CDM")
    return(res)			
			}
# ---
# OUTPUT:
#     return List::create(  
#     		_["LL"] = LL ,  
#     		_["fqkyi"] = fqkyi ,  
#     		_["nik"] = nik ,  
#           _["count_pik"] = pik1  
#    			) ;  			
#############################################################################

  
