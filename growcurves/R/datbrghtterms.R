#' BRIGHT BDI depressive symptom data with (G = 4) session groups divided into separate MM terms.
#' 
#' The Beck Depression Inventory - II scores for the set of de-identified clients who participated
#' in the Building Recovery by Improving Goals, Habits and Thoughts (BRIGHT) study, a community-based
#' effectiveness trial of group cognitive behavioral therapy intervention for treating residential 
#' substance abuse treatment clients experiencing depressive symptoms.  These data include scores for
#' three measurement waves; the first at baseline enrollment to the study, followed by two post-treatment
#' measurements with the aim to test whether clients receiving BRIGHT intervention would experience 
#' sustained improvement.  There 299 participating clients, divided between 159 assigned to the usual care arm,
#' and 140 assigned to CBT.  The data are configured to support model runs using engine function \code{dpgrowmm}.
#' 
#' \itemize{
#'   \item y. Client depressive symptom score responses.  There are \code{N = 815} total measures for \code{P = 299} subjects.  Each entry contains a composite Beck Depression Inventory-II (BDI-II) score.  
#'			The BDI-II score is a sum acoss 21 four-level items (each scored 0 - 3) with a higher score signifying a greater level of depressive symptoms.				  
#'   \item subject. BRIGHT study client identifier \code{(1,2,...,299}.  Note: Participating clients are de-identified in this dataset.
#'   \item trt. Treatment arm identifier of length \code{N} (e.g. \code{(0,0,0,...,1,1,1,...)} , either \code{{0,1}} with usual care (UC) = 0 and group cognitive behavioral therapy (CBT) = 1.
#'			There are 140 clients with CBT = 1 (even though 132 of these actually attend sessions) in order to facilitate an intent-to-treat comparison.	 
#'   \item time. Meaurement times in months for each repeated subject measure of length \code{N}.  There are 3 distinct time points or measurement waves.  e.g. \code{(0,3,6,0,3,6,0,0,3,,,,)}.
#'			The first measure is at baseline when clients enrolled to BRIGHT and at two post-treatment follow-ups at 3 and 6 months with response rates of 86% and 87%, respectively.
#'   \item subj.aff. A list object with each term a vector that indexes \code{n[g]} subjects linked to each of \code{g = 1, ..., (G = 4)} CBT therapy groups.  Each group is specialized
#'			to its own multiple membership (MM) term for employment of engine function \code{dpgrowmmult}.  The number of CBT clients for each group are 
#'			\code{(n[1] = 17, n[2] = 18, n[3] = 19, n[4] = 78)} for a total of \code{Paff = 132} clients that attended sessions, as compared to the 140 assigned to the CBT arm.
#'   \item subj.aff_mat. A matrix object that concatenates the client identifiers in \code{subj.aff} across groups for modeling all groups, together, in a single MM term with engine function \code{dpgrowmm}.
#'   \item W.subj.aff.  A list object containing \code{G, n[g] x S[g]} multiple membership weight matrices that together map the \code{P_aff = 132} affected subjects (in each element of the 	
#'			\code{subj.aff} list) to their particular sessions attended within their assigned group.  There are a total of 245 CBT sessions allocated to the \code{G} CBT therapy groups 
#'			as \code{S[1] = 36, S[2] = 40, S[3] = 40, S[4] = 129}.  The study was designed for clients to attend 16 sessions organized into 4 modules of 4 sessions each.   The 4 modules were offered
#'			on an open enrollment or rotating basis.  Additional make-up sessions added within 2 modules resulted in the possibility for some clients to attend up to 20 sessions.
#'   \item W.subj.aff_mat. An \code{n = 132 x S = 245} matrix object that concatenates the matrix entries of the list object from \code{W.subj.aff} into a block-diagonal matrix 
#'				(with each group of clients and sessions disjoint and non-communicating with the others) for modeling in a single MM term under function \code{dpgrowmm}.
#'   \item group.  A list object of length equal to the number of MM terms under the \code{"mmcar"} or \code{"mmigrp"} prior formulations for session effects.  Each item contains a vector
#'			that specifies sub-group membership for each block of sessions within the \code{G=4} MM terms.   A sub-group would collect sessions that communicate with each other,
#'			but not with the sessions of other sub-groups.  For these data, the sessions within each CBT therapy group all communicate.  Then each vector in \code{group} contains
#'			the single value 1 of length equal to the number of sessions in the applicable therapy group.   This object is input for engine function \code{dpgrowmult} in the
#'			case it is desired to place all of the \code{G = 4} MM terms under prior \code{"mmcar"}.   One may employ the appropriate subset of list entries for those
#'			terms under which it is desired to employ prior \code{"mmcar"}.
#'   \item group_mat. A matrix object denotes the group memberships, from 1 - \code{G = 4} groups, for the \code{S = 245} CBT sessions for modeling under engine function \code{dpgrowmm}.
#'   \item Omega. A list object with each element containing an \code{S[g] x S[g]} CAR adjacency matrix used to model prior association among effects for each MM term under prior \code{"mmcar"}, 
#'			where \code{S[g]} are the number of effects for CBT therapy group \code{g}.  One may employ the approach subset of adjacency matrices for those terms under which one
#'			desires to specify prior \code{"mmcar"}.
#'   \item Omega_mat. An \code{S x S} matrix object encoding the dependence structure among sessions that concatenates the entries of \code{Omega} into a block diagonal structure for use in \code{dpgrowmm}.
#' }
#' 
#' @references
#' 	K. E. Watkins, S. B. Hunter, K. A. Hepner, S. M. Paddock, E. de la Cruz, A. J. Zhou and J. Gilmore (2011) An effectiveness trial of group cognitive behavioral therapy
#'	for patients with persistent depressive symptoms in substance abuse treatment, Archives of General Psychiatry, 68(6), 1- 8.
#' @references
#'	S. M. Paddock and T. D. Savitsky (2012) Bayesian Hierarchical Semiparametric Modeling of Longitudinal Post-treatment Outcomes from Open-enrollment Therapy Groups, invited re-submission to: 
#'      JRSS Series A (Statistics in Society).
#' @docType data
#' @keywords datasets
#' @name datbrghtterms
#' @usage datbrghtterms
#' @format A list object for 815 total observations on 299 subjects
NULL