#'Parallel computing the Likelihood Ratios for the Gene Sets under Scrutiny
#'
#'A parallel version of the function \code{\link{TcGSA.LR}} to be used on a
#'cluster of computing processors.  This function computes the Likelihood
#'Ratios for the gene sets under scrutiny, as well as estimations of genes
#'dynamics inside those gene sets through mixed models.
#'
#'This Time-course Gene Set Analysis aims at identifying gene sets that are not
#'stable over time, either homogeneously or heterogeneously (see \emph{Hejblum
#'et al, 2012}) in terms of their probes.  And when the argument
#'\code{separatePatients} is \code{TRUE}, instead of identifying gene sets that
#'have a significant trend over time (possibly with probes heterogeneity of
#'this trend), \emph{TcGSA} identifies gene sets that have significantly
#'different trends over time depending on the patient.
#'
#'If the \code{monitorfile} argument is a character string naming a file to
#'write into, in the case of a new file that does not exist yet, such a new
#'file will be created. A line is written each time one of the gene sets under
#'scrutiny has been analysed (i.e. the two mixed models have been fitted, see
#'\code{\link{TcGSA.LR}}) by one of the parallelized processors.
#'
#'@aliases TcGSA.LR.parallel
#'
#'@param Ncpus The number of processors available on the cluster.
#'
#'@param type_connec The type of connection between the processors. Supported
#'cluster types are \code{"SOCK"}, \code{"PVM"}, \code{"MPI"}, and
#'\code{"NWS"}. See also \code{\link[parallel:makeCluster]{makeCluster}}.
#'
#'@param expr 
#'a matrix or dataframe of gene expression.  Its dimension are
#'\eqn{n}x\eqn{p}, with the \eqn{p} samples in column and the \eqn{n} genes in
#'row.
#'
#'@param gmt 
#'a \bold{gmt} object containing the gene sets definition.  See
#'\code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and definition on 
#'\href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{www.broadinstitute.org}.
#'
#'@param design
#'a matrix or dataframe containing the experimental variables that used in the model,
#'namely \code{subject_name}, \code{time_name}, and \code{covariates_fixed} 
#'and \code{time_covariates} if applicable.  Its dimension are \eqn{p}x\eqn{m} 
#'and its row are is in the same order as the columns of \code{expr}.
#'
#'@param subject_name
#'the name of the factor variable from \code{design} that contains the information on 
#'the repetition units used in the mixed model, such as the patient identifiers for instance.  
#'Default is \code{'Patient_ID'}.  See Details.
#'
#'@param time_name
#'the name of the numeric or factor variable from \code{design} contains 
#'the information on the time replicates (the time points at which gene 
#'expression was measured).  Default is \code{'TimePoint'}.  See Details.
#'
#'@param crossedRandom
#'logical flag indicating wether the random effects of the subjects and of the time points
#'should be modeled as one crossed random effect or as two separated random effects.  
#'Default is \code{FALSE}. See details.
#'
#'@param covariates_fixed
#'a character vector with the names of numeric or factor variables from the \code{design} 
#'matrix that should appear as fixed effects in the model.  See details.
#'Default is \code{""}, which corresponds to no covariates in the model.
#'
#'@param time_covariates
#'the name of a numeric variable from \code{design} that contains 
#'the information on the time replicates (the time points at which gene 
#'expression was measured).  Default is \code{'TimePoint'}.  See Details.
#'
#'
#'@param time_func 
#'the form of the time trend. Can be either one of \code{"linear"},
#'\code{"cubic"}, \code{"splines"} or specified by the user, or the column name of 
#'a factor variable from \code{design}. If specified by the user, 
#'it must be as an expression using only names of variables from the \code{design} matrix 
#'with only the three following operators: \code{+}, \code{*}, \code{/} . 
#'The \code{"splines"} form corresponds to the natural cubic B-splines 
#'(see also \code{\link[splines:ns]{ns}}).  If there are only a few timepoints, 
#'a \code{"linear"} form should be sufficient. Otherwise, the \code{"cubic"} form is 
#'more parsimonious than the \code{"splines"} form, and should be sufficiently flexible.
#'If the column name of a factor variable from \code{design} is supplied, 
#'then time is considered as discrete in the analysis.
#'If the user specify a formula using column names from design, both factor and numeric
#'variables can be used.
#'
#'@param minGSsize 
#'the minimum number of genes in a gene set.  If there are
#'less genes than this number in one of the gene sets under scrutinity, the
#'Likelihood Ratio of this gene set is not computed (the mixed model are not
#'fitted). Default is \code{10} genes as the minimum.
#'
#'@param maxGSsize 
#'the maximum number of genes in a gene set.  If there are
#'more genes than this number in one of the gene sets under scrutinity, the
#'Likelihood Ratio of this gene set is not computed (the mixed model are not
#'fitted).  This is to avoid very long computation times.  Default is
#'\code{500} genes as the maximum.
#'
#'@param group_name 
#'in the case of several treatment groups, the name of a factor variable 
#'from the \code{design} matrix.  It indicates to which treatment group each sample
#' belongs to.  Default is \code{""}, which means that there is only one 
#' treatment group.  See Details.
#'
#'@param separateSubjects
#'logical flag indicating that the analysis identifies
#'gene sets that discriminates patients rather than gene sets than have a
#'significant trend over time.  Default is \code{FALSE}.  See Details.
#'
#'@param monitorfile
#'a writable \link{connections} or a character string naming a file to write into, 
#'to monitor the progress of the analysis.  
#'Default is \code{""} which is no monitoring.  See Details.
#'
#'@return \code{TcGSA.LR} returns a \code{tcgsa} object, which is a list with
#'the 5 following elements:
#'\itemize{
#'\item fit a data frame that contains the 3 following variables:
#'\itemize{ 
#'\item \code{LR}: the likelihood ratio between the model under the
#'null hypothesis and the model under the alternative hypothesis.  
#'\item
#'\code{CVG_H0}: convergence status of the model under the null hypothesis.
#'\item \code{CVG_H1}: convergence status of the model under the alternative
#'hypothesis.
#'}
#'\item \code{time_func}: a character string passing along the value of the
#'\code{time_func} argument used in the call.
#'\item \code{GeneSets_gmt}: a \code{gmt} object passing along the value of the
#'\code{gmt} argument used in the call.
#'\item \code{group.var}: a factor passing along the \code{group_name} variable
#'from the \code{design} matrix.
#'\item \code{separateSubjects}: a logical flag passing along the value of the
#'\code{separateSubjects} argument used in the call.
#'\item \code{Estimations}: a list of 3 dimensions arrays.  Each element of the
#'list (i.e. each array) corresponds to the estimations of gene expression
#'dynamics for each of the gene sets under scrutiny (obtained from mixed
#'models).  The first dimension of those arrays is the genes included in the
#'concerned gene set, the second dimension is the \code{Patient_ID}, and the
#'third dimension is the \code{TimePoint}.  The values inside those arrays are
#'estimated gene expressions.
#'\item \code{time_DF}: the degree of freedom of the natural splines functions
#'}
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{summary.TcGSA}}, \code{\link{plot.TcGSA}}
#'
#'@references Hejblum BP, Skinner J, Thiebaut R, (2015) 
#'Time-Course Gene Set Analysis for Longitudinal Gene Expression Data. 
#'\emph{PLoS Computat Biol} 11(6): e1004310.
#'doi: 10.1371/journal.pcbi.1004310
#'
#'@importFrom GSA GSA.read.gmt
#'
#'@importFrom lme4 lmer
#'
#'@importFrom stats as.formula deviance fitted
#'
#'@export TcGSA.LR.parallel
#'
#'@examples
#'
#'data(data_simu_TcGSA)
#'
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'                           
#'\dontrun{ 
#'require(doParallel)
#'tcgsa_sim_1grp <- TcGSA.LR.parallel(Ncpus = 2, type_connec = 'SOCK',
#'                             expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                             subject_name="Patient_ID", time_name="TimePoint",
#'                             time_func="linear", crossedRandom=FALSE, 
#'                             separateSubjects=TRUE)
#'}
#'tcgsa_sim_1grp
#'summary(tcgsa_sim_1grp)
#'     
#'
#'

TcGSA.LR.parallel <-
	function(Ncpus, type_connec, 
			 expr, gmt, design, subject_name="Patient_ID", time_name="TimePoint", crossedRandom=FALSE,
			 covariates_fixed="", time_covariates="",
			 time_func = "linear", group_name="", separateSubjects=FALSE,
			 minGSsize=10, maxGSsize=500,
			 monitorfile=""){
		if (requireNamespace("doParallel", quietly = TRUE)) {
			requireNamespace("doParallel")
			if(group_name!="" && separateSubjects){
				stop("'separateSubjects' is TRUE while 'group_name' is not \"\".\n This is an attempt to separate subjects in a multiple group setting.\n This is not handled by the TcGSA.LR function.\n\n")
			}
			
			my_formul <- TcGSA.formula(design=design, subject_name=subject_name, time_name=time_name,  
									   covariates_fixed=covariates_fixed, time_covariates=time_covariates, group_name=group_name,
									   separateSubjects=separateSubjects, crossedRandom=crossedRandom,
									   time_func=time_func)
			time_DF <- my_formul[["time_DF"]]
			
			
			cl <- parallel::makeCluster(Ncpus, type = type_connec)
			doParallel::registerDoParallel(cl)
			
			
			gs <- NULL # This is just to prevent R CMD check to issue a NOTE reanding "TcGSA.LR.parallel: no visible binding for global variable 'gs'". It is extremely ANNOYING !
			res_par <- foreach::"%dopar%"(foreach::foreach(gs=1:length(gmt$genesets), .packages=c("lme4", "reshape2", "splines"), .export=c("TcGSA.dataLME")),
										  {
				probes <- intersect(gmt$genesets[[gs]], rownames(expr))
				
				if(length(probes)>0 && length(probes)<=maxGSsize && length(probes)>=minGSsize){                                                       
					expr_temp <- t(expr[probes, ])
					rownames(expr_temp) <- NULL
					data_lme  <- TcGSA.dataLME(expr=expr_temp, design=design, subject_name=subject_name, time_name=time_name, 
											   covariates_fixed=covariates_fixed, time_covariates=time_covariates,
											   group_name=group_name, time_func=time_func)
					
					if(length(levels(data_lme$probe))>1){
						lmm_H0 <- tryCatch(lmer(formula =my_formul[["H0"]]["reg"], REML=FALSE, data=data_lme),
										   error=function(e){NULL})
						lmm_H1 <- tryCatch(lmer(formula =my_formul[["H1"]]["reg"], REML=FALSE, data=data_lme),
										   error=function(e){NULL})
					}
					else{
						lmm_H0 <- tryCatch(lmer(formula =my_formul[["H0"]]["1probe"], REML=FALSE, data=data_lme),
										   error=function(e){NULL})
						lmm_H1 <- tryCatch(lmer(formula =my_formul[["H1"]]["1probe"], REML=FALSE, data=data_lme),
										   error=function(e){NULL})
					}
					
					if (!is.null(lmm_H0) & !is.null(lmm_H1)) {
						LR <- stats::deviance(lmm_H0, REML=FALSE) - stats::deviance(lmm_H1, REML=FALSE)
						CVG_H0 <- lmm_H0@optinfo[["conv"]]$opt
						CVG_H1 <- lmm_H1@optinfo[["conv"]]$opt
						
						estims <- cbind.data.frame(data_lme, "fitted"=stats::fitted(lmm_H1))
						estims_tab <- reshape2::acast(data=estims, formula = stats::as.formula(paste("probe", subject_name, "t1", sep="~")), value.var="fitted")
						# drop = FALSE by default, which means that missing combination will be kept in the estims_tab and filled with NA
						dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
						estim_expr <- estims_tab
					} 
					else {
						LR <- NA
						CVG_H0 <- NA
						CVG_H1 <- NA
						
						estims <- cbind.data.frame(data_lme, "fitted"=NA)
						estims_tab <- reshape2::acast(data=estims, formula = stats::as.formula(paste("probe", subject_name, "t1", sep="~")), value.var="fitted")
						dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
						estim_expr <- estims_tab
						cat("Unable to fit the mixed models for this gene set\n")
					}
					

					# CONVERGENCE DIAGNOSTICS IN lme4 v1.1-7 (from Nelder Mead optimizer)
					# -3: "nm_forced"
					# -2: "cannot generate a feasible simplex"
					# -1: "initial x is not feasible"
					#  0: "converged"
					
				}
				else{
					LR <- NA
					CVG_H0 <- NA
					CVG_H1 <- NA
					
					estim_expr <- NA
					cat("The size of the gene set",  gmt$geneset.names[[gs]], "is problematic (too many or too few genes)\n")
				}
				
				line_number <- 0
				try(line_number <- length(readLines(monitorfile)), silent=TRUE)
				cat(paste(line_number+1,"/", length(gmt$genesets)," gene sets analyzed (geneset ", gs, ")\n", sep=""), file=monitorfile, append = TRUE)
				
				res <- list("LR"=LR, "CVG_H0"=CVG_H0, "CVG_H1"=CVG_H1, "estim_expr"=estim_expr)
			})
			
			
			cat("Combining the results...")
			
			LR <- sapply(res_par, "[[", "LR") #res_par[[gs]][["LR"]]
			CVG_H0 <- sapply(res_par, "[[", "CVG_H0") #res_par[[gs]][["CVG_H0"]]
			CVG_H1 <- sapply(res_par, "[[", "CVG_H1") #res_par[[gs]][["CVG_H1"]]
			estim_expr <- lapply(res_par, "[[", "estim_expr") #res_par[[gs]][["estim_expr"]]
			
			# LR <- numeric(length(gmt$genesets))
			# CVG_H0 <- numeric(length(gmt$genesets))
			# CVG_H1 <- numeric(length(gmt$genesets))
			# estim_expr <- list()
			# for (gs in 1:length(gmt$genesets)){
			# 	LR[gs] <- res_par[[gs]][["LR"]]
			# 	CVG_H0[gs] <- res_par[[gs]][["CVG_H0"]]
			# 	CVG_H1[gs] <- res_par[[gs]][["CVG_H1"]]
			# 	estim_expr[[gs]] <- res_par[[gs]][["estim_expr"]]
			# }
			
			if(group_name==""){
				gv <- NULL
			} else{
				gv <- design[,group_name]
			}
			
			tcgsa <- list("fit"=as.data.frame(cbind(LR, CVG_H0, CVG_H1)), "time_func"=time_func, "GeneSets_gmt"=gmt, 
						  "group.var"=gv, "separateSubjects"=separateSubjects, "Estimations"=estim_expr, 
						  "time_DF"=time_DF
			)
			class(tcgsa) <- "TcGSA"
			parallel::stopCluster(cl)
			
		} else {
			stop("Package 'doParallel' is not available.\n  -> Try running 'install.packages(\"doParallel\")'\n   or use non parallel version of the function: 'TcGSA.LR'")
		}
		
		return(tcgsa)
	}


