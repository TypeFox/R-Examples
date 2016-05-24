# ParGateRule function is a secondary function which generates decision rules for multistage parallel
# gatekeeping procedures in hypothesis testing problems with multiple families
# of null hypotheses (null hypotheses are assumed to be equally weighted within
# each family)
pargaterule<-function(gateproc,alpha,independence)
# GATEPROC, List of gatekeeping procedure parameters
# ALPHA, Global familywise error rate
# INDEPENDENCE, Boolean indicator (TRUE, Independence condition is imposed; FALSE,
# Independence condition is not imposed)
{

    # Number of families
    nfams<-length(gateproc)

	# Evaluate decision rules for the multistage parallel gatekeeping procedure
	gateproc<-pargateeval(gateproc,alpha,independence)

	cat("Hypothesis testing problem\n\n")
	cat("Global familywise error rate=", alpha, "\n",sep="")
	if (independence==TRUE) cat("Independence condition is imposed (the families are tested from first to last) \n",sep="")
	if (independence==FALSE) cat("Independence condition is not imposed (the families are tested from first to last and then re-tested from last to first) \n",sep="")
	cat("\n\n")

	# Retesting indicator
	retest<-FALSE

    # Sequentially numbered hypotheses labels
    hyplabel<-vector("list",nfams)
    cumtotal<-0

	# Test families from first to last
	for (i in 1:nfams)
	{
		# Local alpha level
		la<-round(gateproc[[i]][[6]],4)
		# Procedure parameter
		pp<-round(gateproc[[i]]$procpar,4)

		# Family and procedure
		cat("Family ", i, " (", gateproc[[i]]$label, ") is tested using", sep="")
		cat(" ", gateproc[[i]]$proc, " procedure (truncation parameter=", pp, ") at alpha", i, "=", la, ".\n\n", sep="")

		# Number of null hypotheses
		nhyps<-length(gateproc[[i]]$rawp)
        # Hypotheses labels
        hyplabel[[i]]<-seq(1:nhyps)+cumtotal
        cumtotal<-cumtotal+nhyps
		# Number of rejected null hypotheses
		rejcount<-sum(gateproc[[i]][[7]])

		for (j in 1:nhyps)
		{
			# Raw p-value
			rp<-round(gateproc[[i]][[2]][j],4)
			# Adjusted p-value
			ap<-round(gateproc[[i]][[5]][j],4)

			cat("Null hypothesis ", hyplabel[[i]][j], " (raw p-value=", rp, ")", sep="")
			if (gateproc[[i]][[7]][j]==TRUE)
			{
				cat(" is rejected.\n\n", sep="")
			}
			if (la==0) cat(" is automatically accepted. \n\n", sep="")
			if (la>0 & gateproc[[i]][[7]][j]==FALSE) cat(" is accepted.\n\n", sep="")

		}

        # Details
        cat("Details on the decision rule for this family can be obtained by running the PValAdjP function for ",gateproc[[i]]$proc," procedure with gamma=",pp," and alpha=",la,".\n\n", sep="")

		# Consclusion
		if (rejcount==0 & i<nfams)
		{
			cat("No null hypotheses are rejected in Family ", i, " and the parallel gatekeeping procedure cannot pass this family.",sep="")
			cat(" Testing stops and all remaining null hypotheses are automatically accepted.\n\n\n")
		}

		if (rejcount>0 & i<nfams)
		{
			cat("One or more null hypotheses are rejected in Family ", i, " and the parallel gatekeeping procedure passes this family.",sep="")
			cat(" Based on the error rate function of ", gateproc[[i]]$proc, " procedure (truncation parameter=", pp, "),", sep="")
			cat(" alpha",i+1,"=",round(gateproc[[i+1]][[6]],4)," is carried over to Family ", i+1, ".\n\n\n",sep="")
		}

		if (i==nfams & independence==FALSE)
		{
			if (rejcount==nhyps)
			{
				retest<-TRUE
				cat("All null hypotheses are rejected in Family ", i, " and the parallel gatekeeping procedure passes this family.", sep="")
				cat(" Retesting begins and alpha",i+1,"=",round(alpha,4)," is carried over to Family ", i-1, ".\n\n\n",sep="")
			}
			if (rejcount<nhyps)
			{
				cat("Some null hypotheses are accepted in Family ", i, " and the parallel gatekeeping procedure cannot pass this family.", sep="")
				cat(" Retesting will not be performed.\n\n\n",sep="")
			}
		}
	}

	# Retest families from last to first if the independence condition is not imposed
	if (independence==FALSE & retest==TRUE)
	{
		for(i in (nfams+1):(2*nfams-1))
		{
			# Family index
			k<-2*nfams-i
			# Local alpha level
			la<-round(gateproc[[i]][[6]],4)
			# Procedure parameter
			pp<-round(gateproc[[i]]$procpar,4)

			cat("Family ", k, " (", gateproc[[k]]$label, ") is retested using", sep="")
			cat(" ", gateproc[[i]]$proc, " procedure (truncation parameter=", pp, ") at alpha", i, "=", la, ".\n\n", sep="")

			# Number of null hypotheses
			nhyps<-length(gateproc[[i]]$rawp)
			# Number of rejected null hypotheses
			rejcount<-sum(gateproc[[i]][[7]])

			for (j in 1:nhyps)
			{
				# Raw p-value
				rp<-round(gateproc[[i]][[2]][j],4)
				# Adjusted p-value
				ap<-round(gateproc[[i]][[5]][j],4)

				cat("Null hypothesis ", hyplabel[[k]][j], " (raw p-value=", rp, ")", sep="")
				if (gateproc[[i]][[7]][j]==TRUE)
				{
					cat(" is rejected.\n\n", sep="")
				}
				if (la==0) cat(" is automatically accepted.\n\n", sep="")
				if (la>0 & gateproc[[i]][[7]][j]==FALSE) cat(" is accepted.\n\n", sep="")
			}

            # Details
            cat("Details on the decision rule for this family can be obtained by running the PValAdjP function for ",gateproc[[i]]$proc," procedure with gamma=",pp," and alpha=",la,".\n\n", sep="")

            # Conclusions
			if (k>1)
			{
				if (rejcount==nhyps)
				{
					retest<-TRUE
					cat("All null hypotheses are rejected in Family ", k, " and the parallel gatekeeping procedure passes this family and", sep="")
					cat(" alpha",i+1,"=",round(alpha,4)," is carried over to Family ", k-1, ".\n\n\n",sep="")
				}
				if (rejcount<nhyps)
				{
					cat("Some null hypotheses are accepted in Family ", k, " and the parallel gatekeeping procedure cannot pass this family.", sep="")
					cat(" Retesting stops.\n\n\n",sep="")
				}
			}
		}
	}

}
# End of pargaterule

