#'
#' Determine the starting dose and the dose-searching space for next subtrial in waterfall design
#'
#' Determine the starting dose and the dose-searching space for next subtrial after
#' the current subtrial is completed when using the waterfall design
#'
#'
#' @param target the target toxicity rate
#' @param npts a \code{J*K} matrix \code{(J<=K)} containing the number of patients treated at each dose combination
#' @param ntox a \code{J*K} matrix \code{(J<=K)} containing the number of patients who experienced dose-limiting
#'             toxicities at each dose combination
#' @param p.saf the highest toxicity probability that is deemed subtherapeutic (i.e. below
#'              the MTD) such that dose escalation should be undertaken. The default value
#'              is \code{p.saf=0.6*target}.
#' @param p.tox the lowest toxicity probability that is deemed overly toxic such that
#'              deescalation is required. The default value is \code{p.tox=1.4*target}.
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety. We recommend
#'                   the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset a small positive number (between 0 and 0.5) to control how strict the
#'               stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a
#'               more strict stopping rule. The default value \code{offset=0.05} generally
#'               works well.
#'
#' @details For the waterfall design, this function is used to obtain the starting dose and
#'          dose-searching space for the next subtrial when the current subtrial is completed.
#'          The input data include: the number of patients treated at each dose combination
#'           (i.e., \code{npts}), the number of patients who experienced dose-limiting
#'           toxicities at each dose combination (i.e., \code{ntox}).
#'
#'
#' @return \code{next.subtrial()} returns the starting dose and the dose-searching space for
#'         the next subtrial.
#'
#' @export
#'
#' @author Suyu Liu and Ying Yuan
#'
#' @references Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'             Trials, Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#'            Lin R. and Yin, G. (2016). Bayesian Optimal Interval Designs for Dose Finding in
#'            Drug-combination Trials, Statistical Methods in Medical Research, to appear.
#'
#'            Zhang L. and Yuan, Y. (2016). A Simple Bayesian Design to Identify the Maximum
#'            Tolerated Dose Contour for Drug Combination Trials, under review.
#'
#' @seealso  Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.2_tutorial.pdf}
#'
#'           Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
#'
#' @examples
#' n<-matrix(c(6, 0, 0, 0,
#'            6, 10, 12, 0,
#'            9, 12, 0, 0), ncol=4, byrow=TRUE)
#' y<-matrix(c(0, 0, 0, 0,
#'             1, 1, 4, 0,
#'             2, 3, 0, 0), ncol=4, byrow=TRUE)
#'
#' next.subtrial(target=0.3, npts=n, ntox=y)
#'
#'
next.subtrial <- function(target, npts, ntox, p.saf="default", p.tox="default",
                          cutoff.eli=0.95, extrasafe=FALSE, offset=0.05){

  waterfall.subtrial.mtd <- function(target, npts, ntox, cutoff.eli=0.95,
                                     extrasafe=FALSE, offset=0.05){
  ## obtain dose escalation and deescalation boundaries
  temp=get.boundary(target, ncohort=150, cohortsize=1, n.earlystop=100,
                    p.saf="default", p.tox="default", cutoff.eli, extrasafe, print=FALSE);
  b.e=temp[2,];   # escalation boundary

  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function (x, wt = rep(1, length(x))){
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  ## determine whether the dose has been eliminated during the trial
  y=ntox; n=npts; ndose=length(n); elimi=rep(0, ndose); is.escalation=0
  for(i in 1:ndose) {
    if(n[i]>=3) {if(1-pbeta(target, y[i]+1, n[i]-y[i]+1)>cutoff.eli) {elimi[i:ndose]=1; break;}}
  }
  if(extrasafe){if(n[1]>=3) {if(1-pbeta(target, y[1]+1, n[1]-y[1]+1)>cutoff.eli-offset) {elimi[1:ndose]=1;}} }

  ## no dose should be selected (i.e., selectdose=99) if the first dose is already very toxic or
  ## all uneliminated doses are never used to treat patients
  if(elimi[1]==1 || sum(n[elimi==0])==0) { selectdose=99; }
  else {
    adm.set = (n!=0) & (elimi==0); adm.index = which(adm.set==T);
    y.adm = y[adm.set]; n.adm = n[adm.set];

    ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
    phat = (y.adm+0.05)/(n.adm+0.1);
    phat.var = (y.adm+0.05)*(n.adm-y.adm+0.05)/((n.adm+0.1)^2*(n.adm+0.1+1));

    ## perform the isotonic transformation using PAVA
    phat = pava(phat, wt=1/phat.var)
    phat = phat + (1:length(phat))*1E-10 ## break ties by adding an increasingly small number
    selectd = sort(abs(phat-target), index.return=T)$ix[1]  ## select dose closest to the target as the MTD
    selectdose = adm.index[selectd];

    if(y[selectdose]<=b.e[n[selectdose]]) { is.escalation=1 }
  }

  list(selectdose=selectdose, is.escalation=is.escalation)
}


##############################################################################################################
## main code for next.subtrial starts here
    n=npts; y=ntox;

	if(sum(y>n)>0) { cat("Error: The data entry may be wrong. Please check it. \n"); return(1); }

    if(nrow(n)>ncol(n) | nrow(y)>ncol(y) ) {cat("Error: npts and ntox should be arranged in a way (i.e., rotated) such that for each of them, the number of rows is less than or equal to the number of columns."); return()}

## get the current and next subtrial indices
	subtrial.space = list()
	subtrial.space[[nrow(n)]] = c(1:(dim(n)[1]-1), (1:dim(n)[2])*dim(n)[1])
	for(j in (dim(n)[1]-1):1) subtrial.space[[j]] = (2:ncol(n))*nrow(n)-(nrow(n)-j)#  dim(n)[1]-j+1 + ((2:dim(n)[2]) -1) * dim(n)[1]

	cur.subtrial=0; nxt.subtrial=0
	for(k in dim(n)[1]:1){
	  if(sum(n[subtrial.space[[k]]])==0) {nxt.subtrial = k; break}
	}
	cur.subtrial = nxt.subtrial + 1
	#if(sum(n[nrow(n),2:dim(n)[2]])!=0) cur.subtrial = nrow(n)
	if(cur.subtrial == 1) { cat("No additional next subtrials are needed!!"); return(); }
	else{
		cur.dosespace = subtrial.space[[cur.subtrial]]
		nxt.dosespace = subtrial.space[[nxt.subtrial]]

	## determine the starting dose for current subtrial
	#	if(length(dose.curr)==0){
	sds= cur.dosespace[which(n[cur.dosespace]>0)[1]] #cur.dosespace[1]
	dj = ifelse(sds%%dim(n)[1]==0, sds%/%dim(n)[1], sds%/%dim(n)[1]+1)
	di = sds - (dj-1) * dim(n)[1]
	dose.curr=c(di, dj)
	#	}

	## if the user does not provide p.saf and p.tox, set them to the default values
		if(p.saf=="default") p.saf=0.6*target;
		if(p.tox=="default") p.tox=1.4*target;

	## simple error checking
		if(npts[dose.curr[1], dose.curr[2]]==0)  {cat("Error: dose entered is not the current dose \n"); return(1);}
		if(target<0.05) {cat("Error: the target is too low! \n"); return(1);}
		if(target>0.6)  {cat("Error: the target is too high! \n"); return(1);}
		if((target-p.saf)<(0.1*target)) {cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return(1);}
		if((p.tox-target)<(0.1*target)) {cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return(1);}
		if(offset>=0.5) {cat("Error: the offset is too large! \n"); return();}
		##if(n.earlystop<=6) {cat("Warning: the value of n.earlystop is too low to ensure
		##                        good operating characteristics. Recommend n.earlystop = 9 to 18 \n");
		##                    return();}

	## obtain dose escalation and de-escalation boundaries
		temp=get.boundary(target, ncohort=150, cohortsize=1, n.earlystop=100, p.saf, p.tox, cutoff.eli, extrasafe, offset, print=FALSE);
		b.e=temp[2,];   # escalation boundary
		b.d=temp[3,];   # deescalation boundary
		b.elim=temp[4,];  # elimination boundary
		lambda1  = log((1-p.saf)/(1-target))/log(target*(1-p.saf)/(p.saf*(1-target)));
		lambda2  = log((1-target)/(1-p.tox))/log(p.tox*(1-target)/(target*(1-p.tox)));

	## determine the starting dose for current subtrial
		earlystop=0;
		d=dose.curr;   # n=npts; y=ntox;
		nc = n[d[1],d[2]];
		ndose=length(npts);
		elimi = matrix(rep(0, ndose),dim(n)[1],dim(n)[2]);  ## indicate whether doses are eliminated

	## determine if early termination is needed
	##	if(n[d[1],d[2]]>=n.earlystop) {
	##		cat("Current subtrial is terminated early because the number of patients treated at the lowest dose (", d[1], ", ", d[2], ") has reached", n.earlystop,  "\n");
	##		d=c(99, 99); #earlystop=1;
	##	}

		if(!is.na(b.elim[nc])) {
			if(d[1]==1 && d[2]==1 && y[d[1],d[2]]>=b.elim[nc]) {
				d=c(99, 99); earlystop=1;
				cat("Current subtrial is terminated because the lowest dose is overly toxic \n");
			}

	## implement the extra safe rule by decreasing the elimination cutoff for the lowest dose
			if(extrasafe) {
				if(d[1]==1 && d[2]==1 && y[1,1]>=3) {
					if(1-pbeta(target, y[1,1]+1, n[1,1]-y[1,1]+1)>cutoff.eli-offset) {
						d=c(99, 99); earlystop=1;
						cat("Current subtrial is terminated because the lowest dose is overly toxic \n");
					}
				}
			}
		}


	## determine elimination status for combinations
		for(i in 1:dim(n)[1])
		{
			for(j in 1:dim(n)[2])
			{
				if(n[i,j]>0 && (!is.na(b.elim[n[i,j]])))
				{
					if(y[i,j]>=b.elim[n[i,j]])
					{
						elimi[i:dim(n)[1], j:dim(n)[2]]=1;
					}
				}
			}
		}

		if(earlystop==0)
		{
			wsmtd = waterfall.subtrial.mtd(target, n[cur.dosespace], y[cur.dosespace], cutoff.eli, extrasafe, offset)
			seldose = cur.dosespace[wsmtd$selectdose]
			if(is.na(seldose)==TRUE){ cat("Current subtrial is terminated early and no MTD is suggested for current subtrial. \n\n")
			}else if(seldose==99){ d=c(99,99)
			}else{
				dj = ifelse(seldose%%dim(n)[1]==0, seldose%/%dim(n)[1], seldose%/%dim(n)[1]+1)
				di = seldose - (dj-1) * dim(n)[1]
				d=c(di, dj)
				dnext = c(di-1, ifelse(dj==dim(n)[2], dj, dj+1))

				FUNC = function(x) paste('(', dnext[1], ', ', x,')', sep='')
				dnextspace = paste(unlist(lapply(2:ncol(n), FUNC)), collapse=', ')

				cat("Next subtrial includes doses: ","\n")
				cat("\t\t", dnextspace, "\n\n")
				cat("The starting dose for this subtrial is:\n",
					"\t\t", paste("(", dnext[1], ", ", dnext[2], ")", sep=''), "\n");
			}
		}
		invisible(d);
	}
}

