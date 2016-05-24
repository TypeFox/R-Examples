#'
#' Get the operating characteristics for phase I/II waterfall design
#'
#' Obtain the operating characteristics of phase I/II waterfall design, which aims to find
#' the optimal dose combination (ODC), defined as the combination that has the highest
#' efficacy among the doses in the MTD contour.
#'
#' @param target the target toxicity rate
#' @param eff.lb the lower bound for efficacy
#' @param p.truetox a \code{J*K} matrix \code{(J<=K)} containing the true toxicity probabilities of
#'               combinations with \code{J} dose levels of agent A and \code{K} dose levels of agent B
#' @param p.trueeff a \code{J*K} matrix \code{(J<=K)} containing the true efficacy probability of
#'              combinations with \code{J} dose levels of agent A and \code{K} dose levels of agent B
#' @param ncohort1 the total number of cohorts for phase I
#' @param cohortsize1 the cohort size for phase I
#' @param n.earlystop the early stopping parameter for phase I. If the number of patients treated
#'  at the current dose reaches \code{n.earlystop}, stop the trial and select the MTD based
#'  on the observed data. The default value \code{n.earlystop=100} essentially turns off this
#'  type of early stopping.
#' @param Nmax1 the maximum number of patients for each subtrial in phase I
#' @param ncohort2 the total number of cohorts for phase II
#' @param cohortsize2 the cohort size for phase II
#' @param cutoff.eli the cutoff for dose elimination rule
#' @param cutoff.eff the cutoff for futility stopping
#' @param p.saf the highest toxicity probability that is deemed subtherapeutic (i.e. below the
#'  MTD) such that dose escalation should be undertaken. The default value is
#'  \code{p.saf=0.6*target}.
#' @param p.tox the lowest toxicity probability that is deemed overly toxic such that
#' deescalation is required. The default value is \code{p.tox=1.4*target}.
#' @param offset  a small positive number (between 0 and 0.5) to control how strict the stopping
#'               rule is when \code{extrasafe=TRUE}. A larger value leads to a more strict stopping
#'               rule. The default value \code{offset=0.05} generally works well.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param ntrial the total number of trials to be simulated
#'
#' @return This function returns the operating characteristics of the waterfall design as a list,
#' including (1) selection percentage at each dose level (\code{selpercent}), (2) the number of
#' patients treated at each dose (\code{npts}), (3) the number of toxicities observed at each
#' dose (\code{ntox}), (4) the number of efficacy/response observed at each dose (\code{neff})
#' (5) the total sample size (\code{totaln}).
#'
#' @export
#'
#' @details \code{get.oc.comb.phase12()} is consisted of two parts.
#' In the phase I part, the waterfall design is used to find the MTD contour
#' on the basis of only toxicity. Once the MTD contour is identified,
#' these MTDs are seamlessly moved to phase II to evaluate efficacy.
#' Each of the MTDs forms a treatment arm. Patients are eqally randomized into these arms to
#' evaluate efficacy. Toxicity and efficacy monitoring will be conducted after every
#' \code{cohortsize2} patients enrolled into each of treatment arms.
#'
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
#' p.truetox<-matrix(c(0.03,0.10,0.28, 0.10,0.30,0.50), nrow=2, byrow=TRUE)
#' p.trueeff<-matrix(c(0.20,0.30,0.50, 0.25,0.35,0.55), nrow=2, byrow=TRUE)
#' get.oc.comb.phase12(p.truetox, p.trueeff, target=0.30, eff.lb=0.2,
#'          ncohort1=12, cohortsize1=3, Nmax1=21, n.earlystop=12,
#'          ncohort2=12, cohortsize2=3, cutoff.eff =0.9, ntrial=10)
#'
#'
get.oc.comb.phase12=function(p.truetox, p.trueeff, target, eff.lb=0.2, ncohort1, cohortsize1,
                             Nmax1, n.earlystop=10, ncohort2, cohortsize2, cutoff.eli=0.95,
                             cutoff.eff, p.saf="default", p.tox="default", extrasafe=FALSE,
                             offset=0.05, ntrial=1000) {
  ############################### function ###############################
  aa=function(x) as.numeric(as.character(x))
  NPTS.all = 0

  waterfall.subtrial.mtd <- function(target, npts,ntox,cutoff.eli=0.95,extrasafe=FALSE,offset=0.05) {
      ## obtain dose escalation and deescalation boundaries
      temp=get.boundary(target,ncohort=150,cohortsize=1,n.earlystop=100,p.saf,p.tox,cutoff.eli,extrasafe,print=FALSE); b.e=temp[2,];
      ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
      pava <- function (x,wt=rep(1,length(x))) {
        n <- length(x); if (n<=1)  return(x)
        if (any(is.na(x)) || any(is.na(wt))) stop("Missing values in 'x' or 'wt' not allowed")
        lvlsets <- (1:n)
        repeat { viol <- (as.vector(diff(x)) < 0)
          if (!(any(viol)))  break
          i <- min((1:(n-1))[viol]); lvl1 <- lvlsets[i]; lvl2 <- lvlsets[i + 1];
          ilvl <- (lvlsets == lvl1 | lvlsets == lvl2); x[ilvl] <- sum(x[ilvl]*wt[ilvl])/sum(wt[ilvl]); lvlsets[ilvl] <- lvl1 }
        x
      }
      ## determine whether the dose has been eliminated during the trial
      y=ntox; n=npts; ndose=length(n); elimi=rep(0,ndose); is.escalation=0
      for (i in 1:ndose) if (n[i] >= 3 & 1-pbeta(target,y[i] + 1,n[i]-y[i] + 1)>cutoff.eli) { elimi[i:ndose]=1; break;}
      if (extrasafe) if (n[1] >= 3 & 1-pbeta(target,y[1] + 1,n[1]-y[1] + 1)>cutoff.eli-offset) {elimi[1:ndose]=1;}
      ## no dose should be selected (i.e.,selectdose=99) if the first dose is already
      ## very toxic or all uneliminated doses are never used to treat patients
      if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) { selectdose=99; }
      else {
        adm.set=(n!=0) & (elimi==0); adm.index=which(adm.set==T); y.adm=y[adm.set]; n.adm=n[adm.set];
        ## poster mean and variance of toxicity probabilities using beta(0.05,0.05) as the prior
        phat=(y.adm+0.05)/(n.adm+0.1); phat.var=(y.adm+0.05)*(n.adm-y.adm+0.05)/((n.adm + 0.1)^2*(n.adm+0.1+1));
        ## perform the isotonic transformation using PAVA
        phat=pava(phat,wt=1/phat.var); phat=phat + (1:length(phat)) * 1E-10 ## break ties by adding an increasingly small number
        selectd=sort(abs(phat-target),index.return=T)$ix[1]  ## select dose closest to the target as the MTD
        selectdose=adm.index[selectd]; if (y[selectdose]<=b.e[n[selectdose]]) { is.escalation=1 }
      }
      list(selectdose=selectdose,is.escalation=is.escalation)
  }

  waterfall.subtrial=function(target,p.truetox,p.trueeff,dosespace,npts,ntox,elimi,ncohort,cohortsize,n.earlystop=20,startdose=1,Nmax=1000,
             n.e,y.e,p.saf="default",p.tox="default",cutoff.eli= 0.95,extrasafe=FALSE,offset=0.05,totaln) {
      ndoses1=nrow(p.truetox); ndoses2=ncol(p.truetox)
      p.true=p.truetox[dosespace]; p.true.eff=p.trueeff[dosespace]
      npts=npts; ntox=ntox; elimi=elimi; n.e=n.e; y.e=y.e
      ## if the user does not provide p.saf and p.tox,set them to the default values
      if (p.saf == "default") p.saf=0.6*target;
      if (p.tox == "default") p.tox=1.4*target;
      ## simple error checking
      if (target < 0.05) { cat("Error: the target is too low! \n"); return();}
      if (target>0.6)  { cat("Error: the target is too high! \n"); return(); }
      if ((target-p.saf) < (0.1*target)) { cat( "Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return();}
      if ((p.tox-target) < (0.1*target)) { cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return(); }
      if (offset >= 0.5) { cat("Error: the offset is too large! \n"); return(); }
      ndose=length(p.true);
      selectdose=0; is.escalation=0# store the selected dose level
      ## obtain dose escalation and deescalation boundaries
      temp=get.boundary(target,ncohort=150,cohortsize=1,n.earlystop=100,p.saf,p.tox,cutoff.eli,extrasafe,print=FALSE)
      b.e=temp[2,]; b.d=temp[3,]; b.elim=temp[4,];
      lambda1 =log((1-p.saf)/(1-target))/log(target*(1-p.saf)/(p.saf*(1-target)));
      lambda2 =log((1-target)/(1-p.tox))/log(p.tox*(1-target)/(target*(1-p.tox)));
      #set.seed(6)
      ################## simulate trials ###################
      y <- rep(0,ndose);  n <- rep(0,ndose);
      ye=rep(0,ndose); ne=rep(0,ndose); earlystop=0;         ## indiate whether the trial terminates early
      d=startdose;  elm=rep(0,ndose);       ## starting dose level, and indicate whether doses are eliminated
      for (icohort in 1:ncohort) {
        ### generate toxicity outcome
        y[d]=y[d] + sum(runif(cohortsize) < p.true[d]); n[d]=n[d] + cohortsize;
        ye[d]=ye[d] + sum(runif(cohortsize) < p.true.eff[d]); ne[d]=ne[d] + cohortsize;
        ## determine if the current dose should be eliminated
        if (!is.na(b.elim[n[d]])) { if (y[d] >= b.elim[n[d]]) { elm[d:ndose]=1;
            if (d == 1) { earlystop=1; break; }}
            ## implement the extra safe rule by decreasing the elimination cutoff for the lowest dose
            if (extrasafe) { if (d == 1 && y[1] >= 3) {
             if(1-pbeta(target,y[1] + .05,n[1]-y[1] +0.05)>cutoff.eli-offset) { earlystop=1; break;}
            }}
        }
        ## dose escalation/de-escalation
        if (y[d]<=b.e[n[d]] && d != ndose) {
          if (elm[d + 1] == 0) d=d + 1;
        } else if (y[d] >= b.d[n[d]] && d != 1) {
          d=d-1;
        }else { d=d }
        if (n[d] == n.earlystop) break;
        if (sum(n) >= Nmax) break;
        if ((totaln + sum(n)) >= (ncohort*cohortsize)) break;
      }

      if (earlystop == 1) { selectdose=99; elm=rep(1,ndose)
      }else{
        wsmtd=waterfall.subtrial.mtd(target,n,y,cutoff.eli,extrasafe,offset)
        selectdose=wsmtd$selectdose; is.escalation=wsmtd$is.escalation
      }
      ## output results
      npts[dosespace]=n;  ntox[dosespace]=y;  elimi[dosespace]=elm; n.e[dosespace]=ne; y.e[dosespace]=ye
      list(ncohort=icohort,ntotal=icohort*cohortsize,startdose=startdose,npts=npts,ntox=ntox,
           totaltox=sum(ntox),totaln=sum(npts),n.e=n.e,y.e=y.e,pctearlystop=sum(selectdose==99)*100,
           selectdose=selectdose,is.escalation=is.escalation,elimi=elimi)
    }

  waterfall.phase1=function(p.truetox,p.trueeff,target,ncohort1,cohortsize1,Nmax1=NULL,n.earlystop=10,cutoff.eli=0.95,
                            epsilon=0.03,p.saf="default",p.tox="default",extrasafe=FALSE,offset=0.05,ntrial=1000) {
      # for comparison purpose,we need to record the number of true MTDs and get the nselpercent
      nMTDs=paste(sort(which(abs(p.truetox-target)<=epsilon)),collapse=',')
      if (is.null(Nmax1)) Nmax1=1000
      aa=function(x) as.numeric(as.character(x))
      ###############################  phase I ###############################
      ## if the user does not provide p.saf and p.tox,set them to the default values
      if (p.saf == "default") p.saf=0.6*target; if (p.tox == "default") p.tox=1.4*target;
      if (target < 0.05) { cat("Error: the target is too low! \n"); return(1); }
      if (target>0.6){ cat("Error: the target is too high! \n"); return(1); }
      if ((target-p.saf)<(0.1*target)){ cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return(1);}
      if ((p.tox-target)<(0.1*target)){ cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return(1);}

      # dose levels for agent 1 and 2
      ndoses1 <- nrow(p.truetox); ndoses2 <- ncol(p.truetox)
      ntrial.phase1=NULL; ntrial.mtd=NULL; ntrial.nt=NULL; ntrial.yt=NULL; ntrial.ne=NULL; ntrial.ye=NULL
      ntrial.ONEMTD=NULL;selonemtdpercent=0# report the selection percentage if only one MTD is recalled
      for (trial in 1:ntrial) {
        # run the subtrials sequentially
        trial.result=NULL
        # store toxicity outcome in y and the number of patients in n
        ntox=matrix(rep(0,(ndoses2)*ndoses1),ncol=ndoses2); colnames(ntox)=paste('ntoxDoseB',1:ndoses2,sep='')
        npts=matrix(rep(0,(ndoses2)*ndoses1),ncol=ndoses2); colnames(npts)=paste('nptsDoseB',1:ndoses2,sep='')
        elimi=matrix(0,nrow=ndoses1,ncol=ndoses2); colnames(elimi)=paste('elimiDoseB',1:ndoses2,sep='')
        mtd=cbind('selectdoseA'=1:ndoses1,'selectdoseB'=rep(NA,ndoses1))
        n.e=matrix(rep(0,(ndoses2)*ndoses1),ncol=ndoses2); colnames(n.e)=paste('nefficacy',1:ndoses2,sep='')
        y.e=matrix(rep(0,(ndoses2)*ndoses1),ncol=ndoses2); colnames(y.e)=paste('yefficacy',1:ndoses2,sep='')
        trial.result=data.frame(cbind('trial'=rep(trial,ndoses1),mtd,npts,ntox,elimi,n.e,y.e))
        totaln=0; startdose=1; dosespace=c(1:(ndoses1-1),(1:ndoses2)*ndoses1); subtriali=1
        while (totaln <= ncohort1*cohortsize1) {
          Nmax=Nmax1
          subtrial=waterfall.subtrial(target,p.truetox=p.truetox,p.trueeff=p.trueeff,dosespace=dosespace,npts=npts,ntox=ntox,elimi=elimi,
            ncohort=ncohort1,cohortsize=cohortsize1,n.earlystop=n.earlystop,startdose=startdose,Nmax=Nmax,
            n.e=n.e,y.e=y.e,p.saf='default',p.tox='default',cutoff.eli,extrasafe,offset,totaln=totaln)
          # update dosespace for next subtrial if further subtrials are needed
          selectdose=ifelse(subtrial$selectdose==99,99,dosespace[subtrial$selectdose])
          if(selectdose==99) break
          dj=ifelse(selectdose %% ndoses1 == 0,selectdose %/% ndoses1,selectdose %/% ndoses1 + 1); di=selectdose-(dj-1)*ndoses1
          totaln=aa(subtrial$totaln); npts=subtrial$npts; ntox=subtrial$ntox; n.e=subtrial$n.e; y.e=subtrial$y.e
          elimi=subtrial$elimi
          if ((subtriali == 1) & (selectdose < ndoses1)) {
            for (a in (di + 1):ndoses1) for (b in 1:ndoses2) elimi[a,b]=1
            if (subtrial$is.escalation == 1) {
                startdose=1;dosespace1=c(di + ((2:ndoses2)-1)*ndoses1); Nmax=Nmax1
                subtrial1=waterfall.subtrial(target,p.truetox=p.truetox, p.trueeff=p.trueeff,dosespace=dosespace1,npts=npts,ntox=ntox,elimi=elimi,
                  ncohort=ncohort1,cohortsize=cohortsize1,n.earlystop=n.earlystop,startdose=startdose,Nmax=Nmax,n.e=n.e,y.e=y.e,p.saf='default',
                  p.tox='default',cutoff.eli,extrasafe,offset,totaln=totaln)
                selectdose1=ifelse(subtrial1$selectdose == 99,selectdose,dosespace1[subtrial1$selectdose])
                if (selectdose1 == 99) break
                dj=ifelse(selectdose1 %% ndoses1 == 0,selectdose1 %/% ndoses1,selectdose1 %/% ndoses1+1); di=selectdose1-(dj-1)*ndoses1
                totaln=aa(subtrial1$totaln); npts=subtrial1$npts; ntox=subtrial1$ntox;
                n.e=subtrial1$n.e; y.e=subtrial1$y.e; elimi=subtrial1$elimi
            }
          }
          subtriali=subtriali + 1
          if (dj < ndoses2) elimi[di,(dj + 1):ndoses2]=1
          if (dj == ndoses2) break
          startdose=dj ; dosespace=di-1 + ((2:ndoses2)-1)*ndoses1
          if (di-1 == 0) break;
        }
        npts=t(apply(npts,1,aa)); ntox=t(apply(ntox,1,aa)); elimi=t(apply(elimi,1,aa))
        n.e =t(apply(n.e, 1,aa)); y.e =t(apply(y.e, 1,aa))

        phat=(ntox + 0.05)/(npts + 0.1); phat=t(apply(phat,1,aa)); colnames(phat)=paste('phat',1:ndoses2,sep='')
		phat[elimi==1] = 1.1
        phat=Iso::biviso(phat,npts + 0.1,warn=TRUE)[,]; phat=phat + (1E-5)*(matrix(rep(1:dim(npts)[1],
             each=dim(npts)[2],len=length(npts)),dim(npts)[1],byrow=T)+matrix(rep(1:dim(npts)[2],each=dim(npts)[1],len=length(npts)),dim(npts)[1]))
        colnames(phat)=paste('phat',1:ndoses2,sep='')
        # save mtd information for one simulation based on elimination & phat information
        for (k in ndoses1:1) {
          kn=npts[k,]; ky=ntox[k,]; kelimi=elimi[k,]; kphat=phat[k,]
          if (kelimi[1] == 1 || sum(npts[kelimi == 0]) == 0) { kseldose=99;
          }else{
            adm.set=(kn != 0) & (kelimi == 0); adm.index=which(adm.set == T);
            y.adm=ky[adm.set]; n.adm=kn[adm.set];
            selectd=sort(abs(kphat[adm.set]-target),index.return=T)$ix[1]
            kseldose=adm.index[selectd];
          }
          mtd[k,2]=kseldose
          if(k<ndoses1) if (mtd[k + 1,2]==ndoses2) mtd[k,2]=ndoses2
          #if(k<ndoses1) if(aa(mtd[k+1,2])==aa(mtd[k,2])) mtd[k,2]=99
        }
        # report ONE MTD for comparison
        onephat=aa(as.matrix(phat))  # ONLY report one MTD
        mtdpos=aa(apply(mtd,1,function(x) (x[2]-1)*ndoses1 + x[1]));
        onemtd=99
        if (sum(mtdpos<=ndoses1*ndoses2)>0) {
          onemtdpos=mtdpos[mtdpos<=ndoses1*ndoses2]
          onemtd=onemtdpos[which.min(abs(target-onephat[onemtdpos]))]
        }
        ntrial.ONEMTD=c(ntrial.ONEMTD,onemtd)
        selonemtdpercent=ifelse(is.element(onemtd,which((abs(p.truetox-target)-epsilon)<=1e-4)) == T,selonemtdpercent + 1,selonemtdpercent)
        trial.result[1:ndoses1,grep('nptsDoseB',colnames(trial.result))]=npts; trial.result[1:ndoses1,grep('ntoxDoseB',colnames(trial.result))]=ntox
        trial.result[1:ndoses1,grep('selectdose',colnames(trial.result))]=mtd; trial.result[1:ndoses1,grep('elimiDoseB',colnames(trial.result))]=elimi
        trial.result[1:ndoses1,grep('nefficacy',colnames(trial.result))]=n.e; trial.result[1:ndoses1,grep('yefficacy',colnames(trial.result))]=y.e
        ntrial.mtd=rbind(ntrial.mtd,cbind('trial'=rep(trial,nrow(mtd)),mtd))
        ## save npts and ntox
        ntrial.nt=rbind(ntrial.nt,cbind('trial'=rep(trial,nrow(npts)),npts)); ntrial.yt=rbind(ntrial.yt,cbind('trial'=rep(trial,nrow(ntox)),ntox))
        ntrial.ne=rbind(ntrial.ne,cbind('trial'=rep(trial,nrow(n.e)),n.e));   ntrial.ye=rbind(ntrial.ye,cbind('trial'=rep(trial,nrow(y.e)),y.e))
        trial.result=cbind(trial.result,phat); ntrial.phase1=rbind(ntrial.phase1,trial.result)
      } # end: for ntrial
      selonemtdpercent=round(selonemtdpercent*100/ntrial,1)
      ntrial.mtd=data.frame(ntrial.mtd); colnames(ntrial.mtd)=c('trial','doseA','doseB')
      alltoxpercent=round(sum(sapply(1:ntrial,function(x)ifelse(sum(ntrial.mtd$doseB[ntrial.mtd$trial == x]==99)==ndoses1,1,0)),na.rm=TRUE)*100/ntrial,3)
      stoppercent=round(sum(sapply(1:ntrial,function(x) sum(ntrial.mtd$doseB[ntrial.mtd$trial == x] == 99)))*100/ntrial,1)
      selpercent=matrix(0,nrow=ndoses1,ncol=ndoses2)
      nselpercent=0 # for comparison purpose,compute the accuracy for returning all the true MTDs
      mtdtable=NULL; mtdlist=list()
      for (triali in 1:ntrial) {
        mtddata=unique(as.matrix(ntrial.mtd[ntrial.mtd$trial == triali,2:3]))
        mtddata=mtddata[!is.na(mtddata[,2]),]
        if (length(mtddata)>0) {
          if (length(mtddata) == 2)
            mtddata=matrix(mtddata,ncol=2); mtdlevel=aa(t(apply(mtddata,1,function(x) (x[2]-1)*ndoses1 + x[1])))
            mtdlevel[mtdlevel>ndoses1*ndoses2]=99
          if (sum(mtdlevel<=ndoses1*ndoses2)>0) {
            selpercent[mtdlevel[mtdlevel<=ndoses1*ndoses2]]=selpercent[mtdlevel[mtdlevel<=ndoses1*ndoses2]] + 1
            mtdlist[[triali]]=mtdlevel[mtdlevel<=ndoses1*ndoses2]
            mtdtable=c(mtdtable,paste(sort(aa(mtdlevel[mtdlevel<=ndoses1*ndoses2])),collapse=','))
            if (paste(sort(aa(mtdlevel[mtdlevel<=ndoses1*ndoses2])),collapse=',')==nMTDs)
              nselpercent=nselpercent + 1
          }else{
            mtdlist[[triali]]=99; mtdtable=c(mtdtable,99)
          }
        }
      } # end of triali in ntrial
      selpercent=round(selpercent*100/ntrial,2)
      rownames(selpercent)=paste('DoseA',1:ndoses1,sep=''); colnames(selpercent)=paste('DoseB',1:ndoses2,sep='')

      list(ntrial.phase1=ntrial.phase1,selpercent=selpercent,alltoxpercent=alltoxpercent,stoppercent=stoppercent,mtdlist=mtdlist,
        ntrial.nt=ntrial.nt,ntrial.yt=ntrial.yt,ntrial.mtd=ntrial.mtd,ntrial.ne=ntrial.ne,ntrial.ye=ntrial.ye,
        nselpercent=round(nselpercent*100/ntrial,2), mtdtable=table(mtdtable),ntrial.ONEMTD=ntrial.ONEMTD,
        selonemtdpercent=selonemtdpercent)
    }

  waterfall.phase2=function(target.t,eff.lb,p.true.t,n.t,y.t,p.true.e,k.arm.spaces,n.e,y.e,ncohort2,cohortsize2,
                            cutoff.eli,cutoff.eff,p.saf='default',p.tox='default',extrasafe=FALSE) {
    NPTS2=0
    aa=function(x) as.numeric(as.character(x))
    ## if the user does not provide p.saf and p.tox,set them to the default values
    if(p.saf == "default") p.saf=0.6*target.t; if(p.tox == "default") p.tox=1.4*target.t;
    ## obtain dose escalation and de-escalation boundaries
    temp=get.boundary(target.t,ncohort=150,cohortsize=1,n.earlystop=100,p.saf='default',p.tox='default',print=FALSE);
    b.e=temp[2,]; b.d=temp[3,]; b.elim=temp[4,];  # elimination boundary
    lambda1 =log((1-p.saf)/(1-target.t))/log(target.t*(1-p.saf)/(p.saf*(1-target.t)));
    lambda2 =log((1-target.t)/(1-p.tox))/log(p.tox*(1-target.t)/(target.t*(1-p.tox)));
    # recommended phase II (RPII) doses
    RPII=NULL # RPII searching space is a matrix with ndoses1 x ndoses2
    ndoses1=nrow(p.true.t); ndoses2=ncol(p.true.t)
    # number of DLT and #.of npts treated at each dose level,number of efficacy/total treated patients
    n.t=n.t; y.t=y.t; n.e=n.e; y.e=y.e; n.e.copy=n.e
    if (sum(k.arm.spaces < 99,na.rm=T) == 0) {
        #cat("Error: no arms are left safe for efficacy evaluation! \n");
        return(list( y.t=y.t,n.t=n.t,y.e=y.e,n.e=n.e,RPII=c(99,99)))
    }else{
## k arm, index, dose space,
      skk = cbind('k'=rep(1, ndoses1+ndoses2-1), 'idx'=1:(ndoses1+ndoses2-1), 'dspc'=c(1:(ndoses1-1),(1:ndoses2)*ndoses1))
      if(ndoses1>=2) {
        for(i in 1:(ndoses1-1)) {
          tmpSKK = cbind('k'=rep(i+1, ndoses2-1), 'idx'=1:(ndoses2-1),  'dspc'=c((1:ndoses2)*ndoses1)[-1] - i )
          skk = rbind(skk, tmpSKK)
        }
      }
      skk=data.frame(skk); skk=skk[order(aa(skk$k), aa(skk$idx)),]
      j=ifelse(skk$dspc %% ndoses1==0,skk$dspc%/% ndoses1,skk$dspc%/% ndoses1 + 1)
      i=skk$dspc-(j-1)*ndoses1;
      skk$i = i; skk$j = j; skk=skk[order(aa(skk$k), aa(skk$idx)),]

      k.arm.spaces=k.arm.spaces[k.arm.spaces<99 & !is.na(k.arm.spaces)]; k.arms=length(k.arm.spaces)
      k.arms.j=ifelse(k.arm.spaces %% ndoses1==0,k.arm.spaces %/% ndoses1,k.arm.spaces %/% ndoses1 + 1)
      k.arms.i=k.arm.spaces-(k.arms.j-1)*ndoses1;
      k.arms.ok=rep(1,length(k.arms.j)); k=sample(k.arms,1,1/k.arms);
      # indiate if the trial terminates early,and if doses are eliminated
      for (i in 1:ncohort2) {
        dA=k.arms.i[k]; d=k.arms.j[k];
        selected.arm = aa(skk$k[skk$i==dA & skk$j==d]); arm.dspaces=skk$dspc[which(aa(skk$k)==selected.arm)]

        # generate toxicity outcome
        y.t[dA,d]=y.t[dA,d]+rbinom(1,cohortsize2,p.true.t[dA,d]); n.t[dA,d]=n.t[dA,d]+cohortsize2;
        y.e[dA,d]=y.e[dA,d]+rbinom(1,cohortsize2,p.true.e[dA,d]); n.e[dA,d]=n.e[dA,d] + cohortsize2;
        npts=n.t[dA,d];  ntox=y.t[dA,d]; npts2=n.e[dA,d]; ntox2=y.e[dA,d]
  NPTS2=NPTS2+cohortsize2
        # if fultility of dose (dA,d): Pr(efficacy>eff.lb | D) is too small,then remove arm k
        if(pbeta(eff.lb,ntox2+.05,npts2-ntox2+0.05)>cutoff.eff) k.arms.ok[k]=0;
        # if too toxic, dose de-escalation
        if(ntox>=b.d[npts]){  #pbeta(target.t,ntox+.05,npts-ntox+0.05)<c.a){ #de-escalation
            if(skk$idx[skk$i==dA & skk$j==d]>1){
              k.arms.i[k]=skk$i[which(skk$i==dA & skk$j==d)-1]; k.arms.j[k]=skk$j[which(skk$i==dA & skk$j==d)-1]
            }else{k.arms.ok[k]=0;}
        }
        if(ntox<=b.e[npts]){#pbeta(target.t,ntox+.05,npts-ntox+0.05)>0.8){ #escalation
            if(skk$idx[skk$i==dA & skk$j==d]!=max(skk$idx[skk$i==dA & skk$j==d])){
              k.arms.i[k]=skk$i[which(skk$i==dA & skk$j==d)+1]; k.arms.j[k]=skk$j[which(skk$i==dA & skk$j==d)+1]
            }
        }
        if (sum(k.arms.ok==1)>0) { k=sample(which(k.arms.ok==1),1,1/sum(k.arms.ok==1)) }else{ break; }
      }# end: for ncohort
      # once the maximum sample size is reached,the dose combination that has the highest posterior mean of efficacy is RPII
      tmpk=cbind(k.arms.i,k.arms.j,k.arms.ok)
      tmpk1=apply(tmpk,1,function(x) ifelse(n.e[[x[1],x[2]]] != 0 & x[3] != 0,(y.e[[x[1],x[2]]] + 0.05)/(n.e[[x[1],x[2]]]+0.1),0))
      if (sum(tmpk1) == 0) { RPII=c(99,99)
      }else{ RPII=tmpk[which.max(tmpk1),1:2] }
      return(list(y.t=y.t,n.t=n.t,y.e=y.e,n.e=n.e,RPII=RPII,NPTS2=NPTS2))
    }
  }

  ##########################################################################################################
  JJ=nrow(p.truetox); KK=ncol(p.truetox)
  if(JJ>KK | (nrow(p.trueeff)>ncol(p.trueeff))) {cat("Error: p.truetox and p.trueeff should be arranged in a way (i.e., rotated) such that the number of rows is less than or equal to the number of columns.");
                  return();}
  if(JJ*KK<=4 & ncohort1<=6) {cat("Warning: the sample size is too small, which may lead to poor operating characteristics. Suggest to increase the number of cohort."); }

  if(JJ*KK>4  & ncohort1<=8) {cat("Warning: the sample size is too small, which may lead to poor operating characteristics. Suggest to increase the number of cohort."); }


  if(JJ>KK){ p.truetox = t(p.truetox); p.trueeff = t(p.trueeff)}
  if(is.null(Nmax1)==TRUE) Nmax1 = round(1.5*ncohort1*cohortsize1/min(dim(p.truetox)),0)

  ###############################  phase I ###############################
  ## if the user does not provide p.saf and p.tox,set them to the default values
  if (p.saf == "default") p.saf=0.6*target;
  if (p.tox == "default") p.tox=1.4*target;

  ## simple error checking
  #if(npts[dose.curr[1],dose.curr[2]]==0)  {cat("Error: dose entered is not the current dose \n"); return(1);}
  if (target < 0.05) { cat("Error: the target is too low! \n"); return(1); }
  if (target>0.6)  { cat("Error: the target is too high! \n"); return(1); }
  if ((target-p.saf)<(0.1*target)) { cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return(1)}
  if ((p.tox-target)<(0.1*target)) { cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return(1)}

  # dose levels for agent 1 and 2
  ndoses1 <- nrow(p.truetox); ndoses2 <- ncol(p.truetox)
  aa=function(x) as.numeric(as.character(x))
  ntrial.phase1=cbind('trial'=rep(1:ntrial,each=ndoses1),matrix(0,nrow=ndoses1*ntrial,ncol=ndoses2*6+2));
  colnames(ntrial.phase1)=c("trial","selectdoseA","selectdoseB",paste("nptsDoseB",1:ndoses2,sep=''),
      paste("ntoxDoseB",1:ndoses2,sep=''),paste("elimiDoseB",1:ndoses2,sep=''),paste("nefficacy",1:ndoses2,sep=''),
      paste("yefficacy",1:ndoses2,sep=''), paste("phat",1:ndoses2,sep=''))
  ntrial.phase2=cbind('trial'=1:ntrial,'RPII.i'=rep(99,ntrial),'RPII.j'=rep(99,ntrial));
  ntrial.nt=cbind('trial'=rep(1:ntrial,each=ndoses1),matrix(0,nrow=ntrial*ndoses1,ncol=ndoses2));
  colnames(ntrial.nt)=c("trial",paste("nptsDoseB",1:ndoses2,sep=''))
  ntrial.yt=cbind('trial'=rep(1:ntrial,each=ndoses1),matrix(0,nrow=ntrial*ndoses1,ncol=ndoses2));
  colnames(ntrial.yt)=c("trial",paste("ntoxDoseB",1:ndoses2,sep=''))
  ntrial.ne=cbind('trial'=rep(1:ntrial,each=ndoses1),matrix(0,nrow=ntrial*ndoses1,ncol=ndoses2));
  colnames(ntrial.ne)=c("trial",paste("nptsDoseB",1:ndoses2,sep=''))
  ntrial.ye=cbind('trial'=rep(1:ntrial,each=ndoses1),matrix(0,nrow=ntrial*ndoses1,ncol=ndoses2));
  colnames(ntrial.ye)=c("trial",paste("nptsDoseB",1:ndoses2,sep=''))

  for (trial in 1:ntrial) {
    ###############################  phase I  ###############################
    trial.result=waterfall.phase1(p.truetox,p.trueeff=p.trueeff,target,ncohort1,cohortsize1,Nmax1,n.earlystop,
                                   cutoff.eli=0.95,epsilon=0.03,p.saf,p.tox,extrasafe,offset,ntrial=1)$ntrial.phase1
    trial.result$trial=rep(trial,nrow(trial.result))
    ntrial.phase1[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=as.matrix(trial.result)
    npts=as.matrix(trial.result[,grep('nptsDose',colnames(trial.result))]); ntox=as.matrix(trial.result[,grep('ntoxDose',colnames(trial.result))])
    ne=as.matrix(trial.result[,grep('nefficacy',colnames(trial.result))]); ye=as.matrix(trial.result[,grep('yefficacy',colnames(trial.result))])
    phase1.elimi=trial.result[,grep('elimiDose',colnames(trial.result))]; phase1.mtd=trial.result[,grep('selectdose',colnames(trial.result))]
    phase1.mtd=phase1.mtd[!is.na(phase1.mtd[,2]),]
    ntrial.nt[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind('trial'=rep(trial,nrow(npts)),npts)
    ntrial.yt[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind('trial'=rep(trial,nrow(ntox)),ntox)
    ntrial.ne[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind('trial'=rep(trial,nrow(ne)),ne)
    ntrial.ye[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind('trial'=rep(trial,nrow(ye)),ye)
## K admissible doses
    k.arm.spaces=apply(trial.result[,2:3],1,function(x) (x[2]-1)*ndoses1+x[1])  #which(npts!=0 & pbeta(target,ntox+.05,npts-ntox+0.05)>c.a)
    if (length(k.arm.spaces)>0) {
      phase2=waterfall.phase2(target.t=target,eff.lb=eff.lb,p.true.t=p.truetox,n.t=npts,y.t=ntox,p.trueeff, k.arm.spaces,
        n.e=ne, y.e=ye, ncohort2, cohortsize2, cutoff.eli, cutoff.eff,
        p.saf='default', p.tox='default', extrasafe=FALSE)
      ntrial.phase2[trial,]=c(trial,aa(phase2$RPII))
      ntrial.nt[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind("trial"=rep(trial,nrow(phase2$n.t)),phase2$n.t)
      ntrial.yt[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind("trial"=rep(trial,nrow(phase2$y.t)), phase2$y.t)
      ntrial.ne[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind(rep(trial, ndoses1), phase2$n.e)
      ntrial.ye[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind(rep(trial, ndoses1), phase2$y.e)
NPTS.all = c(NPTS.all, sum(npts) + sum(phase2$NPTS2))
    }else{
      ntrial.phase2[trial,]=c(trial,99,99)
      ntrial.nt[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind("trial"=rep(trial,nrow(npts)),npts);
      ntrial.yt[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind("trial"=rep(trial,nrow(ntox)),ntox);
      ntrial.ne[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind(rep(trial, ndoses1), matrix(0,nrow=ndoses1,ncol=ndoses2))
      ntrial.ye[((trial-1)*ndoses1 + 1):(trial*ndoses1),]=cbind(rep(trial, ndoses1), matrix(0,nrow=ndoses1,ncol=ndoses2))
NPTS.all = c(NPTS.all, sum(npts) + 0)
    }
  } # end: for ntrial
  ntrial.phase2=data.frame(ntrial.phase2)
  colnames(ntrial.phase2)=c('trial', 'RPII.i', 'RPII.j')

  #########################################################################################################################################
  selpercent=matrix(0, nrow=ndoses1, ncol=ndoses2)
  for (i in 1:ntrial) {
    odcdata=as.matrix(ntrial.phase2[ntrial.phase2$trial == i,2:3]);
    odclevel=aa(t(apply(odcdata,1, function(x) (x[2]-1)*ndoses1 + x[1])))
    odclevel[odclevel>ndoses1*ndoses2]=99
    if (sum(odclevel<=ndoses1*ndoses2)>0)
      selpercent[odclevel[odclevel<=ndoses1*ndoses2]]=selpercent[odclevel[odclevel<=ndoses1*ndoses2]] + 1
  }
  selpercent=round(100*selpercent/ntrial,2)
  rownames(selpercent)=paste('DoseA',1:ndoses1,sep=''); colnames(selpercent)=paste('DoseB',1:ndoses2,sep='')
  #########################################################################################################################################
  ### summary of npts, ntox, neff, totaln
  totaln = round(sum(NPTS.all)/ntrial,1)
  npts = matrix(0, nrow=nrow(p.trueeff),ncol=ncol(p.trueeff))
  neff=matrix(0, nrow=nrow(p.trueeff),ncol=ncol(p.trueeff))
  ntox=matrix(0, nrow=nrow(p.truetox),ncol=ncol(p.truetox))
  ph1 = ntrial.phase1
  ph2ne = ntrial.ne; ph2ye=ntrial.ye
  for(i in 1:ntrial) {
    nt = ph1[ph1[,1]==i, grep('nptsDose', colnames(ph1))]
    yt = ph1[ph1[,1]==i, grep('ntoxDose', colnames(ph1))]
    ne1 = ph1[ph1[,1]==i, grep('nefficacy', colnames(ph1))]
    ye1 = ph1[ph1[,1]==i, grep('yefficacy', colnames(ph1))]
    ne2 = ph2ne[ph2ne[,1]==i, 2:4]
    ye2 = ph2ye[ph2ye[,1]==i, 2:4]
    npts = npts + nt + (ne2 - ne1)
    neff = neff+ ye2-ye1
    ntox = ntox + yt
  }
  npts = round(npts/ntrial,1)
  neff = round(neff/ntrial,1)
  ntox = round(ntox/ntrial,1)


## print summary stats: selpercent
  if(JJ<=KK){
  cat("True toxicity rate of dose combinations:\n");
  for (i in 1:dim(p.truetox)[1]) cat(formatC(p.truetox[i,], digits=2, format="f", width=5), sep="  ", "\n");
  cat("\n");
  cat("True efficacy rate of dose combinations:\n");
  for (i in 1:dim(p.trueeff)[1]) cat(formatC(p.trueeff[i,], digits=2, format="f", width=5), sep="  ", "\n");
  cat("\n");

  cat("selection percentage at each dose combination (%):\n");
  for (i in 1:dim(p.truetox)[1]) cat(formatC(selpercent[i,], digits=2, format="f", width=5), sep="  ", "\n");
  cat("\n");

## print summary stats: npts
  cat("number of patients treated at each dose combination:\n");
  for (i in 1:dim(p.truetox)[1]) cat(formatC(npts[i,], digits=2, format="f", width=5), sep ="  ", "\n");
  cat("\n");

## print summary stats: ntoxdose
  cat("number of toxicity observed at each dose combination:\n");
  for (i in 1:dim(p.truetox)[1]) cat(formatC(ntox[i,], digits=2, format="f", width=5), sep ="  ", "\n");
  cat("\n");
  ## print summary stats: neffdose
  cat("number of efficacy observed at each dose combination:\n");
  for (i in 1:dim(p.truetox)[1]) cat(formatC(neff[i,], digits=2, format="f", width=5), sep ="  ", "\n");
  cat("\n");
  cat("total number of patients:", totaln, "\n");
  }else{
  cat("True toxicity rate of dose combinations:\n");
  for (i in 1:dim(p.truetox)[2]) cat(formatC(p.truetox[,i], digits=2, format="f", width=5), sep="  ", "\n");
  cat("\n");
  cat("True efficacy rate of dose combinations:\n");
  for (i in 1:dim(p.trueeff)[2]) cat(formatC(p.trueeff[,i], digits=2, format="f", width=5), sep="  ", "\n");
  cat("\n");

  cat("selection percentage at each dose combination (%):\n");
  for (i in 1:dim(p.truetox)[2]) cat(formatC(selpercent[,i], digits=2, format="f", width=5), sep="  ", "\n");
  cat("\n");

## print summary stats: npts
  cat("number of patients treated at each dose combination:\n");
  for (i in 1:dim(p.truetox)[2]) cat(formatC(npts[,i], digits=2, format="f", width=5), sep ="  ", "\n");
  cat("\n");

## print summary stats: ntoxdose
  cat("number of toxicity observed at each dose combination:\n");
  for (i in 1:dim(p.truetox)[2]) cat(formatC(ntox[,i], digits=2, format="f", width=5), sep ="  ", "\n");
  cat("\n");
  ## print summary stats: neffdose
  cat("number of efficacy observed at each dose combination:\n");
  for (i in 1:dim(p.truetox)[2]) cat(formatC(neff[,i], digits=2, format="f", width=5), sep ="  ", "\n");
  cat("\n");
  cat("total number of patients:", totaln, "\n");
  npts = t(npts); ntox=t(ntox); neff=t(neff); selpercent=t(selpercent); totaln=totaln
  }

  invisible(list(npts=npts, ntox=ntox, neff=neff, selpercent=selpercent, totaln=totaln))
}



