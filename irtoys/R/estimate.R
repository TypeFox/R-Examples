# prepare and run an ICL setup, return parameter estimates
est.icl = function(resp, model, nqp, est.distr, 
  nch, a.prior, b.prior, c.prior, bilog.defaults, run.name) {
  nit = ncol(resp)
  f = paste(run.name,".tcl", sep="")
  d = paste(run.name,".dat", sep="")
  p = paste(run.name,".iclp",sep="")
  write.table(resp,file=d, append=FALSE, sep="", row.names=FALSE, col.names=FALSE, na=".")
  cat("output -log_file",paste(run.name,".iclo",sep=""),"\n",file=f)
  m = switch(model, "1PL"=1, "2PL"=2, "3PL"=3, 4)
  if (m>3) {
    warning(paste("unknown model",model,"using 2PL instead"))
    m = 2
    model = "2PL"
  } 
  cat("set_default_model_dichtomous",model,"\n",file=f,append=TRUE)
  if (!b.prior) 
    cat("options -default_prior_b none\n",file=f,append=TRUE) else
      if (bilog.defaults)
          cat("options -default_prior_b {normal 0.0 2.0}\n",file=f,append=TRUE)
  if (m > 1) {
    if (!a.prior) 
      cat("options -default_prior_a none\n",file=f,append=TRUE) else
        if (bilog.defaults)
            cat("options -default_prior_a {lognormal 0.0 0.5}\n",file=f,append=TRUE)
  }
  if (m > 2) {
    if (!c.prior) 
      cat("options -default_prior_c none\n",file=f,append=TRUE) else
        if (bilog.defaults) {
            prb = 1 /nch  
            cat("options -default_prior_c {beta",20*prb+1,20*(1-prb)+1,"0.0 1.0}\n",file=f,append=TRUE)
        }
  }
  cat("options -D 1.0\n",file=f,append=TRUE)
  cat("allocate_items_dist",nit,"-num_latent_dist_points",nqp,"\n",file=f,append=TRUE)  
  cat("read_examinees",d,paste(nit,"i1",sep=""),"\n",file=f,append=TRUE)
  cat("starting_values_dichotomous\n",file=f,append=TRUE)
  cat("EM_steps -max_iter 2000",file=f,append=TRUE) 
  if (est.distr) cat(" -estim_distr\n",file=f,append=TRUE) else cat("\n",file=f,append=TRUE) 
  cat("write_item_param",p,"\n",file=f,append=TRUE)
  cat("write_latent_dist",paste(run.name,".icld",sep=""),"\n",file=f,append=TRUE)
  cat("release_items_dist\n",file=f,append=TRUE)
  system(paste("icl",f))
  parms = read.ip.icl(p)
  return(parms) 
}

# prepare and run a BILOG setup, return parameter estimates
est.blm = function(resp, model, nqp, est.distr,
  nch, a.prior, b.prior, c.prior, bilog.defaults, run.name, rasch) {
  nit = ncol(resp)
  if (nit>9999) stop("cannot have more than 9999 items")
  f = paste(run.name,".blm", sep="")
  d = paste(run.name,".dat", sep="")
  p = paste(toupper(run.name),".BLMP",sep="")
  np = nrow(resp)
  if (np>999999) stop("cannot have more than 999999 observations")
  resp = cbind(sprintf("%06d",1:np),resp)
  write.table(resp, file=d, append=FALSE, sep="", row.names=FALSE, col.names=FALSE, na=".", quote=FALSE)
  cat("Running Bilog from R\n\n",file=f)
  m = switch(model, "1PL"=1, "2PL"=2, "3PL"=3, 4)
  if (m>3) {
    warning(paste("unknown model",model,"using 2PL instead"))
    m = 2
    model = "2PL"
  }
  cat(">GLOBAL DFName = '",d,"',\n",sep="",file=f,append=TRUE)
  cat("   NPArm = ",m,",\n",sep="",file=f,append=TRUE) 
  cat("   LOGistic\n",file=f,append=TRUE)
  cat("   SAVE;\n",file=f,append=TRUE) 
  cat(">SAVE PARm = '",p,"';\n",sep="",file=f,append=TRUE)
  cat(">LENGTH NITems = (",nit,");\n",sep="",file=f,append=TRUE)   
  cat(">INPUT NTOtal = ",nit,",\n",sep="",file=f,append=TRUE)
  if (m>2) cat("   NALT = ",nch,",\n",sep="",file=f,append=TRUE)
  if (any(is.na(resp))) {
    cat("   NFName = 'myNF.ile'\n",file=f,append=TRUE)
    cat(paste(rep(".",ncol(resp)),collapse=""),"\n",file="myNF.ile")
  }
  cat("   SAMple = ",np,",\n",sep="",file=f,append=TRUE) 
  cat("   NIDchar = 6;\n",file=f,append=TRUE) 
  ifoo = paste("(ITEM",sprintf("%04d",1),"(1)ITEM",sprintf("%04d",nit),")",sep="")
  cat(">ITEMS INAmes = ",ifoo,";\n",sep="",file=f,append=TRUE)
  cat(">TEST1 TNAme = 'TEST',\n",file=f,append=TRUE)  
  ifoo = paste("(",1,"(1)",nit,")",sep="")
  cat("   INUmber = ",ifoo,";\n",sep="",file=f,append=TRUE) 
  cat("(6A1,",nit,"A1)\n",sep="",file=f,append=TRUE) 
  cat(">CALIB NQPt = ",nqp,",\n",sep="",file=f,append=TRUE)
  if (est.distr) cat("   EMPirical,\n",file=f,append=TRUE) 
  if (m==1 && rasch) cat("   RASCH,\n",file=f,append=TRUE) 
  if (b.prior) cat("   TPRior,\n",file=f,append=TRUE) 
  if (m>1 && !a.prior) cat("   NOSprior,\n",file=f,append=TRUE) 
  if (m>2 && !c.prior) cat("   NOGprior,\n",file=f,append=TRUE) 
  cat("   CYCles = 3000,\n",sep="",file=f,append=TRUE)
  cat("   NEWton = 0;\n",file=f,append=TRUE)
  if (Sys.info()["sysname"]=="Linux") {
  	system(paste("wine","BLM1.EXE",run.name))
  	system(paste("wine","BLM2.EXE",run.name))
  	system(paste("wine","BLM3.EXE",run.name))
  } else {
  	system(paste("blm1",run.name))
  	system(paste("blm2",run.name))
  	system(paste("blm3",run.name))
  }
  parms = read.ip.bilog(p)
  return(parms)  
}

# prepare and run an LTM setup, return parameter estimates
est.ltm = function(resp, model, nqp, rasch) {
  nit = ncol(resp)
  switch(model,
    "1PL" = {
      constr = if (rasch) rbind(c(nit+1, 1)) else NULL
      m = rasch(resp, constraint=constr, control = list(GHk = nqp))},
    "2PL" = {m = ltm(resp ~ z1, control = list(GHk = nqp))},
    "3PL" = {m = tpm(resp, control = list(GHk = nqp), max.guessing=1)}, 
    stop(paste("model",model,"not supported in ltm"))
  )
  p = coef(m)
  p = if(model=="3PL") cbind(p[,3], p[,2], p[,1]) else
      cbind(p[,2],p[,1],0)    
  q = summary(m)$coefficients[,2]
  lq = length(q)
  switch(model, `1PL` = se <- cbind(q[lq],q[-lq],0),
     `2PL` = se <- cbind(matrix(q, ncol=2)[,c(2,1)],0),
     `3PL` = se <- matrix(q,ncol=3)[,c(3,2,1)])
  attr(se, "dimnames") = NULL
  return(list(est=p, se=se)) 
}


#' Estimate item parameters
#' 
#' Estimate IRT item parameters using either ICL, BILOG, or \code{ltm}.
#' Provides access to the most widely used options in these programs.
#' 
#' Estimate the parameters of an IRT model defined in the most general case
#' ("3PL") as
#' \deqn{P(U_{ij}=1|\theta_i,a_j,b_j,c_j)=c_j+(1-c_j)\frac{\displaystyle\exp(a_j(\theta_i-b_j))}{1+\displaystyle\exp(a_j(\theta_i-b_j))}}
#' where \eqn{U_{ij}} is a binary response given by person \eqn{i} to item
#' \eqn{j}, \eqn{\theta_i} is the value of the latent variable ("ability") for
#' person \eqn{i}, \eqn{a_j} is the discrimination parameter for item \eqn{j},
#' \eqn{b_j} is the difficulty parameter for item \eqn{j}, \eqn{c_j} is the
#' asymptote for item \eqn{j}.
#'
#' Some authors prefer to represent the model with a logit \eqn{1.7a^*_j(\theta_i-b_j)}
#' rather than \eqn{a_j(\theta_i-b_j)}. This option has been removed from \code{irtoys}
#' as it is not supported by the remaining functions of the package.  
#' 
#' In the 2PL model (\code{model="2PL"}), all asymptotes \eqn{c_j} are 0. In
#' the 1PL model (\code{model="1PL"}), all asymptotes \eqn{c_j} are 0 and the
#' discriminations \eqn{a_j} are equal for all items (and sometimes to 1).
#' 
#' Package \code{irtoys} provides a simple common interface to the estimation
#' of item parameters with three different programs. It only accesses the most
#' basic and widely used options in these programs. Each of the three programs
#' has a much wider choice of options and cababilities, and serious users must
#' still learn the corresponding syntax in order to access the advanced
#' features.  Even when models are fit "by hand", \code{irtoys} may be useful
#' in plotting results, doing comparisons across programs etc.
#' 
#' Estimation of the more complex IRT models (2PL and 3PL) for some "difficult"
#' data sets often has to use prior distributions for the item parameters.
#' \code{irtoys} adopts the default behaviour of BILOG: no priors for \eqn{b}
#' in any model, priors for \eqn{a} in the 2PL and 3PL models, priors for
#' \eqn{c} in the 3PL model. This can be overriden by changing the values of
#' \code{a.prior}, \code{b.prior}, and \code{c.prior}.
#' 
#' If priors are used at all, they will be the same for all items. Note that
#' both ICL and BILOG can, at some additional effort, set different priors for
#' any individual item. At default, the common priors are the BILOG defaults:
#' \code{normal(0,2)} for \eqn{b}, \code{lognormal (0, 0.5)} for \eqn{a}, and
#' \code{beta(20*p+1, 20(1-p)+1)} for \eqn{c}; \eqn{p} is 1 over the number of
#' choices in the original item formulations, which can be set with the
#' parameter \code{nch}, and is again assumed the same for all items.
#' 
#' When \code{engine="icl"} and \code{bilog.defaults=F}, any priors used will
#' be the ICL default ones, and based on the 4-parameter beta distribution:
#' \code{beta(1.01, 1.01, -6, 6)} for \eqn{b}, \code{beta(1.75, 3, 0, 3)} for
#' \eqn{a}, and \code{beta(3.5, 4, 0, 0.5)} for \eqn{c}.  When
#' \code{engine="ltm"}, all commands involving priors are ignored.
#' 
#' \code{est} only works when some IRT software is installed.  Package
#' \code{ltm} is automatically loaded. ICL can be downloaded from
#' \url{www.b-a-h.com}. BILOG is commercial software sold by SSI --- see
#' \url{www.ssicentral.com} for further detail.  On Windows, make sure that the
#' executable files (\code{icl.exe} for ICL, \code{blm1.exe}, \code{blm2.exe},
#' and \code{blm3.exe}, for BILOG) are located in directories that are included
#' in the PATH variable.
#' 
#' @param resp A matrix of responses: persons as rows, items as columns,
#' entries are either 0 or 1, no missing data
#' @param model The IRT model: "1PL", "2PL", or "3PL". Default is "2PL".
#' @param engine One of "icl", "bilog", or "ltm". Default is "icl".
#' @param nqp Number of quadrature points. Default is 20.
#' @param est.distr T if the probabilities of the latent distribution are to be
#' estimated, F if a normal distribution is assumed. Default is F. Ignored when
#' \code{engine="ltm"}.
#' @param nch Number of choices in the original item formulation. Used to
#' determine the prior for the asymptote when \code{engine="bilog"},
#' \code{model="3PL"}, and \code{c.prior=T}. Default is 5.
#' @param a.prior Whether a prior for the item discriminations is used. Ignored
#' when \code{model="1PL"} or \code{engine="ltm"}. Default is T.
#' @param b.prior Whether a prior for the item difficulties is used. Ignored
#' when \code{engine="ltm"}. Default is F.
#' @param c.prior Whether a prior for the asymptotes is used. Ignored when
#' \code{model="1PL"} or \code{model="2PL"} or \code{engine="ltm"}. Default is
#' T.
#' @param bilog.defaults When \code{engine="icl"} and a prior is used, use the
#' default priors in BILOG rather than the default priors in ICL. Ignored when
#' \code{engine="ltm"}. Default is T.
#' @param rasch When \code{engine="bilog"} and \code{model="1PL"} and
#' \code{"rasch"=T}, the common value for discriminations is forced to 1, and
#' the sum of the difficulties is 0. When \code{engine="ltm"} and
#' \code{model="1PL"} and \code{"rasch"=T}, the common value for
#' discriminations is forced to 1. Ignored in all other cases. Default is F.
#' @param run.name A (short) string used in the names of all files read or
#' written by ICL or BILOG. Default is \code{"mymodel"}. Change to something
#' else to keep the outputs of ICL of BILOG for further use. Ignored when
#' \code{engine="ltm"}
#' @return A list with two elements, \code{est} and \code{se}, for the parameter 
#' estimates and their standard errors, correspondingly. Each element is a  
#' matrix with one row per item, and three columns: [,1] item
#' discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and [,3] asymptote
#' \eqn{c}. For the 1PL and 2PL models, all asymptotes are equal to 0; for the
#' 1PL, the discriminations are all equal but not necessarily equal to 1.
#' When ICL is used as estimation engine, \code{se} is NULL as ICL does not
#' compute standard errors for the item parameter estimates.
#' @author Ivailo Partchev
#' @references Bradley A. Hanson (2002), ICL: IRT Command Language.
#' \url{www.b-a-h.com}
#' 
#' Dimitris Rizopoulos (2006). ltm: Latent Trait Models under IRT.
#' \url{cran.r-project.org}
#' 
#' M. F. Zimowski, E. Muraki, R. J. Mislevy and R. D. Bock (1996), BILOG--MG.
#' Multiple-Group IRT Analysis and Test Maintenance for Binary Items, SSI
#' Scientific Software International, Chicago, IL.  \url{www.ssicentral.com}
#' @keywords models
#' @export
#' @examples
#' 
#' p.1pl <- est(Scored, model="1PL", engine="ltm")
#' p.2pl <- est(Scored, model="2PL", engine="ltm")
#' 
est = function(resp, model="2PL", engine="icl", nqp=20, est.distr=FALSE,
  nch=5, a.prior=TRUE, b.prior=FALSE, c.prior=TRUE, 
  bilog.defaults=TRUE, rasch=FALSE, run.name="mymodel") {
  res = switch(engine,
    "icl"=  est.icl(resp, model, nqp, est.distr, nch, a.prior, b.prior, c.prior, bilog.defaults, run.name),
    "bilog"=est.blm(resp, model, nqp, est.distr, nch, a.prior, b.prior, c.prior, bilog.defaults, run.name, rasch),
    "ltm"=  est.ltm(resp, model, nqp, rasch),
    {
      warning(paste("unknown engine",engine,"using icl instead"))
      est.icl(resp,model, nqp, est.distr, nch, a.prior, b.prior, c.prior, bilog.defaults, run.name)      
    }
  )
  return(res)
}   
