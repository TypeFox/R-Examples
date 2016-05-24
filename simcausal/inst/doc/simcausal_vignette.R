## ----include=FALSE, results='hide'------------------------------------------------------
  exportTexOnly=FALSE # do not include any chunks in the final result (R code is still evaluated)

  require(knitr)
  require(simcausal)
  cache_opt <- TRUE
  opts_chunk$set(fig.path='figure/beamer-',fig.align='center',fig.show='hold',size='footnotesize')
  # to crop white space in output figures:
  knit_hooks$set(pdfcrop = hook_pdfcrop)
  # To disable syntax color highlighing of R code in the entire document
  # opts_chunk$set(highlight=FALSE)
  # To change the background color on knitR output from default gray to white:
  # opts_chunk$set(background='white')
  opts_chunk$set(include=!exportTexOnly)

  # options(width=80)  # make the printing fit on the page
  options(width=90)  # make the printing fit on the page
  set.seed(1121)   # make the results repeatable

## ----CreateRefs, include=FALSE, results='hide'------------------------------------------
  # Need to first load that packages that need to be cited:
  # require(simcausal)
  # require(knitr)
  # require(data.table)
  # require(ggplot2)
  # require(ltmle)
  # require(igraph)
  # require(lavaan); require(lavaan.survey); require(sem); require(semPLS); 
  # # require(OpenMx)
  # require(gems); require(aftgee); require(survsim)
  # # # Select packages to cite:
  # citPkgs <- names(sessionInfo()$otherPkgs)
  # citPkgs <- c(citPkgs, "base")
  # # Write the bibtex file:
  # write_bib(citPkgs, file = "R-Pckgs.bib")

## ----chunkMSMsurvplot, eval=TRUE, echo=FALSE--------------------------------------------
    # get MSM survival predictions from the full data.table in long format (melted) by time (t_vec), and by the MSM term (MSMtermName)
    # predictions from the estimated msm model (based on observational data) can be obtained by passing estimated msm model, est.msm
    # given vector of t (t_vec), results of MSM target.eval and an MSM term get survival table by action
   survbyMSMterm <- function(MSMres, t_vec, MSMtermName, use_actions=NULL, est.msm=NULL) {
      library("data.table")
      # look up MSMtermName in MSMterm map, if exists -> use the new name, if doesn't exist use MSMtermName
      if (!is.null(MSMres$S.msm.map)) {
        mapS_exprs <- as.character(MSMres$S.msm.map[,"S_exprs_vec"])
        XMSMterms <- as.character(MSMres$S.msm.map[,"XMSMterms"])
        map_idx <- which(mapS_exprs%in%MSMtermName)
        XMSMtermName <- XMSMterms[map_idx]
        if (!is.null(XMSMtermName)&&length(XMSMtermName)>0) {
          MSMtermName <- XMSMtermName
        }
      }
      print("MSMtermName used"); print(MSMtermName)
      t_dt <- data.table(t=as.integer(t_vec)); setkey(t_dt, t)
      get_predict <- function(actname) {
        # setkey(MSMres$df_long[[actname]], t)
        setkeyv(MSMres$df_long[[actname]], c("t", MSMtermName))
        # MSMterm_vals <- as.numeric(MSMres$df_long[[actname]][t_dt, mult="first"][[MSMtermName]])
        # print("MSMterm_vals"); print(MSMterm_vals)
        MSMterm_vals <- as.numeric(MSMres$df_long[[actname]][t_dt, mult="last"][[MSMtermName]])
        # print("MSMterm_vals last"); print(MSMterm_vals)
        newdata=data.frame(t=t_vec, MSMterm_vals=MSMterm_vals)
        colnames(newdata) <- c("t", MSMtermName)
        # print("newdata"); print(newdata)
        if (!is.null(est.msm)) {
          pred <- predict(est.msm, newdata=newdata, type="response")
        } else {
          pred <- predict(MSMres$m, newdata=newdata, type="response")
        }
        return(data.frame(t=t_vec, pred=pred))
      }
      action_names <- names(MSMres$df_long)
      if (!is.null(use_actions)) {
        action_names <- action_names[action_names%in%use_actions]
      }
      surv <- lapply(action_names, function(actname) {
                                    res <- get_predict(actname)
                                    if (MSMres$hazard) {
                                      res$surv <- cumprod(1-res$pred)  
                                    } else {
                                      res$surv <- 1-res$pred
                                    }
                                    res$pred <- NULL
                                    res$action <- actname
                                    res
                                  })
      names(surv) <- names(MSMres$df_long)
      surv_melt <- do.call('rbind', surv)
      surv_melt$action <- factor(surv_melt$action, levels=unique(surv_melt$action), ordered=TRUE)
      surv_melt
    }
    plotsurvbyMSMterm <- function(surv_melt_dat) {
      library("ggplot2")
      f_ggplot_surv_wS <- ggplot(data= surv_melt_dat, aes(x=t, y=surv)) +
                          geom_line(aes(group = action, color = action), size=.4, linetype="dashed") + 
                          theme_bw()

    }
    plotsurvbyMSMterm_facet <- function(surv_melt_dat1, surv_melt_dat2, msm_names=NULL) {
      library("ggplot2")
      if (is.null(msm_names)) {
        msm_names <- c("MSM1", "MSM2")
      }
      surv_melt_dat1$MSM <- msm_names[1]
      surv_melt_dat2$MSM <- msm_names[2]

      surv_melt_dat <- rbind(surv_melt_dat1,surv_melt_dat2)
      f_ggplot_surv_wS <- ggplot(data= surv_melt_dat, aes(x=t, y=surv)) +
                          geom_line(aes(group = action, color = action), size=.4, linetype="dashed") + 
                          theme_bw() + 
                          facet_wrap( ~ MSM)
    }   

## ----chunk.est.targetMSM, eval=TRUE, echo=FALSE-----------------------------------------
  # @param DAG Object specifying the directed acyclic graph for the observed data, 
    # must have a well-defined MSM target parameter (\code{set.target.MSM()})
  # @param obs_df Simulated observational data
  # @param Aname Generic names of the treatment nodes (can be time-varying)
  # @param Cname Generic names of the censoring nodes (can be time-varying)
  # @param Lnames Generic names of the time-varying covariates (can be time-varying)
  # @param tvec Vector of time points for Y nodes
  # @param actions Which actions (regimens) should be used in estimation from the observed simulated data. 
    # If NULL then all actions that were defined in DAG will be considered.
  # @param package Character vector for R package name to use for estimation. Currently only "ltmle" is implemented.
  # @param fun Character name for R function name to employ for estimation. Currently only "ltmleMSM" is implemented.
  # @param ... Additional named arguments that will be passed on to ltmleMSM function
  est.targetMSM <- function(DAG, obs_df, Aname="A", Cname="C", Lnames, Ytvec, ACLtvec, actions=NULL, package="ltmle", fun="ltmleMSM", ...) {
    outnodes <- attr(DAG, "target")$outnodes
    param_name <- attr(DAG, "target")$param_name
    if (is.null(outnodes$t)) stop("estimation is only implemented for longitudinal data with t defined")
    if (!param_name%in%"MSM") stop("estimation is only implemented for MSM target parameters")
    if (is.null(actions)) {
      message("actions argument underfined, using all available actions")
      actions <- A(DAG)
    }
    # all time points actually used in the observed data
    t_all <- attr(obs_df, "tvals")  
    tvec <- outnodes$t
    t_sel <- ACLtvec
    # ltmle allows for pooling Y's over smaller subset of t's (for example t=(2:8))
    # in this case summary measures HAVE TO MATCH the dimension of finYnodes, not t_sel
    # currently this is not supported, thus, if tvec is a subset of t_sel this will cause an error
    finYnodes <- outnodes$gen_name%+%"_"%+%Ytvec
    Ynodes <- outnodes$gen_name%+%"_"%+%Ytvec
    Anodes <- Aname%+%"_"%+%ACLtvec
    Cnodes <- Cname%+%"_"%+%ACLtvec
    Lnodes <- t(sapply(Lnames, function(Lname) Lname%+%"_"%+%ACLtvec[-1]))
    Lnodes <- as.vector(matrix(Lnodes, nrow=1, ncol=ncol(Lnodes)*length(Lnames), byrow=FALSE))
    Nobs <- nrow(obs_df)
    #------------------------------------------------------------
    # getting MSM params
    #------------------------------------------------------------
    params.MSM <- attr(DAG, "target")$params.MSM
    working.msm <- params.MSM$form
    msm.family <- params.MSM$family
    if (params.MSM$hazard) stop("ltmleMSM cannot estimate hazard MSMs...")
    # the number of attributes and their dimensionality have to match between different actions 
    n_attrs <- length(attr(actions[[1]], "attnames"))
    #------------------------------------------------------------
    # define the final ltmle arrays
    regimens_arr <- array(dim = c(Nobs, length(ACLtvec), length(actions)))
    summeas_arr <- array(dim = c(length(actions), (n_attrs+1), length(Ytvec)))
    # loop over actions (regimes) creating counterfactual mtx of A's for each action:
    for (action_idx in seq(actions)) {
        # I) CREATE COUNTERFACTUAL TREATMENTS & 
        # II) CREATE summary.measure that describes each attribute by time + regimen
        #------------------------------------------------------------    
        # needs to assign observed treatments and replace the action timepoints with counterfactuals
        A_mtx_act <- as.matrix(obs_df[,Anodes]) 
        #------------------------------------------------------------    
        action <- actions[[action_idx]]
         # action-spec. time-points
        t_act <- as.integer(attr(action, "acttimes"))   
        # action-spec. attribute names
        attnames <- attr(action, "attnames")  
        # list of action-spec. attributes
        attrs <- attr(action, "attrs")        
        # time points for which we need to evaluate the counterfactual treatment assignment as determined by action:
        #------------------------------------------------------------
        # Action t's need always be the same subset of t_sel (outcome-based times), otherwise we are in big trouble
        t_act_idx <- which(t_sel%in%t_act)
        t_chg <- t_sel[t_act_idx]
        # modify only A's which are defined in this action out of all Anodes
        As_chg <- Anodes[t_act_idx] 
        #------------------------------------------------------------
        # creates summary measure array that is of dimension (length(t_chg)) - time-points only defined for this action
        # which t's are in the final pooled MSM => need to save the summary measures only for these ts
        t_vec_idx <- which(t_act%in%tvec) 
        summeas_attr <- matrix(nrow=length(attnames), ncol=length(tvec))
        #------------------------------------------------------------
        # extract values of terms in MSM formula: get all attribute values from +action(...,attrs)
        #------------------------------------------------------------
        obs_df_attr <- obs_df # add all action attributes to the observed data
        for (attr_idx in seq(attnames)) { # self-contained loop # grab values of the attributes, # loop over all attributes
          if (length(attrs[[attnames[attr_idx]]])>1) {
            attr_i <- attnames[attr_idx]%+%"_"%+%t_chg
            val_attr_i <- attrs[[attnames[attr_idx]]][t_act_idx]
          } else {
            attr_i <- attnames[attr_idx]
            val_attr_i <- attrs[[attnames[attr_idx]]]
          } 
          # summary measures, for each action/measure
          summeas_attr[attr_idx,] <- matrix(val_attr_i, nrow=1, ncol=length(t_chg))[,t_vec_idx]
          # observed data values of the attribute
          df_attr_i <- matrix(val_attr_i, nrow=Nobs, ncol=length(val_attr_i), byrow=TRUE)
          # create the combined data.frame (attrs + O.dat)
          colnames(df_attr_i) <- attr_i; obs_df_attr <- cbind(data.frame(df_attr_i), obs_df_attr)
        } # end of loop
        summeas_attr <- rbind(summeas_attr, t_chg[t_vec_idx])
        rownames(summeas_attr) <- c(attnames, "t")
        #------------------------------------------------------------
        # add action specific summary measures to the full array
        summeas_arr[action_idx, , ] <- summeas_attr 
        dimnames(summeas_arr)[[2]] <- rownames(summeas_attr)
        #------------------------------------------------------------
        # GENERATING A MATRIX OF COUNTERFACTUAL TREATMENTS:
        for (Achange in As_chg) { # for each A defined in the action, evaluate its value applied to the observed data
          cur.node <- action[As_chg][[Achange]]
          t <- cur.node$t
          newAval <- with(obs_df_attr, {  # for static no need to sample from a distr
            ANCHOR_VARS_OBSDF <- TRUE
            simcausal:::eval_nodeform(as.character(cur.node$dist_params$prob), cur.node)$evaled_expr
          })

          if (length(newAval)==1) {
            newA <- rep(newAval, Nobs)
          } else {
            newA <- newAval
          }
          A_mtx_act[,which(Anodes%in%Achange)] <- newA
        }
        # Result matrix A_mtx_act has all treatments that were defined in that action replaced with 
        # their counterfactual values
        #------------------------------------------------------------
        # add action specific summary measures to the full array
        regimens_arr[, , action_idx] <- A_mtx_act
        #------------------------------------------------------------
    }
    list(regimens_arr=regimens_arr, summeas_arr=summeas_arr)
  }

## ----message=FALSE----------------------------------------------------------------------
library("simcausal")
D <- DAG.empty()
D <- D + 
  node("race",
        distr = "rcategor",
        probs = c(0.5, 0.25, 0.25)) + 
  node("W1",
        distr = "rnorm",
        mean = ifelse(race == 1, 0, ifelse(race == 2, 3, 10)),
        sd = 1) + 
  node("W2",
        distr = "runif",
        min = 0, max = 1) + 
  node("W3",
        distr = "rbern",
        prob = plogis(-0.5 + 0.7 * W1 + 0.3 * W2)) + 
  node("Anode",
        distr = "rbern",
        prob = plogis(-0.5 - 0.3 * W1 - 0.3 * W2 - 0.2 * W3)) + 
  node("Y",
        distr = "rbern",
        prob = plogis(-0.1 + 1.2 * Anode + 0.1 * W1 + 0.3 * W2 + 0.2 * W3))
Dset <- set.DAG(D)

## ----message=FALSE----------------------------------------------------------------------
str(Dset[1])

## ----eval=FALSE-------------------------------------------------------------------------
#  plotDAG(Dset, xjitter = 0.3, yjitter = 0.04,
#          edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
#          vertex_attrs = list(size = 12, label.cex = 0.8))

## ----DAG1t, fig.pos='H', fig.width=10,fig.height=8,out.width='0.8\\linewidth', echo=FALSE, message=FALSE, fig.cap='Graphical representation of the structural equation model using a DAG.', pdfcrop=TRUE----
plotDAG(Dset, xjitter = 0.3, yjitter = 0.03, 
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
        vertex_attrs = list(size = 12, label.cex = 0.8))

## ----message=FALSE----------------------------------------------------------------------
Odat <- sim(DAG = Dset, n = 100, rndseed = 123)

## ----message=FALSE----------------------------------------------------------------------
Odat[1,]

## ----message=FALSE----------------------------------------------------------------------
A1 <- node("Anode", distr = "rbern", prob = 1)
Dset <- Dset + action("A1", nodes = A1)
A0 <- node("Anode", distr = "rbern", prob = 0)
Dset <- Dset + action("A0", nodes = A0)

## ----message=FALSE----------------------------------------------------------------------
names(A(Dset))
class(A(Dset)[["A0"]])

## ---------------------------------------------------------------------------------------
A(Dset)[["A0"]]$Anode

## ----eval=FALSE-------------------------------------------------------------------------
#  str(A(Dset)[["A0"]])

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
Xdat1 <- sim(DAG = Dset, actions = c("A1", "A0"), n = 1000, rndseed = 123)
names(Xdat1)
nrow(Xdat1[["A1"]])
nrow(Xdat1[["A0"]])  

## ----message=FALSE----------------------------------------------------------------------
Xdat1[["A1"]][1, ]
Xdat1[["A0"]][1, ]

## ----message=FALSE----------------------------------------------------------------------
Dset <- set.targetE(Dset, outcome = "Y", param = "A1")

## ----message=FALSE----------------------------------------------------------------------
eval.target(Dset, data = Xdat1)$res

## ----message=FALSE----------------------------------------------------------------------
eval.target(Dset, n = 1000, rndseed = 123)$res

## ----message=FALSE----------------------------------------------------------------------
Dset <- set.targetE(Dset, outcome = "Y", param = "A1-A0")
eval.target(Dset, data = Xdat1)$res

## ----message=FALSE----------------------------------------------------------------------
Dset <- set.targetE(Dset, outcome = "Y", param = "A1/A0")
eval.target(Dset, data = Xdat1)$res

## ----message=FALSE----------------------------------------------------------------------
A1 <- node("Anode", distr = "rbern", prob = d)
Dset <- Dset + action("A1", nodes = A1, d = 1)
A0 <- node("Anode",distr = "rbern", prob = d)
Dset <- Dset + action("A0", nodes = A0, d = 0)

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
msm.form <- "Y ~ d"
Dset <- set.targetMSM(Dset, outcome = "Y", form = msm.form, family = "gaussian")
msm.res <- eval.target(Dset, n = 1000, rndseed = 123)
msm.res$coef

## ----message=FALSE----------------------------------------------------------------------
distr.list()

## ----message=FALSE----------------------------------------------------------------------
rbern

## ----message=FALSE, results='hide'------------------------------------------------------
rnorm_trunc <- function(n, mean, sd, minval = 0) {
  out <- rnorm(n = n, mean = mean, sd = sd)
  minval <- minval[1]
  out[out < minval] <- minval
  out
}

## ----message=FALSE----------------------------------------------------------------------
Dmin0 <- DAG.empty() 

Dmin0 <- Dmin0 + 
  node("W",    distr = "rbern", 
        prob = plogis(-0.5)) +
  node("Anode", distr = "rbern", 
        prob = plogis(-0.5 - 0.3 * W)) + 
  node("Y", distr = "rnorm_trunc", 
        mean = -0.1 + 1.2 * Anode + 0.3 * W, 
        sd = 10)

Dmin0set <- set.DAG(Dmin0)

## ----message=FALSE----------------------------------------------------------------------
Dmin0 <- Dmin0 + 
  node("Y", distr = "rnorm_trunc", 
        mean = -0.1 + 1.2 * Anode + 0.3 * W, 
        sd = 10, 
        minval = 10)

Dmin10set <- set.DAG(Dmin0)

## ----message=FALSE----------------------------------------------------------------------
Dmin0 <- Dmin0 + 
  node("Y", distr = "rnorm_trunc", 
        mean = -0.1 + 1.2 * Anode + 0.3 * W, 
        sd = 10, 
        minval = ifelse(Anode == 0, 5, 10))

Dminset <- set.DAG(Dmin0)

## ----message=FALSE----------------------------------------------------------------------
vecfun.all.print()

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
power2 <- function(arg) arg^2

D <- DAG.empty()
D <- D + 
  node("W1", distr = "rnorm", 
        mean = 0, sd = 1) + 
  node("W2", distr = "rnorm", 
        mean = 0, sd = 1) + 
  node("W3", distr = "rnorm", 
        mean = power2(W1), sd = power2(W2))

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
D1 <- set.DAG(D)
(tnonvec <- system.time(sim1nonvec <- simobs(D1, n = 1000000, rndseed = 123)))

## ----message=FALSE----------------------------------------------------------------------
vecfun.add(c("power2"))
D1vec <- set.DAG(D)
(tvec <- system.time(sim1vec <- simobs(D1vec, n = 1000000, rndseed = 123)))
all.equal(sim1nonvec,sim1vec)

## ----message=FALSE----------------------------------------------------------------------
vecfun.print()
vecfun.reset()
vecfun.print()

## ----message=FALSE----------------------------------------------------------------------
vecfun.reset()  
ifelse1 <- function(arg) {
  ifelse(arg[1], arg[2], arg[3])
}

D <- DAG.empty()
D <- D + 
  node("W1", distr = "rbern", 
        prob = 0.05) +
  node("W2", distr = "rbern", 
        prob = ifelse1(c(W1, 0.5, 0.1)))

D2nonvec <- set.DAG(D)

## ----message=FALSE----------------------------------------------------------------------
ifelse2 <- function(arg, val1, val2) {
  ifelse(arg, val1, val2)
}

vecfun.add(c("ifelse2"))

D <- DAG.empty()
D <- D + 
  node("W1", distr = "rbern", 
        prob = 0.05) +
  node("W2", distr = "rbern", 
        prob = ifelse2(W1, 0.5, 0.1))
D2vec <- set.DAG(D)

## ----message=FALSE----------------------------------------------------------------------
  (t2nonvec <- system.time(sim2nonvec <- simobs(D2nonvec, n = 100000, rndseed = 123)))
  (t2vec <- system.time(sim2vec <- simobs(D2vec, n = 100000, rndseed = 123)))
  all(unlist(lapply(seq(ncol(sim2nonvec)), 
              function(coli) all.equal(sim2nonvec[, coli], sim2vec[, coli]))))

## ----message=FALSE----------------------------------------------------------------------
library("simcausal")
options(simcausal.verbose=FALSE)

D <- DAG.empty()
D <- D + 
  node("L2", t = 0, distr = "rbern", 
      prob = 0.05) +
  node("L1", t = 0, distr = "rbern", 
      prob = ifelse(L2[0] == 1, 0.5, 0.1)) +
  node("A1", t = 0, distr = "rbern", 
      prob = ifelse(L1[0] == 1 & L2[0] == 0, 0.5,
                ifelse(L1[0] == 0 & L2[0] == 0, 0.1,
                  ifelse(L1[0] == 1 & L2[0] == 1, 0.9, 0.5)))) +
  node("A2", t = 0, distr = "rbern", 
      prob = 0, EFU = TRUE)

## ----message=FALSE----------------------------------------------------------------------
t.end <- 16
D <- D + 
  node("Y", t = 1:t.end, distr = "rbern",
      prob = 
        plogis(-6.5 + L1[0] + 4 * L2[t-1] + 0.05 * sum(I(L2[0:(t-1)] == rep(0, t)))), 
      EFU = TRUE) +
  node("L2", t = 1:t.end, distr = "rbern",
      prob = 
        ifelse(A1[t-1] == 1, 0.1, 
          ifelse(L2[t-1] == 1, 0.9, min(1, 0.1 + t / 16)))) +
  node("A1", t = 1:t.end, distr = "rbern", 
      prob = ifelse(A1[t-1] == 1, 1,
                ifelse(L1[0] == 1 & L2[t] == 0, 0.3,
                  ifelse(L1[0] == 0 & L2[t] == 0, 0.1,
                    ifelse(L1[0] == 1 & L2[t] == 1, 0.7, 0.5))))) +
  node("A2", t = 1:t.end, distr = "rbern", 
      prob = {if(t == 16) {1} else {0}},
      EFU = TRUE)
lDAG <- set.DAG(D)

## ----eval=FALSE-------------------------------------------------------------------------
#  plotDAG(lDAG, xjitter = 0.3, yjitter = 0.01)

## ----eval=FALSE-------------------------------------------------------------------------
#  plotDAG(lDAG, tmax = 3, xjitter = 0.3, yjitter = 0.03,
#        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
#        vertex_attrs = list(size = 12, label.cex = 0.8))

## ----DAGlong, echo=FALSE, message=FALSE, fig.pos='H', out.width='0.8\\linewidth', fig.cap='Graphical representation of the structural equation model using a DAG', pdfcrop=TRUE----
plotDAG(lDAG, xjitter = 0.3, yjitter = 0.01)

## ----DAGlongtmax3, echo=FALSE, message=FALSE, fig.pos='H', out.width='0.8\\linewidth', fig.cap='Graphical representation of a portion of the structural equation model using a DAG. Only the nodes indexed by time points lower than or equal to 3 are represented.', pdfcrop=TRUE----
plotDAG(lDAG, tmax = 3, xjitter = 0.3, yjitter = 0.03,
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8), 
        vertex_attrs = list(size = 12, label.cex = 0.8))

## ----message=FALSE----------------------------------------------------------------------
Odat <- sim(DAG = lDAG, n = 100, rndseed = 123)

## ---------------------------------------------------------------------------------------
Odat[1,1:10]

## ---------------------------------------------------------------------------------------
act_theta <-c(node("A1", t = 0, distr = "rbern", 
                    prob = ifelse(L2[0] >= theta , 1, 0)),
              node("A1", t = 1:(t.end), distr = "rbern", 
                    prob = ifelse(A1[t-1] == 1, 1, ifelse(L2[t] >= theta, 1, 0))))

## ----message=FALSE----------------------------------------------------------------------
Ddyn <- lDAG
Ddyn <- Ddyn + action("A1_th0", nodes = act_theta, theta = 0)
Ddyn <- Ddyn + action("A1_th1", nodes = act_theta, theta = 1)

## ---------------------------------------------------------------------------------------
class(A(Ddyn)[["A1_th0"]])
A(Ddyn)[["A1_th0"]]

## ----eval=FALSE-------------------------------------------------------------------------
#  plotDAG(A(Ddyn)[["A1_th0"]], tmax = 3, xjitter = 0.3, yjitter = 0.03,
#        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
#        vertex_attrs = list(size = 15, label.cex = 0.7))

## ----actDAGlongtmax3, echo=FALSE, fig.pos='H', message=FALSE, out.width='0.8\\linewidth', fig.cap='Graphical representation of the modified structural equation model resulting from a dynamic intervention. Only the nodes indexed by time points lower than or equal to 3 are represented.', pdfcrop=TRUE----
plotDAG(A(Ddyn)[["A1_th0"]], tmax = 3, xjitter = 0.3, yjitter = 0.03,
      edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
      vertex_attrs = list(size = 15, label.cex = 0.7))

## ----message=FALSE----------------------------------------------------------------------
A(Ddyn)[["A1_th0"]]$A1_0
Ddyntry <- Ddyn + 
  action("A1_th0", nodes = node("A1", t = 0, distr = "rbern", prob = 0))
A(Ddyntry)[["A1_th0"]]$A1_0

## ----message=FALSE----------------------------------------------------------------------
A(Ddyntry)[["A1_th0"]]
Ddyntry <- Ddyntry + 
  action("A1_th0", nodes = act_theta, theta = 1, newparam = 100)
A(Ddyntry)[["A1_th0"]]

## ----message=FALSE, warning=FALSE-------------------------------------------------------
act_theta_t <-c(node("A1",t = 0, distr = "rbern", 
                    prob = ifelse(L2[0] >=  theta[t], 1, 0)),
                node("A1",t = 1:t.end, distr = "rbern", 
                    prob = ifelse(A1[t-1]==1, 1, ifelse(L2[t] >= theta[t], 1, 0)))
                )
Ddyntry <- Ddyntry + action("A1_th0", nodes = act_theta_t, theta = rep(0,(t.end)+1))
A(Ddyntry)[["A1_th0"]]

## ----message=FALSE----------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
Dstat <- lDAG
act_A1_tswitch <- node("A1",t = 0:(t.end), distr = "rbern", 
                      prob = ifelse(t >= tswitch, 1, 0))

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
tswitch_vec <- (0:t.end)
for (tswitch_i in tswitch_vec) {
  abar <- rep(0, length(tswitch_vec))
  abar[which(tswitch_vec >= tswitch_i)] <- 1
  Dstat <- Dstat + action("A1_ts"%+%tswitch_i, 
                          nodes = act_A1_tswitch,
                          tswitch = tswitch_i,
                          abar = abar)
}

## ----message=FALSE----------------------------------------------------------------------
A(Dstat)[["A1_ts3"]]

## ----eval=FALSE-------------------------------------------------------------------------
#  plotDAG(A(Dstat)[["A1_ts3"]], tmax = 3, xjitter = 0.3, yjitter = 0.03,
#        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
#        vertex_attrs = list(size = 15, label.cex = 0.7), excludeattrs = "abar")

## ----act2DAGlongtmax3, echo=FALSE,fig.pos='H', message=FALSE, out.width='0.8\\linewidth', fig.cap='Graphical representation of the modified structural equation model resulting from a static intervention. Only the nodes indexed by time points lower than or equal to 3 are represented. The action attribute \\code{abar} is also not represented.', pdfcrop=TRUE----
plotDAG(A(Dstat)[["A1_ts3"]], tmax = 3, xjitter = 0.3, yjitter = 0.03,
      edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
      vertex_attrs = list(size = 15, label.cex = 0.7), excludeattrs = "abar")

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
Xdyn <- sim(Ddyn, actions = c("A1_th0", "A1_th1"), n = 100, rndseed = 123)
nrow(Xdyn[["A1_th0"]])
nrow(Xdyn[["A1_th1"]])
names(Xdyn)

## ----message=FALSE----------------------------------------------------------------------
Xdyn[["A1_th0"]][1, 1:15]
Xdyn[["A1_th1"]][1, 1:15]

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
Xstat <- sim(Dstat, actions = names(A(Dstat)), n = 100, rndseed = 123)
length(Xstat)
nrow(Xstat[["A1_ts3"]])

## ----message=FALSE----------------------------------------------------------------------
Xstat[["A1_ts3"]][1, ]

## ----message=FALSE----------------------------------------------------------------------
Odat.wide <- sim(DAG = lDAG, n = 100, wide = TRUE, rndseed = 123)
Odat.wide[1:2, 1:18]
Odat.long <- sim(DAG = lDAG, n = 100, wide = FALSE, rndseed = 123)
Odat.long[1:5, ]

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
lXdyn <- sim(Ddyn, actions = c("A1_th0", "A1_th1"), n = 1000, wide = FALSE, rndseed = 123)
head(lXdyn[["A1_th0"]], 5)

## ----message=FALSE----------------------------------------------------------------------
Odat.long2 <- DF.to.longDT(Odat.wide)
Odat.long2[1:5, ]

## ----message=FALSE,warning=FALSE--------------------------------------------------------
Odat.wide <- sim(DAG = lDAG, n = 1000, rndseed = 123)
Odat.wide[c(11,76), 1:18]
Odat.wideLTCF <- sim(DAG = lDAG, n = 1000, LTCF = "Y", rndseed = 123)
Odat.wideLTCF[c(11,76), 1:18]

## ----message=FALSE,warning=FALSE--------------------------------------------------------
Odat.wideLTCF2 <- doLTCF(data = Odat.wide, LTCF = "Y")
Odat.wideLTCF2[c(11,76), 1:18]

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
Ddyn <- set.targetE(Ddyn, outcome = "Y", t = 1:16, param = "A1_th1")
surv_th1 <- 1 - eval.target(Ddyn, data = Xdyn)$res
Ddyn <- set.targetE(Ddyn, outcome = "Y", t = 1:16, param = "A1_th0");
surv_th0 <- 1 - eval.target(Ddyn, data = Xdyn)$res

## ----survfig1, fig.pos='htb!', fig.width=8,fig.height=7,out.width='.6\\linewidth', message=FALSE, fig.cap='Estimates of the true survival curves under the two dynamic interventions'----
plotSurvEst(surv = list(d_theta1 = surv_th1, d_theta0 = surv_th0),
            xindx = 1:17,
            ylab = "Counterfactual survival for each intervention",
            ylim = c(0.75,1.0))

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
Ddyn <- set.targetE(Ddyn, outcome = "Y", t = 12, param = "A1_th1-A1_th0")
(psi <- eval.target(Ddyn, data = Xdyn)$res)

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
Ddyn <- set.targetE(Ddyn, outcome = "Y", t = 12, param = "A1_th0/A1_th1")
eval.target(Ddyn, data = Xdyn)$res

## ----message=FALSE----------------------------------------------------------------------
msm.form <- "Y ~ theta + t + I(theta*t)"
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, form = msm.form, 
                      family = "binomial", hazard = FALSE)

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
MSMres <- eval.target(Ddyn, n = 1000, rndseed = 123)
MSMres$coef

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
XdynLTCF <- lapply(Xdyn, doLTCF, LTCF = "Y")
eval.target(Ddyn, data = XdynLTCF)$coef

## ----survfig2, fig.pos='htb!', fig.width=8,fig.height=7,out.width='.6\\linewidth', message=FALSE, fig.cap='Survival curve estimates evaluated based on working MSM 1'----
surv_th0 <- 1 - predict(MSMres$m, newdata = data.frame(theta = rep(0, 16), t = 1:16), 
                        type = "response")
surv_th1 <- 1 - predict(MSMres$m, newdata = data.frame(theta = rep(1, 16), t = 1:16), 
                        type = "response")
plotSurvEst(surv = list(MSM_theta1 = surv_th1, MSM_theta0 = surv_th0),
            xindx = 1:16,
            ylab = "MSM Survival, P(T>t)",
            ylim = c(0.75, 1.0))

## ----message=FALSE, cache=cache_opt-----------------------------------------------------

msm.form <- "Y ~ theta + t + I(t^2) + I(t^3) + I(t^4) + I(t^5) + I(t*theta) + I(t^2*theta) + 
              I(t^3*theta) + I(t^4*theta) + I(t^5*theta)"

Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, formula = msm.form, 
                      family = "binomial", hazard = FALSE)
MSMres2 <- eval.target(Ddyn, n = 1000, rndseed = 123)
MSMres2$coef

## ----survfig3, fig.pos='htb!', fig.width=8,fig.height=7,out.width='.6\\linewidth', message=FALSE, fig.cap='Survival curve estimates evaluated based on working MSM 2'----
surv_th0 <- 1 - predict(MSMres2$m, newdata = data.frame(theta = rep(0, 16), t = 1:16),
                        type = "response")
surv_th1 <- 1 - predict(MSMres2$m, newdata = data.frame(theta = rep(1, 16), t = 1:16),
                        type = "response")
plotSurvEst(surv = list(MSM_theta1 = surv_th1, MSM_theta0 = surv_th0),
            xindx = 1:16,
            ylab = "MSM Survival, P(T>t)",
            ylim = c(0.75, 1.0))  

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
msm.form <- "Y ~ theta + as.factor(t) + as.factor(t):theta "
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, formula = msm.form, 
                      family = "binomial", hazard = FALSE)
MSMres3 <- eval.target(Ddyn, n = 1000, rndseed = 123)
MSMres3$coef

## ----survfig4, fig.pos='htb!', fig.width=8,fig.height=7,out.width='.6\\linewidth', message=FALSE, fig.cap='Survival curve estimates evaluated based on working MSM 3'----
surv_th0 <- 1 - predict(MSMres3$m, newdata = data.frame(theta = rep(0, 16), t = 1:16),
                        type = "response")
surv_th1 <- 1 - predict(MSMres3$m, newdata = data.frame(theta = rep(1, 16), t = 1:16),
                        type = "response")
plotSurvEst(surv = list(MSM_theta1 = surv_th1, MSM_theta0 = surv_th0),
            xindx = 1:16,
            ylab = "MSM Survival, P(T>t)",
            ylim = c(0.75, 1.0))

## ---------------------------------------------------------------------------------------
msm.form <- "Y ~ theta + t + I(theta*t)"
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, form = msm.form, 
                      family = "binomial", hazard = TRUE)

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
MSMres <- eval.target(Ddyn, n = 1000, rndseed = 123)
MSMres$coef

## ----survfig5, fig.pos='htb!', fig.width=8,fig.height=7,out.width='.6\\linewidth', message=FALSE, fig.cap='Hazard estimates evaluated based on working MSM 4'----
h_th0 <- predict(MSMres$m, newdata = data.frame(theta = rep(0, 16), t = 1:16), 
                type = "response")
h_th1 <- predict(MSMres$m, newdata = data.frame(theta = rep(1, 16), t = 1:16), 
                type = "response")
plotSurvEst(surv = list(MSM_theta1 = h_th1, MSM_theta0 = h_th0),
            xindx = 1:16,
            ylab = "MSM hazard function",
            ylim = c(0.0, 0.03))

## ----survfig6, fig.pos='htb!', fig.width=8,fig.height=7,out.width='.6\\linewidth', message=FALSE, fig.cap='Survival curve estimates evaluated based on working MSM 4'----
Surv_h_th0 <- cumprod(1 - h_th0)
Surv_h_th1 <- cumprod(1 - h_th1)
plotSurvEst(surv = list(MSM_theta1 = Surv_h_th1, MSM_theta0 = Surv_h_th0),
            xindx = 1:16,
            ylab = "Survival P(T>t) derived from MSM hazard",
            ylim = c(0.75, 1.0))

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
msm.form_sum <- "Y ~ theta + t + I(theta*t) + I(theta*L1)"
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, form = msm.form_sum, 
                      family = "binomial", hazard = TRUE)
MSMres <- eval.target(Ddyn, n = 1000, rndseed = 123)
MSMres$coef

## ---------------------------------------------------------------------------------------
get_haz <- function(thetaval, L1val) {
  predict(MSMres$m, newdata = data.frame(theta = rep(thetaval, 16), t = 1:16, L1 = L1val),
          type = "response")
}
Sth0L1_0 <- cumprod(1 - get_haz(thetaval = 0, L1val = 0))
Sth1L1_0 <- cumprod(1 - get_haz(thetaval = 1, L1val = 0))
Sth0L1_1 <- cumprod(1 - get_haz(thetaval = 0, L1val = 1))
Sth1L1_1 <- cumprod(1 - get_haz(thetaval = 1, L1val = 1))

## ----survfig7, message=FALSE, fig.pos='htb!', fig.width=16,fig.height=7,out.width='1\\linewidth', fig.cap='Survival curve estimates evaluated based on working MSM 5'----
par(mfrow = c(1,2))
plotSurvEst(surv = list(MSM_theta1 = Sth1L1_0, MSM_theta0 = Sth0L1_0),
            xindx = 1:16,
            ylab = "Survival P(T>t), for L1=0",
            ylim = c(0.5, 1.0))
plotSurvEst(surv = list(MSM_theta1 = Sth1L1_1, MSM_theta0 = Sth0L1_1),
            xindx = 1:16,
            ylab = "Survival P(T>t), for L1=1",
            ylim = c(0.5, 1.0))

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
msm.form.correct <- "Y ~ theta + t + I(theta*t) + I(theta * S(L2[0]))"
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = 1:16, form = msm.form.correct, 
                      family = "binomial", hazard = TRUE)
MSMres.correct <- eval.target(Ddyn, n = 1000, rndseed = 123)
MSMres.correct$coef

## ----message=FALSE, results='hide', cache=cache_opt-------------------------------------
Xts <- sim(Dstat, actions = names(A(Dstat)), n = 1000, rndseed = 123)

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
msm.form_1 <- "Y ~  t + S(mean(abar[0:(t-1)])) + I(t*S(mean(abar[0:(t-1)])))"
Dstat <- set.targetMSM(Dstat, outcome = "Y", t = 1:16, form = msm.form_1, 
                        family = "binomial", hazard = TRUE)
MSMres <- eval.target(Dstat, data = Xts)
MSMres$coef

## ----message=FALSE, cache=cache_opt-----------------------------------------------------
names(MSMres)
MSMres$S.msm.map
names(MSMres$df_long)
MSMres$df_long[["A1_ts2"]]

## ----survfig8, fig.pos='htb!', fig.width=10, fig.height=5, message=FALSE, fig.cap='Survival curve estimates evaluated based on working MSM 7'----
survMSMh_wS <- survbyMSMterm(MSMres = MSMres, t_vec = 1:16, 
                              MSMtermName = "mean(abar[0:(t - 1)])")
print(plotsurvbyMSMterm(survMSMh_wS))

## ----message=FALSE----------------------------------------------------------------------
t.end <- 12
D <- DAG.empty()
D <- D + 
        node("L2", t = 0, distr = "rbern", 
              prob = 0.05) +
        node("m1L2", t = 0, distr = "rconst", 
              const = 1 - L2[0]) +
        node("L1", t = 0, distr = "rbern", 
              prob = ifelse(L2[0] == 1, 0.8, 0.3)) +
        node("A1", t = 0, distr = "rbern",
              prob = ifelse(L1[0] == 1 & L2[0] == 0, 0.5,
                      ifelse(L1[0] == 0 & L2[0] == 0, 0.2,
                        ifelse(L1[0] == 1 & L2[0] == 1, 0.8, 0.5)))) +
        node("A2", t = 0, distr = "rbern", prob = 0, EFU = TRUE)

D <- D + 
        node("Y", t = 1:t.end, distr = "rbern",
              prob = plogis(-7 + 3 * L1[0] + 5 * L2[t-1] + 
                            0.1 * sum(I(L2[0:(t-1)] == rep(0, t)))),
              EFU = TRUE) +
        node("L2", t = 1:t.end, distr = "rbern",
              prob = ifelse(A1[t-1] == 1, 0.1,
                      ifelse(L2[t-1] == 1, 0.9, min(1,0.1 + t/16)))) +
        node("m1L2", t = 1:t.end, distr = "rconst", const = 1-L2[t]) +
        node("A1", t = 1:t.end, distr = "rbern",
              prob = ifelse(A1[t-1] == 1, 1,
                      ifelse(L1[0] == 1 & L2[t] == 0, 0.4,
                        ifelse(L1[0] == 0 & L2[t] == 0, 0.2,
                          ifelse(L1[0] == 1 & L2[t] == 1, 0.8, 0.6))))) +
        node("A2", t = 1:t.end, distr = "rbern",
              prob = {if(t == .(t.end)) {1} else {0}},
              EFU = TRUE)

lDAG <- set.DAG(D)
Ddyn <- lDAG

## ----message=FALSE, echo=FALSE----------------------------------------------------------
anodes <- node("A1", t = (0:t.end), distr = "rbern", 
                prob = {if (t == 0) {ifelse(L2[0] >= theta, 1, 0)} else
                        {ifelse(A1[t-1] == 1, 1, ifelse(L2[t] >= theta, 1, 0))}})
Ddyn <- Ddyn + 
              action("A1_th0", nodes = anodes, theta = 0) + 
              action("A1_th1", nodes = anodes, theta = 1)

## ----message=FALSE, eval=TRUE-----------------------------------------------------------
t0 <- 12
Ddyn <- set.targetE(Ddyn, outcome = "Y", t = (1:t0), param = "A1_th1-A1_th0")

getNP.truetarget <- function() {
  resNP <- eval.target(Ddyn, n = 150000, rndseed = 123)$res
  return(as.vector(resNP[paste0("Diff_Y_", t0)]))
}

f1name <- "vignette_dat/repstudy1_psi0.t0.NP.Rdata"
if (file.exists(f1name)) {
  load(f1name)
} else {
  psi0.t0.NP <- getNP.truetarget()
  save(list = "psi0.t0.NP", file = f1name)
}
psi0.t0.NP

## ----message=FALSE, eval=TRUE-----------------------------------------------------------
MSM_RD_t <- function(resMSM, t) {
  invlogit <- function(x) 1 / (1 + exp(-x))
  Riskth0 <- invlogit(resMSM["(Intercept)"] + resMSM[paste0("as.factor(t)",t)])
  Riskth1 <- invlogit(resMSM["(Intercept)"] + resMSM[paste0("as.factor(t)",t)] +
                          resMSM["theta"] + resMSM[paste0("theta:as.factor(t)",t)])
  return(as.vector(Riskth1-Riskth0))
}

msm.form <- "Y ~ theta + as.factor(t) + as.factor(t):theta "
Ddyn <- set.targetMSM(Ddyn, outcome = "Y", t = (1:t0), formula = msm.form, 
                      family = "binomial", hazard = FALSE)

getMSM.truetarget <- function() {
  resMSM <- eval.target(Ddyn, n = 150000, rndseed = 123)$coef
  return(as.vector(MSM_RD_t(resMSM = resMSM, t = t0)))
}

f2name <- "vignette_dat/repstudy1_psi0.t0.MSM.Rdata"
if (file.exists(f2name)) {
  load(f2name)
} else {
  psi0.t0.MSM <- getMSM.truetarget()
  save(list = "psi0.t0.MSM", file = f2name)
}

psi0.t0.MSM
all.equal(psi0.t0.NP, psi0.t0.MSM)

## ----chunk.simrunltmleMSM1, eval=TRUE, echo=FALSE---------------------------------------
times <- c(0:(t.end-1))
gforms <- c("A1_0 ~ L1_0 + L2_0","A2_0 ~ L1_0")
timesm0 <- times[which(times > 0)]
# correctly specified g:
gforms <- c("A1_0 ~ L1_0 + L2_0","A2_0 ~ L1_0")
gformm0 <- as.vector(sapply(timesm0, function(t) 
              c("A1_"%+%t%+%" ~ A1_"%+%(t-1)%+%" + L1_0"%+%" + L2_"%+%t%+%" + I(L2_"%+%t%+%"*A1_"%+%(t-1)%+%")+ I(L1_0*A1_"%+%(t-1)%+%")",
                "A2_"%+%t%+%" ~ L1_0")))
gforms <- c(gforms, gformm0)
# mis-specified g (no TV covar L2):
gforms_miss <- c("A1_0 ~ L1_0","A2_0 ~ L1_0")
gformm0_miss <- as.vector(sapply(timesm0, function(t) 
                              c("A1_"%+%t%+%" ~ L1_0*A1_"%+%(t-1),
                                "A2_"%+%t%+%" ~ L1_0")))
gforms_miss <- c(gforms_miss, gformm0_miss)
Qformallt <- "Q.kplus1 ~ L1_0"
Lterms <- function(var, tlast){
  tstr <- c(0:tlast)
  strout <- paste(var%+%"_"%+%tstr, collapse = " + ")
  return(strout)
}
tY <- (0:11)
Ynames <- paste("Y_"%+%c(tY+1))
Qforms <- unlist(lapply(tY, function(t)  {
                  a <- Qformallt%+%" + I("%+%Lterms("m1L2",t)%+%") + "%+%Lterms("L2",t)
                  return(a)
  }))
names(Qforms) <- Ynames
survivalOutcome <- TRUE
stratify_Qg <- TRUE
mhte.iptw <- TRUE
Anodesnew <- "A1_"%+%(0:(t.end-1))
Cnodesnew <- "A2_"%+%(0:(t.end-1))
L2nodesnew <- "L2_"%+%(1:(t.end-1))
mL2nodesnew <- "m1L2_"%+%(1:(t.end-1))
Lnodesnew <- as.vector(rbind(L2nodesnew, mL2nodesnew))
Ynodesnew <- "Y_"%+%(1:t.end)
finYnodesnew <- Ynodesnew
dropnms <- c("ID","L2_"%+%t.end,"m1L2_"%+%t.end, "A1_"%+%t.end, "A2_"%+%t.end)
pooledMSM <- FALSE
weight.msm <- FALSE

## ----chunk.simrunltmleMSM2, eval=TRUE, echo=FALSE---------------------------------------
simrun_ltmleMSM <- function(sim, DAG, N, t0,gbounds, gforms) {
  library("ltmle")
  O_datnew <- sim(DAG = DAG, n = N)
  ltmleMSMparams <- est.targetMSM(DAG, O_datnew, Aname = "A1", Cname = "A2", Lnames = "L2", 
                                  Ytvec = (1:t.end), ACLtvec = (0:t.end), package = "ltmle")
  summeas_arr <- ltmleMSMparams$summeas_arr
  regimens_arr <- ltmleMSMparams$regimens_arr[,c(1:t.end),]
  O_datnewLTCF <- doLTCF(data = O_datnew, LTCF = "Y")
  O_dat_selCnew <- O_datnewLTCF[,-which(names(O_datnewLTCF)%in%dropnms)]
  O_dat_selCnew[,Cnodesnew] <- 1-O_dat_selCnew[,Cnodesnew]
  reslTMLE.MSM <- ltmleMSM(data = O_dat_selCnew, Anodes = Anodesnew, Cnodes = Cnodesnew, 
                              Lnodes = Lnodesnew, Ynodes = Ynodesnew, 
                              survivalOutcome = survivalOutcome,
                              gform = gforms, Qform = Qforms,
                              stratify = stratify_Qg, mhte.iptw = mhte.iptw, 
                              iptw.only = FALSE, 
                              working.msm = msm.form, pooledMSM = pooledMSM, 
                              final.Ynodes = finYnodesnew, regimes = regimens_arr, 
                              summary.measures = summeas_arr, weight.msm=weight.msm,
                              estimate.time = FALSE, gbounds = gbounds)
  iptwMSMcoef <- summary(reslTMLE.MSM, estimator = "iptw")$cmat[,1]
  iptwRD <- MSM_RD_t(resMSM = iptwMSMcoef, t = t0)
  tmleMSMcoef <- summary(reslTMLE.MSM, estimator = "tmle")$cmat[,1]
  tmleRD <- MSM_RD_t(resMSM = tmleMSMcoef, t = t0)
  return(c(simN = sim, iptwRD = iptwRD, tmleRD = tmleRD))
}

## ----chunk.simrunltmleMSM3, message=FALSE, eval=FALSE, echo=FALSE-----------------------
#  t0 <- 12
#  Nltmle <- 50000
#  Nsims <- 1000
#  source("./determineParallelBackend.R")
#  sim50K.stratQg.notrunc.g <- foreach(sim = seq(Nsims), .combine = 'rbind')%dopar%{
#                        simrun_ltmleMSM(sim = sim,DAG = Ddyn, N = Nltmle, t0 = t0,
#                          gbounds = c(0.0000001, 1), gforms = gforms)
#              }
#  save(list = "sim50K.stratQg.notrunc.g", file = "vignette_dat/sim50K.stratQg.notrunc.g.Rdata")
#  sim50K.stratQg.notrunc.missg <- foreach(sim = seq(Nsims), .combine = 'rbind')%dopar%{
#                        simrun_ltmleMSM(sim = sim,DAG = Ddyn, N = Nltmle, t0 = t0,
#                          gbounds = c(0.0000001, 1), gforms = gforms_miss)
#              }
#  save(list = "sim50K.stratQg.notrunc.missg", file = "vignette_dat/sim50K.stratQg.notrunc.missg.Rdata")

## ----message=FALSE, warning=FALSE, echo=FALSE-------------------------------------------
library("Hmisc")
load(file = "vignette_dat/repstudy1_psi0.t0.MSM.Rdata")
load(file = "vignette_dat/sim50K.stratQg.notrunc.g.Rdata")
iptw_sims <- sim50K.stratQg.notrunc.g[,"iptwRD"]
tmle_sims <- sim50K.stratQg.notrunc.g[,"tmleRD"]
getrestab <- function(iptw_sims, tmle_sims, origres, modelnm) {
  resrow <- function(sims, estname) {
    data.frame(estname, 
                sprintf("%.3f",mean(sims)), # round(mean(sims),3),
                sprintf("%.4f",mean(sims)-psi0.t0.MSM), # round(((mean(sims)-psi0.t0.MSM)/psi0.t0.MSM)*100,3),
                sprintf("%.4f",sd(sims)), # round(sd(sims),4)
                stringsAsFactors=FALSE
                )
  }
  tabcolnames <- c("Estimator", "$\\psi_{n}$", "Bias", "$\\sigma_{emp}$", "$\\sigma_{emp}^{IPTW}/\\sigma_{emp}^{TMLE}$")    
  restab <- rbind(resrow(tmle_sims, "TMLE "%+%modelnm), resrow(iptw_sims, "IPTW "%+%modelnm))
  restab <- cbind(restab, c(round(sd(iptw_sims)/sd(tmle_sims),2), ""))
  restab <- rbind(c("\\emph{Results from this replication study}", "", "", "", ""), restab)

  restab <- rbind(restab, c("\\emph{Results reported by Neugebauer et al. (2014)}", "", "", "", ""))
  colnames(restab) <- tabcolnames
  colnames(origres) <- tabcolnames
  restab <- rbind(restab, origres)
  restab
}
# Results reported by Neugebauer for correct model:
tmleorigres <- c("TMLE (correct model)", "0.108", "0.0010", "0.0060", "1.392")
iptworigres <- c("IPTW (correct model)", "0.108", "0.0011", "0.0078", "")
# TMLE (correct model), 0.108, 0.001, 0.006, 1.392
# IPTW (correct model), 0.108, 0.0011, 0.0078
restab_corrg <- getrestab(iptw_sims, tmle_sims, origres=rbind(tmleorigres,iptworigres), modelnm="(correct model)")

## ----message=FALSE, echo=FALSE----------------------------------------------------------
load(file = "vignette_dat/sim50K.stratQg.notrunc.missg.Rdata")
iptw_sims <- sim50K.stratQg.notrunc.missg[,"iptwRD"]
tmle_sims <- sim50K.stratQg.notrunc.missg[,"tmleRD"]
# Results reported by Neugebauer for incorrect model 1:
tmleorigres <- c("TMLE (incorrect model 1)", "0.109", "0.0019", "0.0074", "")
iptworigres <- c("IPTW (incorrect model 1)", "0.182", "0.0750", "0.0107", "")
# TMLE (incorrect model 1), 0.109, 0.0019, 0.0074,
# TMLE (incorrect model 1), 0.182, 0.075, 0.0107
restab_missg <- getrestab(iptw_sims, tmle_sims, origres=rbind(tmleorigres,iptworigres), modelnm="(incorrect model 1)")
restab <- rbind(restab_corrg, restab_missg)

## ----message=FALSE, echo=FALSE, results="asis"------------------------------------------
cat("\n")
latex(restab,file = "",where = "!htpb", caption.loc = 'bottom', caption = "Replication of the results from simulation protocol 3 reported in Table 6 of \\citet{neugebauer2014} based on two models for estimating the treatment mechanism: 1) a correctly specified model and 2) a misspecified model missing a term for the time-dependent variable. $\\psi_{n}$ - mean point estimates over 1,000 simulated data sets; $\\sigma_{emp}$ - empirical standard deviation (SD) of point estimates over 1,000 simulated data sets; $\\sigma_{emp}^{IPTW}/\\sigma_{emp}^{TMLE}$ - the relative efficiency measured by the ratio of the empirical SDs associated with the IPW and TMLE point estimates.", label = 'tab6aNeugebauer', booktabs = TRUE, rowname = NULL, landscape = FALSE, col.just = c("l", rep("r", 4)), size = "small")

## ----message=FALSE, results='hide'------------------------------------------------------
rbivNorm <- function(n, whichbiv, norms, mu, var1 = 1, var2 = 1, rho = 0.7) {
  whichbiv <- whichbiv[1]; var1 <- var1[1]; var2 <- var2[1]; rho <- rho[1]
  sigma <- matrix(c(var1, rho, rho, var2), nrow = 2)
  Scol <- chol(sigma)[, whichbiv]
  bivX <- (Scol[1] * norms[, 1] + Scol[2] * norms[, 2]) + mu
  bivX
}

## ----message=FALSE, results='hide'------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
Lnames <- c("LO1", "LO2", "LO3", "LC1")
D <- DAG.empty()

for (Lname in Lnames) {
  D <- D + 
    node(Lname%+%".norm1", distr = "rnorm", mean = 0, sd = 1) + 
    node(Lname%+%".norm2", distr = "rnorm", mean = 0, sd = 1)
}

D <- D +
  node("LO1", t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
        norms = c(LO1.norm1, LO1.norm2), 
        mu = 0) + 
  node("LO2", t = 0:1, distr = "rbivNorm", whichbiv = t + 1,
        norms = c(LO2.norm1, LO2.norm2), 
        mu = 0) + 
  node("LO3", t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
        norms = c(LO3.norm1, LO3.norm2), 
        mu = 0) + 
  node("LC1", t = 0:1, distr = "rbivNorm", whichbiv = t + 1, 
        norms = c(LC1.norm1, LC1.norm2), 
        mu = {if (t == 0) {0} else {-0.30 * A[t-1]}}) + 
  node("alpha", t = 0:1, distr = "rconst", 
        const = {if(t == 0) {log(0.6)} else {log(1.0)}}) + 
  node("A", t = 0:1, distr = "rbern", 
        prob = plogis(alpha[t] + 
                      log(5)*LC1[t] + {if(t == 0) {0} else {log(5)*A[t-1]}})) + 
  node("Y", t = 1, distr = "rnorm", 
        mean = (0.98 * LO1[t] + 0.58 * LO2[t] + 0.33 * LO3[t] + 
                0.98 * LC1[t] - 0.37 * A[t]),
        sd = 1)

DAGO.sc1 <- set.DAG(D)

## ----message=FALSE, results='hide'------------------------------------------------------
defAct <- function (Dact) {
  act.At <- node("A", t = 0:1, distr = "rbern", prob = abar[t])
  Dact <- Dact + 
    action("A00", nodes = act.At, abar = c(0, 0)) + 
    action("A10", nodes = act.At, abar = c(1, 0)) + 
    action("A01", nodes = act.At, abar = c(0, 1)) + 
    action("A11", nodes = act.At, abar = c(1, 1))
  return(Dact)
}

Dact.sc1 <- defAct(DAGO.sc1)
msm.form <- "Y ~ S(abar[0]) + S(abar[1])"
Dact.sc1 <- set.targetMSM(Dact.sc1, outcome = "Y", t = 1, 
                          form = msm.form, family = "gaussian")

## ----message=FALSE, eval=TRUE-----------------------------------------------------------
repstudy2.sc1.truetarget <- function() {
  trueMSMreps.sc1 <- NULL
  reptrue <- 50
  for (i in (1:reptrue)) {
    res.sc1.i <- eval.target(Dact.sc1, n = 500000)$coef
    trueMSMreps.sc1 <- rbind(trueMSMreps.sc1, res.sc1.i)
  }
  return(trueMSMreps.sc1)
}

f1name <- "vignette_dat/trueMSMreps.sc1.Rdata"
if (file.exists(f1name)) {
  load(f1name)
} else {
  trueMSMreps.sc1 <- repstudy2.sc3.truetarget()
  save(list = "trueMSMreps.sc1", file = f1name)
}
(trueMSM.sc1 <- apply(trueMSMreps.sc1, 2, mean))

## ----chunk.runMSMsw, eval=TRUE, echo=FALSE----------------------------------------------
runMSMsw <- function(DAGO, Lnames, trueA, nsamp, nsims) {
  Lnames_0 <- Lnames%+%"_0"
  Lnames_1 <- Lnames%+%"_1"
  gforms <- c("A_0 ~ "%+%paste(Lnames_0, collapse = " + "), "A_1 ~ A_0 + "%+%paste(Lnames_1, collapse = " + "))
  res_sw <- NULL
  for (sims in (1:nsims)) {
    datO <- sim(DAGO, n = nsamp)
    glmA_0 <- glm(datO[,c("A_0",Lnames_0)], formula = gforms[1], family = "binomial")
    glmA_1 <- glm(datO[,c("A_1","A_0",Lnames_0,Lnames_1)], formula = gforms[2], family = "binomial")
    probA0_1 <- predict(glmA_0,  type = "response")
    weight_t0 <- 1 / (probA0_1^(datO$A_0) * (1-probA0_1)^(1-datO$A_0))
    probA1_1 <- predict(glmA_1,  type = "response")
    weight_t1 <- 1 / (probA1_1^(datO$A_1) * (1-probA1_1)^(1-datO$A_1))
    sw1 <- weight_t0*weight_t1
    emp.pA1cA0 <- table(datO$A_1,datO$A_0)/nrow(datO)
    empPA1 <- data.frame(A_0 = c(0,0,1,1),A_1 = c(0,1,1,0))
    empPA1$empPA_1_cA_0 <- apply(empPA1, 1, function(rowA) emp.pA1cA0[as.character(rowA["A_1"]), as.character(rowA["A_0"])])
    empPA1 <- merge(datO[, c("ID","A_0","A_1")],empPA1, sort = FALSE)
    empPA1 <- empPA1[order(empPA1$ID),]
    swts <- empPA1$empPA_1_cA_0*(weight_t0*weight_t1)
    datO$swts <- swts
    MSMres_sw <- glm(datO, formula = "Y_1 ~ A_0 + A_1", weights = swts, family = "gaussian")
    res_sw <- rbind(res_sw, coef(MSMres_sw))
  }
  meanres <- apply(res_sw, 2, mean)
  Varres <- apply(res_sw, 2, var)
  bias <- c(meanres["A_0"]-trueA["A_0"], meanres["A_1"]-trueA["A_1"])
  MSE <- c(bias^2+Varres[c("A_0","A_1")])
  bias10 <- sprintf("%.3f",bias*10)
  MSE10 <- sprintf("%.3f",MSE*10)
  resrow <- c(bias10[1], MSE10[1], bias10[2], MSE10[2])
  col36names <- c("\\specialcell[t]{A(0)\\\\ Bias*10}", 
                  "\\specialcell[t]{A(0)\\\\ MSE*10}", 
                  "\\specialcell[t]{A(1)\\\\ Bias*10}",
                  "\\specialcell[t]{A(1)\\\\ MSE*10}")
  names(resrow) <- col36names
  return(resrow)
}

## ----chunk.lefebrepT2andT4, message=FALSE, eval=TRUE, results='hide', echo=FALSE--------
# recreating Tables 2 & 4 reported in Lefebvre et al.
nsamp <- c(300,1000,10000)
# Lefebvre et al. Tab 2:
covnmT2 <- c(c("\\emph{Lefebvre et al.}: Confounder(s) only", rep("",2)), 
            c("\\emph{Lefebvre et al.}: Confounder(s) &", "risk factors", rep("",1)))
lefebvreT2 <- data.frame(
  covnm = covnmT2,
  N = rep(nsamp,2),
  A0Bias10 = sprintf("%.3f",c(0.768, 0.265, 0.057, 0.757, 0.283, 0.056)),
  A0MSE10 = sprintf("%.3f",c(1.761, 0.761, 0.146, 1.642, 0.718, 0.139)),
  A1Bias10 = sprintf("%.3f",c(0.889, 0.312, 0.086, 0.836, 0.330, 0.081)),
  A1MSE10 = sprintf("%.3f",c(1.728, 0.723, 0.120, 1.505, 0.638, 0.114)),stringsAsFactors = FALSE)
# Lefebvre et al. Tab 4:
covnmT4 <- c(c("\\emph{Lefebvre et al.}: Confounder(s) only", rep("",2)),
              c("\\emph{Lefebvre et al.}: Confounder(s) &", "risk factors", ""),
              c("\\emph{Lefebvre et al.}: Confounder(s) &", "IVs", ""),
              c("\\emph{Lefebvre et al.}: Confounder(s),", "IVs & risk factors",""),
              c("\\emph{Lefebvre et al.}: Mis-specified", rep("",2)),
              c("\\emph{Lefebvre et al.}: Full Model", rep("",2)))
lefebvreT4 <- data.frame(
  covnm = covnmT4,
  N = rep(nsamp,6),
  A0Bias10 = sprintf("%.3f",c(-0.080, -0.371, -0.368, -0.110, -0.330, -0.378, 1.611, 
                              0.824, 0.241, 1.600, 0.867, 0.235, 3.146, 2.460, 2.364, 
                              1.524, 0.878, 0.240)),
  A0MSE10 = sprintf("%.3f",c(1.170, 0.385, 0.056, 1.092, 0.340, 0.051, 3.538, 2.063, 
                              0.684, 3.477, 2.053, 0.676, 3.326, 1.700, 0.832, 3.648, 
                              2.099, 0.679)),
  A1Bias10 = sprintf("%.3f",c(0.099, -0.035, -0.203, 0.112, -0.108, -0.207, 2.069, 1.245, 
                              0.379, 2.143, 1.170, 0.372, 5.591, 5.258, 4.943, 2.221, 1.185,
                              0.377)),
  A1MSE10 = sprintf("%.3f",c(1.155, 0.331, 0.043, 0.865, 0.245, 0.037, 3.841, 2.188, 0.622, 
                              3.598, 2.043, 0.625, 5.494, 3.851, 2.705, 3.907, 2.099, 0.630)),
                  stringsAsFactors = FALSE)
col1name <- "Covariates in $P(A|L)$"
colnames(lefebvreT2)[1] <- colnames(lefebvreT4)[1] <- col1name
col36names <- c("\\specialcell[t]{A(0)\\\\ Bias*10}", 
                "\\specialcell[t]{A(0)\\\\ MSE*10}", 
                "\\specialcell[t]{A(1)\\\\ Bias*10}",
                "\\specialcell[t]{A(1)\\\\ MSE*10}")
colnames(lefebvreT2)[3:6] <- colnames(lefebvreT4)[3:6] <- col36names

## ----chunk.lefebSc1a, message=FALSE, eval=FALSE, results='hide', echo=FALSE-------------
#  trueA <- c(A_0 = -0.294, A_1 = -0.370)
#  nsims <- 10000; restab <- NULL
#  runsim <- function(Lnames, DAGO) {
#    for (nsamp in c(300,1000,10000)) {
#      resSc <- runMSMsw(DAGO = DAGO, Lnames = Lnames, trueA = trueA, nsamp = nsamp, nsims = nsims)
#      restab <- rbind(restab, c(N = nsamp, resSc))
#    }
#    restab
#  }
#  Lnames <- c("LC1")
#  covnm <- c("Confounder(s) only", rep("",2))
#  restab_1 <- cbind(covnm, runsim(Lnames, DAGO.sc1))
#  # restab_1 <- rbind(restab_1, as.matrix(lefebvreT2[1:3,]))
#  Lnames <- c("LC1", "LO1", "LO2", "LO3")
#  covnm <- c("Confounder(s) &", "risk factors", rep("",1))
#  restab_2 <- cbind(covnm, runsim(Lnames, DAGO.sc1))
#  # restab_2 <- rbind(restab_2, as.matrix(lefebvreT2[4:6,]))
#  restab <- rbind(restab_1, restab_2)
#  col1name <- "Covariates in $P(A|L)$"
#  colnames(restab)[1] <- col1name
#  # restabwLef <- restab
#  # save(list = "restabwLef", file = "vignette_dat/restabwLefSc1_all_1Ksims.Rdata");
#  # restab <- restab[c(1:3, 7:9),]
#  save(list = "restab", file = "vignette_dat/restabSc1_all_1Ksims.Rdata");

## ----chunk.lefebSc1b, message=FALSE, echo=FALSE, results="asis"-------------------------
library("Hmisc")
load(file = "vignette_dat/restabSc1_all_1Ksims.Rdata");
cat("\n")
latex(restab, file = "", where = "!htpb", caption.loc = 'bottom', 
    caption = "Replication of the simulation results from \\citet{lefebvre2008} for Scenario 1.", 
    label = 'tab2Lefebvre',booktabs = TRUE,rowname = NULL,landscape = FALSE, 
    col.just = c("l", rep("r", 5)), size = "small")

## ----chunk.lefebSc1c, message=FALSE, echo=FALSE, results="asis"-------------------------
cat("\n")
latex(lefebvreT2, file = "", where = "!htpb", caption.loc = 'bottom', 
    caption = "Simulation results for Scenario 1 as reported in Table II of \\citet{lefebvre2008}.", 
    label = 'origtab2Lefebvre', booktabs = TRUE, rowname = NULL, landscape = FALSE, 
    col.just = c("l", rep("r", 5)), size = "small")

## ----message=FALSE, warning=FALSE, results='hide'---------------------------------------
`%+%` <- function(a, b) paste0(a, b)
Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
D <- DAG.empty()

for (Lname in Lnames) {
  D <- D +
          node(Lname%+%".norm1", distr = "rnorm") +
          node(Lname%+%".norm2", distr = "rnorm")
}

## ----message=FALSE, warning=FALSE, results='hide'---------------------------------------
coefAi <- c(-0.10, -0.20, -0.30)
sdLNi <- c(sqrt(1), sqrt(5), sqrt(10))

for (i in (1:3)) {
  D <- D +
    node("LO"%+%i, t = 0:1, distr = "rbivNorm", whichbiv = t + 1,
          mu = 0,
          params = list(norms = "c(LO"%+%i%+%".norm1, LO"%+%i%+%".norm2)")) +
    node("LE"%+%i, t = 0:1, distr = "rbivNorm", whichbiv = t + 1,
          mu = 0, var1 = 1, var2 = 1, rho = 0.7,
          params = list(norms = "c(LE"%+%i%+%".norm1, LE"%+%i%+%".norm2)")) +
    node("LC"%+%i, t = 0:1, distr = "rbivNorm", whichbiv = t + 1,
          mu = {if (t == 0) {0} else {.(coefAi[i]) * A[t-1]}},
          params = list(norms = "c(LC"%+%i%+%".norm1, LC"%+%i%+%".norm2)")) +
    node("LN"%+%i, t = 0:1, distr = "rnorm",
          mean = 0, sd = .(sdLNi[i]))
}

D <- D + 
  node("alpha", t = 0:1, distr = "rconst",
        const = {if(t == 0) {log(0.6)} else {log(1.0)}}) +
  node("A", t = 0:1, distr = "rbern", 
        prob = plogis(alpha[t] + 
                      log(5) * LC1[t] + log(2) * LC2[t] + log(1.5) * LC3[t] +
                      log(5) * LE1[t] + log(2) * LE2[t] + log(1.5) * LE3[t] +
                      {if (t == 0) {0} else {log(5) * A[t-1]}})) +
  node("Y", t = 1, distr = "rnorm",
        mean = 0.98 * LO1[t] + 0.58 * LO2[t] + 0.33 * LO3[t] +
                0.98 * LC1[t] + 0.58 * LC2[t] + 0.33 * LC3[t] - 0.39 * A[t],
        sd = 1)

DAGO.sc3 <- set.DAG(D)

## ----message=FALSE, eval=TRUE-----------------------------------------------------------
Dact.sc3 <- defAct(DAGO.sc3)
msm.form <- "Y ~ S(abar[0]) + S(abar[1])"
Dact.sc3 <- set.targetMSM(Dact.sc3, outcome = "Y", t = 1, 
                          form = msm.form, family = "gaussian")

repstudy2.sc3.truetarget <- function() {
  trueMSMreps.sc3 <- NULL
  reptrue <- 50
  for (i in (1:reptrue)) {
    res.sc3.i <- eval.target(Dact.sc3, n = 500000)$coef
    trueMSMreps.sc3 <- rbind(trueMSMreps.sc3, res.sc3.i)
  }
  return(trueMSMreps.sc3)
}

f2name <- "vignette_dat/trueMSMreps.sc3.Rdata"
if (file.exists(f2name)) {
  load(f2name)
} else {
  trueMSMreps.sc3 <- repstudy2.sc3.truetarget()
  save(list = "trueMSMreps.sc3", file = f2name)
}
(trueMSM.sc3 <- apply(trueMSMreps.sc3, 2, mean))

## ----chunk.lefebSc3a, message=FALSE, eval=FALSE, echo=FALSE-----------------------------
#  trueA <- c(A_0 = -0.316, A_1 = -0.390)
#  nsims <- 10000; restab <- NULL
#  runsim <- function(Lnames, DAGO) {
#    for (nsamp in c(300,1000,10000)) {
#      resSc <- runMSMsw(DAGO = DAGO, Lnames = Lnames, trueA = trueA, nsamp = nsamp, nsims = nsims)
#      restab <- rbind(restab, c(N = nsamp, resSc))
#    }
#    restab
#  }
#  Lnames <- c("LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s) only", rep("",2))
#  restab_1 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_1 <- rbind(restab_1, as.matrix(lefebvreT4[1:3,]))
#  Lnames <- c("LO1", "LO2", "LO3", "LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s) &", "risk factors", "")
#  restab_2 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_2 <- rbind(restab_2, as.matrix(lefebvreT4[4:6,]))
#  Lnames <- c("LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s) &", "IVs", "")
#  restab_3 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_3 <- rbind(restab_3, as.matrix(lefebvreT4[7:9,]))
#  Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s),", "IVs & risk factors","")
#  restab_4 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_4 <- rbind(restab_4, as.matrix(lefebvreT4[10:12,]))
#  Lnames <- c("LE1", "LE2", "LE3", "LC1")
#  covnm <- c("Mis-specified", rep("",2))
#  restab_5 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_5 <- rbind(restab_5, as.matrix(lefebvreT4[13:15,]))
#  Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3", "LN1", "LN2", "LN3")
#  covnm <- c("Full Model", rep("",2))
#  restab_6 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_6 <- rbind(restab_6, as.matrix(lefebvreT4[16:18,]))
#  restab <- rbind(restab_1, restab_2, restab_3, restab_4, restab_5, restab_6)
#  col1name <- "Covariates in $P(A|L)$"
#  colnames(restab)[1] <- col1name
#  # restabwLef <- restab
#  # save(list = "restabwLef", file = "vignette_dat/restabwLefSc3_all_1Ksims.Rdata");
#  # restab <- restab[c(1:3, 7:9, 13:15, 19:21, 25:27, 31:33),]
#  save(list = "restab", file = "vignette_dat/restabSc3_all_1Ksims.Rdata");

## ----chunk.lefebSc3b, message=FALSE, echo=FALSE, results="asis"-------------------------
library("Hmisc")
load(file = "vignette_dat/restabSc3_all_1Ksims.Rdata");
cat("\n")
latex(restab,file = "",where = "!htpb", caption.loc = 'bottom', 
    caption = "Replication of the simulation results from \\citet{lefebvre2008} for Scenario 3.", 
    label = 'tab4Lefebvre',booktabs = TRUE,rowname = NULL,landscape = FALSE, 
    col.just = c("l", rep("r", 5)), size = "small")

## ----chunk.lefebSc3c, message=FALSE, echo=FALSE, results="asis"-------------------------
cat("\n")
latex(lefebvreT4,file = "",where = "!htpb", caption.loc = 'bottom', 
    caption = "Simulation results for Scenario 3 as reported in Table IV of \\citet{lefebvre2008}.", 
    label = 'origtab4Lefebvre',booktabs = TRUE,rowname = NULL,landscape = FALSE, 
    col.just = c("l", rep("r", 5)), size = "small")

## ----appendixcode1, ref.label = 'chunkMSMsurvplot', eval=FALSE, echo=TRUE, size='tiny'----
#      # get MSM survival predictions from the full data.table in long format (melted) by time (t_vec), and by the MSM term (MSMtermName)
#      # predictions from the estimated msm model (based on observational data) can be obtained by passing estimated msm model, est.msm
#      # given vector of t (t_vec), results of MSM target.eval and an MSM term get survival table by action
#     survbyMSMterm <- function(MSMres, t_vec, MSMtermName, use_actions=NULL, est.msm=NULL) {
#        library("data.table")
#        # look up MSMtermName in MSMterm map, if exists -> use the new name, if doesn't exist use MSMtermName
#        if (!is.null(MSMres$S.msm.map)) {
#          mapS_exprs <- as.character(MSMres$S.msm.map[,"S_exprs_vec"])
#          XMSMterms <- as.character(MSMres$S.msm.map[,"XMSMterms"])
#          map_idx <- which(mapS_exprs%in%MSMtermName)
#          XMSMtermName <- XMSMterms[map_idx]
#          if (!is.null(XMSMtermName)&&length(XMSMtermName)>0) {
#            MSMtermName <- XMSMtermName
#          }
#        }
#        print("MSMtermName used"); print(MSMtermName)
#        t_dt <- data.table(t=as.integer(t_vec)); setkey(t_dt, t)
#        get_predict <- function(actname) {
#          # setkey(MSMres$df_long[[actname]], t)
#          setkeyv(MSMres$df_long[[actname]], c("t", MSMtermName))
#          # MSMterm_vals <- as.numeric(MSMres$df_long[[actname]][t_dt, mult="first"][[MSMtermName]])
#          # print("MSMterm_vals"); print(MSMterm_vals)
#          MSMterm_vals <- as.numeric(MSMres$df_long[[actname]][t_dt, mult="last"][[MSMtermName]])
#          # print("MSMterm_vals last"); print(MSMterm_vals)
#          newdata=data.frame(t=t_vec, MSMterm_vals=MSMterm_vals)
#          colnames(newdata) <- c("t", MSMtermName)
#          # print("newdata"); print(newdata)
#          if (!is.null(est.msm)) {
#            pred <- predict(est.msm, newdata=newdata, type="response")
#          } else {
#            pred <- predict(MSMres$m, newdata=newdata, type="response")
#          }
#          return(data.frame(t=t_vec, pred=pred))
#        }
#        action_names <- names(MSMres$df_long)
#        if (!is.null(use_actions)) {
#          action_names <- action_names[action_names%in%use_actions]
#        }
#        surv <- lapply(action_names, function(actname) {
#                                      res <- get_predict(actname)
#                                      if (MSMres$hazard) {
#                                        res$surv <- cumprod(1-res$pred)
#                                      } else {
#                                        res$surv <- 1-res$pred
#                                      }
#                                      res$pred <- NULL
#                                      res$action <- actname
#                                      res
#                                    })
#        names(surv) <- names(MSMres$df_long)
#        surv_melt <- do.call('rbind', surv)
#        surv_melt$action <- factor(surv_melt$action, levels=unique(surv_melt$action), ordered=TRUE)
#        surv_melt
#      }
#      plotsurvbyMSMterm <- function(surv_melt_dat) {
#        library("ggplot2")
#        f_ggplot_surv_wS <- ggplot(data= surv_melt_dat, aes(x=t, y=surv)) +
#                            geom_line(aes(group = action, color = action), size=.4, linetype="dashed") +
#                            theme_bw()
#  
#      }
#      plotsurvbyMSMterm_facet <- function(surv_melt_dat1, surv_melt_dat2, msm_names=NULL) {
#        library("ggplot2")
#        if (is.null(msm_names)) {
#          msm_names <- c("MSM1", "MSM2")
#        }
#        surv_melt_dat1$MSM <- msm_names[1]
#        surv_melt_dat2$MSM <- msm_names[2]
#  
#        surv_melt_dat <- rbind(surv_melt_dat1,surv_melt_dat2)
#        f_ggplot_surv_wS <- ggplot(data= surv_melt_dat, aes(x=t, y=surv)) +
#                            geom_line(aes(group = action, color = action), size=.4, linetype="dashed") +
#                            theme_bw() +
#                            facet_wrap( ~ MSM)
#      }

## ----appendixcode2, ref.label='chunk.est.targetMSM', eval=FALSE, echo=TRUE, size='tiny'----
#    # @param DAG Object specifying the directed acyclic graph for the observed data,
#      # must have a well-defined MSM target parameter (\code{set.target.MSM()})
#    # @param obs_df Simulated observational data
#    # @param Aname Generic names of the treatment nodes (can be time-varying)
#    # @param Cname Generic names of the censoring nodes (can be time-varying)
#    # @param Lnames Generic names of the time-varying covariates (can be time-varying)
#    # @param tvec Vector of time points for Y nodes
#    # @param actions Which actions (regimens) should be used in estimation from the observed simulated data.
#      # If NULL then all actions that were defined in DAG will be considered.
#    # @param package Character vector for R package name to use for estimation. Currently only "ltmle" is implemented.
#    # @param fun Character name for R function name to employ for estimation. Currently only "ltmleMSM" is implemented.
#    # @param ... Additional named arguments that will be passed on to ltmleMSM function
#    est.targetMSM <- function(DAG, obs_df, Aname="A", Cname="C", Lnames, Ytvec, ACLtvec, actions=NULL, package="ltmle", fun="ltmleMSM", ...) {
#      outnodes <- attr(DAG, "target")$outnodes
#      param_name <- attr(DAG, "target")$param_name
#      if (is.null(outnodes$t)) stop("estimation is only implemented for longitudinal data with t defined")
#      if (!param_name%in%"MSM") stop("estimation is only implemented for MSM target parameters")
#      if (is.null(actions)) {
#        message("actions argument underfined, using all available actions")
#        actions <- A(DAG)
#      }
#      # all time points actually used in the observed data
#      t_all <- attr(obs_df, "tvals")
#      tvec <- outnodes$t
#      t_sel <- ACLtvec
#      # ltmle allows for pooling Y's over smaller subset of t's (for example t=(2:8))
#      # in this case summary measures HAVE TO MATCH the dimension of finYnodes, not t_sel
#      # currently this is not supported, thus, if tvec is a subset of t_sel this will cause an error
#      finYnodes <- outnodes$gen_name%+%"_"%+%Ytvec
#      Ynodes <- outnodes$gen_name%+%"_"%+%Ytvec
#      Anodes <- Aname%+%"_"%+%ACLtvec
#      Cnodes <- Cname%+%"_"%+%ACLtvec
#      Lnodes <- t(sapply(Lnames, function(Lname) Lname%+%"_"%+%ACLtvec[-1]))
#      Lnodes <- as.vector(matrix(Lnodes, nrow=1, ncol=ncol(Lnodes)*length(Lnames), byrow=FALSE))
#      Nobs <- nrow(obs_df)
#      #------------------------------------------------------------
#      # getting MSM params
#      #------------------------------------------------------------
#      params.MSM <- attr(DAG, "target")$params.MSM
#      working.msm <- params.MSM$form
#      msm.family <- params.MSM$family
#      if (params.MSM$hazard) stop("ltmleMSM cannot estimate hazard MSMs...")
#      # the number of attributes and their dimensionality have to match between different actions
#      n_attrs <- length(attr(actions[[1]], "attnames"))
#      #------------------------------------------------------------
#      # define the final ltmle arrays
#      regimens_arr <- array(dim = c(Nobs, length(ACLtvec), length(actions)))
#      summeas_arr <- array(dim = c(length(actions), (n_attrs+1), length(Ytvec)))
#      # loop over actions (regimes) creating counterfactual mtx of A's for each action:
#      for (action_idx in seq(actions)) {
#          # I) CREATE COUNTERFACTUAL TREATMENTS &
#          # II) CREATE summary.measure that describes each attribute by time + regimen
#          #------------------------------------------------------------
#          # needs to assign observed treatments and replace the action timepoints with counterfactuals
#          A_mtx_act <- as.matrix(obs_df[,Anodes])
#          #------------------------------------------------------------
#          action <- actions[[action_idx]]
#           # action-spec. time-points
#          t_act <- as.integer(attr(action, "acttimes"))
#          # action-spec. attribute names
#          attnames <- attr(action, "attnames")
#          # list of action-spec. attributes
#          attrs <- attr(action, "attrs")
#          # time points for which we need to evaluate the counterfactual treatment assignment as determined by action:
#          #------------------------------------------------------------
#          # Action t's need always be the same subset of t_sel (outcome-based times), otherwise we are in big trouble
#          t_act_idx <- which(t_sel%in%t_act)
#          t_chg <- t_sel[t_act_idx]
#          # modify only A's which are defined in this action out of all Anodes
#          As_chg <- Anodes[t_act_idx]
#          #------------------------------------------------------------
#          # creates summary measure array that is of dimension (length(t_chg)) - time-points only defined for this action
#          # which t's are in the final pooled MSM => need to save the summary measures only for these ts
#          t_vec_idx <- which(t_act%in%tvec)
#          summeas_attr <- matrix(nrow=length(attnames), ncol=length(tvec))
#          #------------------------------------------------------------
#          # extract values of terms in MSM formula: get all attribute values from +action(...,attrs)
#          #------------------------------------------------------------
#          obs_df_attr <- obs_df # add all action attributes to the observed data
#          for (attr_idx in seq(attnames)) { # self-contained loop # grab values of the attributes, # loop over all attributes
#            if (length(attrs[[attnames[attr_idx]]])>1) {
#              attr_i <- attnames[attr_idx]%+%"_"%+%t_chg
#              val_attr_i <- attrs[[attnames[attr_idx]]][t_act_idx]
#            } else {
#              attr_i <- attnames[attr_idx]
#              val_attr_i <- attrs[[attnames[attr_idx]]]
#            }
#            # summary measures, for each action/measure
#            summeas_attr[attr_idx,] <- matrix(val_attr_i, nrow=1, ncol=length(t_chg))[,t_vec_idx]
#            # observed data values of the attribute
#            df_attr_i <- matrix(val_attr_i, nrow=Nobs, ncol=length(val_attr_i), byrow=TRUE)
#            # create the combined data.frame (attrs + O.dat)
#            colnames(df_attr_i) <- attr_i; obs_df_attr <- cbind(data.frame(df_attr_i), obs_df_attr)
#          } # end of loop
#          summeas_attr <- rbind(summeas_attr, t_chg[t_vec_idx])
#          rownames(summeas_attr) <- c(attnames, "t")
#          #------------------------------------------------------------
#          # add action specific summary measures to the full array
#          summeas_arr[action_idx, , ] <- summeas_attr
#          dimnames(summeas_arr)[[2]] <- rownames(summeas_attr)
#          #------------------------------------------------------------
#          # GENERATING A MATRIX OF COUNTERFACTUAL TREATMENTS:
#          for (Achange in As_chg) { # for each A defined in the action, evaluate its value applied to the observed data
#            cur.node <- action[As_chg][[Achange]]
#            t <- cur.node$t
#            newAval <- with(obs_df_attr, {  # for static no need to sample from a distr
#              ANCHOR_VARS_OBSDF <- TRUE
#              simcausal:::eval_nodeform(as.character(cur.node$dist_params$prob), cur.node)$evaled_expr
#            })
#  
#            if (length(newAval)==1) {
#              newA <- rep(newAval, Nobs)
#            } else {
#              newA <- newAval
#            }
#            A_mtx_act[,which(Anodes%in%Achange)] <- newA
#          }
#          # Result matrix A_mtx_act has all treatments that were defined in that action replaced with
#          # their counterfactual values
#          #------------------------------------------------------------
#          # add action specific summary measures to the full array
#          regimens_arr[, , action_idx] <- A_mtx_act
#          #------------------------------------------------------------
#      }
#      list(regimens_arr=regimens_arr, summeas_arr=summeas_arr)
#    }

## ----appendixcodesimMSM1, ref.label='chunk.simrunltmleMSM1', eval=FALSE, echo=TRUE, size='tiny'----
#  times <- c(0:(t.end-1))
#  gforms <- c("A1_0 ~ L1_0 + L2_0","A2_0 ~ L1_0")
#  timesm0 <- times[which(times > 0)]
#  # correctly specified g:
#  gforms <- c("A1_0 ~ L1_0 + L2_0","A2_0 ~ L1_0")
#  gformm0 <- as.vector(sapply(timesm0, function(t)
#                c("A1_"%+%t%+%" ~ A1_"%+%(t-1)%+%" + L1_0"%+%" + L2_"%+%t%+%" + I(L2_"%+%t%+%"*A1_"%+%(t-1)%+%")+ I(L1_0*A1_"%+%(t-1)%+%")",
#                  "A2_"%+%t%+%" ~ L1_0")))
#  gforms <- c(gforms, gformm0)
#  # mis-specified g (no TV covar L2):
#  gforms_miss <- c("A1_0 ~ L1_0","A2_0 ~ L1_0")
#  gformm0_miss <- as.vector(sapply(timesm0, function(t)
#                                c("A1_"%+%t%+%" ~ L1_0*A1_"%+%(t-1),
#                                  "A2_"%+%t%+%" ~ L1_0")))
#  gforms_miss <- c(gforms_miss, gformm0_miss)
#  Qformallt <- "Q.kplus1 ~ L1_0"
#  Lterms <- function(var, tlast){
#    tstr <- c(0:tlast)
#    strout <- paste(var%+%"_"%+%tstr, collapse = " + ")
#    return(strout)
#  }
#  tY <- (0:11)
#  Ynames <- paste("Y_"%+%c(tY+1))
#  Qforms <- unlist(lapply(tY, function(t)  {
#                    a <- Qformallt%+%" + I("%+%Lterms("m1L2",t)%+%") + "%+%Lterms("L2",t)
#                    return(a)
#    }))
#  names(Qforms) <- Ynames
#  survivalOutcome <- TRUE
#  stratify_Qg <- TRUE
#  mhte.iptw <- TRUE
#  Anodesnew <- "A1_"%+%(0:(t.end-1))
#  Cnodesnew <- "A2_"%+%(0:(t.end-1))
#  L2nodesnew <- "L2_"%+%(1:(t.end-1))
#  mL2nodesnew <- "m1L2_"%+%(1:(t.end-1))
#  Lnodesnew <- as.vector(rbind(L2nodesnew, mL2nodesnew))
#  Ynodesnew <- "Y_"%+%(1:t.end)
#  finYnodesnew <- Ynodesnew
#  dropnms <- c("ID","L2_"%+%t.end,"m1L2_"%+%t.end, "A1_"%+%t.end, "A2_"%+%t.end)
#  pooledMSM <- FALSE
#  weight.msm <- FALSE

## ----appendixcodesimMSM2, ref.label='chunk.simrunltmleMSM2', eval=FALSE, echo=TRUE, size='tiny'----
#  simrun_ltmleMSM <- function(sim, DAG, N, t0,gbounds, gforms) {
#    library("ltmle")
#    O_datnew <- sim(DAG = DAG, n = N)
#    ltmleMSMparams <- est.targetMSM(DAG, O_datnew, Aname = "A1", Cname = "A2", Lnames = "L2",
#                                    Ytvec = (1:t.end), ACLtvec = (0:t.end), package = "ltmle")
#    summeas_arr <- ltmleMSMparams$summeas_arr
#    regimens_arr <- ltmleMSMparams$regimens_arr[,c(1:t.end),]
#    O_datnewLTCF <- doLTCF(data = O_datnew, LTCF = "Y")
#    O_dat_selCnew <- O_datnewLTCF[,-which(names(O_datnewLTCF)%in%dropnms)]
#    O_dat_selCnew[,Cnodesnew] <- 1-O_dat_selCnew[,Cnodesnew]
#    reslTMLE.MSM <- ltmleMSM(data = O_dat_selCnew, Anodes = Anodesnew, Cnodes = Cnodesnew,
#                                Lnodes = Lnodesnew, Ynodes = Ynodesnew,
#                                survivalOutcome = survivalOutcome,
#                                gform = gforms, Qform = Qforms,
#                                stratify = stratify_Qg, mhte.iptw = mhte.iptw,
#                                iptw.only = FALSE,
#                                working.msm = msm.form, pooledMSM = pooledMSM,
#                                final.Ynodes = finYnodesnew, regimes = regimens_arr,
#                                summary.measures = summeas_arr, weight.msm=weight.msm,
#                                estimate.time = FALSE, gbounds = gbounds)
#    iptwMSMcoef <- summary(reslTMLE.MSM, estimator = "iptw")$cmat[,1]
#    iptwRD <- MSM_RD_t(resMSM = iptwMSMcoef, t = t0)
#    tmleMSMcoef <- summary(reslTMLE.MSM, estimator = "tmle")$cmat[,1]
#    tmleRD <- MSM_RD_t(resMSM = tmleMSMcoef, t = t0)
#    return(c(simN = sim, iptwRD = iptwRD, tmleRD = tmleRD))
#  }

## ----appendixcodesimMSM3, ref.label='chunk.simrunltmleMSM3', eval=FALSE, echo=TRUE, size='tiny'----
#  t0 <- 12
#  Nltmle <- 50000
#  Nsims <- 1000
#  source("./determineParallelBackend.R")
#  sim50K.stratQg.notrunc.g <- foreach(sim = seq(Nsims), .combine = 'rbind')%dopar%{
#                        simrun_ltmleMSM(sim = sim,DAG = Ddyn, N = Nltmle, t0 = t0,
#                          gbounds = c(0.0000001, 1), gforms = gforms)
#              }
#  save(list = "sim50K.stratQg.notrunc.g", file = "vignette_dat/sim50K.stratQg.notrunc.g.Rdata")
#  sim50K.stratQg.notrunc.missg <- foreach(sim = seq(Nsims), .combine = 'rbind')%dopar%{
#                        simrun_ltmleMSM(sim = sim,DAG = Ddyn, N = Nltmle, t0 = t0,
#                          gbounds = c(0.0000001, 1), gforms = gforms_miss)
#              }
#  save(list = "sim50K.stratQg.notrunc.missg", file = "vignette_dat/sim50K.stratQg.notrunc.missg.Rdata")

## ----appendixcode5, ref.label='chunk.runMSMsw', eval=FALSE, echo=TRUE, size='tiny'------
#  runMSMsw <- function(DAGO, Lnames, trueA, nsamp, nsims) {
#    Lnames_0 <- Lnames%+%"_0"
#    Lnames_1 <- Lnames%+%"_1"
#    gforms <- c("A_0 ~ "%+%paste(Lnames_0, collapse = " + "), "A_1 ~ A_0 + "%+%paste(Lnames_1, collapse = " + "))
#    res_sw <- NULL
#    for (sims in (1:nsims)) {
#      datO <- sim(DAGO, n = nsamp)
#      glmA_0 <- glm(datO[,c("A_0",Lnames_0)], formula = gforms[1], family = "binomial")
#      glmA_1 <- glm(datO[,c("A_1","A_0",Lnames_0,Lnames_1)], formula = gforms[2], family = "binomial")
#      probA0_1 <- predict(glmA_0,  type = "response")
#      weight_t0 <- 1 / (probA0_1^(datO$A_0) * (1-probA0_1)^(1-datO$A_0))
#      probA1_1 <- predict(glmA_1,  type = "response")
#      weight_t1 <- 1 / (probA1_1^(datO$A_1) * (1-probA1_1)^(1-datO$A_1))
#      sw1 <- weight_t0*weight_t1
#      emp.pA1cA0 <- table(datO$A_1,datO$A_0)/nrow(datO)
#      empPA1 <- data.frame(A_0 = c(0,0,1,1),A_1 = c(0,1,1,0))
#      empPA1$empPA_1_cA_0 <- apply(empPA1, 1, function(rowA) emp.pA1cA0[as.character(rowA["A_1"]), as.character(rowA["A_0"])])
#      empPA1 <- merge(datO[, c("ID","A_0","A_1")],empPA1, sort = FALSE)
#      empPA1 <- empPA1[order(empPA1$ID),]
#      swts <- empPA1$empPA_1_cA_0*(weight_t0*weight_t1)
#      datO$swts <- swts
#      MSMres_sw <- glm(datO, formula = "Y_1 ~ A_0 + A_1", weights = swts, family = "gaussian")
#      res_sw <- rbind(res_sw, coef(MSMres_sw))
#    }
#    meanres <- apply(res_sw, 2, mean)
#    Varres <- apply(res_sw, 2, var)
#    bias <- c(meanres["A_0"]-trueA["A_0"], meanres["A_1"]-trueA["A_1"])
#    MSE <- c(bias^2+Varres[c("A_0","A_1")])
#    bias10 <- sprintf("%.3f",bias*10)
#    MSE10 <- sprintf("%.3f",MSE*10)
#    resrow <- c(bias10[1], MSE10[1], bias10[2], MSE10[2])
#    col36names <- c("\\specialcell[t]{A(0)\\\\ Bias*10}",
#                    "\\specialcell[t]{A(0)\\\\ MSE*10}",
#                    "\\specialcell[t]{A(1)\\\\ Bias*10}",
#                    "\\specialcell[t]{A(1)\\\\ MSE*10}")
#    names(resrow) <- col36names
#    return(resrow)
#  }

## ----appendixcode6, ref.label='chunk.lefebrepT2andT4', eval=FALSE, echo=TRUE, size='tiny'----
#  # recreating Tables 2 & 4 reported in Lefebvre et al.
#  nsamp <- c(300,1000,10000)
#  # Lefebvre et al. Tab 2:
#  covnmT2 <- c(c("\\emph{Lefebvre et al.}: Confounder(s) only", rep("",2)),
#              c("\\emph{Lefebvre et al.}: Confounder(s) &", "risk factors", rep("",1)))
#  lefebvreT2 <- data.frame(
#    covnm = covnmT2,
#    N = rep(nsamp,2),
#    A0Bias10 = sprintf("%.3f",c(0.768, 0.265, 0.057, 0.757, 0.283, 0.056)),
#    A0MSE10 = sprintf("%.3f",c(1.761, 0.761, 0.146, 1.642, 0.718, 0.139)),
#    A1Bias10 = sprintf("%.3f",c(0.889, 0.312, 0.086, 0.836, 0.330, 0.081)),
#    A1MSE10 = sprintf("%.3f",c(1.728, 0.723, 0.120, 1.505, 0.638, 0.114)),stringsAsFactors = FALSE)
#  # Lefebvre et al. Tab 4:
#  covnmT4 <- c(c("\\emph{Lefebvre et al.}: Confounder(s) only", rep("",2)),
#                c("\\emph{Lefebvre et al.}: Confounder(s) &", "risk factors", ""),
#                c("\\emph{Lefebvre et al.}: Confounder(s) &", "IVs", ""),
#                c("\\emph{Lefebvre et al.}: Confounder(s),", "IVs & risk factors",""),
#                c("\\emph{Lefebvre et al.}: Mis-specified", rep("",2)),
#                c("\\emph{Lefebvre et al.}: Full Model", rep("",2)))
#  lefebvreT4 <- data.frame(
#    covnm = covnmT4,
#    N = rep(nsamp,6),
#    A0Bias10 = sprintf("%.3f",c(-0.080, -0.371, -0.368, -0.110, -0.330, -0.378, 1.611,
#                                0.824, 0.241, 1.600, 0.867, 0.235, 3.146, 2.460, 2.364,
#                                1.524, 0.878, 0.240)),
#    A0MSE10 = sprintf("%.3f",c(1.170, 0.385, 0.056, 1.092, 0.340, 0.051, 3.538, 2.063,
#                                0.684, 3.477, 2.053, 0.676, 3.326, 1.700, 0.832, 3.648,
#                                2.099, 0.679)),
#    A1Bias10 = sprintf("%.3f",c(0.099, -0.035, -0.203, 0.112, -0.108, -0.207, 2.069, 1.245,
#                                0.379, 2.143, 1.170, 0.372, 5.591, 5.258, 4.943, 2.221, 1.185,
#                                0.377)),
#    A1MSE10 = sprintf("%.3f",c(1.155, 0.331, 0.043, 0.865, 0.245, 0.037, 3.841, 2.188, 0.622,
#                                3.598, 2.043, 0.625, 5.494, 3.851, 2.705, 3.907, 2.099, 0.630)),
#                    stringsAsFactors = FALSE)
#  col1name <- "Covariates in $P(A|L)$"
#  colnames(lefebvreT2)[1] <- colnames(lefebvreT4)[1] <- col1name
#  col36names <- c("\\specialcell[t]{A(0)\\\\ Bias*10}",
#                  "\\specialcell[t]{A(0)\\\\ MSE*10}",
#                  "\\specialcell[t]{A(1)\\\\ Bias*10}",
#                  "\\specialcell[t]{A(1)\\\\ MSE*10}")
#  colnames(lefebvreT2)[3:6] <- colnames(lefebvreT4)[3:6] <- col36names

## ----appendixcodeSc1a, ref.label='chunk.lefebSc1a', eval=FALSE, echo=TRUE, size='tiny'----
#  trueA <- c(A_0 = -0.294, A_1 = -0.370)
#  nsims <- 10000; restab <- NULL
#  runsim <- function(Lnames, DAGO) {
#    for (nsamp in c(300,1000,10000)) {
#      resSc <- runMSMsw(DAGO = DAGO, Lnames = Lnames, trueA = trueA, nsamp = nsamp, nsims = nsims)
#      restab <- rbind(restab, c(N = nsamp, resSc))
#    }
#    restab
#  }
#  Lnames <- c("LC1")
#  covnm <- c("Confounder(s) only", rep("",2))
#  restab_1 <- cbind(covnm, runsim(Lnames, DAGO.sc1))
#  # restab_1 <- rbind(restab_1, as.matrix(lefebvreT2[1:3,]))
#  Lnames <- c("LC1", "LO1", "LO2", "LO3")
#  covnm <- c("Confounder(s) &", "risk factors", rep("",1))
#  restab_2 <- cbind(covnm, runsim(Lnames, DAGO.sc1))
#  # restab_2 <- rbind(restab_2, as.matrix(lefebvreT2[4:6,]))
#  restab <- rbind(restab_1, restab_2)
#  col1name <- "Covariates in $P(A|L)$"
#  colnames(restab)[1] <- col1name
#  # restabwLef <- restab
#  # save(list = "restabwLef", file = "vignette_dat/restabwLefSc1_all_1Ksims.Rdata");
#  # restab <- restab[c(1:3, 7:9),]
#  save(list = "restab", file = "vignette_dat/restabSc1_all_1Ksims.Rdata");

## ----appendixcodeSc1b, ref.label='chunk.lefebSc1b', eval=FALSE, echo=TRUE, size='tiny'----
#  library("Hmisc")
#  load(file = "vignette_dat/restabSc1_all_1Ksims.Rdata");
#  cat("\n")
#  latex(restab, file = "", where = "!htpb", caption.loc = 'bottom',
#      caption = "Replication of the simulation results from \\citet{lefebvre2008} for Scenario 1.",
#      label = 'tab2Lefebvre',booktabs = TRUE,rowname = NULL,landscape = FALSE,
#      col.just = c("l", rep("r", 5)), size = "small")

## ----appendixcodeSc1c, ref.label='chunk.lefebSc1c', eval=FALSE, echo=TRUE, size='tiny'----
#  cat("\n")
#  latex(lefebvreT2, file = "", where = "!htpb", caption.loc = 'bottom',
#      caption = "Simulation results for Scenario 1 as reported in Table II of \\citet{lefebvre2008}.",
#      label = 'origtab2Lefebvre', booktabs = TRUE, rowname = NULL, landscape = FALSE,
#      col.just = c("l", rep("r", 5)), size = "small")

## ----appendixcodeSc3a, ref.label='chunk.lefebSc3a', eval=FALSE, echo=TRUE, size='tiny'----
#  trueA <- c(A_0 = -0.316, A_1 = -0.390)
#  nsims <- 10000; restab <- NULL
#  runsim <- function(Lnames, DAGO) {
#    for (nsamp in c(300,1000,10000)) {
#      resSc <- runMSMsw(DAGO = DAGO, Lnames = Lnames, trueA = trueA, nsamp = nsamp, nsims = nsims)
#      restab <- rbind(restab, c(N = nsamp, resSc))
#    }
#    restab
#  }
#  Lnames <- c("LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s) only", rep("",2))
#  restab_1 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_1 <- rbind(restab_1, as.matrix(lefebvreT4[1:3,]))
#  Lnames <- c("LO1", "LO2", "LO3", "LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s) &", "risk factors", "")
#  restab_2 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_2 <- rbind(restab_2, as.matrix(lefebvreT4[4:6,]))
#  Lnames <- c("LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s) &", "IVs", "")
#  restab_3 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_3 <- rbind(restab_3, as.matrix(lefebvreT4[7:9,]))
#  Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3")
#  covnm <- c("Confounder(s),", "IVs & risk factors","")
#  restab_4 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_4 <- rbind(restab_4, as.matrix(lefebvreT4[10:12,]))
#  Lnames <- c("LE1", "LE2", "LE3", "LC1")
#  covnm <- c("Mis-specified", rep("",2))
#  restab_5 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_5 <- rbind(restab_5, as.matrix(lefebvreT4[13:15,]))
#  Lnames <- c("LO1", "LO2", "LO3", "LE1", "LE2", "LE3", "LC1", "LC2", "LC3", "LN1", "LN2", "LN3")
#  covnm <- c("Full Model", rep("",2))
#  restab_6 <- cbind(covnm, runsim(Lnames, DAGO.sc3))
#  # restab_6 <- rbind(restab_6, as.matrix(lefebvreT4[16:18,]))
#  restab <- rbind(restab_1, restab_2, restab_3, restab_4, restab_5, restab_6)
#  col1name <- "Covariates in $P(A|L)$"
#  colnames(restab)[1] <- col1name
#  # restabwLef <- restab
#  # save(list = "restabwLef", file = "vignette_dat/restabwLefSc3_all_1Ksims.Rdata");
#  # restab <- restab[c(1:3, 7:9, 13:15, 19:21, 25:27, 31:33),]
#  save(list = "restab", file = "vignette_dat/restabSc3_all_1Ksims.Rdata");

## ----appendixcodeSc3b, ref.label='chunk.lefebSc3b', eval=FALSE, echo=TRUE, size='tiny'----
#  library("Hmisc")
#  load(file = "vignette_dat/restabSc3_all_1Ksims.Rdata");
#  cat("\n")
#  latex(restab,file = "",where = "!htpb", caption.loc = 'bottom',
#      caption = "Replication of the simulation results from \\citet{lefebvre2008} for Scenario 3.",
#      label = 'tab4Lefebvre',booktabs = TRUE,rowname = NULL,landscape = FALSE,
#      col.just = c("l", rep("r", 5)), size = "small")

## ----appendixcodeSc3c, ref.label='chunk.lefebSc3c', eval=FALSE, echo=TRUE, size='tiny'----
#  cat("\n")
#  latex(lefebvreT4,file = "",where = "!htpb", caption.loc = 'bottom',
#      caption = "Simulation results for Scenario 3 as reported in Table IV of \\citet{lefebvre2008}.",
#      label = 'origtab4Lefebvre',booktabs = TRUE,rowname = NULL,landscape = FALSE,
#      col.just = c("l", rep("r", 5)), size = "small")

