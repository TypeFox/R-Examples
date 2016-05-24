  
####################################################################
###########  function 6: check basic numeric properties  ############
####################################################################
# small function
CheckIt = function(vector_in) {
   tmp1 = any(is.na(vector_in) | is.nan(vector_in) | is.infinite(vector_in))
   cat("^^any na or nan:",tmp1,"...")
   
   if(is.list(vector_in)) stop("groupwise data is list","\n")
   
   if (is.vector(vector_in)) {
   cat("Min is:",min(vector_in)," Max is:", max(vector_in),"...")
   cat("length is",length(vector_in),"\n") }
   
   if (is.matrix(vector_in)) {
      nr = nrow(vector_in); nc = ncol(vector_in)
      cat("nrow:",nr,"; ncol:",nc,"\n")
     if (nc == 2) {
          rs1 =  rowSums(vector_in)
         cat("^^min rowSum:",min(rs1),"^^Max rowSum:",max(rs1),"\n") 
         }
    if(nc > 2) {
       rs1 =  rowSums(vector_in[,1:floor((nc/2))])
       rs2 =  rowSums(vector_in[,((1+floor((nc/2))):nc)]) 
     cat("^^Min rowSum gp1:",min(rs1),"; min of rowSum gp2:",min(rs2),"; Max rowTotal",max(rs1+rs2),"\n")  }
     }
   }  
 

###########################################################################
#### function: updated main grouping with specified minimal group size 
#######################################################################
### function
# FOR PRACTICAL PURPOSED to_merge is always true.

eNetBuilder = function(div_mat = NULL,ngp = NULL, merge_size = NULL,rad = NULL)
  {
  if (!is.matrix(div_mat))
    stop("^^ A matrix with entries as divergences is needed ...","\n")
  
  m_lt = nrow(div_mat)

  ngrptmp = 1; cnt_idx = double(0);  div_mat_in = div_mat
  all_idx = 1:m_lt;  # caution, no longer all_idx = 1:m 
  
  idx_left = all_idx ; itn = length(idx_left)
  centers_id = double(0) ;    grp_idx = vector("list",0)

  while(itn & (ngrptmp <= ngp)) {
     cat("^^Forming group", ngrptmp, "....","\n")
   ball = vector("list",itn);  ball_size = double(itn)

   # compare itn balls
   for (i in 1:itn) {
        ball_id_tmp = which(abs(div_mat_in[i,]) <= rad)
        ball[[i]] = idx_left[ball_id_tmp]    # ids in ball from original 1:nrow(div_mat) indices
        ball_size[i] = length(ball[[i]])    }

    smax_id = which.max(ball_size);  id_mxball = idx_left[smax_id]
    ball_max = ball[[smax_id]];  size_mxball = ball_size[smax_id]

    cat("^^Ids: e-ball center", id_mxball,". Ball size:", size_mxball, "\n") # "ids in ball:",ball_max,"\n")

    centers_id = c(centers_id,id_mxball);     grp_idx[[ngrptmp]] = ball_max
    cnt_idx = c(cnt_idx,ball_max)
    # added check
    cat("^^# of Id used",length(cnt_idx),"\n")
    
    if (any(cnt_idx <0))     stop("^^negative id in current ids","\n") 
    
    idx_left = all_idx[-cnt_idx]  # take current ids from 1:m

    itn = length(idx_left)
    cat("^^# of ids left currently",itn,"\n") # ". Idx left:", idx_left,"\n")

    # Follwing is the best partition, case 1:  ngrptmp <= ngp-1
    if (ngrptmp <= ngp-1) {
       if (itn <= merge_size) {  
         # the above condition should be itn <= merge_size instead of itn < merge_size since it is when ngrptmp <= ngp-1
         cat("^^Can not form",ngp,"groups each with cardinality no less than",merge_size,"\n")
         cat("^^---- Regroup with smaller e-net size ---- ","\n")
         rad = rad/2  # adjust this
         idx_left = 1:m_lt; itn = length(idx_left)
         ngrptmp = 0;  cnt_idx = double(0);  div_mat_in = div_mat
         }

       if (itn == merge_size & ngrptmp == ngp-1) {
         cat("^^Merging",itn,"ids.",ngp,"e-balls reached as requested","\n")
         grp_idx[[ngp]] =  idx_left
         itn = 0
        }
       
       if (itn > merge_size) {
          cat("^# of id's left to group is", itn,",continue grouping","\n")
         div_mat_in = div_mat[idx_left,][,idx_left]
         }
    }  # end  if (ngrptmp <= ngp-1)

    # case 2:  ngrptmp = ngp
    if (ngrptmp == ngp) {
      # nothing left 
      if (itn == 0)
        cat("^^",ngp,"e-balls reached as requested","\n")

      if (itn > merge_size) {
       cat("^^---Some ids left for more groups. Regroup with larger e-net size --- ","\n")
        rad = 1.5*rad  # adjust this
        idx_left = 1:m_lt; itn = length(idx_left)
         ngrptmp = 0;  cnt_idx = double(0);  div_mat_in = div_mat
       }  
       
       if (itn <= merge_size) {
           grp_idx[[ngp]] = c(grp_idx[[ngp]],idx_left)
          cat("^^Merging",itn,"ids.",ngp,"e-balls reached as requested","\n")
           }
      } # end if (ngrptmp == ngp)

    ngrptmp = ngrptmp + 1

  } # end of while
      return (grp_idx)
} # end of fcuntion

### ###########################
eNetFull = function(metrics_in = NULL, ngrp_in = NULL, merge_size = NULL, rad_in = NULL, mgpsize_in = NULL) {

    if(is.null(mgpsize_in))
     stop("^^Specify minimal group size (MGS)")
    rad_tmp = rad_in;  mgps_tmp = 1  # mininal number of groups and group size

    while (mgps_tmp < mgpsize_in) {

      groups_tmp = eNetBuilder(metrics_in,ngrp_in,merge_size,rad_tmp)
      ngp_tmp = length(groups_tmp)
  #    cat("^^Current # of groups",ngp_tmp,"\n")
      mgps_tmp = length(groups_tmp[[ngp_tmp]])

    } # end of while
    return(groups_tmp)
}   # end of function
 


########################## ########################## ##########################  
################## function:  divergence  matrix   ########################## 
########################## ########################## ##########################

GetDivergenceMatrix = function(scalfac=NULL,pvSpList=NULL)  {   
  
  lgA1 = length(pvSpList)     # caution: no longer m but lgA1
  lgt_div = (lgA1-1)*lgA1/2
  cat("^^Computing", lgt_div, "pairwise divergences. Takes time ...","\n")
  
  chi_mat = matrix(0,lgA1,lgA1); infNorm_mat = matrix(0,lgA1,lgA1); div_mat = matrix(0,lgA1,lgA1)
  
  # start of computing pairwise quantities
  for (i in 2:lgA1) { 
       sp_pv1 = pvSpList[[i]];  lg1 = length(sp_pv1)     #pvSpList has been formated
           
     for (j in 1:(i-1)) {       
       sp_pv2 = pvSpList[[j]]; lg2 = length(sp_pv2)    #pvSpList has been formated 
       
       chi_mat[i,j] = abs(lg1 - lg2)
     
       teva = union(sp_pv1,sp_pv2)  # evaluate cdf at union of pvalue supports 
       cdfn1 = pvalueDist(teva,sp_pv1);  cdfn2 = pvalueDist(teva,sp_pv2)
       
       infNorm_mat[i,j] = max(abs(cdfn1-cdfn2)) 
        }
    }   # end of computing pairwise quantities
    
      ## compute divergences
      chiWgt = max(chi_mat); infNormWgt = max(infNorm_mat)
      
      chi_mat = chi_mat + t(chi_mat) ;  infNorm_mat = infNorm_mat + t(infNorm_mat)
      
      if (chiWgt ==0 | infNormWgt == 0)   cat("^^ Homogeneous null distributions ...","\n")
  
      if (chiWgt ==0 & infNormWgt != 0)  {
         cat("^^ Homogeneous supports ...","\n") 
         div_mat = scalfac*(infNorm_mat/infNormWgt)
      }
      
      if (chiWgt != 0 & infNormWgt == 0)  {
         cat("^^ Homogeneous masses ...","\n") 
         div_mat = scalfac*(chi_mat/chiWgt)
      }
      
      if (chiWgt != 0 & infNormWgt != 0) {
        cat("^^ Heterogeneous supports and masses ...","\n")
        div_mat = scalfac*(chi_mat + infNorm_mat)/(chiWgt + infNormWgt)
      }
      
      div_mat = div_mat + t(div_mat) # since diagnal is zero, this is correct
      cat("^^Finished computing matrix of pairwise divergences...","\n")
       
      return(list(div_mat,chi_mat/chiWgt,infNorm_mat,chiWgt))
      
  }
  
######################### ########################## ##########################   
  ### Identify Approximate Uniform, this function also formats different psupport
  ######################### ########################## ########################## 
  
  # caution: pvalueByBinoSupport   for Binomial test     support<- c(meantmp,pvalue, psupport)
# pvalueSupport  for FET    support<- c(meantmp,psupport)
# pvalueByNegativeBinoSupportApp for ENT   support<- c(meantmp,pvalue,psupport)
 
  Div_Appr_Unif = function(pvSpList=NULL,test_used=NULL,appr_unif_tol = NULL) {
      lg3 = length(pvSpList)
      pv_supp_formatted = vector("list",lg3)
      id_appr_unif = double(0)
       
     # caution: psupport has different elements 
       if (test_used == "Binomial Test" | test_used == "Exact Negative Binomial Test") 
         {
           for (i in 1:lg3) { 
             sp_pv = pvSpList[[i]][-1][-1] 
             pv_supp_formatted[[i]] = sp_pv
             
             sp_pv_tmp = c(0,sp_pv)
             if (max(abs(diff(sp_pv_tmp))) < appr_unif_tol)      
              id_appr_unif = c(id_appr_unif,i)
          } 
       }  
        
       if (test_used == "Fisher's Exact Test") {
          for (i in 1:lg3) { 
             sp_pv = pvSpList[[i]][-1]
             pv_supp_formatted[[i]] = sp_pv
             
             sp_pv_tmp = c(0,sp_pv)
            if (max(abs(diff(sp_pv_tmp))) < appr_unif_tol)      
               id_appr_unif = c(id_appr_unif,i)
            } 
          }                
    return(list(pv_supp_formatted,id_appr_unif))       
  } # end of func
 
 
########################################################################
### Compute supremum norm to reference Unif
########################################################################
 Div_Ref_Unif = function(pvSpList=NULL,test_used=NULL) {

      lg3 = length(pvSpList)
      d_ref_unif = double(lg3)
      s_supp = double(lg3)

     # caution: psupport has different elements
       for (i in 1:lg3) {
         if (test_used == "Binomial Test" | test_used == "Exact Negative Binomial Test")
           sp_pv = pvSpList[[i]][-1][-1]

         if (test_used == "Fisher's Exact Test")
           sp_pv = pvSpList[[i]][-1]

         # using unif as reference distribution
         sp_pv_tmp = c(0,sp_pv)
         d_ref_unif[i] = max(abs(diff(sp_pv_tmp)))
         s_supp[i] = length(sp_pv)
       }
    return(list(d_ref_unif,s_supp))
  } # end of func
   
  

################## edgeR extracted function: estimateCommonDisp
#  added: if one library size is computed zero, stop since normalization is not possible
# changed: number of times optimize is applied from 2 to 5
# changed: end point of search interval in optimize
# added: cat("^^EdgeR: estimate dispersion in process. delta = ")

estimateCommonDisp_adjusted = function (object, nopsinestdisp = NULL, tol = 1e-06, rowsum.filter = 5, verbose = FALSE)
{
    if (!is(object, "DGEList"))
        stop("Currently supports DGEList objects")
    group <- object$samples$group <- as.factor(object$samples$group)
    if (all(tabulate(group) <= 1)) {
        warning("There is no replication, setting dispersion to NA.")
        object$common.dispersion <- NA
        return(object)
    }
    
    # added 3/20/2013
    if (any(object$samples$lib.size == 0)) 
          stop("^^edgeR: zero library size found in DGEList object ...")   #
    
    tags.used <- rowSums(object$counts) > rowsum.filter
    pseudo.obj <- object[tags.used, ]
        
    disp <- 0.01    # default in edgeR
    
    # edegR defult  for (i in 1:2)
    for (i in 1:nopsinestdisp) {
        out <- equalizeLibSizes_adjusted(object, dispersion = disp)   # changed from equalizeLibSizes into equalizeLibSizes_adjusted
        pseudo.obj$counts <- out$pseudo.counts[tags.used, , drop = FALSE]
        y <- splitIntoGroups(pseudo.obj)  
            
        # interval need to be changed, optimize default in edgeR interval = c(1e-04,100/(100 + 1))
          delta <- optimize(commonCondLogLikDerDelta, interval = c(1e-04,
            1000/(1000 + 1)), tol = tol, maximum = TRUE, y = y,
            der = 0)
                
       # cat("^^EdgeR: estimate dispersion in process. delta = ")  
                      
        delta <- delta$maximum
        print(delta)
        
        disp <- delta/(1 - delta)
    }
    if (verbose)
        cat("EdgeR: Disp =", round(disp, 5), ", BCV =", round(sqrt(disp),
            4), "\n")
    object$common.dispersion <- disp
    object$pseudo.counts <- out$pseudo.counts
    effective.lib.size <- object$samples$lib.size * object$samples$norm.factors
    abundance <- mglmOneGroup(object$counts, dispersion = disp,
        offset = log(effective.lib.size))
    object$logCPM <- log1p(exp(abundance + log(1e+06)))/log(2)
    object$pseudo.lib.size <- out$common.lib.size
    object
}


## keyfunction DGEList. 
# only the message  "Calculating library sizes from column totals" changed to 
# "EdgeR: Calculating library sizes from column totals."

DGEList_adjusted = function (counts = matrix(0, 0, 0), lib.size = NULL, norm.factors = NULL,
    group = rep.int(1, ncol(counts)), genes = NULL, remove.zeros = FALSE)
{
    counts <- as.matrix(counts)
    nlib <- ncol(counts)
    ntags <- nrow(counts)
    if (nlib > 0 && is.null(colnames(counts)))
        colnames(counts) <- paste("Sample", 1:ncol(counts), sep = "")
    rn <- rownames(counts)
    if (is.null(lib.size)) {
        lib.size <- colSums(counts)
        message("EdgeR: Calculating library sizes from column totals.")
    }
    else {
        if (nlib != length(lib.size))
            stop("Length of 'lib.size' must equal number of columns in 'counts'")
    }
    if (is.null(norm.factors)) {
        norm.factors <- rep(1, nlib)
    }
    else {
        if (nlib != length(norm.factors))
            stop("Length of 'norm.factors' must equal number of columns in 'counts'")
    }
    group <- as.factor(group)
    if (nlib != length(group))
        stop("Length of 'group' must equal number of columns in 'counts'")
    samples <- data.frame(group = group, lib.size = lib.size,
        norm.factors = norm.factors)
    row.names(samples) <- colnames(counts)
    x <- new("DGEList", list(counts = counts, samples = samples))
    if (!is.null(genes)) {
        genes <- as.data.frame(genes, stringsAsFactors = FALSE)
        if (nrow(genes) != ntags)
            stop("counts and genes have different nrows")
        if (!is.null(rn))
            rownames(genes) <- rn
        x$genes <- genes
    }
    if (remove.zeros) {
        all.zeros <- rowSums(counts, na.rm = TRUE) == 0
        if (any(all.zeros)) {
            x <- x[!all.zeros, ]
            message("Removing ", sum(all.zeros), " rows with all zero counts.")
        }
    }
    x
}

### in the original, argument dispersion =0, now it is set dispersion = NULL. otherwise running
## on server error: unsued arguments(0): dispersion = disp when this function is 
# called and executed by estimateCommonDisp_adjusted
equalizeLibSizes_adjusted = function (object, dispersion = NULL, common.lib.size = NULL) 
{
    d <- dim(object)
    ntags <- d[1]
    nlibs <- d[2]
    lib.size <- object$samples$lib.size * object$samples$norm.factors
    if (is.null(common.lib.size)) 
        common.lib.size <- exp(mean(log(lib.size)))
    levs.group <- unique(object$samples$group)
    if (length(dispersion) == 1) 
        dispersion <- rep(dispersion, ntags)
    input.mean <- output.mean <- matrix(0, ntags, nlibs)
    for (i in 1:length(levs.group)) {
        j <- object$samples$group == levs.group[i]
        beta <- mglmOneGroup(object$counts[, j, drop = FALSE], 
            dispersion = dispersion, offset = log(lib.size[j]))
        lambda <- exp(beta)
        input.mean[, j] <- matrix(lambda, ncol = 1) %*% matrix(lib.size[j], 
            nrow = 1)
        output.mean[, j] <- matrix(lambda, ncol = 1) %*% matrix(common.lib.size, 
            nrow = 1, ncol = sum(j))
    }
    pseudo <- q2qnbinom(object$counts, input.mean = input.mean, 
        output.mean = output.mean, dispersion = dispersion)
    pseudo[pseudo < 0] <- 0
    list(pseudo.counts = pseudo, common.lib.size = common.lib.size)
}


## edgeR: non adjusted q2qnbinom
q2qnbinom_notadjusted = function (x, input.mean, output.mean, dispersion = 0) 
{
    if (any(x < 0)) 
        stop("x must be non-negative")
    if (any(input.mean < 0)) 
        stop("input.mean must be non-negative")
    if (any(output.mean < 0)) 
        stop("output.mean must be non-negative")
    if (any(dispersion < 0)) 
        stop("dispersion must be non-negative")
    eps <- 1e-14
    zero <- input.mean < eps | output.mean < eps
    input.mean[zero] <- input.mean[zero] + 0.25
    output.mean[zero] <- output.mean[zero] + 0.25
    ri <- 1 + dispersion * input.mean
    vi <- input.mean * ri
    ro <- 1 + dispersion * output.mean
    vo <- output.mean * ro
    i <- (x >= input.mean)
    j <- !i
    p1 <- p2 <- q1 <- q2 <- x
    if (any(i)) {
        p1[i] <- pnorm(x[i], mean = input.mean[i], sd = sqrt(vi[i]), 
            lower.tail = FALSE, log.p = TRUE)
        p2[i] <- pgamma(x[i], shape = input.mean[i]/ri[i], scale = ri[i], 
            lower.tail = FALSE, log.p = TRUE)
        q1[i] <- qnorm(p1[i], mean = output.mean[i], sd = sqrt(vo[i]), 
            lower.tail = FALSE, log.p = TRUE)
        q2[i] <- qgamma(p2[i], shape = output.mean[i]/ro[i], 
            scale = ro[i], lower.tail = FALSE, log.p = TRUE)
    }
    if (any(j)) {
        p1[j] <- pnorm(x[j], mean = input.mean[j], sd = sqrt(vi[j]), 
            lower.tail = TRUE, log.p = TRUE)
        p2[j] <- pgamma(x[j], shape = input.mean[j]/ri[j], scale = ri[j], 
            lower.tail = TRUE, log.p = TRUE)
        q1[j] <- qnorm(p1[j], mean = output.mean[j], sd = sqrt(vo[j]), 
            lower.tail = TRUE, log.p = TRUE)
        q2[j] <- qgamma(p2[j], shape = output.mean[j]/ro[j], 
            scale = ro[j], lower.tail = TRUE, log.p = TRUE)
    }
    (q1 + q2)/2
}

### key function in edgeR :  mglmOneGroup
# a dispersion is estimated from estimateCommonDisp, then its fed into  mglmOneGroup
# to get abundance, then abundance is used for normalization

mglmOneGroup_notadjusted = function (y, dispersion = 0, offset = 0, maxit = 50, trace = FALSE, 
    tol = 1e-06) 
{
    y <- as.matrix(y)
    if (any(y < 0)) 
        stop("y must be non-negative")
    ntags <- nrow(y)
    nlibs <- ncol(y)
    beta <- rep(-Inf, ntags)
    names(beta) <- rownames(y)
    if (any(dispersion < 0)) 
        stop("dispersion must be non-negative")
    N <- exp(offset)
    if (all(dispersion == 0)) {
        if (is.null(dim(N))) 
            m <- mean(N)
        else m <- .rowMeans(N, ntags, nlibs)
        return(log(.rowMeans(y/m, ntags, nlibs)))
    }
    dispersion <- rep(dispersion, length = ntags)
    offset <- expandAsMatrix(offset, dim(y))
    N <- expandAsMatrix(N, dim(y))
    beta <- log(.rowMeans(y/N, ntags, nlibs))
    if (nlibs == 1) 
        return(beta)
    iter <- 0
    i <- is.finite(beta)
    while (any(i)) {
        iter <- iter + 1
        if (iter > maxit) {
            warning("max iterations exceeded")
            return(beta)
        }
        if (trace) 
            cat("Iter=", iter, "Still converging=", sum(i), "\n")
        mu <- exp(beta[i] + offset[i, , drop = FALSE])
        var.div.mu <- 1 + dispersion[i] * mu
        m <- nrow(mu)
        dl <- .rowSums((y[i, , drop = FALSE] - mu)/var.div.mu, 
            m, nlibs)
        info <- .rowSums(mu/var.div.mu, m, nlibs)
        step <- dl/info
        beta[i] <- beta[i] + step
        i[i] <- abs(step) > tol
    }
    beta
}     
 