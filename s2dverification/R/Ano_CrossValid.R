Ano_CrossValid <- function(var_exp, var_obs, memb = TRUE) {
  #
  #  Enlarge the number of dimensions of var_exp and var_obs to 7 if necessary
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimexp <- dim(var_exp)
  dimobs <- dim(var_obs)
  if (length(dimexp) < 4 | length(dimobs) < 4) {
    stop("At least 4 dim needed : c(nexp/nobs, nmemb, nsdates, nltime)")
  }
  for (jn in 3:max(length(dimexp), length(dimobs))) {
    if (dimexp[jn] != dimobs[jn]) { 
      stop("Wrong input dimensions")
    }
  }
  var_exp <- Enlarge(var_exp, 7)
  var_obs <- Enlarge(var_obs, 7)
  #
  #  Computation of anomalies in cross-validation mode
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enl_ano_exp <- array(dim = dim(var_exp))
  enl_ano_obs <- array(dim = dim(var_obs))
  for (jsdat in 1:dimexp[3]) {
    tmp1 <- array(var_exp[, , jsdat, , , , ], dim = dim(var_exp)[-3])
    tmp2 <- array(var_obs[, , jsdat, , , , ], dim = dim(var_obs)[-3])
    tmp3 <- array(var_exp[, , -jsdat, , , , ], dim = c(dimexp[1:2], 
                                                       dimexp[3] - 1, 
                                                       dim(var_exp)[4:7]))
    tmp4 <- array(var_obs[, , -jsdat, , , , ], dim = c(dimobs[1:2], 
                                                       dimobs[3] - 1, 
                                                       dim(var_obs)[4:7]))
    tmp <- Clim(tmp3, tmp4, memb)
    if (memb) { 
      clim_exp <- tmp$clim_exp
      clim_obs <- tmp$clim_obs
    } else {
      clim_exp <- InsertDim(tmp$clim_exp, 2, dimexp[2]) 
      clim_obs <- InsertDim(tmp$clim_obs, 2, dimobs[2]) 
    }
    enl_ano_exp[, , jsdat, , , , ] <- tmp1 - clim_exp
    enl_ano_obs[, , jsdat, , , , ] <- tmp2 - clim_obs  
  } 
  #
  #  Reduce the number of dimensions to the original one 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  ano_exp <- array(dim = dimexp)
  ano_obs <- array(dim = dimobs)
  ano_exp[] <- enl_ano_exp
  ano_obs[] <- enl_ano_obs
  #
  #  Outputs
  # ~~~~~~~~~
  #
  invisible(list(ano_exp = ano_exp, ano_obs = ano_obs))
}
