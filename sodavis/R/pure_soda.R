library(nnet)

# create predictor matrix from terms
create_pmatrix_from_terms = function(xx, terms)
{
  nt = length(terms);
  nr = nrow(xx);
  pmatrix = matrix(0, nr, 0);
  
  if (nt > 0)
  for(it in 1:nt)
  {
    term = terms[it];
    if (grepl("*",term,fixed=TRUE))
    {
      splits = strsplit(term,"*",fixed=TRUE)[[1]];
      id1 = as.numeric(splits[1]);
      id2 = as.numeric(splits[2]);
      pmatrix = cbind(pmatrix, term=xx[,id1]*xx[,id2]);
    }
    else
    {
      id  = as.numeric(term);
      pmatrix = cbind(pmatrix, term=xx[,id]);
    }
  }
  return(pmatrix);
}


calc_lda_BIC = function(xx, yy, cur_set, D, K, debug=FALSE, gam=0)
{
  N = nrow(xx);
  D = ncol(xx);
  K = max(yy);
  d = length(cur_set);  
  
  ll = 0;
  if (d == 0)
  {
    p = numeric(K);
    for(i in 1:N)
    {
      p[yy[i]] = p[yy[i]] + 1.0/N;
    }
    for(k in 1:K)
    {
      ll = ll + sum(yy==k)*log(p[k]);
    }
    BIC = -2*ll + (K-1)*(log(N) + 2*gam*log(D));
    return(BIC);
  } 
  else 
  {
    pmatrix = as.matrix(xx[,cur_set]);
    colnames(pmatrix) = NULL;
    lgt = multinom(yy ~ pmatrix, family = "binomial", trace=FALSE);
    BIC = lgt$deviance + (K-1)*(1+d)*(log(N) + 2*gam*log(D));
    return(BIC);
  }
}


#
#      xx: explanatory variables
#      yy: response variable
# cur_set: current set of selected variables
#   debug: if shows debug information
#     gam: gamma in EBIC
#   terms: selected linear and interaction terms
#
calc_BIC = function(xx, yy, terms, D, K, debug=FALSE, gam=0)
{
  N = length(yy);
  D = ncol(xx);
  K = max(yy);
  d = length(terms);  
  
  ll = 0;
  if (d == 0)
  {
    p = numeric(K);
    for(i in 1:N)
    {
      p[yy[i]] = p[yy[i]] + 1.0/N;
    }
    for(k in 1:K)
    {
      ll = ll + sum(yy==k)*log(p[k]);
    }
    BIC = (K-1) * (log(N) + 2*gam*log(D));
    BIC = BIC - 2*ll;
    return(BIC);
  } 
  else
  {
    pmatrix = create_pmatrix_from_terms(xx, terms);
    colnames(pmatrix) = NULL;
    lgt = multinom(yy ~ pmatrix, family = "multinomial", trace=FALSE);
    BIC = lgt$deviance;
    BIC = BIC + (K-1)*(1+ncol(pmatrix))*(log(N) + 2*gam*log(D)); # BIC with quadratic penalty
    
    #     cat(sprintf("Degrees of freedom  = %d\n", (1+ncol(pmatrix))));
    #     print(pmatrix[1:2,])
    #     if (readline() != "")
    #       stop("Stop!");
    return(BIC);
  }
}

get_term_name = function(x_names, cur_set, term)
{
  if (length(cur_set) == 0)
  {
    return("(Empty)");
  }
  str = "(Empty)";
  splits = strsplit(term, ".", fixed=TRUE)[[1]]
  first = T;
  for(i in 1:length(cur_set))
  {
    split = splits[i];
    if (split=="1")
    {
      if (first)
      {
        first = F;
        str = x_names[cur_set[i]];
      }
      else
      {
        str = paste0(str, "*", x_names[cur_set[i]]);
      }
    }
    if (split=="2")
    {
      if (first)
      {
        first = F;
        str = paste0(x_names[cur_set[i]], "^2");
      }
      else
      {
        str = paste0(str, "*", paste0(x_names[cur_set[i]], "^2"));
      }
    }      
  }
  return (str);
}

get_term_name_2 = function(x_names, term)
{
  if (grepl("*",term,fixed=TRUE))
  {
    splits = strsplit(term,"*",fixed=TRUE)[[1]];
    id1 = as.numeric(splits[1]);
    id2 = as.numeric(splits[2]);
    return(paste(x_names[id1], x_names[id2], sep="*"))
  }
  else
  {
    id  = as.numeric(term);
    return(x_names[id])
  }
}

get_lin_terms = function(n_terms)
{
  terms = c();
  for(i in 1:n_terms)
  {
    arr = numeric(n_terms);
    arr[i] = 1;
    terms = c(terms, paste0(arr, collapse = "."));
  }
  return(terms);
}

get_lin_terms_vec = function(c_set)
{
  terms = as.character(c_set);
  return(terms);
}

get_quad_terms_vec = function(c_set)
{
  terms = c();
  if (length(c_set) <= 0)
    return(terms);
  
  terms = get_lin_terms_vec(c_set);  
  for(i in 1:length(c_set))
    for(j in i:length(c_set))
      terms = c(terms, paste(c_set[i],c_set[j],sep="*"));
  
  return(terms);
}

get_inter_terms_vec = function(c_set)
{
  terms = c();
  if (length(c_set) < 2)
    return(terms);
  
  for(i in 1:(length(c_set)-1))
  {
    for(j in (i+1):length(c_set))
      terms = c(terms, paste(c_set[i],c_set[j],sep="*"));
  }
  
  return(terms);
}

get_set_from_terms = function(terms)
{
  nt = length(terms);
  c_set = c();
  if (nt > 0)
  for(it in 1:nt)
  {
    term = terms[it];
    if (grepl("*",term,fixed=TRUE))
    {
      splits = strsplit(term,"*",fixed=TRUE)[[1]];
      id1 = as.numeric(splits[1]);
      id2 = as.numeric(splits[2]);
      c_set = union(c_set, id1);
      c_set = union(c_set, id2);
    }
    else
    {
      id  = as.numeric(term);
      c_set = union(c_set, id);
    }
  }
  return(c_set)
}

trim_terms = function(terms)
{
  if (length(terms) == 0)
    return(c());
  
  splits = strsplit(terms[1], ".", fixed=TRUE)[[1]]
  Nvar = length(splits);
  var_set = rep(F, Nvar);
  
  if (Nvar == 0)
  {
    return(c());
  }
  
  for(term in terms)
  {
    splits = strsplit(term, ".", fixed=TRUE)[[1]]
    Nvar = length(splits);
    for(i in 1:Nvar)
    {
      split = splits[i];
      if (split != "0")
        var_set[i] = TRUE;
    }
  }
  
  new_terms = c();
  for(term in terms)
  {
    splits = strsplit(term, ".", fixed=TRUE)[[1]]
    splits = splits[which(var_set)];
    new_term = paste0(splits, collapse = ".");
    new_terms = c(new_terms, new_term)
  }

  return(new_terms);
}

nqnorm = function(data)
{
  if (class(data)=="numeric")
  {
    N = length(data);
    qs = rank(data)/N - 0.5/N;
    data = qnorm(qs);
    return(data);
  }
  else
  {
    N = dim(data)[1];
    D = dim(data)[2];
    for(d in 1:D)
    {
      qs = rank(data[,d])/N - 0.5/N;
      data[,d] = qnorm(qs);
    }
    return(data);
  }
}

#
#    xx: explanatory variables
#    yy: response variable
#  norm: if TRUE, xx are quantile normalized to normal
# debug: if shows debug information
#   gam: gamma in EBIC
#  minF: minimum number of forward interaction screening steps
#
soda = function(xx, yy, norm=FALSE, debug=FALSE, gam=0, minF = 3)
{
  if (min(yy) == 0)
    yy = yy + 1;
  
  K = max(yy);
  N = dim(xx)[1];
  D = dim(xx)[2];
  
  minF = min(D, minF);
  
  if (norm)
  {
    for(k in 1:K)
      xx[yy=k,] = nqnorm(xx[yy=k,]);
  }
  
  x_names = colnames(xx);
  if (is.null(x_names))
    x_names = paste0("X",1:D);
  
  set_all = 1:D;
  cur_set = c();
  
  BIC  = c();
  Type = c();
  Var  = list();
  Term = list();
  
  BIC[1]    = calc_BIC(xx, yy, c(), D, K, debug, gam=gam);
  Type[1]   = "Init";
  Var[[1]]  = cur_set;
  Term[[1]] = c();
  
  cur_score = BIC[1];
  
  cat(paste0("Initialization: empty set, BIC = ", sprintf("%.3f", BIC[1]), "\n\n"));
  
  tt = 1;
  ########################
  # Linear Forward Stage #
  ########################
  cat(paste0("Forward Stage - Main effects:\n"));
  while(T) 
  {
    ops = list();
    n_ops = 0;
    
    ######################
    # Forward Operations #
    ######################
    not_set = setdiff(set_all, cur_set);
    Nnset = length(not_set);
    if (Nnset > 0)
    {
      for(j in 1:Nnset)
      {
        jj = not_set[j];
        new_set   = sort(c(jj, cur_set));
        new_score = calc_lda_BIC(xx, yy, new_set, D, K, debug, gam=gam);
        if (debug)
          cat(paste0("  Trying to add variable ", jj , ": ", x_names[jj], " into main effect set...  D_Score: ", cur_score-new_score, "\n\n"));
        if (new_score < cur_score)
        {
          n_ops = n_ops + 1;
          ops[[n_ops]] = list();
          ops[[n_ops]]$new_set   = new_set;
          ops[[n_ops]]$new_score = new_score;
          ops[[n_ops]]$print     = paste0("  Main effects: add variable ", jj , ": ", x_names[jj], " ...  df = ", length(new_set)+1, ",  BIC = ", sprintf("%.3f",new_score));
        }
      }
    }
    
    #######################
    # The Best Operations #
    #######################
    if (n_ops == 0)
    {
      break;
    }
    
    toprint = "";
    for(i in 1:n_ops)
    {
      if (ops[[i]]$new_score < cur_score)
      {
        cur_score = ops[[i]]$new_score;
        cur_set   = ops[[i]]$new_set;
        toprint   = ops[[i]]$print;
      }
    }
    
    tt = tt + 1;
    BIC[tt]    = cur_score;
    Type[[tt]] = "Forward (Main)";
    Var[[tt]]  = cur_set;
    Term[[tt]] = get_lin_terms_vec(cur_set);
    
    cat(paste0(toprint,"\n"));
  }
  
  linear_set = cur_set;
  cur_terms  = get_lin_terms_vec(linear_set);
  cur_set    = c();
  
  ###########################
  # Quadratic Forward Stage #
  ###########################
  cat(paste0("\nForward Stage - Interactions: \n"));
  while(T) 
  {
    ops = list();
    n_ops = 0;
    
    ######################
    # Forward Operations #
    ######################
    not_set = setdiff(set_all, cur_set);
    Nnset = length(not_set);
    if (Nnset > 0)
    {
      for(j in 1:Nnset)
      {
        jj = not_set[j];
        new_set   = sort(c(jj, cur_set));
        new_terms = union(cur_terms, get_quad_terms_vec(new_set));
        
        new_score = calc_BIC(xx, yy, new_terms, D, K, debug, gam=gam);        
        if (debug)
          cat(paste0("  Trying to add variable ", jj , ": ", x_names[jj], " into interaction set...  D_Score: ", cur_score-new_score, "\n"));
        if (new_score < cur_score || length(cur_set) < minF)
        {
          n_ops = n_ops + 1;
          ops[[n_ops]] = list();
          ops[[n_ops]]$new_set   = new_set;
          ops[[n_ops]]$new_score = new_score;
          ops[[n_ops]]$new_terms = new_terms;
          ops[[n_ops]]$print     = paste0("  Interactions: add variable ", jj , ": ", x_names[jj], " ...  df = ", length(new_terms)+1, ",  BIC = ", sprintf("%.3f",new_score));
        }
      }
    }
    
    ######################
    # The Best Operation #
    ######################
    if (n_ops == 0)
    {
      break;
    }
    
    toprint = "";
    
    if (length(cur_set) < minF)
      cur_score = 1e6;
    
    for(i in 1:n_ops)
    {
      if (ops[[i]]$new_score < cur_score)
      {
        cur_score = ops[[i]]$new_score;
        cur_set   = ops[[i]]$new_set;
        toprint   = ops[[i]]$print;
        cur_terms = ops[[i]]$new_terms;
      }
    }
    
    tt = tt + 1;
    BIC[tt]    = cur_score;
    Type[[tt]] = "Forward (Int)";
    Var[[tt]]  = c(setdiff(linear_set, cur_set), cur_set);
    Term[[tt]] = cur_terms;
    
    cat(paste0(toprint,"\n"));
  }
  
  # set of variables at end of forward stage
  cur_set = c(setdiff(linear_set, cur_set), cur_set);
  int_terms = get_inter_terms_vec(cur_set);
  cur_terms = union(cur_terms, get_lin_terms_vec(linear_set))
  cur_score = calc_BIC(xx, yy, cur_terms, D, K, debug, gam=gam);
  
  cat(paste0("\nBackward stage: \n"));
  ##################
  # Backward Stage #
  ##################
  if (length(cur_set) > 0)
  {
    while(T) 
    {
      ops = list();
      n_ops = 0;
      
      #######################
      # Backward Operations #
      #######################
      Nterms = length(cur_terms);
      if(Nterms > 0)
      {
        for(j in 1:Nterms)
        {
          term = cur_terms[j];  
          new_terms = setdiff(cur_terms, term);
          new_score = calc_BIC(xx, yy, new_terms, D, K, debug, gam=gam);
          if (debug)
          {
            term_name = get_term_name_2(x_names, term);
            cat(paste0("  Trying to remove term ", term_name, " ...  Score: ", cur_score - new_score, "\n\n"));
          }
          if (new_score < cur_score)
          {
            n_ops = n_ops + 1;
            term_name = get_term_name_2(x_names, term);
            ops[[n_ops]] = list();
            ops[[n_ops]]$new_terms = new_terms;
            ops[[n_ops]]$new_score = new_score;
            ops[[n_ops]]$print     = paste0("  Remove term ", term_name, " ...  df = ", length(new_terms)+1, ",  BIC = ", sprintf("%.3f", new_score));
          }
        }
      }
      
      #######################
      # The Best Operations #
      #######################
      if (n_ops == 0)
      {
        break;
      }
      
      toprint = "";
      for(i in 1:n_ops)
      {
        if (ops[[i]]$new_score < cur_score)
        {
          cur_score = ops[[i]]$new_score;
          cur_terms = ops[[i]]$new_terms;
          cur_set   = get_set_from_terms(ops[[i]]$new_terms);
          toprint   = ops[[i]]$print;
        }
      }
      
      tt = tt + 1;
      BIC[tt]    = cur_score;
      Type[[tt]] = "Backward";
      Var[[tt]]  = cur_set;
      Term[[tt]] = cur_terms;
      
      cat(paste0(toprint,"\n"));
    }
  }
  
  result = list();
  result$EBIC = BIC;
  result$Type = Type;
  result$Var  = Var;
  result$Term = Term;
  if (tt > 0)
  {
    result$final_EBIC = BIC[tt];
    result$final_Var  = Var[[tt]];
    result$final_Term = Term[[tt]];
  } 
  else 
  {
    result$final_EBIC = BIC[1];
    result$final_Var  = Var[[1]];
    result$final_Term = Term[[1]];
  }
  cat(paste("\nFinal selected variables: ", paste0(x_names[result$final_Var], collapse=", ")));
  term_names = c();
  for(term in result$final_Term)
  {
    term_name = get_term_name_2(x_names, term);
    term_names = c(term_names, term_name);
  }
  cat(paste("\n                   terms: ", paste0(term_names, collapse=", "), "\n"));

  return(result)
}


#
# Calculate a trace of cross-validation error rate for SODA forward-backward procedure
#
soda_trace_CV = function(xx, yy, res_SODA)
{
  if (min(yy) == 0)
    yy = yy + 1;
  
  N_CV = 20;  
  Np = length(res_SODA$Var);
  errors_ss = matrix(0, Np, N_CV);
  ss_V    = numeric(Np);
  ss_MT   = numeric(Np);
  ss_IT   = numeric(Np);
  ss_EBIC = numeric(Np);
  ss_Typ  = character(Np);
  for(i in 1:Np)
  {
    cat(paste0("Calculating CV error for step ", i, " ...\n"));
    SS = res_SODA$Var[[i]];
    TT = res_SODA$Term[[i]];
    for(icv in 1:N_CV)
    {
      errors_ss[i,icv] = 1 - logistic_terms_CV(xx, yy, TT, 10)$m_sp;
    }
    ss_V[i]    = length(SS);
    ss_MT[i]   = length(TT) - sum(grepl("*",TT,fixed=T));
    ss_IT[i]   = sum(grepl("*",TT,fixed=T));
    ss_EBIC[i] = res_SODA$EBIC[i]
    ss_Typ[i]  = res_SODA$Type[i]
  }
  ss_mean = apply(errors_ss, 1, mean);
  tab = data.frame(ss_Typ, ss_EBIC, ss_V, ss_MT, ss_IT, ss_mean);
  colnames(tab) = c("Step Type", "EBIC", "# Variables", "# Main terms", "# Interactions", "CV Error");
  return(tab);
}


logistic_terms_CV = function(xx, yy, terms, KK, Debug=FALSE)
{
  if (min(yy) == 0)
    yy = yy + 1;
  
  N = length(yy);
  K = max(yy);
  D = dim(xx)[2];
  o = sample(1:N);
  
  n = floor(N/KK);
  
  if (is.null(terms))
  {
    xx = matrix(0, N, 0);
  }
  else
  {
    xx = create_pmatrix_from_terms(as.matrix(xx), terms);
  }
  
  xx = as.matrix(xx[o,]);
  yy = yy[o];
  
  m_succ = 0;
  c_succ = 0;
  
  for(kk in 1:KK)
  {
    set_tr = setdiff(1:N,((kk-1)*n+1):((kk)*n));
    set_te = ((kk-1)*n+1):((kk)*n);
    xx_tr = as.matrix(xx[set_tr,]);
    yy_tr = yy[set_tr];
    xx_te = as.matrix(xx[set_te,]);
    yy_te = yy[set_te];
    
    pmatrix = rep(1,length(yy_tr));
    if (length(xx_tr) > 0)
      pmatrix = xx_tr;
    colnames(pmatrix) = NULL;
    fit = multinom(yy_tr ~ pmatrix, family = "multinomial", trace=F);
    cef = coef(fit);
    if (K==2)
      cef = t(as.matrix(coef(fit)));
    zz  = matrix(0, K, length(yy_te))
    for(k in 2:K)
    {
      pmatrix = rep(1,length(yy_te));
      if (length(xx_te) > 0)
        pmatrix = xx_te;
      zz[k,] = cbind(1,pmatrix) %*% cef[k-1,];
    }
    pp  = apply(zz,2,which.max);
    
    m_succ = m_succ + sum(pp == yy_te);
  }
  
  res = list();
  res$m_sp = m_succ/N;
  
  return(res);
}
