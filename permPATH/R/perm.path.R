perm.path <-
function(expr, y, local.test, global.test="wilcoxon", B, gset, min.num=2, max.num, imputeval=NULL, transfun=function(x){x}, sort="pval", anno=NULL)
  {
    if ((local.test=="ttest" | local.test=="wilcoxon" | local.test=="pearson" | local.test=="spearman" | local.test=="jt") & (global.test=="wilcoxon" | global.test=="maxmean" | global.test=="mean" | global.test=="meanabs")){
      # remove samples with missing outcome
      expr = expr[,!is.na(y)]
      y = y[!is.na(y)]

      if (!is.null(imputeval))
        expr[is.na(expr)] = imputeval

      expr = transfun(expr)

      # subset gset
      gset = lapply(gset, function(x){intersect(x,rownames(expr))})
   
      # remove sets with number of probes < min.num and > max.num
      gset = gset[unlist(lapply(gset, function(x){ifelse(length(x)<min.num | length(x)>max.num, FALSE, TRUE)}))]

      # subset expr
      expr = expr[unique(unlist(gset)),]

      probeindex = sapply(unlist(gset), function(x){which(rownames(expr)==x)})

      probeid = rownames(expr)

      K = nrow(expr)
      n = ncol(expr)

      numpath = unlist(lapply(gset, length))
      paths = length(numpath)

      scores_T = rep(0,(B+1)*paths+K)

      if (local.test == "wilcoxon"){
        naind = which(is.na(expr), arr.ind=TRUE)
        rdat = apply(expr, 1, function(x){rank(x, ties.method = "average")})
        rdat[naind] = NA
        scores_T = .C("perm_path_w", as.double(rdat), as.integer(y), as.integer(probeindex), as.integer(numpath), as.integer(paths), as.integer(n), as.integer(K), as.integer(B), as.character(global.test), as.double(scores_T), NAOK=TRUE, PACKAGE="permPATH")[[10]]
      }
      else if (local.test == "ttest")
        scores_T = .C("perm_path_t", as.double(t(expr)), as.integer(y), as.integer(probeindex), as.integer(numpath), as.integer(paths), as.integer(n), as.integer(K), as.integer(B), as.character(global.test), as.double(scores_T), NAOK=TRUE, PACKAGE="permPATH")[[10]]
      else if (local.test == "pearson")
        scores_T = .C("perm_path_pearson", as.double(t(expr)), as.double(y), as.integer(probeindex), as.integer(numpath), as.integer(paths), as.integer(n), as.integer(K), as.integer(B), as.character(global.test), as.double(scores_T), NAOK=TRUE, PACKAGE="permPATH")[[10]]
      else if (local.test == "spearman"){
        naind = which(is.na(expr), arr.ind=TRUE)
        rdat = apply(expr, 1, function(x){rank(x, ties.method = "average")})
        rdat[naind] = NA
        ry = rank(y, ties.method="average")
        scores_T = .C("perm_path_spearman", as.double(rdat), as.double(ry), as.integer(probeindex), as.integer(numpath), as.integer(paths), as.integer(n), as.integer(K), as.integer(B), as.character(global.test), as.double(scores_T), NAOK=TRUE, PACKAGE="permPATH")[[10]]
      }
      else if (local.test == "jt"){
        #naind = which(is.na(expr), arr.ind=TRUE)
        #expr[naind] = NA
        y = y - 1
        cl = length(unique(y))
        m = 0
        v = 0
        scores_T = .C("perm_path_jt", as.double(t(expr)), as.integer(y), as.integer(cl), as.integer(probeindex), as.integer(numpath), as.integer(paths), as.integer(n), as.integer(K), as.double(m), as.double(v), as.integer(B), as.character(global.test), as.double(scores_T), NAOK=TRUE, PACKAGE="permPATH")[[13]]
      }
      else
        print("local test not implemented")

      T = scores_T[((B+1)*paths+1):((B+1)*paths+K)]
      names(T) = rownames(expr)
      scores = matrix(scores_T[1:((B+1)*paths)], B+1, paths, byrow=TRUE)

      scores_max = apply(scores[-1,], 1, function(x){max(abs(x))})

      # compute p-values
      pval = apply(scores, 2, function(x){mean(abs(x[-1])>=abs(x[1]))})
      pfwer = sapply(scores[1,], function(x){mean(abs(x)<=scores_max)})

      # compute FDR
      fdr = p.adjust(pval, "fdr")

      # compute Bonferroni
      bonferroni = p.adjust(pval, "bonferroni")

      # round p-values
      pval = round(pval, 3)
      pfwer = round(pfwer, 3)
      fdr = round(fdr, 3)
      bonferroni = round(bonferroni, 3)
    
      res = data.frame(Score=scores[1,], pval, pfwer, fdr, bonferroni)

      gnames = unlist(lapply(gset, function(x){paste(x, collapse=";")}))[names(gset)]
   
      res = data.frame(Pathway=names(gset), Genes=gnames, Size=numpath, res)

      if (!is.null(anno))        
        res = data.frame(res, anno=unlist(anno[as.character(res[["Pathway"]])]))

      if (sort ==  "pval")
        ind = sort(abs(res$pval), index.return=TRUE)$ix
      if (sort == "score")
        ind = sort(abs(res$Score), index.return=TRUE, decreasing=TRUE)$ix
      res = res[ind,]

      return(list(res=res, stats=T, scores=scores)) 
    }
    else
      print("local and /or global test statistic not implemented")
  }
