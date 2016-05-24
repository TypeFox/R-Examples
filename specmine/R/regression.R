"linregression_onevar" = function(dataset, x.val, metadata.vars, combination) {
  values = get_data_values(dataset, x.val)
  sub.df = dataset$metadata[,metadata.vars]
  sub.df = cbind(values,sub.df)
  terms = names(sub.df)[2:ncol(sub.df)]
  terms = cbind(terms, combination)
  reslm = lm(reformulate(terms, "values"), data = sub.df)
  lm.summary = summary(reslm)
  lm.summary
}

"linreg_all_vars" = function(dataset, metadata.vars, combination)
{
  m = vector("list",nrow(dataset$data))
  for(i in 1:nrow(dataset$data))
    m[[i]] = linregression_onevar(dataset, rownames(dataset$data)[i], metadata.vars, combination)
  names(m) = rownames(dataset$data)  
  m
}

linreg_coef_table = function(linreg.results, write.file = F, file.out = "linreg-coefs.csv"){
  num_vars = dim(linreg.results[[1]]$coefficients)[1]
  m = matrix(NA, length(linreg.results), num_vars)
  rownames(m) = names(linreg.results)
  for (i in 1:length(linreg.results)){
    m[i,] = linreg.results[[i]]$coefficients[,1]
  }
  coef.table = as.data.frame(m)
  colnames(coef.table) = trim(rownames(linreg.results[[1]]$coefficients))
  if (write.file) write.csv(coef.table, file = file.out)
  coef.table  
}

linreg_pvalue_table = function(linreg.results, write.file = F, file.out = "linreg-pvalues.csv"){
  num_vars = dim(linreg.results[[1]]$coefficients)[1]
  m = matrix(NA, length(linreg.results), num_vars)
  rownames(m) = names(linreg.results)
  for (i in 1:length(linreg.results)){
    m[i,] = linreg.results[[i]]$coefficients[,4]
  }
  pv.table = as.data.frame(m)
  colnames(pv.table) = trim(rownames(linreg.results[[1]]$coefficients))
  if (write.file) write.csv(pv.table, file = file.out)
  pv.table  
}

linreg_rsquared = function(linreg.results, write.file = F, file.out = "linreg-rsquared.csv"){
  m = matrix(NA, length(linreg.results), 2)
  rownames(m) = names(linreg.results)
  for (i in 1:length(linreg.results)){
    m[i,1] = linreg.results[[i]]$r.sq
    m[i,2] = linreg.results[[i]]$adj.r.sq
  }
  rsq.table = as.data.frame(m)
  colnames(rsq.table) = c("r.squared", "adj.r.squared")
  if (write.file) write.csv(rsq.table, file = file.out)
  rsq.table  
}

plot_regression_coefs_pvalues = function(linreg.results, bar.col = NULL, coef.size = 5, ...){
	coefs = linreg_coef_table(linreg.results)
	pvalues = linreg_pvalue_table(linreg.results)
	df.coefs = data.frame(t(coefs))
	df.coefs = data.frame(round(df.coefs,2))
	df.pvalues = data.frame(t(pvalues))
	df.pvalues = data.frame(-log10(df.pvalues))
	df.pvalues.coefs = data.frame(cbind(df.pvalues, df.coefs))
	num.variables = length(rownames(pvalues))
	colnames(df.pvalues.coefs) = paste("X",1:(num.variables*2), sep="")
	df.pvalues.coefs$names = colnames(pvalues)
	plots = list()
	if (is.null(bar.col)){
		bar.col = rep("steelblue", num.variables)
	}
	
	for (count in 1:num.variables){
			p = ggplot2::ggplot(data = df.pvalues.coefs, ggplot2::aes_string(x="names", y=paste("X",count,sep="") )) + 
            ggplot2::geom_bar(stat="identity", position = ggplot2::position_dodge(), fill = bar.col[count], ...) + 
            ggplot2::geom_text(ggplot2::aes_string(label=paste("X",count+num.variables,sep="")), vjust=-0.3, size=coef.size, ...) + 
            ggplot2::ylab("-log10(pvalue)") + ggplot2::xlab(rownames(pvalues)[count])
		plots[[count]] = p
	}
	
	num.cols = num.variables %/% 2
	
	multiplot(plots, cols = num.cols)
}