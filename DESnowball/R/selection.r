fs.selection <- function(full.list,
			 cutoff.p=0.05,
			 p.adjust.method="BH"
			 )
  ## select the significant features based on p values
{
    selected.list <- subset(full.list,
			    subset=p.adjust(full.list$pval,
					    method=p.adjust.method) < cutoff.p)
    selected.list
}
