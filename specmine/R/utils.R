
# creates names of samples given two lists of "groups" for two "factors" 
"create_names_samples" = function(names1, names2, sep= "_", exclude.list = c())
{
  rep1 = rep(names1, each = length(names2))
  rep2 = rep(names2, length(names1))
  full.res = paste(rep1, rep2, sep=sep)
  full.res[!full.res %in% exclude.list]
}


"rebuild_factors_df" = function(df)
{
  for (i in 1:ncol(df))
    if (is.factor(df[[i]]) )
        df[[i]] = factor(df[[i]])
  df
}
# attention because order of factors is lost ...