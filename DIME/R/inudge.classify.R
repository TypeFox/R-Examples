inudge.classify <-
function(data, obj, obj.cutoff = 0.1, obj.sigma.diff.cutoff = NULL,
  obj.mu.diff.cutoff=NULL)
{
obj <- DIME.classify(data,obj,obj.cutoff,obj.sigma.diff.cutoff,
  obj.mu.diff.cutoff);
return(obj);
}