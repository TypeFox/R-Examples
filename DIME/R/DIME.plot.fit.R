DIME.plot.fit <-
function(data, obj, ...)
{
obj <- obj$best;
if (obj$name =="NUDGE"){
  nudge.plot.fit(data, obj, ...);
  }else if (obj$name=="iNUDGE"){
  inudge.plot.fit(data, obj, ...);
  }else{
  gng.plot.fit(data, obj, ...);
  }
}

