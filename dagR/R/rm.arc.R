rm.arc <-
function(dag, arc)
{ # this conveniently removes an arc;
  # it removes the evaluated paths etc;
  dag$arc<-dag$arc[-arc,];
  dag$arc.type<-dag$arc.type[-arc];
  dag$curve.x<-dag$curve.x[-arc];
  dag$curve.y<-dag$curve.y[-arc];
  dag$pathsN<-NULL;
  dag$paths<-NULL;
  dag$path.status<-NULL;
  dag$searchType <- NULL;
  dag$searchRes <- NULL;
  return(dag);
}

