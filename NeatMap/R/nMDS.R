nMDS<-function(data,embed.dim=2,n.iters=300,metric="pearson",random.seed=NULL)
{
  data<-as.matrix(data);
  n.points<-dim(data)[1];
  length.profile<-dim(data)[2];
  if(is.null(random.seed))
  {
        random.seed=sample(1:10000,1);
  }
  METRICS=c("pearson","euclidean");
  meth<-pmatch(metric,METRICS);
  if(is.na(meth)) stop("Invalid Method");
  if(meth == -1) stop("Ambiguous Method");
  metr=switch(meth,0,1);
  pos=matrix(.C("nMDS_R",as.integer(n.points),as.integer(length.profile),
	as.double(data),as.integer(embed.dim),as.integer(n.iters),
	positions.out=double(n.points*embed.dim),
	as.integer(random.seed),as.integer(metr),
	NAOK=T)$positions.out,nrow=n.points,ncol=embed.dim);
  rownames(pos)<-rownames(data);
  colnames(pos)<-paste("nMDS_",1:embed.dim,sep="");
  res<-list(x=pos);
  class(res)<-"nMDS";
  res
}
