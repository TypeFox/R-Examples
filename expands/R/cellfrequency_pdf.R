cellfrequency_pdf <-function(af,cnv,pnb,freq, max_PM=6, snv_cnv_flag=3, SP_cnv=NA, PM_cnv=NA){
  .jaddClassPath("ExPANdS.jar")
  .jinit(classpath="ExPANdS.jar")
  javaImport(packages = "core.analysis.ngs.algorithms.*")
  ##Compute all possible solutions for variable PM_cnv, PM_B and e
  if(snv_cnv_flag==4){
    expands <-.jnew("ExPANdS", as.double(af),as.double(cnv),as.double(SP_cnv),as.integer(PM_cnv),as.integer(pnb));
  }else{
    expands <-.jnew("ExPANdS", as.double(af),as.double(cnv),as.integer(pnb),as.integer(max_PM));
  }
  .jcall(expands,,"run",as.integer(snv_cnv_flag))
  bestF<-.jcall(expands,"D","getF");
  fits<-.jcall(expands,"Ljava/util/Collection;","solutions");
  fits<-.jcall(fits,"Ljava/util/Iterator;","iterator");
  results=c();
  while (.jcall(fits,"Z","hasNext")){
    solutions=.jcall(fits,"Ljava/lang/Object;","next")
    solutions=.jcall(solutions,"Ljava/util/Iterator;","iterator");
    while (.jcall(solutions,"Z","hasNext")){
      rawR<-.jcall(solutions,"Ljava/lang/Object;","next")
      rawR<-.jcall(rawR,"[D","toDouble");
      results=rbind(results,rawR);
    }
  }
  ##Keep only valid solutions
  fit=matrix(results, nrow = nrow(results), ncol = ncol(results), 
             dimnames = list(paste(1:nrow(results)), .jfield(expands,,"SOLUTION_ENTITIES")))
  fit=fit[fit[,"f"]>=freq[1] & fit[,"f"]<=freq[length(freq)],];
  
  ##Normalize deviation by proximity of SP size solution to other SP sizes (i.e. centrality of SP size)
  kstest=ks.test(freq,"punif",freq[1],freq[length(freq)]); 
  ##If freq is uniformely distributed than this function has been called 
  ##for clustering purposes and so probabilities need to be adjusted to 
  ##take into account SP size centrality. This is not necessary for SNV assignment.
  if (kstest$p.value==1){
    fit=.addColumn(fit,"corrFactor",NA);
    dist=c(); 
    for (j in 1:length( freq)){ 
      dist[j]=1/mean(abs(freq[j]-freq))
    }
    dist=dist/max(dist)
    for (i in 1:nrow(fit)){
      ij=which.min(abs(freq-fit[i,"f"]));
      fit[i,"corrFactor"]=dist[ij]*dist[ij];
    }
    fit[,"dev"]=fit[,"dev"]*fit[,"corrFactor"]  
  }
  
  z=round(fit[,"f"]*100);
  z1=sort(unique(z));
  dm=matrix(nrow=length(z1), ncol=ncol(fit), 
            dimnames=list(paste(1:length(z1)),colnames(fit)));
  tfit=t(fit);
  for (i in 1:length(z1)){
    f=z1[i];
    similarFrequencies=t(tfit[,z==f]);
    ia=which.min(similarFrequencies[,"dev"]);
    dm[i,]=similarFrequencies[ia,];
  }
  
  ##create frequency array weghted by deviation
  normdev=dm[,"dev"];
  normdev=round(-100*.sigmoid(normdev*50, 2, 4)+101);
  # normdev=round(-1*log10(normdev/max(normdev)))+1;
  #normdev=1+round(100*(1-normdev));
  # normdev=round(100/double(dm(:,"dev")));
  f=matrix(NA,sum(normdev),1);
  for (i in 1:nrow(dm)){
    if (i==1){
      start = 1;
    }else{
      start=sum(normdev[1:i-1])+1;
    }
    idx=start:sum(normdev[1:i]);
    f[idx]=dm[i,"f"];
  }
  
  nComponents=ceil(fit[1,"CN_Estimate"])+1; ##number of components in gaussian mixture model
  p=matrix(NA,1,length(freq));
  
  if(length(unique(f))>1){
    obj=suppressWarnings(densityMclust(f[!is.na(f)],G=max(1,nComponents-1):nComponents));
    stowarn<-warnings();
    for (w in names(stowarn)){
      if ( isempty(agrep("occurs at min or max",w,max.distance=0.3)) ){
        warning(as.character(stowarn[[w]]),": ",w)
      }
    }
    p=predict(obj, t(freq));  
    #p=p/sum(p); comparability across mutations maintained if implemented? Should behaviour be distinct for clustering vs. mutation assignment?
  }else{
    p=dnorm(freq, mean = unique(f), sd = 1/length(f));
    p=p/sum(p);
    p=p*exp((length(f)-1)/10)/exp(10);#from relative (within this SNV's solution) to absolute probabilities (across multiple SNV solutions)
    # p=p*1/(102-length(f))^10
    
#     p=matrix(0,1,length(freq));
#     p[which.min(abs(unique(f)-freq))]=1;
  }
  output=list("p"=p,"bestF"=bestF,"fit"=fit,"errors"=NULL);
  return(output)
}

.sigmoid<-function(x, t1, t2){
  res = 1/(1 + exp(-t1*(x-t2)));
  return(res);
}

