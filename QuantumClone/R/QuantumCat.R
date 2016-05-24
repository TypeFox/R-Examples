
#' Data generation
#'
#' Creates plausible data as would be oserved by genome sequencing
#' @param number_of_clones The wanted number of observable clones (meaning bearing at least 1 mutation)
#' @param number_of_mutations The total observed number of mutations (across all clones)
#' @param number_of_samples The number of samples on which the data should be simulated. Default is 2.
#' @param depth The depth of sequencing (does not account for contamination). Default is 100x
#' @param ploidy The general ploidy of the tumor. Default is 2. If "disomic" : only AB regions will be generated.
#' @param Random_clones Should the number of clones be generated randomly (sample(1:10))
#' @param contamination A numeric vector indicating the fraction of normal cells in each sample.
#' @param Subclonal.CNA.fraction Cell fraction of the subclone that has subclonal CNA
#' @keywords Data generation phylogeny
#' @export
#' @examples
#' print("Generate small set of mutations from 2 differents clones...")
#' print("...in 1 sample, contaminated at 10% by normal cells")
#'  
#' QuantumCat(number_of_clones=2,number_of_mutations=50,number_of_samples=1,contamination=0.1)
QuantumCat<-function(number_of_clones,number_of_mutations,ploidy=2,depth=100,number_of_samples=2,Random_clones=F,contamination=NULL,
                     Subclonal.CNA.fraction=NULL){
  if(Random_clones){
    number_of_clones<-sample(1:10)
    print(paste("Number of clones",number_of_clones))
  }
  #Tree generation
  if(number_of_clones<=0){
    warning("Invalid number of clones")
    return(NA)
  }
  else{#creates phylogenetic tree of clones
    Tree<-phylo_tree_generation(number_of_clones = number_of_clones,number_of_samples)
  }
  #print(Tree) 
  Cellularities<-cbind(as.matrix(Tree[Tree[,'mutated']==1,3:dim(Tree)[2]]))
  Clonal_attribution<-(rep(1,times = number_of_mutations))
  while(length(unique(Clonal_attribution))<number_of_clones){
    Clonal_attribution<-sample(x=(1:number_of_clones),size=number_of_mutations,replace = T)
  }
  if(is.null(contamination)){
    conta<-rep(0,times = number_of_samples)
  }
  else{
    conta<-contamination
  }
  A<-list()
  B<-list()
  number_of_copies<-list()
  Genotype<-list()
  i<-0
  test<-T
  if(ploidy=="disomic"){
   for(i in 1:number_of_samples){
     A[[i]]<-rep(1,times = number_of_mutations)
     B[[i]]<-rep(1,times = number_of_mutations)
     Genotype[[i]]<-rep("AB",times = number_of_mutations)
     number_of_copies[[i]]<-rep(1,times = number_of_mutations)
   }
  }
  else if(class(ploidy)=="character"){
    if(length(ploidy)==1){
      for(i in 1:number_of_samples){
        A[[i]]<-rep(strcount(x = ploidy,pattern = 'A'),times = number_of_mutations)
        B[[i]]<-rep(strcount(x = ploidy,pattern = 'B'),times = number_of_mutations)
        Genotype[[i]]<-rep(ploidy,times = number_of_mutations)
        number_of_copies[[i]]<-sample(x = 1:(A[[i]][1]),size = number_of_mutations,replace = TRUE)
      }
    }
    else{
      for(i in 1:number_of_samples){
        A[[i]]<-rep(strcount(x = ploidy[i],pattern = 'A'),times = number_of_mutations)
        B[[i]]<-rep(strcount(x = ploidy[i],pattern = 'B'),times = number_of_mutations)
        Genotype[[i]]<-rep(ploidy[i],times = number_of_mutations)
        number_of_copies[[i]]<-sample(x = 1:(A[[i]][1]),size = number_of_mutations,replace = TRUE)
      }
    }
  }
  else{
    for(n in 1:number_of_mutations){
      c<-rpois(number_of_samples,ploidy/2)
      while(sum(c==0)>0){
        for(s in which(c==0)){
          c[s]<-rpois(1,ploidy/2)
        }
      }
      for(i in 1:number_of_samples){
        if(length(A)<number_of_samples){
          A[[i]]<-c[i]
          B[[i]]<-sample(0:c[i],size=1)
          Genotype[[i]]<-paste(paste(rep('A',times=A[[i]]),collapse=''),paste(rep('B',times=B[[i]]),collapse=''),sep='')
        }
        else{
          A[[i]]<-c(A[[i]],c[i])
          B[[i]]<-c(B[[i]],sample(0:c[i],size=1))
          Genotype[[i]]<-c(Genotype[[i]],paste(paste(rep('A',times=tail(A[[i]],1)),collapse=''),paste(rep('B',times=tail(B[[i]],1)),collapse=""),sep=''))
        }
      }
    }
    
    for(j in 1:number_of_samples){
      number_of_copies[[j]]<-sample(x=1:A[[j]][1],size=1)
      
      for(i in 2:number_of_mutations){
        number_of_copies[[j]]<-c(number_of_copies[[j]],sample(x=1:A[[j]][i],size=1))
      }
    }
  }
  if(number_of_samples>1){
  Cell<-Cellularities[Clonal_attribution[1:number_of_mutations],]
  }
  else{
    Cell<-Cellularities[Clonal_attribution[1:number_of_mutations]]
  }
  recap<-list()
  if(number_of_samples>1){
    for(i in 1:number_of_samples){
      frequency<-Cell[,i]*number_of_copies[[i]]/(A[[i]]+B[[i]])*(1-conta[i])
      SNP_depth<-rnbinom(n=number_of_mutations,size = 4.331601,mu = depth)
      Alt_depth<-numeric()
      for(j in 1:number_of_mutations){
        Alt_depth[j]<-rbinom(n = 1,size=SNP_depth[j],prob= frequency[j]/100)
      }
      recap[[i]]<-data.frame(Clonal_attribution,1:length(Clonal_attribution),Genotype[[i]],Cell[,i],number_of_copies[[i]],frequency,SNP_depth,Alt_depth)
      colnames(recap[[i]])<-c("Chr","Start",'Genotype','Cellularity','number_of_copies',"Frequency",'Depth',"Alt")
    }
  }
  else if(number_of_samples==1){
    frequency<-Cell*number_of_copies[[1]]/(A[[1]]+B[[1]])*(1-conta)
    SNP_depth<-rnbinom(n=number_of_mutations,size = 4.331601,mu = depth)
    Alt_depth<-numeric()
    for(j in 1:number_of_mutations){
      Alt_depth[j]<-rbinom(n = 1,size=SNP_depth[j],prob= frequency[j]/100)
    }
    Observed_freq<-Alt_depth/SNP_depth*100
    recap[[1]]<-data.frame(Clonal_attribution,1:length(Clonal_attribution),Genotype[[1]],Cell,number_of_copies[[1]],frequency,SNP_depth,Alt_depth)
    colnames(recap[[1]])<-c("Chr","Start",'Genotype','Cellularity','number_of_copies',"Frequency",'Depth',"Alt")
  }
  return(recap)
}

#' Data generation
#'
#' WARNING: Tree_generation recreates a tree from data while phylo_tree_generation randomly creates a phylogeny
#' @param number_of_clones The wanted number of observable clones (meaning bearing at least 1 mutation)
#' @param number_of_samples The number of samples on which the data should be simulated
#' @keywords Data creation phylogeny
phylo_tree_generation<-function(number_of_clones,number_of_samples){
  Localization<-as.character('Ancestral')
  Proportion<-list()
  for(i in 1:number_of_samples){
    Proportion[[i]]<-100
  }
  j<-1
  mutated<-1
  has_progeny<-as.numeric(0)
  while(sum(mutated)<number_of_clones){
    j=j+1
    if(length(which(has_progeny==0))==1){
      selection<-which(has_progeny==0)
    }
    else{
      selection<-sample(which(has_progeny==0),size = 1)
    }
    has_progeny[selection]<-1
    
    if(j==2){
      Localization<-c(Localization,'L','R')
      
      Lfreq<-sample(x = 1:99,size=number_of_samples,replace = T)
      for(i in 1:number_of_samples){
        Proportion[[i]]<-c(Proportion[[i]],Lfreq[i],100-Lfreq[i])
      }
      mutated<-c(mutated,1,0)
    }else{
      for(i in 1:number_of_samples){
        Lfreq[i]<-sample(x=1:Proportion[[i]][selection],size = 1)
        Proportion[[i]]<-c(Proportion[[i]],Lfreq[i],Proportion[[i]][selection]-Lfreq[i])
        
      }
      Localization<-c(Localization,paste(Localization[selection],'L',sep=""),paste(Localization[selection],'R',sep=""))
      mutated<-c(mutated,1,0)
    }
    has_progeny<-c(has_progeny,rep(0,times=2))
  }
  points<-data.frame(mutated,Localization)
  for(i in 1:number_of_samples){
    points<-cbind(points,Proportion[[i]])
  }
  colnames(points)<-c('mutated','Localization',1:number_of_samples)
  return(points)
}
#' Multiple testing
#'
#' Generates multiple data with the QuantumCat function, and assesses the clustering algorithm (output)
#' @param number_of_tests The number of times the QuantumCat function should be called
#' @param number_of_clones The wanted number of observable clones (meaning bearing at least 1 mutation)
#' @param number_of_mutations The total number of observed number of mutations (across all clones)
#' @param number_of_samples The number of samples on which the data should be simulated. Default is 2.
#' @param depth The depth of sequencing (does not account for contamination). Default is 100x
#' @param ploidy The general ploidy of the tumor. Default is 2.
#' @param Sample_name Name of the folder which will be created to save the data
#' @param plot_results Should the intermediate results be plotted and saved? Default is TRUE
#' @param save_Cell Should the dataframe with the cellularities be saved? Default is FALSE
#' @param preclustering Should the data be pre-clustered by k-means?
#' @param contamination A numeric vector indicating the fraction of normal cells in each sample.
#' @param nclone_range Range of values tested for the maximum likelihood
#' @param maxit Number of trials different initial conditions to be used to evaluate the mximum likelihood
#' @param epsilon The stop value of the EM algorithm
#' @param ncores Number of cores to be used for the optimization research
#' @keywords Data generation phylogeny
#' @export
#' @examples
#' \dontrun{
#' print("Tests multiple conditions of data generation...")
#' print("... and gives information on the clustering quality")
#' Multitest(number_of_tests=1,Sample_name="Example",ploidy="AB",
#' number_of_clones=2,plot_results=FALSE,depth=500,number_of_mutations=10)
#' }
Multitest<-function(number_of_tests,number_of_samples=2,ploidy=2,Sample_name='Multitest',depth=100,number_of_mutations=150,
                    number_of_clones=NULL,plot_results=T,save_Cell=F,preclustering = T ,
                    contamination = NULL,nclone_range = 2:5 ,maxit = 1, epsilon = 5*(10**(-3)), ncores = 1){
  set.seed(1)
  if(is.null(contamination)){
    contamination<-rep(0,times=number_of_samples)
  }
  if(is.null(number_of_clones)){
    Random_clones<-T
  }
  else{
    Random_clones<-F
  }
  for(i in 1:number_of_tests){
    if(Random_clones==T){
      number_of_clones<-sample(2:5,size = 1)
    }
    if(plot_results | save_Cell){
      dir.create(path = paste(Sample_name,i,sep=''), showWarnings = FALSE)
    }
    data_generated<-QuantumCat(number_of_clones = number_of_clones,number_of_mutations = number_of_mutations,contamination = contamination,
                               depth = depth,ploidy = ploidy,number_of_samples = number_of_samples,Random_clones=Random_clones)
    
    tempo_col<-colnames(data_generated[[1]])
    data_generated[[1]]<-cbind(paste(Sample_name,i,sep=''),data_generated[[1]])
    colnames(data_generated[[1]])<-c('Samples',tempo_col)
    Cluster<-as.character(data_generated[[1]][,'Chr'])
    if(plot_results){
      if(number_of_samples>1){
        U<-expand.grid(1:length(data_generated),1:length(data_generated))
        U<-U[U[,1]<U[,2],]
        for(p in 1:dim(U)[1]){
          q<-ggplot2::qplot(x=data_generated[[U[p,1]]][,'Frequency'],y=data_generated[[U[p,2]]][,'Frequency'],asp = 1,main=paste('Freq plot multitest',i),
                            xlim=c(0,100),ylim=c(0,100),
                            xlab=paste('Exact freq',U[p,1],sep=''),ylab=paste('Exact freq',U[p,2],sep=''),colour=Cluster)
          ggplot2::ggsave(plot = q, filename = paste(Sample_name,i,'/', 'Exact_freq',U[p,1],'_',U[p,2],'.png',sep=''),width = 6.04,height = 6.04)
        
        }
      }
    }
    t<-One_step_clustering(SNV_list = data_generated,simulated = T,save_plot = plot_results,epsilon = epsilon,ncores = ncores,
                           preclustering = ,contamination = contamination,nclone_range = nclone_range,maxit = maxit)
    percentage_misclustered<-0
    if(number_of_samples>1){
      #combining results
      Cell<-data.frame(data_generated[[1]][,'Cellularity'])
      for(j in 2:number_of_samples){
        Cell<-cbind(Cell,data_generated[[j]][,'Cellularity'])
      }
      Cell<-unique(Cell[order(data_generated[[1]][,'Chr']),])/100
    }
    else{
      Cell<-data.frame(data_generated[[1]][,'Cellularity'])
    }
    if(i==1){
      statistics<-data.frame(number_of_clones,max(t$cluster))
      P_cluster<-table(t$cluster)/length(t$cluster)
      P_clone<-table(t$filtered.data[[1]]$Chr)/length(t$filtered.data[[1]]$Chr)
      H_clone<--sum(P_clone*log(P_clone))
      H_cluster<--sum(P_cluster*log(P_cluster))
      A<-aggregate(rep(1, times = length(t$cluster)), by = list(x= t$cluster,y=t$filtered.data[[1]]$Chr ), sum)
      L<-log(A[,3]/(length(t$cluster)*P_cluster[A[,1]]*P_clone[A[,2]]))
      NMI<-2*sum(A[,3]/length(t$cluster)*L)/(H_clone+H_cluster)      
#       for(j in 1:number_of_clones){
#         cluster<-t$cluster[t$filtered.data[[1]]$Chr==j]
#         corresponding_cluster<-as.numeric(names(table(cluster))[which.max(table(cluster))])
#         #percentage_misclustered<-percentage_misclustered+sum(cluster!=corresponding_cluster)
#         if(length(corresponding_cluster)==0){
#           warning(paste("One unobserved clone:",cluster))
#         }
#         else{
#           if(j==1){
#             distance<-sqrt(sum((t$pamobject$medoids[corresponding_cluster,]-Cell[j,])**2))
#           }
#           else{
#             distance<-c(sqrt(sum((t$pamobject$medoids[corresponding_cluster,]-Cell[j,])**2)))
#           }
#         }
#       }
      
      #percentage_misclustered<-percentage_misclustered/length(t$pamobject$clustering)*100
#       average_distance<-mean(distance)
      statistics<-cbind(statistics,NMI)
      colnames(statistics)<-colnames(statistics)<-c("number_of_clones",'clones_observed','NMI')
    }
    else{
      P_cluster<-table(t$cluster)/length(t$cluster)
      P_clone<-table(t$filtered.data[[1]]$Chr)/length(t$filtered.data[[1]]$Chr)
      H_clone<--sum(P_clone*log(P_clone))
      H_cluster<--sum(P_cluster*log(P_cluster))
      A<-aggregate(rep(1, length(t$cluster)), by = list(x= t$cluster,y=t$filtered.data[[1]]$Chr ), sum)
      L<-log(A[,3]/(length(t$cluster)*P_cluster[A[,1]]*P_clone[A[,2]]))
      NMI<-2*sum(A[,3]/length(t$cluster)*L)/(H_clone+H_cluster)
#       for(j in 1:number_of_clones){
#         #percentage_misclustered<-percentage_misclustered+sum(cluster!=corresponding_cluster)
#         cluster<-t$pamobject$clustering[t$filtered.data[[1]]$Chr]==j]
#         corresponding_cluster<-as.numeric(names(table(cluster))[which.max(table(cluster))])
#         if(length(corresponding_cluster)==0){
#           warning(paste("One unobserved clone:",cluster))
#         }
#         else{
#           if(j==1){
#             distance<-sqrt(sum((t$pamobject$medoids[corresponding_cluster,]-Cell[j,])**2))
#           }
#           else{
#             distance<-c(sqrt(sum((t$pamobject$medoids[corresponding_cluster,]-Cell[j,])**2)))
#           }
#         }
#       }
      #percentage_misclustered<-percentage_misclustered/length(t$pamobject$clustering)*100
      #average_distance<-mean(distance)
      statistics<-rbind(statistics,c(number_of_clones,max(t$cluster),NMI))
    }
#     if(plot_results){
#       if(number_of_samples>1){
#         for(p in 1:dim(U)[1]){
#           q<-ggplot2::qplot(x=c(t$pamobject$medoids[,U[p,1]],Cell[,U[p,1]]),y=c(t$pamobject$medoids[,U[p,2]],Cell[,U[p,2]]),
#                             xlab=paste('Cellularity of centers, sample',U[p,1]),
#                             ylab=paste('Cellularity of centers, sample',U[p,2]),asp=1,main='Comparison of centers',
#                             colour=factor(c(rep(1,times=dim(t$pamobject$medoids)[1]),rep(2,times=dim(Cell)[1]))),
#                             shape=factor(c(rep(1,times=dim(t$pamobject$medoids)[1]),rep(2,times=dim(Cell)[1]))))+scale_shape_discrete(name="Centers",breaks=c(1,2),labels=c("Estimated","Real"))+scale_colour_discrete(name="Centers",breaks=c(1,2),labels=c("Estimated","Real"))+coord_cartesian(xlim=c(0,1),ylim=c(0,1))
#           ggplot2::ggsave(plot = q, filename = paste(Sample_name,i,'/', 'Comparison_simulation', U[p,1],'_',U[p,2], '.png',sep=''),width = 6.04,height = 6.04)
#         }
#       }
#       else{
#         q<-ggplot2::qplot(x=c(t$pamobject$medoids[,1],unique(Cell[,1])/100),y=0.5,
#                           xlab=paste('Cellularity of centers'),
#                           ylab='',asp=1,main='Comparison of centers',
#                           colour=factor(c(rep(1,times=dim(t$pamobject$medoids)[1]),rep(2,times=dim(unique(Cell))[1]))),
#                           shape=factor(c(rep(1,times=dim(t$pamobject$medoids)[1]),rep(2,times=dim(unique(Cell))[1]))))+scale_shape_discrete(name="Centers",breaks=c(1,2),labels=c("Estimated","Real"))+scale_colour_discrete(name="Centers",breaks=c(1,2),labels=c("Estimated","Real"))+coord_cartesian(xlim=c(0,1),ylim=c(0,1))
#         ggplot2::ggsave(plot = q, filename = paste(Sample_name,i,'/', 'Comparison_simulation', '.png',sep=''),width = 6.04,height = 6.04)
#         
#       }
     
  }
  return(statistics)
}
#' Multiple testing
#'
#' Call the Multitest function on different parameters to evaluate the influence of each according to the model. Does not save plot results.
#' @param number_of_tests_per_condition The number of times the QuantumCat function should be called
#' @param range_clones The range of observable clones (meaning bearing at least 1 mutation)
#' @param range_mutations The total number of observed number of mutations (across all clones)
#' @param range_samples The number of samples on which the data should be simulated. Default is 2.
#' @param range_depth The depth of sequencing (does not account for contamination). Default is 100x
#' @param range_ploidy The general ploidy of the tumor. Default is 2.
#' @param folder Name of the folder in which the statistics will be saved
#' @param plot_results Should the intermediate results be plotted and saved? Default is TRUE
#' @param preclustering Should the data be pre-clustered by k-means?
#' @param save_Cell Should the dataframe with the cellularities be saved? Default is FALSE
#' @param range_contamination A numeric vector indicating the fraction of normal cells in samples.
#' @param nclone_range Range of values tested for the maximum likelihood
#' @param maxit Number of trials different initial conditions to be used to evaluate the mximum likelihood
#' @param epsilon The stop value of the EM algorithm
#' @param ncores Number of cores to be used for the optimization research
#' @keywords Data generation phylogeny
#' @export
#' @examples 
#' \dontrun{
#' statistics_on_Multitest(number_of_tests_per_condition=2,range_clones=2:3,range_mutations=c(20),
#'                                  range_ploidy="AB",range_samples=2,range_depth=1000,
#'                                  range_contamination = 0,folder='stat',
#'                                  plot_results=FALSE,save_Cell=FALSE,preclustering = TRUE ,
#'                                  nclone_range = 2 ,maxit = 1, epsilon = 5*(10**(-3)),ncores = 1)
#'}

statistics_on_Multitest<-function(number_of_tests_per_condition=100,range_clones=2:5,range_mutations=c(50,80,100,250,500),
                                  range_ploidy=c(2,3),range_samples=c(2,3,4,5,7,10),range_depth=c(50,80,100,250,500,50000),range_contamination = c(0,0.25,0.5,0.75),folder='stat',
                                  plot_results=T,save_Cell=F,preclustering = T ,
                                  nclone_range = 2:5 ,maxit = 4, epsilon = 5*(10**(-3)),ncores = 4){
  dir.create(path = folder, showWarnings = FALSE)
  Test_code<-NULL
  for(i in range_clones){
    for(j in range_mutations){
      for(k in range_ploidy){
        for(l in range_samples){
          for(m in range_depth){
            for(cont in range_contamination){
              for(eps in epsilon){
                contamination<-rep(x = cont,times=l)
                statistics<-Multitest(number_of_tests = number_of_tests_per_condition,preclustering = preclustering,epsilon = eps,
                                      contamination = contamination,nclone_range = nclone_range,maxit = maxit,
                                      number_of_samples = l,ploidy = k,Sample_name = "Test",number_of_mutations = j,
                                      number_of_clones = i,plot_results = plot_results,save_Cell = F,depth = m, ncores = ncores )
                write.table(x = statistics,file = paste(folder,'/C',i,'M',j,'P',k,'S',l,'D',m,"eps",eps,'.txt',sep=''))
                if(is.null(Test_code)){
                  Test_code<-paste('C',i,'M',j,'P',k,'S',l,'D',m,"Co",floor(cont*100),"eps",eps,sep="")
                  mean_clone_found<-mean(statistics[,2],na.rm = T)
                  sd_clone_found<-sd(statistics[,2],na.rm = T)
                  mean_percentage_misclustered<-mean(statistics[,3],na.rm = T)
                  sd_percentage_misclustered<-sd(statistics[,3],na.rm = T)
                  #               mean_distance<-mean(statistics[,4])
                  #              sd_distance<-sd(statistics[,4])
                  obs<-sum((!is.na(statistics[,2]) & !is.na(statistics[,3])))
                }
                else{
                  mean_clone_found<-c(mean_clone_found,mean(statistics[,2],na.rm = T))
                  sd_clone_found<-c(sd_clone_found,sd(statistics[,2],na.rm = T))
                  mean_percentage_misclustered<-c(mean_percentage_misclustered,mean(statistics[,3],na.rm = T))
                  sd_percentage_misclustered<-c(sd_percentage_misclustered,sd(statistics[,3],na.rm = T))
                  #               mean_distance<-c(mean_distance,mean(statistics[,4]))
                  #              sd_distance<-c(sd_distance,sd(statistics[,4]))
                  obs<-c(obs,sum((!is.na(statistics[,2]) & !is.na(statistics[,3])))
                  )
               }
              }
            }
          }
        }
      }  
    }
  }
  result<-data.frame(Test_code,mean_clone_found,sd_clone_found,mean_percentage_misclustered,sd_percentage_misclustered,obs)
  colnames(result)<-c("Test_code",'mean_clone_found','sd_clone_found','mean_NMI','sd_NMI',"number_of_obs")
  return(result)
}
