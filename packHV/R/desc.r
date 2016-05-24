#' Making descriptive statistics
#'
#' Makes descriptive statistics of a data frame according to a group covariate or not, can export the results
#'
#' @param data data frame to describe in which we can find \code{vars} and \code{group}
#' @param vars vector of character strings of the covariates to describe
#' @param group character string, statistics created for each levels of this covariate
#' @param whole boolean, \code{TRUE} to add a column with the whole statistics when comparing groups (set to \code{FALSE} if \code{group=NULL})
#' @param vars.labels vector of character string for sweeter names of covariates in the output
#' @param group.labels vector of character string for sweeter column names
#' @param type.quanti character string, \code{"med"} to compute median [Q1;Q3], \code{"mean"} to compute mean (sd), \code{"mean_med"} to compute both mean (sd) and median [Q1;Q3] or \code{"med_mm"}, \code{"mean_mm"} or \code{"mean_med_mm"} to add (min;max)
#' @param test.quanti character string, \code{"param"} to compute parametric tests for quantitative covariates (t-test or ANOVA) or \code{"nonparam"} for non parametric tests (Wilcoxon test or Kruskal-Wallis test)
#' @param test boolean, \code{TRUE} to perform tests (\code{FALSE} if \code{group} is \code{NULL}): Khi-2 or Fisher exact test for categorical covariates, t-test/ANOVA or Wilcoxon/Kruskal-Wallis Rank Sum Test for numerical covariates
#' @param noquote boolean, \code{TRUE} to hide quotes when printing the table
#' @param justify boolean, \code{TRUE} to justify columns on right or left (\code{FALSE} if export)
#' @param digits number of digits of the statistics (mean, sd, median, min, max, Q1, Q3, \%), p-values always have 3 digits
#' @param file.export character string, name of the XLS file exported
#' @param language character string, \code{"french"} or \code{"english"}
#' @return A matrix of the descriptive statistics
#' @author Hugo Varet
#' @examples
#' cgd$steroids=factor(cgd$steroids)
#' cgd$status=factor(cgd$status)
#' desc(cgd,vars=c("center","sex","age","height","weight","steroids","status"),group="treat")

# last updated: august 12, 2013

desc=function(data,vars,group=NULL,whole=TRUE,vars.labels=vars,group.labels=NULL,type.quanti="mean",
                  test.quanti="param",test=TRUE,noquote=TRUE,justify=TRUE,digits=2,file.export=NULL,
                  language="english"){
                  
  name_data=deparse(substitute(data))
  
  ############################### checkings ####################################
  if (!is.null(file.export)){justify=FALSE}
  if (!is.data.frame(data)){stop(paste(name_data," is not a data frame\n",sep=""))}
  if (test & !I(test.quanti %in% c("param","nonparam"))){stop("argument test.quanti must be equal to \"param\" or \"nonparam\"")}
  if (language=="french"){
    whole_pop="Echantillon entier"; pvalue="p-valeur"; parameter="Variable";
  } else{
    whole_pop="Whole sample"; pvalue="p-value"; parameter="Covariate";
  }
  if (is.null(group)){
    test=FALSE; whole=FALSE; groupname=""; sep=""; data.g=data;
    group=factor(rep(whole_pop,nrow(data))); group.g=factor(rep(whole_pop,nrow(data)));
  } else{
    if (!I(group %in% names(data))){stop(paste("The variable of interest ",group," is not in ",name_data,sep=""))}
    groupname=group
    if (any(is.na(data[,groupname]))){cat(paste(sum(is.na(data[,groupname]))," missing values in ",groupname,"\n",sep=""))}  
    sep=" "; data.g=data[!is.na(data[,groupname]),]; 
    group=factor(data[,groupname]); group.g=factor(data.g[,groupname]);
  }
  if (length(vars)!=length(vars.labels)){stop("vars and vars.labels do not have the same length")}
  if (any(!I(vars %in% names(data)))){
    cat(paste("Variable(s) ",paste(vars[!I(vars %in% names(data))],collapse=", ")," not in ",name_data,"\n",sep=""))
    vars.labels=vars.labels[vars %in% names(data)]
    vars=vars[vars %in% names(data)]
  }
  # deleting covariates != numeric, character, logical or factor
  delete=vars[!I(apply(as.data.frame(data[,vars]),2,class) %in% c("integer","numeric","factor","character","logical"))]
  if (length(delete)>0){
    cat(paste("Variable(s) ",paste(delete,collapse=", ")," not numeric, neither factor, neither character, neither logical\n",sep=""))
  }
  vars.labels=vars.labels[!I(vars %in% delete)]; vars=vars[!I(vars %in% delete)];
    
  ########################## functions fun.quanti ##############################
  if (!I(type.quanti %in% c("mean","med","mean_mm","med_mm","mean_med","mean_med_mm"))){
    stop("type.quanti must be equal to \"mean\", \"med\", \"mean_mm\", \"med_mm\", \"mean_med\" or \"mean_med_mm\"")
  }
  type.test.quanti=ifelse(test,ifelse(test.quanti=="nonparam"," and Wilcoxon/Kruskal-Wallis test"," and t-test/anova"),"")
  khi2_fisher=ifelse(test," and khi-2 or Fisher exact test","")
  fun.quanti=function(x,type.quanti){
    mean=round(mean(x,na.rm=TRUE),digits); sd=round(sd(x,na.rm=TRUE),digits);
    q=round(quantile(x,probs=c(0,0.25,0.5,0.75,1),na.rm=TRUE),digits)
    names(q)=c("min","q1","med","q3","max")
    return(switch(type.quanti,   
       "mean" = paste(mean," (",sd,")",sep=""),
       "med" = paste(q["med"]," [",q["q1"],";",q["q3"],"]",sep=""),
       "mean_mm" = paste(mean," (",sd,")"," (",q["min"],";",q["max"],")",sep=""),
       "med_mm" = paste(q["med"]," [",q["q1"],";",q["q3"],"]"," (",q["min"],";",q["max"],")",sep=""),
       "mean_med" = paste(mean," (",sd,") ",q["med"]," [",q["q1"],";",q["q3"],"]",sep=""),
       "mean_med_mm" = paste(mean," (",sd,") ",q["med"]," [",q["q1"],";",q["q3"],"]"," (",q["min"],";",q["max"],")",sep="")))
  }
  if (type.quanti=="med"){cat(paste("med [Q1;Q3]",type.test.quanti," for numeric variables\n",sep=""))} 
  if (type.quanti=="med_mm"){cat(paste("med [Q1;Q3] (min;max)",type.test.quanti," for numeric variables\n",sep=""))} 
  if (type.quanti=="mean"){cat(paste("mean (sd)",type.test.quanti," for numeric variables\n",sep=""))}
  if (type.quanti=="mean_mm"){cat(paste("mean (sd) (min;max)",type.test.quanti," for numeric variables\n",sep=""))}
  if (type.quanti=="mean_med"){cat(paste("mean (sd) med [Q1;Q3]",type.test.quanti," for numeric variables\n",sep=""))}
  if (type.quanti=="mean_med_mm"){cat(paste("mean (sd) med [Q1;Q3] (min;max)",type.test.quanti," for numeric variables\n",sep=""))}    
  cat(paste("N (%)",khi2_fisher," for categorical variables\n",sep=""))
  
  ####################### setting the output matrix ############################
  levgp=levels(group); nb.group=length(levgp);
  ncol=2+nb.group+1*whole+1*test
  out=matrix("",nrow=2,ncol=ncol)
  if (is.null(group.labels)){colnames=paste(groupname,levgp,sep=sep)} else{colnames=group.labels}
  out[1,]=c(parameter,"",if(whole){whole_pop}else{NULL},colnames,if(test){pvalue}else{NULL})
  out[2,]=c("","",if(whole){paste("(N=",nrow(data),")",sep="")}else{NULL},paste("(N=",table(group.g),")",sep=""),if(test){""}else{NULL})    

  ######################## loop for each covariate #############################
  for (i in 1:length(vars)){                   
    name_var=vars.labels[i]
    var=data[,vars[i]]
    
    if (is.numeric(var)){                                                       # for each numeric variable
      var.g=data.g[,vars[i]]
      res=NULL
      for (lev_i in levgp){res=c(res,fun.quanti(var.g[group.g==lev_i],type.quanti))}
      if (whole){res=c(fun.quanti(var,type.quanti),res)}
      p.value=NULL
      if (test){
        p.value=NA
        if (test.quanti=="param"){
          if (nb.group>2){
            p.value=tryCatch(anova(lm(var.g~group.g))[1,5],
                     error=function(e){cat(paste("Error: problem of ANOVA with covariate",vars[i],"\n"));return(NA);},
                     warning=function(w){options(warn=-1);
                                         cat(paste("Warning: problem of test with covariate",vars[i],"\n"));
                                         p=anova(lm(var.g~group.g))[1,5];
                                         options(warn=0);
                                         return(p)})
          } else{
            p.value=tryCatch(t.test(var.g~group.g)$p.value,
                     error=function(e){cat(paste("Error: problem of t-test with covariate",vars[i],"\n"));return(NA);},
                     warning=function(w){options(warn=-1);
                                         cat(paste("Warning: problem of test with covariate",vars[i],"\n"));
                                         p=t.test(var.g~group.g)$p.value;
                                         options(warn=0);
                                         return(p)})
          }
        } else {
          if (nb.group>2){
            p.value=tryCatch(kruskal.test(var.g~group.g)$p.value,
                     error=function(e){cat(paste("Error: problem of Kruskal-Wallis test with covariate",vars[i],"\n"));return(NA);},
                     warning=function(w){options(warn=-1);
                                         cat(paste("Warning: problem of test with covariate",vars[i],"\n"));
                                         p=kruskal.test(var.g~group.g)$p.value;
                                         options(warn=0);
                                         return(p)})
          } else{
            p.value=tryCatch(wilcox.test(var.g~group.g)$p.value,
                     error=function(e){cat(paste("Error: problem of Wilcoxon test with covariate",vars[i],"\n"));return(NA);},
                     warning=function(w){options(warn=-1);
                                         cat(paste("Warning: problem of test with covariate",vars[i],"\n"));
                                         p=wilcox.test(var.g~group.g)$p.value;
                                         options(warn=0);
                                         return(p)})
          }
        }
      }
      if (!is.null(p.value)){p.value=round(p.value,3)}
      out=rbind(out,c(name_var,"",res,p.value))
      # getting the number of NA if any
      if (I(any(is.na(var)) & whole) | any(is.na(var.g))){
        out.na=NULL
        for (lev_i in levgp){out.na=c(out.na,sum(is.na(var.g[group.g==lev_i])))}
        if (whole){out.na=c(sum(is.na(var)),out.na)}
        out=rbind(out,c("","NA",out.na,if(test){""}else{NULL}))
      }
    
    } else{                                                                     # for each factor/character/logical variable
      var=factor(var,levels=if(is.factor(var)){                                 # transforming characters and factors in factors
                              levels(var)
                            }else{
                              sort(unique.default(as.character(var)),na.last=TRUE)
                            })
      var.g=var[!is.na(group)]                                 
      tab=table(var.g,group.g)
      tab.test=table(var.g,group.g)      
      tab.p=round(100*prop.table(tab,margin=2),digits)
      if (whole){
        tab=cbind(table(var),tab)
        tab.p=cbind(round(100*prop.table(table(var)),digits),tab.p)
      }
      tab.p=ifelse(is.nan(tab.p),0,tab.p)
      p.value=NULL
      if (test){
        p.value=NA
        if (nrow(tab.test)>1){
          if (any(tab.test<5)){
            p.value=tryCatch(fisher.test(tab.test)$p.value,
                     error=function(e){cat(paste("Error: problem of Fisher test with covariate",vars[i],"\n"));return(NA);},
                     warning=function(w){options(warn=-1);
                                         cat(paste("Warning: problem of Fisher test with covariate",vars[i],"\n"));
                                         p=fisher.test(tab.test)$p.value;
                                         options(warn=0);
                                         return(p)})
          } else{
            p.value=tryCatch(chisq.test(tab.test,correct=FALSE)$p.value,
                     error=function(e){cat(paste("Error: problem of khi-2 test with covariate",vars[i],"\n"));return(NA);},
                     warning=function(w){options(warn=-1);
                                         cat(paste("Warning: problem of khi-2 test with covariate",vars[i],"\n"));
                                         p=chisq.test(tab.test,correct=FALSE)$p.value;
                                         options(warn=0);
                                         return(p)})
          }
        }
      }
      if (nrow(tab)>0){
        res=matrix(paste(tab," (",tab.p,"%)",sep=""),ncol=ncol(tab))
      } else{
        res=matrix(nrow=0,ncol=ncol(tab))
      }
      res.out=matrix("",nrow=nrow(res)+1*I(I(any(is.na(var)) & whole) | any(is.na(var.g))),ncol=ncol)
      res.out[1,1]=name_var
      if (!is.null(p.value)){res.out[1,ncol]=round(p.value,3)}
      if (nrow(tab)>0){
        res.out[1:nrow(res),3:(2+nb.group+1*whole)]=res
      }
      res.out[1:nrow(res.out),2]=c(levels(var),if(I(I(any(is.na(var)) & whole) | any(is.na(var.g)))){"NA"}else{NULL})
      if (I(any(is.na(var)) & whole) | any(is.na(var.g))){
        res.na=NULL
        for (lev_i in levgp){res.na=c(res.na,sum(is.na(var.g[group.g==lev_i])))}
        if (whole){res.na=c(sum(is.na(var)),res.na)}
        res.out[nrow(res.out),3:(2+nb.group+1*whole)]=res.na
      }      
      out=rbind(out,res.out)
    }
  }
  
  ############################# last checkings #################################
  if (test){
    out[,ncol(out)]=ifelse(out[,ncol(out)]==0,"<0.001",out[,ncol(out)])
    out[,ncol(out)]=ifelse(is.na(out[,ncol(out)]),"---",out[,ncol(out)])
  }
  dimnames(out)=list(rep("",nrow(out)),rep("",ncol(out)))
  if (language=="french"){out[-1,-c(1,2)]=gsub("\\.", ",", out[-1,-c(1,2)])}    # replace . by , (\\. ou [.])  
  if (justify){
    out[,1]=format(out[,1],justify="right")
    out[,2]=format(out[,2],justify="left")
    out[,-c(1:2,if(test){ncol(out)}else{NULL})]=format(out[,-c(1:2,if(test){ncol(out)}else{NULL})],justify="right")
    if (test){out[,ncol(out)]=format(out[,ncol(out)],justify="left")}
  }  
  if (all(out[,2]==rep("",nrow(out)))){out=out[,-2]}
  
  ################################# export #####################################  
  if (!is.null(file.export)){
    out.df=as.data.frame(out,row.names=1:nrow(out))
    WriteXLS("out.df",file.export,"Descriptive statitics",col.names=FALSE,AdjWidth=TRUE)
    cat("Export created: ",file.export,"\n",sep="")
    if (noquote){out=noquote(out)}
    return(invisible(out))
  } else{
    if (noquote){out=noquote(out)}
    return(out)
  }
}
