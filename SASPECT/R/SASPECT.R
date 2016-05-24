########################################################
###  SASPECT spectral counting algorithm
###
###   Author:  Pei Wang, Fred Hutchinson Cancer Research Center 
###
###  Implemented in R for use with CPAS version 8.1, www.labkey.org
###  To call, see RunScoringFunction
###
############ parameters
### peptideData     : list of peptide measurement matrices, contains
###     pep.count       : numeric matrix (px(n1+n2)), peptide spectral counts
###     pep.data        : numeric matrix (px(n1+n2)), peptide prophet score
### proteinInfo     : optional dataframe containing protein-dependent data (not used here)
### pep.set         : character vector (1xp), peptide name for the two data matrix
### pep.pro.name    : numeric matrix (?x2), first column is for protein name and the second column is for peptide name
### run.group.info  : data frame , first column is run group name, second is number of runs in group
### permu.iter      : integer, number of iterations for permutation to infer FDR
### filter.run      : integer, filter criteria to remove peptides observed in too few number of samples.
### filter.score    : integer, minimum peptide prophet score to be counted in filter


SASPECT<-function(peptideData, pep.set,  pep.pro.name, run.group.info, permu.iter=50, filter.run=2, filter.score=.95)
{

    pep.count=peptideData$PeptideCount 
    pep.data=peptideData$PeptideConfidence
    n1=as.integer(run.group.info$count[1])
    n2=as.integer(run.group.info$count[2])
    RunGrp1=as.character(run.group.info$label[1])
    RunGrp2=as.character(run.group.info$label[2])

    ########################################################
    ############ constants   
    NSP.ratio.sm<-c(3.6926159, 2.7337886, 1.8848071, 1.2373227, 0.8057487, 
                0.5403145, 0.3933011, 0.3128286, 0.2695229, 0.2445471, 0.2191232, 
                0.1966275, 0.1693055, 0.1358687, 0.1045338)


    pool<-as.vector(pep.data)
    p0<-mean(pool) #### average probability for one pep is observed
    PP0<-mean(pool[pool>0])
    N=nrow(pep.data)
    M=ncol(pep.data)

    pro.set<-unique(pep.pro.name[,1]) ###### unique protein names



    ################################################################
    ########### shrink factor for different sample
    m.v<-apply(pep.data>0, 2, sum)
    S.scale<-(N-m.v)/m.v*p0/(1-p0)

    S.scale.new<-S.scale
    for(i in 1:length(m.v))
    {
        cur.m<-m.v[i]
        S.scale.new[i]<-dbinom(cur.m, size=N-1, prob=p0)/dbinom(cur.m-1, size=N-1, prob=p0)
    }


    ################################################################
    ################## using information of NSP (number of sibling peptides) to update peptide prophet

    NSP.data<-get.NSP.new(pro.set, pep.pro.name, pep.set, pep.data)
    NSP.d.tr<-NSP.data
    NSP.d.tr[NSP.data>15]<-15



    ########################################################## 
    ###  update the peptide prophet score

    pep.data.new<-apply(as.matrix(1:M), 1, modify.PP, pep.m=pep.data, NSP.m=NSP.d.tr, NSP.ratio=NSP.ratio.sm)


    ##########################################################
    ################# filtering data

    pick.f<-apply(pep.data>filter.score, 1, sum)>filter.run #### only keep those peptides having at least two reliable identifications.

    pep.set.f<-pep.set[pick.f]
    pep.data.new.f<-pep.data.new[pick.f, ]
    pep.count.f<-pep.count[pick.f, ]

    S.scale.f<-S.scale


    #######################################################################################
    ##################################### get concise protein list 

    pep.pro.name.f<-pep.pro.name[is.element(pep.pro.name[,2], pep.set.f),]

    pro.set.temp<-unique(pep.pro.name.f[,1])
    pro.set.detail<-pro.set.temp
    for(i in 1:length(pro.set.temp))
    {
        cur<-pep.pro.name.f[pep.pro.name.f[,1]==pro.set.temp[i],2]
        pro.set.detail[i]<-paste(sort(cur), collapse=".")
    }

    pep.pro.group<-NULL
    temp<-unique(pro.set.detail)
    for(i in 1:length(temp))
    {
        cur.pep<-temp[i]
        cur.pro<-pro.set.temp[pro.set.detail==cur.pep]
        pep.pro.group<-rbind(pep.pro.group, c(cur.pro[1], paste(cur.pro, collapse=".")))
    }

    N.pro<-nrow(pep.pro.group)

    ########################
    temp.pro<-unique(pep.pro.name.f[,1])
    temp.pep<-unique(pep.pro.name.f[,2])
    pro.pep.count<-1:nrow(pep.pro.name.f)
    pep.pro.count<-1:nrow(pep.pro.name.f)

    for(i in 1:length(temp.pro))
    {
        cur.pro<-temp.pro[i]
        cur.pick<-pep.pro.name.f[,1]==cur.pro
        pro.pep.count[cur.pick]<-length(unique(pep.pro.name.f[cur.pick,2]))
    }

    for(i in 1:length(temp.pep))
    {
        cur.pep<-temp.pep[i]
        cur.pick<-pep.pro.name.f[,2]==cur.pep
        pep.pro.count[cur.pick]<-length(unique(pep.pro.name.f[cur.pick,1]))
    }


    ########################################################## 
    ######################################### begin EM
    ##print("Begin EM")
    obs.s.n<-apply(as.matrix(1:N.pro), 1, score.onepro.nnew, pro.set=pep.pro.group[,1], pep.pro.name=pep.pro.name.f, PP0=PP0,
               pep.set=pep.set.f, pepdata=pep.data.new.f, pep.count=pep.count.f, pep.pro.count=pep.pro.count,
               Sscale=S.scale.f, n1=n1, n2=n2)

    ###########################################################
    ######################################### permutation
    permuate.score.new<-array(0, dim=c(3, N.pro, permu.iter))
    for(i in 1:permu.iter)
    {
        ##print(i)
        index<-sample(1:ncol(pep.data.new.f),ncol(pep.data.new.f), replace=F)
        permu.data<-pep.data.new.f[,index]
        permu.count<-pep.count.f[,index]
        permu.scale<-S.scale.f[index]
        cur.score<-apply(as.matrix(1:N.pro), 1, score.onepro.nnew, pro.set=pep.pro.group[,1], pep.pro.name=pep.pro.name.f, PP0=PP0,
               pep.set=pep.set.f, pepdata=permu.data, pep.count=permu.count, pep.pro.count=pep.pro.count, 
               Sscale=permu.scale, n1=n1, n2=n2)
        permuate.score.new[,,i]<-cur.score
    }

    permu.score.sort<-apply(permuate.score.new[3,,1:permu.iter], 2, sort)
    null.p.mean<-apply(permu.score.sort, 1, mean)

    score.order<-order(obs.s.n[3,])
    obs.s.new.sort<-sort(obs.s.n[3,])

    Qvalue.new<-apply(as.matrix(1:ncol(obs.s.n)), 1, qvalue.cut, p.mean=null.p.mean, 
                 p.obs=obs.s.new.sort, p.pool=as.vector(permu.score.sort))



    ################################################################ output

    report.pro.n<-pep.pro.group[score.order,]
    report.pep.count<-1:nrow(report.pro.n)
    for(i in 1:length(report.pep.count))
    {
        cur.pro<-report.pro.n[i,1]
        report.pep.count[i]<-length(unique(pep.pro.name.f[pep.pro.name.f[,1]==cur.pro, 2]))
    }

    report.obs.s<-t(obs.s.n[,score.order])

    result=cbind(report.pro.n[,1:2], report.obs.s, Qvalue.new, report.pep.count)
    
    n=nrow(result)
    result=result[n:1,]
    colnames(result)=as.character(c("Protein", "ProteinsInGroup" , "Score1",  "Score2",   "Score",  "Qvalue",   "PeptideNumber"))

    return(result)
}

#############################################################################
#############################################################################


################################################### basic functions

qvalue.cut<-function(i, p.mean, p.obs, p.pool)
{
  delta<-abs(p.obs[i]-p.mean[i])
  
  up=5
  if(sum(p.obs>0 & p.obs-p.mean>=delta)>0)
      up<-min(p.obs[p.obs>0 & p.obs-p.mean>=delta])
  
  obs.count<-sum(p.obs>=up)
  null.count<-sum(p.pool>=up)
 
  if(obs.count>0)
   {
  result<-null.count/length(p.pool)*length(p.obs)/obs.count
   } else {
  result<-1
   }
  return(result)
}

my.ttest<-function(x.v, y.v)
{
 
   x.m<-mean(x.v)
   y.m<-mean(y.v)

   n<-length(x.v)
   m<-length(y.v)

   t.score<-0
   if( (n+m)>5 )
     {
      #var.all<-(var(x.v)*(n-1)+var(y.v)*(m-1))/(m+n-2)
      var.all<-var(c(x.v, y.v))
   
      if(var.all>0)
       {
        t.score<-(y.m-x.m)/((var.all*(1/n+1/m))^0.5+0.01)
       }
      }       
   return(t.score)
}



################################################################################
################### 5/2007

score.onepro.nnew<-function(i, pro.set, pep.pro.name, pep.set, pepdata, pep.count,pep.pro.count, Sscale, n1, n2, PP0, EM=TRUE)
{ #### input variables ####
   ### pro.set:  unique protein name in the experiment
   ### pep.pro.name: matching table between peptide and protein
   ### pep.set:  peptide name for each row of "pepdata"
   ### pepdata:  data base search confidence for each peptide
   ### pep.count: spectral count for each peptide in each sample
   ### pep.pro.count: number of proteins each peptide associating to. Length is equal to row # of pep.pro.name

  #### take the target protein
     pep.list=pep.pro.name[pep.pro.name[,1]==pro.set[i],2]
     pick<-is.element(pep.set, pep.list)
   
   ##print(i)

   score<-NA

   if(sum(pick)==1)
     {
      xdata<-pepdata[pick,]
      ydata<-pep.count[pick,]
      score<-final.score.simple.new(xdata, ydata, pp0=PP0,s.scale=Sscale, n1=n1, n2=n2)
     } 
     
   if(sum(pick)>1)
    {
      xdata<-pepdata[pick,]
      ydata<-pep.count[pick,]
      
      pep.weight<-1:sum(pick)
      for(k in 1:sum(pick))
       {
         cur.pep<-pep.set[pick][k]
         pep.weight[k]<-1/mean(pep.pro.count[pep.pro.name[,2]==cur.pep])
       }

      if(EM)
       {

        score<-final.score.new(xdata, ydata, pep.weight=pep.weight, pp0=PP0,s.scale=Sscale, n1=n1, n2=n2)
       } else{
        score<-final.score.noEM.new(xdata, ydata, pep.weight=pep.weight, pp0=PP0, s.scale=Sscale, n1=n1, n2=n2)
       }  
    }
   return(score)
}


#######################################################
##################### if protein only has one peptide:

final.score.simple.new<-function(xdata,  ydata, pp0, s.scale, n1, n2)
{
     ################# adjust the pp scores
     xdata.ad<-adjust.pp(xdata, p0=pp0)
   
     ######################## adjust scale
     zfinal<-g.scale(s.scale, xdata.ad)         
     x.indi<-zfinal[1:n1]
     y.indi<-zfinal[-(1:n1)]
     score1<-my.ttest(x.indi, y.indi)

     ######################## t score based on spectral count
     spec.ad<- zfinal*ydata
     score2<-my.score2(spec.ad, n1, n2)
     
     ######################## 
     result<-score1^2+score2^2
     return(c(score1, score2, result))
}

my.score2<-function(spec.ad, n1, n2)
{
        x.raw<-spec.ad[1:n1]
        y.raw<-spec.ad[-(1:n1)]
        
        x.v=0
        y.v=0
          
        if(sum(x.raw>0)>0)
             x.v<-x.raw[x.raw>0]
        if(sum(y.raw>0)>0)  
             y.v<-y.raw[y.raw>0]    
       
        score2<-my.ttest(x.v, y.v)
        return(score2)
}

#################################################################################
##################### if protein only has multiple peptides:
###
###  TODO:  iter is never passed in by callers of this function.  should it be passed down from permu.iter at the top level ?
###

final.score.noEM.new<-function(xdata, ydata, pep.weight=pep.weight, pp0, s.scale, n1, n2, iter=100)
{
    ######################################## part 1: for phat1-phat2

  nulldata<-xdata[,1:n1]
  alterdata<-xdata[,-(1:n1)]
  
  znull.0<-1-apply(nulldata, 2, prob.notE)
  zalter.0<-1-apply(alterdata, 2, prob.notE)
 
  znull.old<-znull.0
  zalter.old<-zalter.0

########################## adjust the pp scores

  nulldata.ad<-adjust.pp(nulldata, p0=pp0)
  alterdata.ad<-adjust.pp(alterdata,p0=pp0)  

  para<-para.esti(znull.old, zalter.old, nulldata.ad, alterdata.ad, s.scale[1:n1], s.scale[-(1:n1)], H0=F)
  znull.new<-z.esti.pp(para[1], para[-(1:2)], nulldata.ad, s.scale[1:n1])
  zalter.new<-z.esti.pp(para[2], para[-(1:2)], alterdata.ad, s.scale[-(1:n1)])

  zfinal<-g.scale(s.scale,c(znull.new, zalter.new)) 
  x.indi<-zfinal[1:n1]
  y.indi<-zfinal[-(1:n1)]
  score1<-my.ttest(x.indi, y.indi)
  
  ####################################### part2:  for u1-u2
   ydata.ad<-ydata*matrix(zfinal, nrow=nrow(ydata), ncol=ncol(ydata), byrow=TRUE)
   score2<-nrow(ydata)^0.5* sum(apply(ydata.ad, 1, my.score2, n1=n1, n2=n2) * pep.weight)/sum(pep.weight)
      
  ####################################### combine part1 and part2
   result<-score1^2+score2^2
   return(c(score1, score2,result))
}

###############################################################################
###  TODO:  iter is never passed in by callers of this function.  should it be passed down from permu.iter at the top level ?
###

final.score.new<-function(xdata, ydata,pep.weight=pep.weight, pp0, s.scale, n1, n2, iter=100)
{
  ######################################## part 1: for phat1-phat2
  ######################## inital the EM

  nulldata<-xdata[,1:n1]
  alterdata<-xdata[,-(1:n1)]
  
  if(n1>1)
  {
  znull.0<-1-apply(nulldata, 2, prob.notE)
  } else{
  znull.0<-1-prob.notE(nulldata)
  }

  if(n2>1)
  {
  zalter.0<-1-apply(alterdata, 2, prob.notE)
  } else{
  zalter.0<-1-prob.notE(alterdata)
  }

  znull.old<-znull.0
  zalter.old<-zalter.0

  ########################## adjust the pp scores

  nulldata.ad<-adjust.pp(nulldata, p0=pp0)
  alterdata.ad<-adjust.pp(alterdata,p0=pp0)  

  psuedo.x<-as.numeric(apply(xdata>0, 2, sum)>0)
  flag<-TRUE

  ########################### do EM
   for(i in 1:iter)
    {
      para<-para.esti(znull.old, zalter.old, nulldata.ad, alterdata.ad, s.scale[1:n1], s.scale[-(1:n1)], H0=F)
      znull.new<-z.esti.pp(para[1], para[-(1:2)], nulldata.ad, s.scale[1:n1])
      zalter.new<-z.esti.pp(para[2], para[-(1:2)], alterdata.ad, s.scale[-(1:n1)])
 
      flag<-all(c(zalter.old, znull.old)==c(zalter.new, znull.new))
      if(flag)
        break;
       
      znull.old<-znull.new
      zalter.old<-zalter.new
 
     }
  ############################# adjust scale
     zfinal<-g.scale(s.scale,c(znull.new, zalter.new)) 
     x.indi<-zfinal[1:n1]
     y.indi<-zfinal[-(1:n1)]
     score1<-my.ttest(x.indi, y.indi)
  
  ####################################### part2:  for u1-u2
     ydata.ad<-ydata*matrix(zfinal, nrow=nrow(ydata), ncol=ncol(ydata), byrow=TRUE)
     score2<-nrow(ydata)^0.5* sum(apply(ydata.ad, 1, my.score2, n1=n1, n2=n2) * pep.weight)/sum(pep.weight)
     
  ####################################### combine part1 and part2
     result<-score1^2+score2^2
     return(c(score1, score2,result))
}


###################################################
z.esti.pp<-function(p0.old, pk, data, ss)
{
  n<-length(ss)

  p0<-f.scale(ss, p0.old)

  temp<-exp(sum(log(1-pk)))

  temp1=(p0*temp+(1-p0))
  temp1[temp1==0]=0.0001
  p.new<-p0*temp/temp1
  
  if(ncol(as.matrix(data))>1)
  {
   p.x<-apply(data, 2, prob.notE)
  } else {
   p.x<-prob.notE(data)
  }

  znew<-1-p.x + p.x*p.new
  znew
}

####################################################
prob.notE<-function(pp.v)
{
  exp(sum(log(1-pp.v)))
}

####################################################
modify.PP<-function(i, pep.m, NSP.m, NSP.ratio) #### maximum value of NSP.v does not exceed the length of NSP.ratio
{
pep.v<-pep.m[,i]
NSP.v<-NSP.m[,i]
result<-pep.v
pick<-(1:length(pep.v))[pep.v!=0]
result[pick]<-pep.v[pick]/(pep.v[pick]+(1-pep.v[pick])*NSP.ratio[(NSP.v[pick])])
return(result)
}

####################################################
########### 2/2007
get.NSP.new<-function(pro.set, pep.pro.name, pep.set, p.m)
{
  n<-length(pro.set)
  result<-matrix(0, nrow=nrow(p.m), ncol=ncol(p.m))

  for(i in 1:n)
    {
     pep.list=pep.pro.name[pep.pro.name[,1]==pro.set[i],2]
     pick<-is.element(pep.set, pep.list)

     if(sum(pick)>1)
      {  
       nsp.v<-apply(p.m[pick, ]>0, 2, sum)
       temp<-((p.m[pick,]>0)+0)*matrix(nsp.v, nrow=sum(pick), ncol=length(nsp.v), byrow=TRUE) 
       result[pick, ]<-pmax(temp, result[pick, ])
      } else {
       nsp.v<-(p.m[pick,]>0)+0
       result[pick,]<-pmax(nsp.v, result[pick,])
      }
    }
   result
}

###########
g.scale<-function(s.scale, z.v)
{
s.scale*z.v/(1-z.v+s.scale*z.v)
}

###########
f.scale<-function(s.scale, p.v)
{
p.v/(p.v+(1-p.v)*s.scale)
}

############
adjust.pp<-function(data.m, p0)
{
  data.ad<-data.m
  data.ad[data.m==0]=1-p0
 
  temp<-data.m[data.m!=0]
  data.ad[data.m!=0]=temp+(1-p0)*(1-temp)
  return(data.ad)
}


para.esti<-function(znull, zalter, nulldata, alterdata,s.scale.null, s.scale.alter, H0)
 {
    znull.old<-g.scale(s.scale.null, znull)
    zalter.old<-g.scale(s.scale.alter, zalter) 

    if(H0)
     {
       temp<-mean(c(znull.old, zalter.old))
       p0=temp
       p1=temp
    } else{     
    p0<-mean(znull.old)
    p1<-mean(zalter.old)
     }

    zv<-c(znull, zalter)
    xdata<-cbind(nulldata, alterdata)
  
    pk<-apply(xdata, 1, pk.esti, zv=zv)
    result<-c(p0, p1, pk)
    return(result)
}

z.esti<-function(p0.old, pk, data, z0, ss)
{
  znew<-z0
  n<-length(z0)
  change<-(1:n)[z0==0]

  p0<-f.scale(ss, p0.old)

  temp<-exp(sum(log(1-pk)))
  p.new<-p0*temp/(p0*temp+(1-p0))
  znew[change]<-p.new[change]
  znew
}


pk.esti<-function(xv, zv)
{
   sum(xv*zv)/sum(zv)
}


#########################################
###
###  TODO:  following two functions are not called.  remove?
###

my.SOnePro<-function(i, pepset, proname, pepdata, Sscale, n1, n2)
{
   pick<-proname==pepset[i]
   if(sum(pick)==1)
     {
      xdata<-pepdata[pick,]
      pep.count<-1
      control.count<-sum(xdata[1:n1]>0)
      case.count<-sum(xdata[-(1:n1)]>0)
   } else {
      xdata<-pepdata[pick,]
      pep.count<-nrow(xdata)
      control.count<-sum(as.vector(xdata[,1:n1]>0))
      case.count<-sum(as.vector(xdata[,-(1:n1)]>0))
   }
   c(
, control.count, case.count)
}
####################################################
countscans<-function(pep.m) 
{
pep.v<-pep.m[,1]
scans.v<-pep.m[,2]
result<-scans.v
pick<-(1:length(pep.v))[pep.v>=.95]

result[]<-
return(result)
}