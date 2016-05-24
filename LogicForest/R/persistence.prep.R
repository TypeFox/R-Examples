persistence.prep <-
function(fit, preds, PI)
{
 per1<-persist.match(fit=fit, preds=preds, match.list=PI)
 PI<-PI
 PIsz<-per1$uniq.sz
 PIsizes<-per1$all.sz
 names(PI)<-PI
 mches<-per1$matches[[1]][[1]]
 frq<-fit$PI.frequency
 PIoI.freq<-ifelse (PI%in%names(frq), frq[PI], 0)
 allPI.list<-names(frq)
 m<-length(unique(unlist(mches)))
 if (m==0 & PIoI.freq==0 | m==1 & PIoI.freq>0) stop("There are no matches for this PI")
 sub.sizes<-c()
 if (PIoI.freq==0) {nms<-unique(unlist(mches))} 
 if (PIoI.freq>0) {nms<-unique(unlist(mches))[2:m]} 
 all.frq<-c(PIoI.freq, frq[c(nms)])
 for(i in 1:length(nms))
   {
   PIset<-setdiff(unlist(strsplit(nms[i], " ")), "&")
   sub.sizes<-append(sub.sizes, length(PIset))
   }
 names(sub.sizes)<-nms
 if(PIoI.freq==0) {lger.PIs<-unique(unlist(mches))}
 if(PIoI.freq>0) {lger.PIs<-unique(unlist(mches))[2:m]}
 match.mat<-matrix(0, nrow=length(lger.PIs), ncol=length(lger.PIs))
 rownames(match.mat)<-colnames(match.mat)<-nms
 if (PIsz%in%PIsizes==0) {comp<-length(mches)-1}
 if (PIsz%in%PIsizes==1) {comp<-length(mches)-2}
 for(i in 1:comp)#loop to fill in matrix
   {
   if (PIsz%in%PIsizes==0) {PI.nms<-as.vector(unique(mches[[i]]))}
   if (PIsz%in%PIsizes==1) {PI.nms<-as.vector(unique(mches[[i+1]]))}
   PI.ids<-which(colnames(match.mat)%in%PI.nms)
   if (length(PI.nms!=0)){
     sb<-subs(fit=fit, allPI.list=allPI.list, preds=preds, match.list=PI.nms)$matches[[1]]
     for(j in 1:length(sb))
       {
       PI.id<-PI.ids[j]
       sb.mtches<-as.vector(unique(unlist(sb[[j]])))
       sb.ids<-which(rownames(match.mat)%in%sb.mtches)
       match.mat[c(sb.ids),PI.id]<-1
       }
     }
   }
 ids1<-which(rowSums(match.mat)==0)
 primes<-names(which(rowSums(match.mat)==0))
 prime.nms<-c()
 for (i in 1:length(primes))
   {
   PIoI<-setdiff(unlist(strsplit(PI, " ")), "&")
   pr.PI<-setdiff(unlist(strsplit(primes[i], " ")), "&")
   diff<-setdiff(pr.PI, PIoI)
   if (length(diff)==1) {newnm<-diff}
   if (length(diff)==2) {newnm<-paste(diff[1], "&", diff[2])}
   if (length(diff)==3) {newnm<-paste(diff[1], "&", diff[2], "&", diff[3])}
   if (length(diff)==4) {newnm<-paste(diff[1], "&", diff[2], "&", diff[3], "&", diff[4])}
   if (length(diff)==5) {newnm<-paste(diff[1], "&", diff[2], "&", diff[3], "&", diff[4], "&", diff[5])}
   prime.nms<-append(prime.nms, newnm)
   }
 len1<-sub.sizes[c(primes)]
 le<-length(len1)
 pos1<-rad(seq(0, 359, by=360/le))
 po<-length(pos1)
 if(le==po) {pos<-pos1; len<-len1}
 if(le!=po) {m<-min(le,po);pos<-pos1[1:m];len<-len1[1:m]}
 rad.mat<-cbind(len, pos)
 names(pos)<-primes
 prime.xy<-matrix(0, nrow=length(len), ncol=2)
 freq1<-c()#frequency occurence primary nodes
 for (i in 1:length(len))
   {
   prime.xy[i,]<-c(cos(pos[i])*len[i], sin(pos[i])*len[i])
   prime.freq<-frq[primes[i]]
   freq1<-append(freq1, prime.freq)
   }
 ids.2plus<-which(rowSums(match.mat)>0)
 mtch.2plus<-names(which(rowSums(match.mat)>0))
 sec.mtch<-c()
 if (length(ids.2plus)>0)
   {
   for (i in 1:length(ids.2plus))
     {
     PInm<-mtch.2plus[i]
     mtch.ids<-which(match.mat[PInm,]==1)
     if(length(mtch.ids)==sum(as.numeric(mtch.ids%in%ids1))) {sec.mtch<-append(sec.mtch, PInm)}
     }
   }
 ids2<-which(colnames(match.mat)%in%sec.mtch)
 ids12<-c(as.vector(ids1),ids2)
 mtch.3plus<-colnames(match.mat[,-c(ids12)])
 tert.mtch<-c()
 if (length(mtch.3plus)>0)
   {
   for(i in 1:length(mtch.3plus))
     {
     PInm<-mtch.3plus[i]
     mtch.ids<-as.vector(which(match.mat[PInm,]==1))
     if(length(mtch.ids)==sum(as.numeric(mtch.ids%in%ids12))) 
       {
       tert.mtch<-append(tert.mtch, PInm)
       new.ids<-sec.mtch[c(which(sec.mtch%in%colnames(match.mat[, mtch.ids])))]
       match.mat[PInm,]<-0
       match.mat[PInm,new.ids]<-1
       }
     }
   }
 ids3<-which(colnames(match.mat)%in%tert.mtch)
 ids123<-c(ids12, ids3)
 mtch.4plus<-colnames(match.mat[,-c(ids123)])
 quat.mtch<-c()
 if (length(mtch.4plus)>0)
   {
   for(i in 1:length(mtch.4plus))
     {
     PInm<-mtch.4plus[i]
     mtch.ids<-which(match.mat[PInm,]==1)
     if(length(mtch.ids)==sum(as.numeric(mtch.ids%in%ids123))) 
       {
       quat.mtch<-append(quat.mtch, PInm)
       new.ids<-tert.mtch[c(which(tert.mtch%in%colnames(match.mat[,mtch.ids])))]
       match.mat[PInm,]<-0
       match.mat[PInm,new.ids]<-1
       }
     }
   }
 ids4<-which(colnames(match.mat)%in%quat.mtch)
 ids1234<-c(ids123, ids4)
 mtch.5plus<-colnames(match.mat[,-c(ids1234)])
 quin.mtch<-c()
 if (length(mtch.5plus)>0)
   {
   for(i in 1:length(mtch.5plus))
     {
     PInm<-mtch.5plus[i]
     mtch.ids<-which(match.mat[PInm,]==1)
     if(length(mtch.ids)==sum(as.numeric(mtch.ids%in%ids1234))) 
       {
       quin.mtch<-append(quin.mtch, PInm)
       new.ids<-quat.mtch[c(which(quat.mtch%in%colnames(match.mat[,mtch.ids])))]
       match.mat[PInm,]<-0
       match.mat[PInm,new.ids]<-1
       }
     }
   }
 add.rad<-rad(c(0,3,3,6,6,9,9,12,12,15,15,18,18,21,21,24,24,27,27,30,30,33,33,36,36,39,39,42,42))
 mtch.list<-list(primes, sec.mtch, tert.mtch, quat.mtch, quin.mtch)
 posit.list<-list(pos)
 freq.sec<-c()
 for (i in 2:length(mtch.list))
   {
   if(length(mtch.list[[i]])>0)
     {
     if(i==2)
       {
       sec.nms<-c()
       names2<-c()
       freq2<-c()
       x1s<-c(); y1s<-c(); x2s<-c(); y2s<-c()
       for (j in 1:length(mtch.list[[i-1]]))
         {
         prime.nm<-mtch.list[[i-1]][j]
         p.nm<-setdiff(unlist(strsplit(prime.nm, " ")), "&")
         mch.nm2<-names(which(match.mat[,prime.nm]==1))
         if (length(mch.nm2)>0)
           {
           if(length(mch.nm2)>29) {
                           wrn<-paste("There are >29 subset matches for PI = ", prime.nm, ", only the first 29 will be included on the plot.", sep="")
                           warning(wrn)
                           }
           for (k in 1:min(length(mch.nm2), 29))
             {
             pos1<-posit.list[[i-1]][j]
             names2<-append(names2, mch.nm2[k])
             freq2<-append(freq2, frq[mch.nm2[k]])
             freq.sec<-append(freq.sec, frq[mch.nm2[k]])
             nm2<-setdiff(setdiff(unlist(strsplit(mch.nm2[k], " ")), "&"), p.nm)
             if (length(nm2)==1) {nm.2<-nm2}
             if (length(nm2)==2) {nm.2<-paste(nm2[1], "&", nm2[2])}
             if (length(nm2)==3) {nm.2<-paste(nm2[1], "&", nm2[2], "&", nm2[3])}
             if (length(nm2)==4) {nm.2<-paste(nm2[1], "&", nm2[2], "&", nm2[3], "&", nm2[4])}
             sec.nms<-append(sec.nms, nm.2)
             if(even(k)==1) {new.pos2<-pos1+add.rad[k]}
             if(odd(k)==1) {new.pos2<-pos1-add.rad[k]}
             sz1<-sub.sizes[prime.nm]
             sz2<-sub.sizes[mch.nm2[k]]
             x1<-cos(pos1)*sz1; y1<-sin(pos1)*sz1
             x2<-cos(new.pos2)*sz2; y2<-sin(new.pos2)*sz2
             x1s<-append(x1s, x1); y1s<-append(y1s, y1)
             x2s<-append(x2s, x2); y2s<-append(y2s, y2)
             } 
           }
         }
       pos2.mat<-cbind(x1s, y1s, x2s, y2s, freq2, sec.nms)
       rownames(pos2.mat)<-names2
       posit.list[[i]]<-pos2.mat
       }
     if(i==3)
       {
       tert.nms<-c()
       names3<-c()
       freq3<-c()
       x1s<-c(); y1s<-c(); x2s<-c(); y2s<-c()
       for (j in 1:length(mtch.list[[i-1]]))
         {
         sec.nm<-mtch.list[[i-1]][j]
         s.nm<-setdiff(unlist(strsplit(sec.nm, " ")), "&")
         mch.nm3<-names(which(match.mat[,sec.nm]==1))
         if (length(mch.nm3)>0)
           {
           for (k in 1:length(mch.nm3))
             {
             pos2<-acos(as.numeric(posit.list[[i-1]][sec.nm, 3])/sub.sizes[sec.nm])
             names3<-append(names3, mch.nm3[k])
             freq3<-append(freq3, frq[mch.nm3[k]])
             freq.sec<-append(freq.sec, frq[mch.nm2[k]])
             nm.3<-setdiff(setdiff(unlist(strsplit(mch.nm3[k], " ")), "&"), s.nm)
             if (length(nm.3)==1) {nm3<-nm.3}
             if (length(nm.3)==2) {nm3<-paste(nm.3[1],"&",nm.3[2])}
             if (length(nm.3)==3) {nm3<-paste(nm.3[1],"&",nm.3[2],"&",nm.3[3])}
             tert.nms<-append(tert.nms, nm3)
             if(even(k)==1) {new.pos3<-pos2+add.rad[k]}
             if(odd(k)==1) {new.pos3<-pos2-add.rad[k]}
             sz2<-sub.sizes[sec.nm]
             sz3<-sub.sizes[mch.nm3[k]]
             x1<-cos(pos2)*sz2; y1<-sin(pos2)*sz2
             x2<-cos(new.pos3)*sz3; y2<-sin(new.pos3)*sz3
             x1s<-append(x1s, x1); y1s<-append(y1s, y1)
             x2s<-append(x2s, x2); y2s<-append(y2s, y2)
             } 
           }
         }
       pos3.mat<-cbind(x1s, y1s, x2s, y2s, freq3, tert.nms)
       rownames(pos3.mat)<-names3
       posit.list[[i]]<-pos3.mat
       }
     if(i==4)
       {
       quat.nms<-c()
       names4<-c()
       freq4<-c()
       x1s<-c(); y1s<-c(); x2s<-c(); y2s<-c()
       for (j in 1:length(mtch.list[[i-1]]))
        {
        tert.nm<-mtch.list[[i-1]][j]
        t.nm<-setdiff(unlist(strsplit(tert.nm, " ")), "&")
        mch.nm4<-names(which(match.mat[,tert.nm]==1))
        if (length(mch.nm4)>0)
          {
          for (k in 1:length(mch.nm4))
            {
            pos3<-acos(as.numeric(posit.list[[i-1]][tert.nm, 3])/sub.sizes[tert.nm])
            names4<-append(names4, mch.nm4[k])
            freq4<-append(freq4, frq[mch.nm4[k]])
            freq.sec<-append(freq.sec, frq[mch.nm4[k]])
            nm.4<-setdiff(setdiff(unlist(strsplit(mch.nm4[k], " ")), "&"), t.nm)
            if (length(nm.4)==1) {nm4<-nm.4}
            if (length(nm.4)==2) {nm4<-paste(nm.4[1],"&",nm.4[2])}
            quat.nms<-append(quat.nms, nm4)
            if(even(k)==1) {new.pos4<-pos3+add.rad[k]}
            if(odd(k)==1) {new.pos4<-pos3-add.rad[k]}
            sz3<-sub.sizes[tert.nm]
            sz4<-sub.sizes[mch.nm4[k]]
            x1<-cos(pos3)*sz3; y1<-sin(pos3)*sz3
            x2<-cos(new.pos4)*sz4; y2<-sin(new.pos4)*sz4
            x1s<-append(x1s, x1); y1s<-append(y1s, y1)
            x2s<-append(x2s, x2); y2s<-append(y2s, y2)
            } 
          }
        }
      pos4.mat<-cbind(x1s, y1s, x2s, y2s, freq4, quat.nms)
      rownames(pos4.mat)<-names4
      posit.list[[i]]<-pos4.mat
      }
    if(i==5)
      {
      quin.nms<-c()
      names5<-c()
      freq5<-c()
      x1s<-c(); y1s<-c(); x2s<-c(); y2s<-c()
      for (j in 1:length(mtch.list[[i-1]]))
        {
        quat.nm<-mtch.list[[i-1]][j]
        q.nm<-setdiff(unlist(strsplit(quat.nm, " ")), "&")
        mch.nm5<-names(which(match.mat[,quat.nm]==1))
        if (length(mch.nm5)>0)
          {
          for (k in 1:length(mch.nm5))
            {
            pos4<-acos(as.numeric(posit.list[[i-1]][quat.nm, 3])/sub.sizes[quat.nm])
            names5<-append(names5, mch.nm5[k])
            freq5<-append(freq5, frq[mch.nm5[k]])
            freq.sec<-append(freq.sec, frq[mch.nm5[k]])
            nm5<-setdiff(setdiff(unlist(strsplit(mch.nm5[k], " ")), "&"), q.nm)
            quin.nms<-append(quin.nms, nm5)
            if(even(k)==1) {new.pos5<-pos4+add.rad[k]}
            if(odd(k)==1) {new.pos5<-pos4-add.rad[k]}
            sz4<-sub.sizes[quat.nm]
            sz5<-sub.sizes[mch.nm5[k]]
            x1<-cos(pos4)*sz4; y1<-sin(pos4)*sz4
            x2<-cos(new.pos5)*sz5; y2<-sin(new.pos5)*sz5
            x1s<-append(x1s, x1); y1s<-append(y1s, y1)
            x2s<-append(x2s, x2); y2s<-append(y2s, y2)
            } 
          }
        }
      pos5.mat<-cbind(x1s,y1s,x2s,y2s,freq5,quin.nms)
      rownames(pos5.mat)<-names5
      posit.list[[i]]<-pos5.mat
      }
     }   
   }
 freqs<-c(PIoI.freq, freq1, freq.sec)
 prime.mat<-cbind(prime.xy, freq1, prime.nms)
 szes<-sort(unique(sub.sizes), decreasing=TRUE)
 ans<-list(all.freqs=all.frq, rad.mat=rad.mat,main.freq=PIoI.freq, primary.info=prime.mat,
           position.list=posit.list, sizes=szes)
}
