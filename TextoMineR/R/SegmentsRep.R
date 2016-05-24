SegmentsRep <-function(base, num.text, lower= TRUE, sep.strong="[()?./:!=+;{}-]", nxlon=10,
               nfreq=10, nfreq2=10, nfreq3=10) 
{

if (!is.null(num.text)) {
 if (is.character(num.text)) 
 num.text<- which(colnames(base) %in% num.text)
  if (is.numeric(num.text)) 
  num.text<- num.text
 if(length(num.text)==1)
  num.text<-num.text 
 if(length(num.text)>1){
   for(i in 1:length(num.text)){
    if(i==1)
    text1<-base[,num.text[1]]
    else text1<-paste(text1,base[,num.text[i]],sep=".") 
   }
    base[,(ncol(base)+1)]<-text1
    num.text<-ncol(base)
  }
}
# auxiliary functions
maj.in.min = lower
REPWEAK<-function(chaine,sep.weak) res<-str_replace_all(chaine,sep.weak, " ")
REPSTRONG<-function(chaine,sep.strong) res<-str_replace_all(chaine,sep.strong, " zzwwxxyystr ")
PROCHE<-function(ideb,ifin,ITEX,ITDR,ITRE,nfreq,nfreq2,nfreq3,ltrou,long,nxlon,nbseg)
      {   
# the function proche detects the first sublist of adresses in ITDR corresponding a a same successor
# if this successor is not "end of answer" or "strong separator", we have located an admissible segment
	list.segment<-list()
       ad.segment<-vector()
       te.segment<-NULL
       rep.segment<-vector()
       mfrec<-0
       ipunt<-ideb-1
       isucc<-ITEX[ITDR[ideb]+long-1] 
       while(  (ITEX[ITDR[ipunt+1]+long-1]==isucc) &  (ipunt < ifin) )
            {
            if (!( (isucc=="zzwwxxyystr") | (isucc=="zzwwxxyyendrep")))
               {
               mfrec<-mfrec+1
               ad.segment[mfrec]<-ITDR[ipunt+1]
               rep.segment[mfrec]<-min(which(ITRE>ad.segment[mfrec]))
               }
            ipunt<-ipunt + 1
            }

       ifin<-ipunt
       nfreq.threshold<-nfreq
       if (long==1)  nfreq.threshold<-999999999999
       if (long==2)  nfreq.threshold<-nfreq2
       if (long==3)  nfreq.threshold<-nfreq3
       ltrou<-( !( (isucc=="zzwwxxyystr") | (isucc=="zzwwxxyyendrep"))  & (mfrec >= nfreq)) 
       ltrouseg<-( ltrou & (mfrec >= nfreq.threshold))      
       if (ltrouseg) 
          { 
            for (i in 1:long) 
              {
	         te.segment<-paste(te.segment,ITEX[ITDR[ideb]+(i-1)],sep=" ")
               te.segment<-str_trim(te.segment)
              }
           
           lo.segment<-long
           fr.segment<-mfrec
           nr.segment<-nbseg+1
           list.segment<-list(te.segment,fr.segment,ad.segment,rep.segment,lo.segment,nr.segment)
           names(list.segment)<-c("text","frequency","adresses","documents","length","nr.seg")
          }       
        
       return(list(ifin=ifin,ltrou=ltrou,ltrouseg=ltrouseg,list.segment=list.segment))
       }                                                                                                                                                                                                                                                                                                                                                                                                                                               



ORD.EXT<-function(ICRIT,ADR,long1)
      {
# Ordering adresses from successors in text
       ICRIT_ord<-order(ICRIT)
       ADR<-ADR[c(ICRIT_ord)]
       return(list(ADR=ADR))
       }

sep.weak ="([\u0027\u02BC\u0060]|[,;'?\n<U+202F><U+2009><U+0028>]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+"
if (nfreq2<nfreq) nfreq2<-nfreq
if (nfreq3<nfreq) nfreq3<-nfreq
if (nfreq3>nfreq2)nfreq3<-nfreq2

text<-as.character(base[,num.text])
text<-as.matrix(text)
text1<-apply(text,1,FUN=REPSTRONG,sep.strong)
text1<-as.matrix(text1)
text2<-apply(text1,1,FUN=REPWEAK,sep.weak)
text2<-as.vector(text2)
text3<-str_c(text2,"zzwwxxyyendrep",sep=" ")

nrep<-NROW(text3)

listrep<-(strsplit(as.character(text3),split=" ")) 

ITEX <- unlist(listrep) 

# ITEX is a vector of words, but with fictitious "empty" words because of the multiple spaces
# these fictitious words have to be eliminated
	sel <- which(ITEX=="")      
	if (length(sel)==0){
          if (maj.in.min == TRUE)  ITEX <- tolower(ITEX)
	}
	if (length(sel)!=0){	 		
    	ITEX <- ITEX[-sel]           
	if (maj.in.min == TRUE)  ITEX <- tolower(ITEX)       
	}

# The text is in the form of a vector ocf occurrences  
ITEX.f<-as.factor(ITEX)
FREQ.mots<-table(ITEX.f)
FREQ.cum<-cumsum(FREQ.mots)
Vplus<-dim(FREQ.mots)

# To conserve the addresses when ordering ITEX
ITDR<-order(ITEX) 

# adress of the answers (=adress of the first word corresponding to the answer) in ITEX
ITRE<-vector()
ITRE<-which(ITEX=="zzwwxxyyendrep")

##############################################################
#the data structures are built:  principal part of the function
##############################################################
Nplus<-length(ITEX)
ITDR<-seq(1,Nplus,1)
lpil<-vector()
list.tot.segment<-list()
# global initialisations
ideb<-1
ifin<-Nplus
long<-0
nbseg<-0
# for all the distinct words, we have to detect the segments beginning with this word
ltrou<-((ifin-ideb+1) >= nfreq)  

while(ltrou)      
     {                
      while (ltrou & (long<=nxlon))       #exploration of the possible segments issued from word_in_course
         {
            long1<-long
	      long<-long+1
            lpil[long]<-ifin
            res.ORD.EXT<-ORD.EXT(ITEX[ITDR[ideb:ifin]+long1],ITDR[ideb:ifin],long1)
            ITDR[ideb:ifin]<-res.ORD.EXT$ADR
            ltrou<-FALSE
            res.proch<-PROCHE(ideb,ifin,ITEX,ITDR,ITRE,nfreq,nfreq2,nfreq3,ltrou,long,nxlon,nbseg)
            ifin<-res.proch$ifin
            ltrou<-res.proch$ltrou
            ltrouseg<-res.proch$ltrouseg
            if (ltrouseg)
              {
               nbseg<-res.proch$list.segment[[6]]
               list.segment<-res.proch$list.segment
               list.tot.segment[[nbseg]]<-list.segment 
               }         
         }

  
          ltrou<-FALSE
          while (!ltrou & (long>=1) & (ifin<Nplus))
          {
            ideb=ifin+1
            ifin=lpil[long]
            while (  ( (ideb+nfreq)>ifin) & (long > 1) )
          {
              ideb=ifin+1
              long<-long-1
              ifin=lpil[long]
           }

        if (long>=1)
            {
              ltrou<-FALSE
              res.proch<-PROCHE(ideb,ifin,ITEX,ITDR,ITRE,nfreq,nfreq2,nfreq3,ltrou,long,nxlon,nbseg)
              ifin<-res.proch$ifin
              ltrou<-res.proch$ltrou
              ltrouseg<-res.proch$ltrouseg
              if (ltrouseg)
              {
                 nbseg<-res.proch$list.segment[[6]]
                 list.segment<-res.proch$list.segment
                 list.tot.segment[[nbseg]]<-list.segment
               }                
             }     

           }


        }   
 # all the segments have been detected and the doc_segments (tab.seg) will be created

tab.seg<-matrix(0,nrow=nrep,ncol=nbseg)
if (nbseg==0) print ("no segments fullfil the conditions")
if (nbseg>0)
     {
      for (iseg in 1:nbseg)
        {
         list.segment<-list.tot.segment[[iseg]]
         mfreq<-list.segment[[2]]
         long.seg<-list.segment[[5]]
         nseg<-list.segment[[6]]
         for (i in 1:mfreq) 
           {
            rep<-list.segment[[4]][i]
            tab.seg[rep,nseg]<-tab.seg[rep,nseg]+1
            }
         }
      row.names(tab.seg)<-row.names(base)
      nom.col<-vector()
      for (iseg in 1:nbseg) nom.col[iseg]<-(list.tot.segment[[iseg]]$text)
      colnames(tab.seg)<-nom.col
     }

impri.segment<-data.frame(ncol=3)
# liste des segments en ordre alphabetique
for (iseg in 1:nbseg) 
    {
      impri.segment[iseg,1]<-list.tot.segment[[iseg]]$text
      impri.segment[iseg,2]<-list.tot.segment[[iseg]]$frequency
      impri.segment[iseg,3]<-list.tot.segment[[iseg]]$length
     }
colnames(impri.segment)<-c("segment","frequency","long")
segOrderFreq<-with(impri.segment,impri.segment[order(frequency,long,decreasing=TRUE),])
segOrderlist<-impri.segment
Glossary.segments<-list(segOrderFreq=segOrderFreq, segOrderlist=segOrderlist)

namesSeg<-colnames(tab.seg)
numSeg<-rep(1:ncol(tab.seg),1)
colnames(tab.seg) = paste(numSeg, namesSeg, sep=":")

 res<-list(list.tot.segment=list.tot.segment,nbseg=nbseg,tab.seg=tab.seg, Glossary.segments=Glossary.segments)
 class(res)<-c("SegmentsRep","list")   
 return(res)
  
 }
