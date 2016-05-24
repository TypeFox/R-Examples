isopattern <-
function(iso_list,compound,limit){
  
    ############################################################################
    # (1) generate isotope list for given compound & convert compound name string to a list:
    isos<-iso_list;
    isos<-isos[isos[,4]!=0,];    
    getit<-seq(1,length(isos[,1]),1);
    getthat<-c();
    element<-c();                          # store element names
    number<-c();                           # store number of atoms of this element
    ende<-nchar(compound);
    i<-c(1);
    while(i<=ende){
      warn<-TRUE;
      # on element names with square brackets:
      if(substr(compound,i,i)==c("[")){
        k<-i;
        while(any(substr(compound,i,i)==c("]"))!=TRUE)
            {i<-c(i+1);};
        while(any(substr(compound,i,i)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE)
            {i<-c(i+1);};
        m<-c(i-1);
        element<-c(element,substr(compound,k,m));
        getthat<-c(getthat,getit[isos[,1]==substr(compound,k,m)]);
        warn<-FALSE;          
      }
      # on element names without square brackets:
      if(any(substr(compound,i,i)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE){
        k<-i;
        while(any(substr(compound,i,i)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE)
          {i<-c(i+1);};
        m<-c(i-1);
        i<-c(i-1);
        element<-c(element,substr(compound,k,m));
        getthat<-c(getthat,getit[isos[,1]==substr(compound,k,m)]);
        warn<-FALSE;
        };
      # on element numbers:
      if(any(substr(compound,i,i)==c("0","1","2","3","4","5","6","7","8","9"))==TRUE){
        k<-i;
        while(any(substr(compound,i,i)==c("0","1","2","3","4","5","6","7","8","9"))==TRUE){i<-c(i+1);}
        m<-c(i-1);
        i<-i-1;
        number<-c(number,as.numeric(substr(compound,k,m)));
        warn<-FALSE;
        };
    if(warn==TRUE){stop("Calculation interrupted: compound chemical formula incorrect!");}
    i<-i+1;
    } # end while loop
    # check if all elements with entry in iso_list
    for(i in 1:length(element)){if(any(isos[,1]==element[i])==FALSE){stop("Calculation interrupted: one element not found in iso_list");};};
    ########## add [13]C to that algorithm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    isos<-t(isos[getthat,]);
    ############################################################################
  
    ############################################################################
    # (2) generate mass and abundance of monoisotopic peak = starting peak:
    monomass<-as.double(0);              # save mass of the monoisotopic peak, one value
    monoabund<-as.double();                   # save abundance of the monoisotopic peak, vector per element block
    for(i in 1:length(element)){
      getit<-isos[,isos[1,]==element[i]];
      if(length(dim(getit))>0){
          monomass<-c(monomass+(as.double(getit[3,as.numeric(getit[4,])==max(as.numeric(getit[4,]))])*as.numeric(number[i])));
          monoabund<-c(monoabund,(as.double(getit[4,as.numeric(getit[4,])==max(as.numeric(getit[4,]))])^as.numeric(number[i])));
        }else{
          monomass<-c(monomass+(as.double(getit[3])*as.numeric(number[i])));
          monoabund<-c(monoabund,(as.double(getit[4])^as.numeric(number[i])));
        }
    }
    if(monomass==0){stop("Calculation interrupted: monoisotopic mass could not be calculated!");};
    if(length(monoabund)==0){stop("Calculation interrupted: monoisotopic abundance incorrect!");};
    if(length(monoabund)!=length(number)){stop("Calculation interrupted: not all elements found in iso_list");};
    ############################################################################
  
    ############################################################################
    # (3) setup matrices & parameters:
    # (a) indicate which entries in compo are monoisotopic peaks:
    getit<-seq(1,length(isos[1,]),1);
    road<-matrix(nrow=length(element),ncol=10,0);
    rownames(road)=element;
    colnames(road)=c("monoiso",rep("nonmono",9));
    for(i in 1:length(element)){    
      getthat<-getit[isos[1,]==element[i]]
      road[i,1]=getthat[as.numeric(isos[4,isos[1,]==element[i]])==max(as.numeric(isos[4,isos[1,]==element[i]]))];
      getthat<-getthat[as.numeric(isos[4,isos[1,]==element[i]])!=max(as.numeric(isos[4,isos[1,]==element[i]]))];
      if(length(getthat)>0){road[i,c(2:(1+length(getthat)))]=getthat};
      };
    ############################################################################    
    # (b) matrix to contain peak masses and abundances & initialize for monoisotopic peak:
    leng1<-length(isos[1,]);
    peaks<-matrix(ncol=(4+leng1),nrow=500000,as.double(0));
    leng2<-length(peaks[1,]);
    leng3<-length(peaks[,1]);
    peaks[1,1]=monomass;                      # mass
    peaks[1,2]=prod(monoabund);               # abundance over all element blocks
    peaks[1,3]=0;                             # generation number
    peaks[1,4]=0;                             # stop further calculation if below threshold?
    peaks[1,4+road[,1]]=number;               # atom number distribution
    colnames(peaks)=c("mass","abundance","generation","stop?",isos[2,]);
    if(length(monoabund)==1 && monoabund==1){peaks[1,4]=1};
    # (c) matrix indices:
    k<-c(1); # matrix index: where to write from
    m<-c(2); # matrix index: where to write to
    b_start<-c(2); # matrix index: where is start of block
    b_end<-c(1); # matrix index: where is end of block
    g<-c(1); # generation number
    ############################################################################
    
    ############################################################################
    # (4) reset limit by maximum permutational number found among all elements & isotopes
    # ...skipped 
    minlimit<-limit;
    isonumber<-c(); # store number of isotopes; needed from step(6) onwards!  
    #isocounts<-c(); # store approx. max combination numbers per element
    getit<-rep(1,length(road[1,]));
    for(i in 1:length(road[,i])){
      isonumber<-c(isonumber,length(getit[road[i,]!=0]));
      #isocounts<-c(isocounts,combi(rep(floor(number[i]/length(getit[road[i,]!=0])),length(getit[road[i,]!=0]))));
    };
    #minlimit<-limit/max(isonumber);
    ############################################################################
    
    ############################################################################
    # (5.a) start algorithm for first generation peaks based on monoisotopic peak
    for(i in 1:length(road[,1])){
      if(isonumber[i]>1){
          for(j in 2:isonumber[i]){                    
            #if(peaks[k,4+road[i,1]]!=0){
              # insert composition of parent peak:
              peaks[m,c(5:(4+leng1))]=peaks[k,c(5:(4+leng1))];
              # assign generation number:
              peaks[m,3]=g;
              # calculate new mass from old one:
              peaks[m,1]=c(peaks[k,1]-(as.numeric(isos[3,road[i,1]]))+(as.numeric(isos[3,road[i,j]])));
              # calculate new abundance from old one; simple for first substitution:
              peaks[m,2]=c(peaks[k,2]*  as.numeric(isos[4,road[i,j]]) / (peaks[m,4+road[i,j]]+1) * peaks[m,4+road[i,1]] / as.numeric(isos[4,road[i,1]])  );
              # update composition
              peaks[m,4+road[i,1]]=c(peaks[m,4+road[i,1]]-1);
              peaks[m,4+road[i,j]]=c(peaks[m,4+road[i,j]]+1);
              # below limit = stop?
              if(peaks[m,2]<minlimit){peaks[m,4]=1};
              # where to write to next entry into peak list?
              m<-m+1;
          #} # subtraction possible? used up? -> not necessary for first; enable later!      
        } # over j along row of road, with road entry != 0
      } # if isonumber[i]>1
    } # over i along columns of road
    ############################################################################
    # (5.b) sort outcomes, check (a) for same masses and (b) abundances, then (c) for same composition
    # = not necessary for first generation !
    g<-c(g+1); # generation number
    k<-c(k+1); # matrix index: where to write from
    b_end<-c(m-1);
    ############################################################################
  
    ############################################################################
    # (6.a) repeat algorithm over all other generations g2 to gx ###############
    # stop criteria: (a) limit & (b) length of peak vector
    while(any(peaks[peaks[,3]==(g-1),4]==0) && m<leng3 && b_end>=b_start){ # && g<=100 
    ############################################################################
      # along all peaks of one generation g:
      while(k<=b_end){
        if(peaks[k,2]>=limit){
          for(i in 1:length(road[,1])){
            if(isonumber[i]>1){
                if(peaks[k,4+road[i,1]]!=0 && m<leng3 ){ # any monoisotopic atoms left to substract for that element?
                  for(j in 2:isonumber[i]){                    
                    # insert composition of parent peak:
                    peaks[m,c(5:(4+leng1))]=peaks[k,c(5:(4+leng1))];
                    # assign generation number:
                    peaks[m,3]=g;                    
                    # calculate new abundance from old one:
                    # keep that order to minimize rounding errors!
                    peaks[m,2]=c(peaks[k,2]*  as.numeric(isos[4,road[i,j]]) / (peaks[m,4+road[i,j]]+1) * peaks[m,4+road[i,1]] / as.numeric(isos[4,road[i,1]])  );
                    # calculate new mass from old one:
                    peaks[m,1]=round(peaks[k,1]-(as.numeric(isos[3,road[i,1]]))+(as.numeric(isos[3,road[i,j]])),digits=9);                    
                    # update composition
                    peaks[m,4+road[i,1]]=c(peaks[m,4+road[i,1]]-1);
                    peaks[m,4+road[i,j]]=c(peaks[m,4+road[i,j]]+1);
                    # below limit = stop?
                    if(peaks[m,2]<minlimit){peaks[m,4]=1};
                    # where to write to next entry into peak list?
                    m<-m+1;
                }; # over j along row of road, with road entry != 0
              }; # subtraction possible? used up? -> not necessary for first; enable later!      
            }; # if isonumber[i]>1
          };  # over i along columns of road 
          k<-c(k+1);
        }else{ # if below limit: no further calculation wanted = step to next entry in peak list
          k<-c(k+1);  
        }; # check if above limit!
      };  # while k within generation g
      ############################################################################
      # (6.b) sort outcomes, check (a) for same masses and (b) abundances, then (c) for same composition
      b_end<-c(m-1); # where new block ends
      b_start<-c(k); # where new block starts
      if(b_end>b_start){ # more than one entry to compare at all? 
        getit<-seq(b_start,b_end,1);
        getthat<-c();
        back<-getit[order(peaks[getit,1],decreasing=FALSE)];
        peaksort<-peaks[back,];
        for(i in 1:(length(peaksort[,1])-1)){                             # over all entries in generation g
          if(peaksort[i,1]==peaksort[i+1,1]){                             # same masses?
            if(round(peaksort[i,2],digits=3)==round(peaksort[i+1,2],digits=3)){                           # same abundances?
              if(all(peaksort[i,c(5:leng2)]==peaksort[i+1,c(5:leng2)])){  # same composition?
                  getthat<-c(getthat,back[i+1]);
              }
            }
          }
        }
        leng4<-length(getthat);
        if(leng4>0){
          peaks<-peaks[-getthat,];
          m<-c(m-leng4);
          b_end<-c(b_end-leng4);
          leng3<-c(leng3-leng4);
          rm(peaksort);
        }
      } # if b_end>b_start
      g<-c(g+1); # update generation counter    
    ############################################################################
    } # while on limit and/or generation number
    ############################################################################
    
    if(m>=leng3){warning("Storage maximum for number of peaks reached!");};
    
    ###########################################################################
    # clean up ...
    rm(road, isonumber, k, b_start, b_end, leng1, leng2, leng3);
    peaks<-peaks[peaks[,1]!=0,];
    if(m>2){ # must be more than one element / isotope!
      peaks<-peaks[peaks[,2]!=0,] # for zero-probabaility isotopes eg. 36-Cl
      peaks<-peaks[order(peaks[,1],decreasing=FALSE),]
      }
    # scale to monositopic peak!
    # peaks[,2]<-c(peaks[,2]/peaks[1,2]);
    return(peaks);
} # function end

