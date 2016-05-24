check_chemform <-
function(isotopes,chemforms,get_sorted=FALSE){

    ############################################################################
    # internal function definitions ############################################
    # (A) Multiplier ###########################################################
    multif<-
    function(formula1,fact,numbers){
        formulas<-c();
        ########################################################################
        # on first chemical formula ############################################
        formula1<-gsub("D","[2]H",formula1);
        ende1<-nchar(formula1);
        element1<-c();
        number1<-c();
        ########################################################################
        # on formula1 
        j<-c(1);
        while(j<=ende1){
          if(substr(formula1,j,j)==c("[")){
                b<-j
                while(any(substr(formula1,j,j)==c("]"))!=TRUE){
                    j<-c(j+1);
                };
                k<-j
                while(any(substr(formula1,j,j)==numbers)!=TRUE){
                    j<-c(j+1);
                };
                m<-c(j-1);
                element1<-c(element1,substr(formula1,b,m));
          }
          if(any(substr(formula1,j,j)==numbers)!=TRUE){
                k<-j;
                while(any(substr(formula1,j,j)==numbers)!=TRUE){
                  j<-c(j+1);
                };
                m<-c(j-1);
                j<-c(j-1);
                element1<-c(element1,substr(formula1,k,m));
          };
          if(any(substr(formula1,j,j)==numbers)==TRUE){
                k<-j;
                while(any(substr(formula1,j,j)==numbers)==TRUE){
                  j<-c(j+1);
                }
                m<-c(j-1);
                j<-c(j-1);
                number1<-c(number1,as.numeric(substr(formula1,k,m)));
          };
        j<-j+1;
        } # end while loop
        ########################################################################
        # multiply ! ###########################################################
        number1<-fact*number1
        formula_fin<-""
        for(p in 1:length(element1)){
          formula_fin<-paste(formula_fin,element1[p],number1[p],sep="")
        }
        formulas<-c(formulas,formula_fin);
      return(formulas)
      ##########################################################################
    }
    ############################################################################
    
    ############################################################################
    capitals <- c("[","A", "B", "C", "D", "E", "F", "G", "H",
        "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
        "T", "U", "V", "W", "X", "Y", "Z")
    numbers <- c("0", "1", "2", "3", "4", "5", "6", "7", "8",
        "9")
    allpossible<-c(capitals,numbers,"(",")","]","a","b","c","d","e","f","g","h","i",
      "j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"
    )
    masses <- c()
    warn <- c()
    elem <- unique(as.character(isotopes[, 1]))
    isotopes2 <- matrix(nrow = length(elem), ncol = 2)
    isotopes2[, 1] <- elem                                     
    for (i in 1:length(elem)) {
        intermed <- isotopes[isotopes[, 1] == elem[i], ]
        if (is.vector(intermed) == TRUE) {
            isotopes2[, 2][isotopes2[, 1] == elem[i]] <- intermed[3]
        }
        else {
            isotopes2[, 2][isotopes2[, 1] == elem[i]] <- intermed[,
                3][as.numeric(intermed[, 4]) == max(as.numeric(intermed[,
                4]))]
        }
    }
    info <- isotopes2           
    for (i in 1:length(chemforms)) {
        mass <- c(0); 
        warnit <- FALSE;
        if(chemforms[i]==""){
          warnit <- TRUE;
        }
        # (0) split string #####################################################
        if(warnit==FALSE){
          formel <- as.character(chemforms[i]);
          formel <- strsplit(formel," ")[[1]];
          m <- strsplit(formel, as.character());
          m <- m[[1]];
        }
        # (1) all characters plausible? ########################################
        if(warnit==FALSE){
          for(j in 1:length(m)){
            if(any(allpossible==m[j])==FALSE){
              warnit<-TRUE
            }
          }
        }
        # (2) do all [(bracket)] types close & [] only contain numbers? ########
        if(warnit==FALSE){        
          if( any(m=="[") || any(m=="]") || any(m=="(") || any(m==")") ){        
            getit1<-c(0);
            getit2<-c(0);
            a<-c(1);
            while((a)<=length(m)){
              if(m[a]=="["){getit1<-getit1+1;};
              if(m[a]=="]"){getit1<-getit1-1;};
              if(m[a]=="("){getit2<-getit2+1;};
              if(m[a]==")"){getit2<-getit2-1;};              
              if( getit1>0 & (any(numbers==m[a])==FALSE & m[a]!="[" & m[a]!="]") ){warnit<-TRUE};
              if( getit1<0 || getit2<0 ){warnit<-TRUE}
              a<-a+1;
            }
            if( getit1!=0 || getit2!=0 ){
              warnit<-TRUE;
            }
          }
        }
        # (3) start correct? ###################################################
        if(warnit==FALSE){
          if(any( any(capitals==m[1])==FALSE & m[1]!="(" & m[1]!=")")){
            warnit<-TRUE;
          }
          if(length(m)==1){
            m<-c(m,"1");
            formel<-paste(formel,"1",sep="");
          }
        }
        # (4) empty brackets? ##################################################
        if(warnit==FALSE){         
            for(k in 2:length(m)){
              if(m[k-1]=="(" & m[k]==")"){
                warnit<-TRUE;
              }
            }
        }
        # (5) insert 1 where missing ###########################################
        if(warnit==FALSE){
          # for closing )-brackets #############################################  
          if(any(m=="(")){
            m2<-c();
            for(j in 1:(length(m)-1)){
              if( 
                m[j]==")" & 
                any(m[j+1]==numbers)==FALSE
              ){
                m2<-c(m2,m[j],"1");
              }else{
                m2<-c(m2,m[j]);
              }
            }
            m2<-c(m2,m[length(m)]);
            if(m[length(m)]==")"){
              m2<-c(m2,"1");
            }
            m<-m2;
          }
          # for all other cases ################################################
          m2<-m[1];
          for(j in 2:length(m)){
            if(
              ( any(m[j]==capitals) || m[j]==")" || m[j]=="(" ) &
                all(m[j-1]!=numbers) & m[j-1]!="(" & m[j-1]!="]"      
            ){
              m2<-c(m2,"1",m[j]);
            }else{
              m2<-c(m2,m[j]);
            } 
          }
          if(all(m[length(m)]!=numbers)){
            m2<-c(m2,"1")
          }
          m<-m2;
          formel<-""
          for(k in 1:length(m)){
            formel<-paste(formel,m[k],sep="")
          }
        }
        # (6) multiply for square brackets, with nesting #######################
        if(warnit==FALSE){
          while(any(m=="(")){
            a<-c(1);
            getit1<-1;
            getit2<-1;
            while( getit1!=0 & getit2!=0 & a<=length(m)){              
              if(m[a]=="("){
                getit1<-2;  
                from<-a
              }
              if(m[a]==")"){
                getit2<-2;  
                to<-a
              }
              if(getit1==2 & getit2==2){
                 b<-a+1;
                 count<-""
                 while(any(m[b]==numbers & b<=length(m))){
                    count<-paste(count,m[b],sep="");
                    b<-b+1;
                  }
                  count<-as.numeric(count);
                  m2<-""
                  for(k in (from+1):(to-1)){
                    m2<-paste(m2,m[k],sep="")
                  }
                  m2<-multif(m2,count,numbers);
                  m2<-strsplit(m2, as.character())[[1]];
                  m3<-c();
                  doneit<-FALSE;
                  for(z in 1:length(m)){
                    if( z<from || z>=b){
                      m3<-c(m3,m[z])
                    }else{
                      if(doneit==FALSE & (z>=from || z<b)){
                        m3<-c(m3,m2)                     
                        doneit<-TRUE
                      }
                    }
                  }
                  m<-m3;
                  getit1<-0;
                  getit2<-0;
              }
              a<-a+1;
            }
          }
          formel<-""
          for(k in 1:length(m)){
            formel<-paste(formel,m[k],sep="")
          }
        }
        # (7) dissassemble #####################################################
        if(warnit==FALSE){
          element1<-c();
          number1<-c();
          ######################################################################
          j<-c(1);
          while(j<=nchar(formel)){
            if(substr(formel,j,j)==c("[")){
                  b<-j;
                  while(
                    any(substr(formel,j,j)==c("]"))!=TRUE &
                    j<=nchar(formel)
                  ){
                      j<-c(j+1);
                  };
                  k<-j;
                  while(any(substr(formel,j,j)==numbers)!=TRUE){
                      j<-c(j+1);
                  };
                  z<-c(j-1);
                  element1<-c(element1,substr(formel,b,z));
            }
            if(any(substr(formel,j,j)==numbers)!=TRUE){
                  k<-j;
                  while(
                    any(substr(formel,j,j)==numbers)!=TRUE &
                    j<=nchar(formel)
                  ){
                    j<-c(j+1);
                  };
                  z<-c(j-1);
                  j<-c(j-1);
                  element1<-c(element1,substr(formel,k,z));
            };
            if(any(substr(formel,j,j)==numbers)==TRUE){
                  k<-j;
                  while(
                    any(substr(formel,j,j)==numbers)==TRUE &
                    j<=nchar(formel)
                  ){
                    j<-c(j+1);
                  }
                  z<-c(j-1);
                  j<-c(j-1);
                  number1<-c(number1,as.numeric(substr(formel,k,z)));
            };
          j<-j+1;
          }# end while
        }
        # (8) check if all elements present in isotopes list ###################
        if(warnit==FALSE){
          for(j in 1:length(element1)){ 
            if(any(element1[j]==as.character(isotopes[,1]))==FALSE){
              warnit<-TRUE;
            }
          }
          if(length(element1)!=length(number1)){
            warnit<-TRUE;
          }
        }
        # (9) merge non-unique elements ########################################
        if(warnit==FALSE){
          element2<-c();
          number2<-c();
          doneit<-rep(FALSE,length(element1));
          for(j in 1:length(element1)){
            if(doneit[j]==FALSE){
              doneit[element1==element1[j]]<-TRUE;
              element2<-c(element2,element1[element1==element1[j]][1])
              number2<-c(number2,as.character(sum(as.numeric(number1[element1==element1[j]]))))
            }
          }
          element1<-element2;rm(element2);
          number1<-number2;rm(number2);
		  if(get_sorted){ # ensure unambiguous order of elements in the formula
				this<-order(match(element1,info))
				number1<-number1[this]
				element1<-element1[this]
		  }
          formel<-""
          for(k in 1:length(element1)){
            formel<-paste(formel,element1[k],number1[k],sep="")
            mass<-mass+(
              as.numeric(info[info[,1]==element1[k],2][1])*as.numeric(number1[k])
            )  
          }
       }        
        ########################################################################
        # make the final entry #################################################
        if(warnit==FALSE){
          warn<-c(warn,FALSE);
          masses<-c(masses,mass);
          chemforms[i]<-formel;
        }else{
          warn<-c(warn,TRUE);
          masses<-c(masses,-9999);  
        }
    }    
    ############################################################################    
    checked <- data.frame(warn, chemforms, masses)
    names(checked) <- c("warning", "new_formula", "monoisotopic_mass")
    checked[,2]<-as.character(checked[,2]);
    return(checked);
  
}        
        
        
        
           
      
        
        
        
        
        
        
        
      