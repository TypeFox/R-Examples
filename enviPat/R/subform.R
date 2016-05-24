subform<-
function(formula1,formula2){

  ##############################################################################
  # on second chemical formula (deduct) ########################################
  formula2<-gsub("D","[2]H",formula2);
  ende2<-nchar(formula2);
  element2<-c();
  number2<-c();
  ##############################################################################
  j<-c(1);
  while(j<=ende2){
    if(substr(formula2,j,j)==c("[")){
          b<-j
          while(any(substr(formula2,j,j)==c("]"))!=TRUE){
              j<-c(j+1);
          };
          k<-j
          while(any(substr(formula2,j,j)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE){
              j<-c(j+1);
          };
          m<-c(j-1);
          element2<-c(element2,substr(formula2,b,m));
    }
    if(any(substr(formula2,j,j)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE){
          k<-j;
          while(any(substr(formula2,j,j)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE){
            j<-c(j+1);
          };
          m<-c(j-1);
          j<-c(j-1);
          element2<-c(element2,substr(formula2,k,m));
    };
    if(any(substr(formula2,j,j)==c("0","1","2","3","4","5","6","7","8","9"))==TRUE){
          k<-j;
          while(any(substr(formula2,j,j)==c("0","1","2","3","4","5","6","7","8","9"))==TRUE){
            j<-c(j+1);
          }
          m<-c(j-1);
          j<-c(j-1);
          number2<-c(number2,as.numeric(substr(formula2,k,m)));
    };
  j<-j+1;
  } # end while loop
  ##############################################################################
  
  ##############################################################################
  # on first chemical formula(s) ###############################################
  formulas<-c()
  for(i in 1:length(formula1)){
    formula1[i]<-gsub("D","[2]H",formula1[i]);
    ende1<-nchar(formula1[i]);
    element1<-c();
    number1<-c();
    # on formula1 = vector of formulas
    j<-c(1);
    while(j<=ende1){
      if(substr(formula1[i],j,j)==c("[")){
            b<-j
            while(any(substr(formula1[i],j,j)==c("]"))!=TRUE){
                j<-c(j+1);
            };
            k<-j
            while(any(substr(formula1[i],j,j)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE){
                j<-c(j+1);
            };
            m<-c(j-1);
            element1<-c(element1,substr(formula1[i],b,m));
      }
      if(any(substr(formula1[i],j,j)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE){
            k<-j;
            while(any(substr(formula1[i],j,j)==c("0","1","2","3","4","5","6","7","8","9"))!=TRUE){
              j<-c(j+1);
            };
            m<-c(j-1);
            j<-c(j-1);
            element1<-c(element1,substr(formula1[i],k,m));
      };
      if(any(substr(formula1[i],j,j)==c("0","1","2","3","4","5","6","7","8","9"))==TRUE){
            k<-j;
            while(any(substr(formula1[i],j,j)==c("0","1","2","3","4","5","6","7","8","9"))==TRUE){
              j<-c(j+1);
            }
            m<-c(j-1);
            j<-c(j-1);
            number1<-c(number1,as.numeric(substr(formula1[i],k,m)));
      };
    j<-j+1;
    } # end while loop
    # check: subtract possible? ##################################################
    formula_all<-TRUE
    for(i in 1:length(element2)){
      if(any(element2[i]==element1)==FALSE){
        formula_all<-paste(element2[i]," from formula 2 not part of formula1",sep="")
      }else{
        if(number2[i]>number1[element2[i]==element1]){
          formula_all<-paste("Atom number of ",element2[i]," from formula 2 not fully subset of formula1 atom number",sep="")
        }
      }
    }
    ##############################################################################
    # subtract ! #################################################################
    if(formula_all==TRUE){  
      formula_all<-""
      counts<-c();
      for(i in 1:length(element1)){
         if(any(element2==element1[i])){
             counts<-c(counts,
              number1[i]-(number2[element2==element1[i]])
             );
         }else{
             counts<-c(counts,number1[i])
         }
      }
      element1<-element1[counts>0];
      counts<-counts[counts>0];
      for(i in 1:length(counts)){
        formula_all<-paste(formula_all,element1[i],counts[i],sep="")
      }
    }
    formulas<-c(formulas,formula_all)
  }
  return(formulas)
  ##############################################################################

}

