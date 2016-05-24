mergeform<-
function(formula1,formula2){

  ##############################################################################
  # on second chemical formula #################################################
  formula2<-gsub("D","[2]H",formula2);
  ende2<-nchar(formula2);
  element2<-c();
  number2<-c();
  ##############################################################################
  # on formula2 (deduct)
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
  formulas<-c();
  for(i in 1:length(formula1)){
      ##############################################################################
      # on first chemical formula ##################################################
      formula1[i]<-gsub("D","[2]H",formula1[i]);
      ende1<-nchar(formula1[i]);
      element1<-c();
      number1<-c();
      ##############################################################################
      # on formula1 (deduct)
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
      # merge ! ####################################################################
      both<-unique(c(element1,element2));
      counts<-c()
      for(i in 1:length(both)){
        if(any(element1==both[i])){
          it1<-c(number1[element1==both[i]]);
        }else{
          it1<-c(0)
        }
        if(any(element2==both[i])){
          it2<-c(number2[element2==both[i]]);
        }else{
          it2<-c(0)
        }
        counts<-c(counts,it1+it2);
      }
      formula_all<-""
      for(i in 1:length(both)){
        formula_all<-paste(formula_all,both[i],counts[i],sep="")
      }
  formulas<-c(formulas,formula_all);
  }
  return(formulas)
  ##############################################################################

}
