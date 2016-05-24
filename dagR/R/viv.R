viv <-
function(v1,v2){
  v1s<-unique(na.omit(v1));
  v2s<-unique(na.omit(v2));
  if(is.numeric(v1s)==FALSE)
  { cat('v1 is empty');
    rv<-TRUE; # empty always contained in v2 !
  } else
  { if(is.numeric(v2s)==FALSE)
    { cat('v2 is empty');
      rv<-FALSE; # only empty can be in empty !
    }
    else
    { score<-0;
      for(i in 1:length(v1s))
      { score<-score+sum(v2s==(v1s[i]), na.rm=TRUE);
      }
      if(score==length(v1s)) rv<-TRUE
       else                  rv<-FALSE;
    }
  }
  return(rv);
}