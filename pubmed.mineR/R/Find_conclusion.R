Find_conclusion = function(y) {
check1=lapply(y@Abstract, function(x) {SentenceToken(x)});
check2 =lapply(check1, function(x) { check3 = regexpr("conclusion", x, ignore.case=T); 
check4 = which(check3 != -1); 
if (length(check4) != 0) 
    { 
     check5 = x[check4:length(x)]; return(check5)
    } 
else 
    {
if(length(x)>3){
 check5 = x[(length(x)-3):length(x)]; return(check5)
}
    }
  })
}
