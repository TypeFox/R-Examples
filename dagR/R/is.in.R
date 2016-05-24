is.in <-
function(x,c=NULL)
{ # used internally by eval.paths();
  # only for numeric arguments!!!
  i<-0;
  rv<-FALSE;
  while(i<length(c))
  {
    i<-i+1;
    if(identical(x,as.numeric(c[i]))==TRUE)
    {
      rv<-TRUE;
      i<-length(c);
    }
  }
  is.in<-rv;
}

