`REFLECT` <-
function(A)
{
    #  /* Reflect to lower hemisphere */
   h = Preflect(A$az, A$dip );
   A$az=h$az; A$dip=h$dip;
   A = tocartL(A);
    return(A);
}

