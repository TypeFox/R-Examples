ht_mal_int <-
function(lower, upper, g, d, smoking)
{
 INT <- HT_mal(upper, g, d, smoking) - HT_mal(lower, g, d, smoking)
 return(INT)
}
