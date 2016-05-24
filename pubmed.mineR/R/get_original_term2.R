get_original_term2 = function(x,y)
{
dummy_cluster = list(list(Words=x,Frequencies=0))
res = get_original_term(y,dummy_cluster)
return(res);
}
