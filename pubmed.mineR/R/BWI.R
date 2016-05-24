BWI=
function(current,previous,n,N){ a = (sum(tdm_for_lsa(current,n))+1)/(sum(tdm_for_lsa(current,N))+1); b = (sum(tdm_for_lsa(previous,n))+1)/(sum(tdm_for_lsa(previous,N))+1); c = sum(tdm_for_lsa(current,n)); if (c != 0) {d = log(c)*(a/b); return(d)} else return(0)}
