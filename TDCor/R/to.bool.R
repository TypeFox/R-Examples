to.bool <-
function(rd,thr)

{rd_bool=norm.data(rd)

 rd_bool[rd_bool<thr]=0

 return(sign(rd_bool))}
