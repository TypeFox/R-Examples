GQD.remove <-
function()
{
   namess=c('G0','G1','G2','Q0','Q1','Q2','a00','a10','a20','a01','a02','a11',
             'b00','b10','b20','b01','b02','b11',
             'c00','c10','c20','c01','c02','c11',
             'd00','d10','d20','d01','d02','d11',
             'e00','e10','e20','e01','e02','e11',
             'f00','f10','f20','f01','f02','f11','priors')
   func.list=rep(0,length(namess))
   obs=objects(pos=1)
   for(i in 1:length(namess))
   {
     if(sum(obs==namess[i])){func.list[i]=1}
   }
   removes=namess[which(func.list==1)]
   Info="Removed : "

   for(i in 1:length(removes))
   {
     Info=paste0(Info,' ',removes[i])
   }
   print(Info,row.names = FALSE,right=F)
   remove(list=removes,pos=1)
}
