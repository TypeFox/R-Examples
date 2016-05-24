"shannonbioboot" <-
function(data1,B=1000){

boot(data1,function(data1,i) shannonbio(data1[i,]),R=B)

}

