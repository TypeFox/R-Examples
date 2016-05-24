hstar = function(abdord,wtord,const,n,y){
    m = length(abdord)
    ic = 0
    icm = 0
    hstar = 0
    while(ic < 1){
       icm = icm + 1
       if(abdord[icm] <= y){
          hstar = hstar + wtord[icm]
       } else {
          ic = 2
       }
       if(icm == m){ic=2}
     }
     hstar = hstar/((n-1)*const)
     hstar
}
