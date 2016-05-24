`part` <-
function(ms1,ms2){

             par1<-list(setdiff(ms1,ms2))

             par2<-union(par1, list(intersect(ms1,ms2)))

             par3<-union(par2, list(setdiff(ms2,ms1)))

             par3

}

