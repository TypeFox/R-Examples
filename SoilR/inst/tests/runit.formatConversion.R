#
# vim:set ff=unix expandtab ts=2 sw=2:
test.BoundFc=function(){
   tstart=0
   tend=10
   f=function(t){2*t}
   # the first two objects can be created because the formats are supported
   obj1=BoundFc(f,tstart,tend,format="Delta14C") 
   obj2=BoundFc(f,tstart,tend,format="AbsoluteFractionModern") 
   
   # but we expect trouble for the next line because "foo-bar" 
   # is not supported as format for atmospheric 14C
   checkException(new(Class="BoundFc",tstart,tend,f,format="foo-bar"))
   
   # now we test the back and forth transformation from one Format to another
   tt=seq(from=tstart,to=tend,by=(tend-tstart)/100)

   obj3=Delta14C(AbsoluteFractionModern(obj1))
   fref=getFunctionDefinition(obj3)
   checkEquals(fref(tt),f(tt))
   
   
   obj4=AbsoluteFractionModern(Delta14C(obj2))
   fref=getFunctionDefinition(obj4)
   checkEquals(fref(tt),f(tt))
   
   # we add some tests to check the conversion of C14 data
   # produce some data
   dat=f(tt)
   # interpret it as Delta14C values and check if it is converted correctly
   checkEquals(AbsoluteFractionModern_from_Delta14C(dat),dat/1000+1)
   # now interpret the same data as Absolute Fraction Modern  values and check if it is converted correctly
   checkEquals(Delta14C_from_AbsoluteFractionModern(dat),(dat-1)*1000)
   
   # now do it with matrices
   m1=matrix(nrow=3,ncol=length(dat),dat)
   checkEquals(AbsoluteFractionModern_from_Delta14C(m1),m1/1000+1)
   checkEquals(Delta14C_from_AbsoluteFractionModern(m1),(m1-1)*1000)
   
}

