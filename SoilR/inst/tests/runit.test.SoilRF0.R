#
# vim:set ff=unix expandtab ts=2 sw=2:
test.ConstFc=function(){
   # the first two objects can be created because the formats are supported
   obj1=ConstFc(c(0,0),"Delta14C")
   obj2=ConstFc(c(0,0),"AbsoluteFractionModern")

   # but we expect trouble for the next line because "foo-bar" 
   # is not supported as format for atmospheric 14C
   checkException( ConstFc(c(0,0),"foo-bar"))
   
   # check conversion back and forth
   obj3=Delta14C(AbsoluteFractionModern(obj1))

   checkEquals(getFormat(obj3),getFormat(obj1))
   checkEquals(getValues(obj3),getValues(obj1))
   
	       
   obj4=AbsoluteFractionModern(Delta14C(obj2))

   checkEquals(getFormat(obj4),getFormat(obj2))
   checkEquals(getValues(obj4),getValues(obj2))

}
