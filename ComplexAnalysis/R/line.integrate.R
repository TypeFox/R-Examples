line.integrate <-
function(f,a,b,complex=TRUE){
result<-NULL
if(Im(a)==0 & Im(b)==0 & complex!=TRUE){result<-integrate(f,Re(a),Re(b))$value}else{
Re.func<-function(t){Re(f(Re(a)+(Re(b)-Re(a))*t+1i*(Im(a)+(Im(b)-Im(a))*t))*(Re(b)-Re(a)+1i*(Im(b)-Im(a))))}
Im.func<-function(t){Im(f(Re(a)+(Re(b)-Re(a))*t+1i*(Im(a)+(Im(b)-Im(a))*t))*(Re(b)-Re(a)+1i*(Im(b)-Im(a))))}
result<-integrate(Re.func,0,1)$value+integrate(Im.func,0,1)$value*1i}
return(result)}
