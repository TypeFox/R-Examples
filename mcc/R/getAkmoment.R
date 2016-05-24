getAkmoment <-
function(x,y){
 m=dim(x)[1]
 n=dim(x)[2]
 sumx=rowSums(x)
 sumx2=rowSums(x^2)
 sumx3=rowSums(x^3)
 sumx4=rowSums(x^4)
 sumy=sum(y)
 sumy2=sum(y^2)
 sumy3=sum(y^3)
 sumy4=sum(y^4)
# do moment1
 final1=sumx*sumy/n
# do moment2
 final2=sumx2*sumy2/n+(sumx^2-sumx2)*(sumy^2-sumy2)/(n^2-n)
# do moment3
 size1=n
 size2=3*n*(n-1)
 size3=n^3-size2-size1
 alt.term1x=sumx3
 alt.term2x=3*(sumx2*sumx-sumx3)
 alt.term3x=sumx^3-alt.term2x-alt.term1x
 alt.term1y=sumy3
 alt.term2y=3*(sumy2*sumy-sumy3)
 alt.term3y=sumy^3-alt.term2y-alt.term1y
 if (n==2){final3=alt.term1x*alt.term1y/size1+alt.term2x*alt.term2y/size2}
 if (n>2){final3=alt.term1x*alt.term1y/size1+alt.term2x*alt.term2y/size2+alt.term3x*alt.term3y/size3}
# do moment4
 size1=n
 size2=4*n*(n-1)
 size3=3*n*(n-1)
 size4=6*n*(n-1)*(n-2)
 size5=n^4-6*n^3+11*n^2-6*n
 sumy=sum(y)
 sumy2=sum(y^2)
 sumy3=sum(y^3)
 sumy4=sum(y^4)
 alt.term1x=sumx4
 alt.term2x=4*(sumx3*sumx-alt.term1x)
 alt.term3x=3*(sumx2^2-alt.term1x)
 alt.term4x=6*sumx2*sumx^2-12*sumx3*sumx-6*sumx2^2+12*sumx4
 alt.term5x=sumx^4-alt.term4x-alt.term3x-alt.term2x-alt.term1x
 alt.term1y=sumy4
 alt.term2y=4*(sumy3*sumy-alt.term1y)
 alt.term3y=3*(sumy2^2-alt.term1y)
 alt.term4y=6*sumy2*sumy^2-12*sumy3*sumy-6*sumy2^2+12*sumy4
 alt.term5y=sumy^4-alt.term4y-alt.term3y-alt.term2y-alt.term1y
 if (n==2){
 final4=alt.term1x*alt.term1y/size1+alt.term2x*alt.term2y/size2+alt.term3x*alt.term3y/size3}
 if (n==3){
 final4=alt.term1x*alt.term1y/size1+alt.term2x*alt.term2y/size2+alt.term3x*alt.term3y/size3+
       alt.term4x*alt.term4y/size4}
 if (n>3){
 final4=alt.term1x*alt.term1y/size1+alt.term2x*alt.term2y/size2+alt.term3x*alt.term3y/size3+
       alt.term4x*alt.term4y/size4+alt.term5x*alt.term5y/size5} 
 sigma2=final2-final1^2
 myskew=(final3-3*final1*sigma2-final1^3)/sigma2^1.5
 mu4=final4-4*final1*final3+6*final1^2*final2-3*final1^4
 mykurt=mu4/sigma2^2-3
 return(list(final1=final1,final2=final2,final3=final3,final4=final4,sigma2=sigma2,myskew=myskew,mykurt=mykurt))
 }

