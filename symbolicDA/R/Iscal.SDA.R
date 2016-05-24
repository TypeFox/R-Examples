iscal.SDA<-function(x,d=2,calculateDist=FALSE){

if(.is.symbolic(x)){
x<-SO2Simple(x)
}
dist_min<-x[,,1]
dist_max<-x[,,2]

z<-ncol(dist_min)
inter<-interscal.SDA(x,d,calculateDist)$xprim
k_dolne<-inter[,,1]
k_gorne<-inter[,,2]
x<-(k_dolne+k_gorne)/2
x<-array(0,c(z,d))
for(i in 1:z){
x[i,1]<-(k_dolne[i,1])/2
x[i,2]<-(k_dolne[i,2])/2
}


losoweR<-sample((-z):(z),z*d)
r<-array(losoweR,c(z,d))
r<-(k_gorne-k_dolne)


losoweX<-1:z*d/max(k_gorne)
x<-array(losoweX,c(z,d))

for(i in 1:z){
x[i,1]<-(k_dolne[i,1]+k_gorne[i,1])/2
x[i,2]<-(k_dolne[i,2]+k_gorne[i,2])/2
r[i,1]<-(k_gorne[i,1]-k_dolne[i,1])
r[i,2]<-(k_gorne[i,2]-k_dolne[i,2])
}


kryteriumStopu<-10^(-6)
maxIteracji<-200
STRESSSym<-0
STRESSSymOld<-0
liczbaIteracji<-1
d_min<-array(0,c(z,z))
d_max<-array(0,c(z,z))
wagi<-1/(z^2)
while ((STRESSSymOld-STRESSSym>kryteriumStopu && liczbaIteracji<=maxIteracji) || (liczbaIteracji<=2))
{
   for (i in 1:z)
   {
      for (j in 1:z)
      {
	
	d_max[i,j]<-d_max[i,j]<- sqrt( ( abs(x[i,1]-x[j,1])+(r[i,1]+r[j,1]) )^2 + ( abs(x[i,2]-x[j,2])+(r[i,2]+r[j,2]) )^2 )

	d_min[i,j]<-d_min[i,j]<- sqrt( (max (0, abs(x[i,1]-x[j,1])+(r[i,1]+r[j,1]) ) )^2 + (max (0, abs(x[i,2]-x[j,2])+(r[i,2]+r[j,2]) ) )^2 )

      }
   }
alfa1_d<-array(1,c(z,z,d))
alfa2_d<-array(1,c(z,z,d))
alfa3_d<-array(1,c(z,z))
alfa4_d<-array(1,c(z,z,d))
alfa5_d<-array(1,c(z,z,d))
beta1_d<-array(1,c(z,z,d))
beta2_d<-array(1,c(z,z,d))
beta3_d<-array(1,c(z,z,d))
beta4_d<-array(1,c(z,z,d))
beta5_d<-array(1,c(z,z,d))

   for (i in 1:z)
   {
      for (j in 1:z)
      {
         for (s in 1:d)
         {
            if(i!=j)
            {
            alfa1_d[i,j,s]<-wagi*( 1+(r[i,s]+r[j,s])/abs(x[i,s]-x[j,s]) )
            alfa2_d[i,j,s]<-wagi*( abs(x[i,s]-x[j,s])+(r[i,s]+r[j,s]) )/r[i,s]
            alfa4_d[i,j,s]<-2*wagi*( 1+r[j,s]/r[i,s] )            
            }
            else
            {
            alfa1_d[i,j,s]<-wagi
            alfa2_d[i,j,s]<-wagi
            alfa4_d[i,j,s]<-wagi
            }

            alfa3_d[i,j]<-2*wagi


            if ( abs(x[i,s]-x[j,s])>=(r[i,s]+r[j,s]) )
            {

               alfa5_d[i,j,s]<-( wagi*as.numeric(dist_min[i,j])*max( 0, (abs(x[i,s]-x[j,s])-r[i,s]-r[j,s]) ) / (r[i,s]*d_min[i,j]))
            }
            else
            {
               alfa5_d[i,j,s]<-0
            }

            if ( abs(x[i,s]-x[j,s])>0 && d_max[i,j]>0 )
            {
               beta1_d[i,j,s]<-( wagi*as.numeric(dist_max[i,j])*( abs(x[i,s]-x[j,s])+r[i,s]+r[j,s] ) ) / ( abs(x[i,s]-x[j,s])*d_max[i,j] )
            }
            else
            {
               beta1_d[i,j,s]<-0
            }

            if (d_max[i,j]>0)
            {
               beta2_d[i,j,s]<-( wagi*as.numeric(dist_max[i,j])*( abs(x[i,s]-x[j,s])+r[i,s]+r[j,s] ) ) / as.numeric(d_max[i,j])
            }
            else
            {
               beta2_d[i,j,s]<-0
            }

            if ( abs(x[i,s]-x[j,s])>=(r[i,s]+r[j,s]) )
            {
               beta4_d[i,j,s]<-wagi*( abs(x[i,s]-x[j,s])+r[i,s]+r[j,s] )
            }
            else
            {
               beta4_d[i,j,s]<-2*wagi*(r[i,s]+r[j,s])
            }

            if ( abs(x[i,s]-x[j,s])>=(r[i,s]+r[j,s]) && abs(x[i,s]-x[j,s])>0 )
            {
               beta3_d[i,j,s]<-( wagi*( abs(x[i,s]-x[j,s])+r[i,s]+r[j,s] ) ) / ( abs(x[i,s]-x[j,s]) )
            }
            
            if ( abs(x[i,s]-x[j,s])<(r[i,s]+r[j,s]) && abs(x[i,s]-x[j,s])>0 )
            {
               beta3_d[i,j,s]<-2*wagi
            }

            if ( abs(x[i,s]-x[j,s])==0 )
            {
               beta3_d[i,j,s]<-0
            }

         }
      }
   }

A1_d<-array(0, c(z,z,d))
A2_d<-array(0, c(z,z,d))
B1_d<-array(0, c(z,z,d))
B2_d<-array(0, c(1,z,d))
A1_plus_d<-array(0, c(z,z,d))
wektor<-array(1, c(1,z))
wektor1<-as.vector(wektor)
      for (i in 1:z)
      {
         for (j in 1:z)
         {
            for (s in 1:d)
            {
            
            if(!(i==j))
            {
               A1_d[i,j,s]<- -(alfa1_d[i,j,s]+alfa3_d[i,j])
               A1_d[i,i,s]<- A1_d[i,i,s]-A1_d[i,j,s]
               
               B1_d[i,j,s]<- -(beta1_d[i,j,s]+beta3_d[i,j,s]+beta5_d[i,j,s])
               B1_d[i,i,s]<- B1_d[i,i,s]-B1_d[i,j,s]

               B2_d[1,i,s]<- B2_d[1,i,s]+beta2_d[i,j,s]+beta4_d[i,j,s]
               
               A2_d[i,i,s]<- A2_d[i,i,s]+alfa2_d[i,j,s]+alfa4_d[i,j,s]+alfa5_d[i,j,s]

               A1_d[is.nan(A1_d)]<-1
               A2_d[is.nan(A2_d)]<-1
               B1_d[is.nan(B1_d)]<-1
               B2_d[is.nan(B2_d)]<-1  

            }

               

               r[i,s]<- B2_d[1,i,s]/A2_d[i,i,s]

            }
         }
      }
     A1_plus_d[,,s]<- solve(A1_d[,,s]+z^(-1)*t(wektor)%*%wektor)-z^(-1)*t(wektor)%*%wektor
     x[,s]<-( A1_plus_d[,,s] %*% B1_d[,,s] %*% x[,s] )

liczbaIteracji<-liczbaIteracji+1
STRESSSymOld<-STRESSSym
STRESSSym<-0
  for (i in 1:(z-1))
   {
      for (j in (i+1):z)
      {
      d_max[i,j]<-sqrt(abs(x[i,1]-x[j,1])+(r[i,1]+r[j,1]))^2+((abs(x[i,2]-x[j,2])+(r[i,2]+r[j,2]))^2)
      d_min[i,j]<-sqrt(max(0,abs(x[i,1]-x[j,1])-r[i,1]-r[j,1])^2+max(0,abs(x[i,2]-x[j,2])-(r[i,2]+r[j,2])^2))



	STRESSSym<-STRESSSym+(wagi*(as.numeric(dist_max[i,j])-d_max[i,j])^2+wagi*(as.numeric(dist_min[i,j])-d_min[i,j])^2)
  
   }

if(sum(is.nan(A2_d))!=0)
{
A2_d[is.nan(A2_d)]<-0
}
}
    mian1<-0
    mian2<-0
    
    for (k in 1:(z-1))
    {
      for (l in (k+1):z)
      {
        mian1<-mian1+(wagi*dist_min[k,l]^2)
        mian2<-mian2+(wagi*dist_max[k,l]^2)
      }
    }
    STRESSSym<-STRESSSym/(mian1+mian2)
}

wymiar1<-array(,c(z,d))
wymiar2<-array(,c(z,d))
	for (i in 1:z)
	{
			
			wymiar1[i,1]<-x[i,1]-0.5*r[i,1]
			wymiar2[i,1]<-x[i,2]-0.5*r[i,2]
			
			wymiar1[i,2]<-x[i,1]+0.5*r[i,1]
			wymiar2[i,2]<-x[i,2]+0.5*r[i,2]
	}

xprim<-array(,c(z,d,2))
for (i in 1:z)
{
	for (j in 1:d){
    xprim[i,j,1]<-wymiar1[i,j]
    xprim[i,j,2]<-wymiar2[i,j]
	}
}

     res<-list(xprim=xprim,stress.sym<-STRESSSym)
      res
 
}
