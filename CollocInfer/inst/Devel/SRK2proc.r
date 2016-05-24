####### Squared Process Error ######################

#  this file contains functions for evaluating the total SRK2 for ODE
#  and its derivative values:

#    fn      = SRK2proc
#    dfdc    = dSRK2proc.dc
#    dfdp    = dSRK2proc.dp
#    d2fdc2  = d2SRK2proc.dc2
#    d2fdcdp = d2SRK2proc.dcdp

#  It uses function make.SRK2lik but substitutes derivative value
#  for data value and quadrature point values for time values

################################################################################

make.SRK2proc <- function()
{

    ##############################################################################

    SRK2proc <- function(coefs,bvals,pars,more)
    {
       names(pars) = more$parnames
       vals = SRK2fns()$SRK2common(coefs,bvals,pars,more)

       f = with(vals,{
         f = sum( w1*( (devals2-devals1) - diag(h/4)%*%fdevals1 - diag(h/2)%*%fdevals1a - diag(h/4)%*%fdevals1b)^2 ) +
            sum( w3*( (devals1a-devals1) - diag(2*h/3)%*%fdevals1)^2 )
       })

       return( f )

    }

    ##############################################################################


    dSRK2proc.dc <- function(coefs,bvals,pars,more)
    {
       names(pars) = more$parnames
       vals1 = SRK2fns()$SRK2bvals(bvals)
       vals2 = SRK2fns()$dSRK2dx(coefs,bvals,pars,more,return.diffs=TRUE)

       g = with(c(vals1,vals2),{
           g11 = array(0,c(nrow(diffs1),ncol(diffs1)))
           g21 = array(0,c(nrow(diffs1),ncol(diffs1)))
           g12 = array(0,c(nrow(diffs2),ncol(diffs2)))
           g22 = array(0,c(nrow(diffs2),ncol(diffs2)))
    
           for(i in 1:ncol(g11)){
              g11[,i] = apply(dfdx11[,,i]*diffs1,1,sum)
              g12[,i] = apply(dfdx12[,,i]*diffs1,1,sum)
              g21[,i] = apply(dfdx21[,,i]*diffs2,1,sum)
              g22[,i] = apply(dfdx22[,,i]*diffs2,1,sum)
           }
    
    
           g = 2*(t(bvals1)%*%(g11+g21)+t(bvals1a)%*%(g12+g22) + t(bvals2)%*%diffs1)
       })

       return(g)
    }



    ##############################################################################

    dSRK2proc.dp <- function(coefs,bvals,pars,more)
    {
      vals = SRK2fns()$dSRK2dp(coefs,bvals,pars,more,return.diffs=TRUE)

      g = with(vals,{
          g = matrix(0,dim(dfdpvals1)[1],dim(dfdpvals1)[3])
    
          for(i in 1:ncol(g)){
            g[,i] = apply(dfdp1[,,i]*diffs1,1,sum) + apply(dfdp2[,,i]*diffs2,1,sum)
          }
    
          g = 2*apply(g,2,sum)
      })

      return(g)
    }




    ##############################################################################

    d2SRK2proc.dc2 <- function(coefs,bvals,pars,more)
    {
       names(pars) = more$parnames
       vals1 = SRK2fns()$SRK2bvals(bvals)
       vals2 = SRK2fns()$d2SRK2dx2(coefs,bvals,pars,more)
  
       H = with(c(vals1,vals2),{  
          d2fdx212 = 0*d2fdx211
          d2fdx222 = 0*d2fdx211
       
        
           for(k in 1:ncol(diffs1)){
             for(i in 1:ncol(diffs1)){
              for(j in 1:ncol(diffs1)){
                  d2fdx111[,i,j,k] = d2fdx111[,i,j,k]*diffs1[,i] + w1*dfdx11[,i,j]*dfdx11[,i,k]
                  d2fdx122[,i,j,k] = d2fdx122[,i,j,k]*diffs1[,i] + w1*dfdx12[,i,j]*dfdx12[,i,k]
                  d2fdx112[,i,j,k] = d2fdx112[,i,j,k]*diffs1[,i] + w1*dfdx11[,i,j]*dfdx12[,i,k]
                  d2fdx211[,i,j,k] = d2fdx211[,i,j,k]*diffs2[,i] + w3*dfdx21[,i,j]*dfdx21[,i,k]
                  d2fdx212[,i,j,k] = w3*dfdx21[,i,j]*dfdx22[,i,k]
                  d2fdx222[,i,j,k] = w3*dfdx22[,i,j]*dfdx22[,i,k]
               }
             }
           }
    
           H = list(len=dim(bvals$bvals)[2])
            for(i in 1:dim(diffs1)[2]){
              H[[i]] = list(len=dim(diffs1)[2])
              for(j in 1:dim(diffs1)[2]){
                  H[[i]][[j]] = t(bvals1)%*%diag(as.vector(apply(d2fdx111[,,i,j]+d2fdx211[,,i,j],1,sum)))%*%bvals1 +
                       t(bvals1a)%*%diag(apply(d2fdx122[,,i,j]+d2fdx222[,,i,j],1,sum))%*%bvals1a +
                       t(bvals1)%*%diag(apply(d2fdx112[,,i,j]+d2fdx212[,,i,j],1,sum))%*%bvals1a +
                       t(bvals1a)%*%diag(apply(d2fdx112[,,j,i]+d2fdx212[,,j,i],1,sum))%*%bvals1 +
                       t(bvals2)%*%diag(w1*dfdx11[,i,j])%*%bvals1 +
                       t(bvals1)%*%diag(w1*dfdx11[,j,i])%*%bvals2 +
                       t(bvals2)%*%diag(w1*dfdx12[,i,j])%*%bvals1a +
                       t(bvals1a)%*%diag(w1*dfdx12[,j,i])%*%bvals2
              }
              H[[i]][[i]] = H[[i]][[i]] + t(bvals2)%*%diag(w1)%*%bvals2
            }

            H = 2*blocks2mat(H)
        })

        return(H)
    }

    ##############################################################################

    d2SRK2proc.dcdp <- function(coefs,bvals,pars,more)
    {
       names(pars) = more$parnames
       vals1 = SRK2fns()$SRK2bvals(bvals)
       vals2 = SRK2fns()$d2SRK2dxdp(coefs,bvals,pars,more)
      
       H = with(c(vals1,vals2),{

           for(k in 1:length(pars)){
             for(i in 1:ncol(diffs1)){
              for(j in 1:ncol(diffs1)){
                  d2fdx11[,i,j,k] = d2fdx11[,i,j,k]*diffs1[,i] + w1*dfdx11[,i,j]*dfdp1[,i,k]
                  d2fdx21[,i,j,k] = d2fdx21[,i,j,k]*diffs2[,i] + w3*dfdx21[,i,j]*dfdp2[,i,k]
                  d2fdx12[,i,j,k] = d2fdx12[,i,j,k]*diffs1[,i] + w1*dfdx12[,i,j]*dfdp1[,i,k]
                  d2fdx22[,i,j,k] = d2fdx22[,i,j,k]*diffs2[,i] + w3*dfdx22[,i,j]*dfdp2[,i,k]
              }
             }
           }
    
          H = c()
 
          for(i in 1:ncol(diffs1)){
              H = rbind(H,as.matrix(t(bvals1)%*%apply(d2fdx11[,,i,]+d2fdx21[,,i,],c(1,3),sum) +
                          t(bvals1a)%*%apply(d2fdx12[,,i,]+d2fdx22[,,i,],c(1,3),sum) +
                          t(bvals2)%*%diag(w1)%*%dfdp1[,i,]))
          }
      H})

      return(2*H)
    }

    ##############################################################################

    return(
        list(
            fn = SRK2proc,
            dfdc = dSRK2proc.dc,
            dfdp = dSRK2proc.dp,
            d2fdc2 = d2SRK2proc.dc2,
            d2fdcdp = d2SRK2proc.dcdp
       )
    )
}


########################################################

SRK2fns = function()
{
    SRK2bvals = function(bvals){
       if(is.null(bvals$I)){ bvals$I = SRK2fns()$SRKindeces(nrow(bvals$bvals)) }

       bvals1 = bvals$bvals[bvals$I[,1],]
       bvals2 = bvals$bvals[bvals$I[,2],]
       bvals1a = bvals$bvals[bvals$I[,3],]

       return(list(bvals1 = bvals1,
                   bvals2 = bvals2,
                   bvals1a = bvals1a)
             )

    }

    #####################################

    SRK2indeces = function(m)
    {
      if( m%%2 == 0){ error('Number of quadrature points must be odd.') }
      I = cbind( seq(1,m-2,2),seq(1,m-2,2)+2, seq(1,m-2,2)+1)
      return(I)
    }

    #####################################

    SRK2common = function(coefs,bvals,pars,more){
       devals = as.matrix(bvals$bvals%*%coefs)  
       colnames(devals) = more$names
       if(is.null(bvals$I)){ bvals$I = SRK2indeces(nrow(devals)) }
       
        h = more$qpts[bvals$I[,2]] - more$qpts[bvals$I[,1]]

       q1 = more$qpts[bvals$I[,1]]
       q2 = more$qpts[bvals$I[,2]]
       q3 = more$qpts[bvals$I[,3]]

       devals1 = devals[bvals$I[,1],]
       devals2 = devals[bvals$I[,2],]
       devals1a = devals[bvals$I[,3],]

       fdevals1 = more$fn(q1,devals1,pars,more$more)
       fdevals1a = more$fn(q3,devals1a,pars,more$more)

       devals1b = devals1 - diag(h/3)%*%fdevals1 + diag(h)%*%fdevals1a    #  -devals1 + devals1a
       
       fdevals1b = more$fn(q3,devals1b,pars,more$more)

       w1 = more$weights[bvals$I[,1]]
       w3 = more$weights[bvals$I[,3]]
       
       diffs1 = diag(w1)%*%( (devals2-devals1) - diag(h/4)%*%fdevals1 - diag(h/2)%*%fdevals1a - diag(h/4)%*%fdevals1b)
       diffs2 = diag(w3)%*%( (devals1a-devals1) - diag(2*h/3)%*%fdevals1)       

       return(list(devals1 = devals1,
                   devals2 = devals2,
                   devals1a = devals1a,
                   devals1b = devals1b,
                   fdevals1 = fdevals1,
                   fdevals1a = fdevals1a,
                   devals1b = devals1b,
                   fdevals1b = fdevals1b,
                   diffs1 = diffs1,
                   diffs2 = diffs2,
                   h = h,
                   w1 = w1,
                   w3 = w3,
                   q1 = q1,
                   q2 = q2,
                   q3 = q3)
             )
    }

    #####################################

    dSRK2dx = function(coefs,bvals,pars,more,return.diffs=FALSE)
    {
       names(pars) = more$parnames
       vals = SRK2fns()$SRK2common(coefs,bvals,pars,more)
 #      attach(vals)

       return.list = with(vals,{
           dfdxvals1 = more$dfdx(q1,devals1,pars,more$more)
           dfdxvals1a = more$dfdx(q3,devals1a,pars,more$more)
           dfdxvals1b = more$dfdx(q3,devals1b,pars,more$more)
    
           dfdx11 = 0*dfdxvals1
           dfdx12 = 0*dfdxvals1
           dfdx21 = 0*dfdxvals1
           dfdx22 = 0*dfdxvals1
            
           for(i in 1:dim(dfdx11)[2]){
    
            for(j in 1:dim(dfdx11)[3]){
              dfdx11[,i,j] =  - h/4*dfdxvals1[,i,j] + h^2/12*apply(dfdxvals1b[,i,]*dfdxvals1[,,j],1,sum) - 
                  h/4*dfdxvals1b[,i,j]
                  
              dfdx21[,i,j] =  - 2*h/3*dfdxvals1[,i,j]
    
              dfdx12[,i,j] = - h/2*dfdxvals1a[,i,j] - h^2/4*apply(dfdxvals1b[,i,]*dfdxvals1a[,,j],1,sum)
              }
    
              dfdx11[,i,i] = dfdx11[,i,i] - 1
              dfdx21[,i,i] = dfdx21[,i,i] - 1
              dfdx22[,i,i] = 1
           }
    
          if(return.diffs){
            return.list = list(dfdx11=dfdx11,dfdx12=dfdx12,dfdx21=dfdx21,dfdx22=dfdx22,
                        dfdxvals1=dfdxvals1,dfdxvals1a=dfdxvals1a,dfdxvals1b=dfdxvals1b,
                        diffs1=diffs1,diffs2=diffs2)
          }
          else{
            return.list = list(dfdx11=dfdx11,dfdx12=dfdx12,dfdx21=dfdx21,dfdx22=dfdx22,
                        dfdxvals1=dfdxvals1,dfdxvals1a=dfdxvals1a,dfdxvals1b=dfdxvals1b)
          }
      return.list})
      
      return(return.list)
    }

    #####################################

    dSRK2dp = function(coefs,bvals,pars,more,return.diffs=FALSE)
    {
       names(pars) = more$parnames
       vals = SRK2fns()$SRK2common(coefs,bvals,pars,more)
      
       return.list=with(vals,{
           fdevals1b = more$fn(q3,devals1b,pars,more$more)
    
           dfdpvals1 = more$dfdp(q1,devals1,pars,more$more)
           dfdpvals1a = more$dfdp(q3,devals1a,pars,more$more)
           dfdpvals1b = more$dfdp(q3,devals1b,pars,more$more)
    
           dfdxvals1b = more$dfdx(q3,devals1b,pars,more$more)
    
           dfdp1 = 0*dfdpvals1
           dfdp2 = 0*dfdpvals1
    
           for(i in 1:dim(dfdpvals1)[2]){
            for(j in 1:dim(dfdpvals1)[3]){
              dfdp1[,i,j] = -h/4*dfdpvals1[,i,j]-h/2*dfdpvals1a[,i,j] - h/4*dfdpvals1b[,i,j] -
                h^2/4*apply(dfdxvals1b[,i,]*(dfdpvals1a[,,j]-dfdpvals1[,,j]/3),1,sum)
              dfdp2[,i,j] = -2*h/3*dfdpvals1[,i,j]
             }
           }
           
           
           
          if(return.diffs){
            return.list = list(dfdp1=dfdp1,dfdp2=dfdp2,
              dfdpvals1=dfdpvals1,dfdpvals1a=dfdpvals1a,dfdpvals1b=dfdpvals1b,
                        diffs1=diffs1,diffs2=diffs2)
          }
          else{
            return.list = list(dfdp1=dfdp1,dfdp2=dfdp2,
              dfdpvals1=dfdpvals1,dfdpvals1a=dfdpvals1a,dfdpvals1b=dfdpvals1b)
          }
          
      return.list})

      return(return.list)
      
    }

    #####################################
    
    d2SRK2dxdp = function(coefs,bvals,pars,more)
    {  
      names(pars) = more$parnames
      vals1 = SRK2fns()$SRK2common(coefs,bvals,pars,more)
      vals2 = SRK2fns()$dSRK2dx(coefs,bvals,pars,more,return.diffs=FALSE)
      vals3 = SRK2fns()$dSRK2dp(coefs,bvals,pars,more)
      
      return.list = with(c(vals1,vals2,vals3),{    
    
           d2fdxdpvals1 = more$d2fdxdp(q1,devals1,pars,more)
           d2fdxdpvals1a = more$d2fdxdp(q3,devals1a,pars,more)
           d2fdxdpvals1b = more$d2fdxdp(q3,devals1b,pars,more)
    
           d2fdx2vals1b =  more$d2fdx2(q3,devals1b,pars,more)
    
           d2fdx11 = 0*d2fdxdpvals1;
           d2fdx21 = 0*d2fdxdpvals1;
           d2fdx12 = 0*d2fdxdpvals1a;
           d2fdx22 = 0*d2fdxdpvals1a;
    
           for(k in 1:length(pars)){
             for(i in 1:ncol(devals1)){
               for(j in 1:ncol(devals1)){
                 d2fdx11[,i,j,k] = -h/4*d2fdxdpvals1[,i,j,k] + h^2/12*apply(d2fdxdpvals1b[,i,,k]*dfdxvals1[,,j],1,sum) +
                             h^2/12*apply(dfdxvals1b[,i,]*d2fdxdpvals1[,,j,k],1,sum) - h/4*d2fdxdpvals1b[,i,j,k] -
                             h^2/4*apply(d2fdx2vals1b[,i,,j]*(dfdpvals1a[,,k] - dfdpvals1[,,k]/3),1,sum)
    
                 d2fdx21[,i,j,k] = -2*h/3*d2fdxdpvals1[,i,j,k]
    
                 d2fdx12[,i,j,k] = -h/2*d2fdxdpvals1a[,i,j,k] - h^2/4*apply(d2fdxdpvals1b[,i,,k]*dfdxvals1a[,,j],1,sum) -
                             h^2/4*apply(dfdxvals1b[,i,]*d2fdxdpvals1a[,,j,k],1,sum)
    
                 for(l in 1:ncol(devals1)){
                    d2fdx11[,i,j,k] = d2fdx11[,i,j,k] + 
                      h^3/12*apply(d2fdx2vals1b[,i,l,]*(dfdpvals1a[,,k] - dfdpvals1[,,k]/3),1,sum)*dfdxvals1[,l,j]
                    d2fdx12[,i,j,k] = d2fdx12[,i,j,k] - 
                      h^3/4*apply(d2fdx2vals1b[,i,l,]*(dfdpvals1a[,,k] - dfdpvals1[,,k]/3),1,sum)*dfdxvals1a[,l,j]
                 }
               }
    
             }
           }
           
           return.list=list(d2fdx11 = d2fdx11,
                            d2fdx12 = d2fdx12,
                            d2fdx21 = d2fdx21,
                            d2fdx22 = d2fdx22,
                            dfdx11 = dfdx11,
                            dfdx12 = dfdx12,
                            dfdx21 = dfdx21,
                            dfdx22 = dfdx22,
                            dfdp1 = dfdp1,
                            dfdp2 = dfdp2,
                            diffs1 = diffs1,
                            diffs2 = diffs2,
                            w1 = w1,
                            w3 = w3)
      return.list})
                            
      return(return.list)     
           
    }    

    #####################################
    
    d2SRK2dx2 = function(coefs,bvals,pars,more)
    {
      names(pars) = more$parnames
      vals1 = SRK2fns()$SRK2common(coefs,bvals,pars,more)
      vals2 = SRK2fns()$dSRK2dx(coefs,bvals,pars,more,return.diffs=FALSE) 
       
      return.list = with(c(vals1,vals2),{        
           d2fdx2vals1 =  more$d2fdx2(q3,devals1,pars,more)
           d2fdx2vals1a =  more$d2fdx2(q3,devals1a,pars,more)
           d2fdx2vals1b =  more$d2fdx2(q3,devals1b,pars,more)
       
           d2fdx111 = array(0,c(nrow(devals1),ncol(devals1),ncol(devals1),ncol(devals1)))
           d2fdx112 = array(0,c(nrow(devals1),ncol(devals1),ncol(devals1),ncol(devals1)))
           d2fdx122 = array(0,c(nrow(devals1),ncol(devals1),ncol(devals1),ncol(devals1)))
    
           d2fdx211 = array(0,c(nrow(devals1),ncol(devals1),ncol(devals1),ncol(devals1)))
    
    
           for(i in 1:ncol(devals1)){
            for(j in 1:ncol(devals1)){
              for(k in 1:ncol(devals1)){
                d2fdx111[,i,j,k] = -h/4*d2fdx2vals1[,i,j,k] + h^2/12*apply(dfdxvals1b[,i,]*d2fdx2vals1[,,j,k],1,sum) +
                    h^2/12*(apply(d2fdx2vals1b[,i,,k]*dfdxvals1[,,j],1,sum) + apply(d2fdx2vals1b[,i,j,]*dfdxvals1[,,k],1,sum))-
                    h/4*d2fdx2vals1b[,i,j,k]
                d2fdx211[,i,j,k] = -2*h/3*d2fdx2vals1[,i,j,k]
                
                d2fdx122[,i,j,k] = -h/2*d2fdx2vals1a[,i,j,k] - h^2/4*apply(dfdxvals1b[,i,]*d2fdx2vals1a[,,j,k],1,sum)
    
                d2fdx112[,i,j,k] = - h^2/4*apply(d2fdx2vals1b[,i,j,]*dfdxvals1a[,,k],1,sum)
    
                for(l in 1:ncol(devals1)){
                  d2fdx111[,i,j,k] = d2fdx111[,i,j,k] - h^3/36*apply(d2fdx2vals1b[,i,,l]*dfdxvals1[,,j],1,sum)*dfdxvals1[,l,k]
                  d2fdx122[,i,j,k] = d2fdx122[,i,j,k] - h^3/4*apply(d2fdx2vals1b[,i,,l]*dfdxvals1a[,,j],1,sum)*dfdxvals1a[,l,k]
                  d2fdx112[,i,j,k] = d2fdx112[,i,j,k] + h^3/12*apply(d2fdx2vals1b[,i,,l]*dfdxvals1[,,j],1,sum)*dfdxvals1a[,l,k]
                }
              }
            }
          }    
           return.list=list(d2fdx111 = d2fdx111,
                            d2fdx112 = d2fdx112,
                            d2fdx122 = d2fdx122,
                            d2fdx211 = d2fdx211,
                            dfdx11 = dfdx11,
                            dfdx12 = dfdx12,
                            dfdx21 = dfdx21,
                            dfdx22 = dfdx22,
                            diffs1 = diffs1,
                            diffs2 = diffs2,
                            w1 = w1,
                            w3 = w3)    
    
      return.list})
      
      return(return.list)
    
    }
    
    return(list(SRK2bvals=SRK2bvals,
                SRK2common=SRK2common,
                SRK2indeces=SRK2indeces,
                dSRK2dx=dSRK2dx,
                dSRK2dp=dSRK2dp,
                d2SRK2dx2=d2SRK2dx2,
                d2SRK2dxdp=d2SRK2dxdp
                )
    )
}