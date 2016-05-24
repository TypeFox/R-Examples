
"sgpls" <-
function( x, y, K, eta, scale.x=TRUE,
        eps=1e-5, denom.eps=1e-20, zero.eps=1e-5, maxstep=100,
        br=TRUE, ftype='iden' )
{
    y.u <- unique(y)
    ngroups <- length(y.u)
    
    if ( ngroups >= 2 )
    {
        if ( ngroups == 2 )
        {
            object <- sgpls.binary( x=x, y=y, K=K, eta=eta, scale.x=scale.x,
                eps=eps, denom.eps=denom.eps, zero.eps=zero.eps, maxstep=maxstep,
                br=br, ftype=ftype )
        }
        if ( ngroups > 2 )
        {
            object <- sgpls.multi( x=x, y=y, K=K, eta=eta, scale.x=scale.x,
                eps=eps, denom.eps=denom.eps, zero.eps=zero.eps, maxstep=maxstep,
                br=br, ftype=ftype )
        }
    } else
    {
        stop( "Check the response vector!" )
    }
    
    return( object )
}
