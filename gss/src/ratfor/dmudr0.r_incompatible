subroutine  dmudr0 (vmu,_
                    s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,_     # inputs
                    tol, init, prec, maxite,_                       # tune para
                    theta, nlaht, score, varht, c, d,_              # outputs
                    wk, info)

integer           vmu
integer           lds, nobs, nnull, ldqr, ldqc, nq, init, maxite,_
                  info

double precision  s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec,_
                  theta(*), nlaht, score, varht, c(*), d(*),_
                  wk(*)

character*1       vmu1

if ( vmu == 1 )  vmu1 = 'v'
if ( vmu == 2 )  vmu1 = 'm'
if ( vmu == 3 )  vmu1 = 'u'

call  dmudr (vmu1, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol,
             init, prec, maxite, theta, nlaht, score, varht, c, d, wk,
             info)

return
end
