subroutine  dsidr0 (vmu,_
                    s, lds, nobs, nnull, y, q, ldq,_       # data
                    tol, job, limnla,_                     # job requests
                    nlaht, score, varht, c, d,_            # output
                    qraux, jpvt, wk,_                      # work arrays
                    info)                                  # error message

integer           vmu
integer           lds, nobs, nnull, ldq, job, jpvt(*), info
double precision  s(lds,*), y(*), q(ldq,*), tol, limnla(2), nlaht, score(*),_
                  varht, c(*), d(*), qraux(*), wk(*)

character*1       vmu1

if ( vmu == 1 )  vmu1 = 'v'
if ( vmu == 2 )  vmu1 = 'm'
if ( vmu == 3 )  vmu1 = 'u'

call  dsidr (vmu1, s, lds, nobs, nnull, y, q, ldq, tol, job, limnla,
             nlaht, score, varht, c, d, qraux, jpvt, wk, info)

return
end
