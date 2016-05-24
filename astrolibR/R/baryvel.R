#source('premat.R')

baryvel = function( dje, deq) {
    dc2pi = 2*pi
    cc2pi = 2*pi
    dc1 = 1.0e0
    dcto = 2415020.0e0
    dcjul = 36525.0e0                     #days in julian year
    dcbes = 0.313e0
    dctrop = 365.24219572e0               #days in tropical year (...572 insig)
    dc1900 = 1900.0e0
    au = 1.4959787e8
    dcfel = c(1.7400353e00, 6.2833195099091e02,  5.2796e-6
             ,6.2565836e00, 6.2830194572674e02, -2.6180e-6
             ,4.7199666e00, 8.3997091449254e03, -1.9780e-5
             ,1.9636505e-1, 8.4334662911720e03, -5.6044e-5
             ,4.1547339e00, 5.2993466764997e01,  5.8845e-6
             ,4.6524223e00, 2.1354275911213e01,  5.6797e-6
             ,4.2620486e00, 7.5025342197656e00,  5.5317e-6
             ,1.4740694e00, 3.8377331909193e00,  5.6093e-6)
    dcfel = matrix(dcfel,3,8)
    dceps = c(4.093198e-1, -2.271110e-4, -2.860401e-8 )
    ccsel = c(1.675104e-2, -4.179579e-5, -1.260516e-7
    ,2.220221e-1,  2.809917e-2,  1.852532e-5
    ,1.589963e00,  3.418075e-2,  1.430200e-5
    ,2.994089e00,  2.590824e-2,  4.155840e-6
    ,8.155457e-1,  2.486352e-2,  6.836840e-6
    ,1.735614e00,  1.763719e-2,  6.370440e-6
    ,1.968564e00,  1.524020e-2, -2.517152e-6
    ,1.282417e00,  8.703393e-3,  2.289292e-5
    ,2.280820e00,  1.918010e-2,  4.484520e-6
    ,4.833473e-2,  1.641773e-4, -4.654200e-7
    ,5.589232e-2, -3.455092e-4, -7.388560e-7
    ,4.634443e-2, -2.658234e-5,  7.757000e-8
    ,8.997041e-3,  6.329728e-6, -1.939256e-9
    ,2.284178e-2, -9.941590e-5,  6.787400e-8
    ,4.350267e-2, -6.839749e-5, -2.714956e-7
    ,1.348204e-2,  1.091504e-5,  6.903760e-7
    ,3.106570e-2, -1.665665e-4, -1.590188e-7 )
    ccsel = matrix(ccsel,3,17)
    dcargs = c(5.0974222e0, -7.8604195454652e2
    ,3.9584962e0, -5.7533848094674e2
    ,1.6338070e0, -1.1506769618935e3
    ,2.5487111e0, -3.9302097727326e2
    ,4.9255514e0, -5.8849265665348e2
    ,1.3363463e0, -5.5076098609303e2
    ,1.6072053e0, -5.2237501616674e2
    ,1.3629480e0, -1.1790629318198e3
    ,5.5657014e0, -1.0977134971135e3
    ,5.0708205e0, -1.5774000881978e2
    ,3.9318944e0,  5.2963464780000e1
    ,4.8989497e0,  3.9809289073258e1
    ,1.3097446e0,  7.7540959633708e1
    ,3.5147141e0,  7.9618578146517e1
    ,3.5413158e0, -5.4868336758022e2)
    dcargs = matrix(dcargs,2,15)
    ccamps =
        c(-2.279594e-5,  1.407414e-5,  8.273188e-6,  1.340565e-5, -2.490817e-7
          ,-3.494537e-5,  2.860401e-7,  1.289448e-7,  1.627237e-5, -1.823138e-7
          , 6.593466e-7,  1.322572e-5,  9.258695e-6, -4.674248e-7, -3.646275e-7
          , 1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7
          , 9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7
          , 7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7
          ,-2.603449e-6,  7.359472e-6,  3.168357e-6,  1.119056e-6, -1.655307e-7
          ,-3.228859e-6,  1.308997e-7,  1.013137e-7,  2.403899e-6, -3.736225e-7
          , 3.442177e-7,  2.671323e-6,  1.832858e-6, -2.394688e-7, -3.478444e-7
          , 8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8
          ,-1.488378e-6, -1.251789e-5,  5.226868e-7, -2.049301e-7,  0.e0
          ,-8.043059e-6, -2.991300e-6,  1.473654e-7, -3.154542e-7,  0.e0
          , 3.699128e-6, -3.316126e-6,  2.901257e-7,  3.407826e-7,  0.e0
          , 2.550120e-6, -1.241123e-6,  9.901116e-8,  2.210482e-7,  0.e0
          ,-6.351059e-7,  2.341650e-6,  1.061492e-6,  2.878231e-7,  0.e0)
    ccamps = matrix(ccamps,5,15)
    ccsec3 = -7.757020e-8
    ccsec = c(1.289600e-6, 5.550147e-1, 2.076942e00
    ,3.102810e-5, 4.035027e00, 3.525565e-1
    ,9.124190e-6, 9.990265e-1, 2.622706e00
    ,9.793240e-7, 5.508259e00, 1.559103e01)
    ccsec = matrix(ccsec,3,4)
    dcsld = 1.990987e-7                   #sidereal rate in longitude
    ccsgd = 1.990969e-7                   #sidereal rate in mean anomaly
    cckm = 3.122140e-5
    ccmld = 2.661699e-6
    ccfdi = 2.399485e-7
    dcargm = c(5.1679830e0,  8.3286911095275e3
    ,5.4913150e0, -7.2140632838100e3
    ,5.9598530e0,  1.5542754389685e4)
    dcargm = matrix(dcargm,2,3)
    ccampm = c(1.097594e-1, 2.896773e-7, 5.450474e-2,  1.438491e-7
    ,-2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8
    , 1.148966e-2, 5.658888e-8, 8.249439e-3,  4.063015e-8)
    ccampm = matrix(ccampm,4,3)
    ccpamv = c(8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12)
    dc1mme = 0.99999696e0
    dt = (dje - dcto) / dcjul
    tvec = c(1e0, dt, dt*dt)
    temp = (tvec %*% dcfel) %% dc2pi
    dml = temp[1]
    forbel = temp[2:8]
    g = forbel[1]                         #old fortran equivalence
    deps = sum(tvec*dceps) %% dc2pi
    sorbel = (tvec %*% ccsel) %% dc2pi
    e = sorbel[1]                         #old fortran equivalence
    dummy=cos(2.0)
    sn = sin((tvec[1:2] %*% ccsec[2:3,]) %% cc2pi)
    pertl = sum(ccsec[1,] * sn) + dt*ccsec3*sn[2]
    pertld = 0.0
    pertr = 0.0
    pertrd = 0.0
    for(k in 1:15) {
        a = (dcargs[1,k]+dt*dcargs[2,k]) %% dc2pi
        cosa = cos(a)
        sina = sin(a)
        pertl = pertl + ccamps[1,k]*cosa + ccamps[2,k]*sina
        pertr = pertr + ccamps[3,k]*cosa + ccamps[4,k]*sina
        if(k<12 ){
            pertld = pertld + (ccamps[2,k]*cosa-ccamps[1,k]*sina)*ccamps[5,k]
            pertrd = pertrd + (ccamps[4,k]*cosa-ccamps[3,k]*sina)*ccamps[5,k]
        }
    }
    phi = (e*e/4e0)*(((8e0/e)-e)*sin(g) +5*sin(2*g) +(13/3e0)*e*sin(3*g))
    f = g + phi
    sinf = sin(f)
    cosf = cos(f)
    dpsi = (dc1 - e*e) / (dc1 + e*cosf)
    phid = 2*e*ccsgd*((1 + 1.5*e*e)*cosf + e*(1.25 - 0.5*sinf*sinf))
    psid = ccsgd*e*sinf / sqrt(dc1 - e*e)
    d1pdro = dc1+pertr
    drd = d1pdro * (psid + dpsi*pertrd)
    drld = d1pdro*dpsi * (dcsld+phid+pertld)
    dtl = (dml + phi + pertl) %% dc2pi
    dsinls = sin(dtl)
    dcosls = cos(dtl)
    dxhd = drd*dcosls - drld*dsinls
    dyhd = drd*dsinls + drld*dcosls
    pertl = 0.0
    pertld = 0.0
    pertp = 0.0
    pertpd = 0.0
    for(k  in 1:3) {
        a = (dcargm[1,k] + dt*dcargm[2,k]) %% dc2pi
        sina = sin(a)
        cosa = cos(a)
        pertl = pertl + ccampm[1,k]*sina
        pertld = pertld + ccampm[2,k]*cosa
        pertp = pertp + ccampm[3,k]*cosa
        pertpd = pertpd - ccampm[4,k]*sina
    }
    tl = forbel[2] + pertl
    sinlm = sin(tl)
    coslm = cos(tl)
    sigma = cckm / (1.0 + pertp)
    a = sigma*(ccmld + pertld)
    b = sigma*pertpd
    dxhd = dxhd + a*sinlm + b*coslm
    dyhd = dyhd - a*coslm + b*sinlm
    dzhd= -sigma*ccfdi*cos(forbel[3])
    dxbd = dxhd*dc1mme
    dybd = dyhd*dc1mme
    dzbd = dzhd*dc1mme
    for(k in 1:4) {
        plon = forbel[k+3]
        pomg = sorbel[k+1]
        pecc = sorbel[k+9]
        tl = (plon + 2.0*pecc*sin(plon-pomg)) %% cc2pi
        dxbd = dxbd + ccpamv[k]*(sin(tl) + pecc*sin(pomg))
        dybd = dybd - ccpamv[k]*(cos(tl) + pecc*cos(pomg))
        dzbd = dzbd - ccpamv[k]*sorbel[k+13]*cos(plon - sorbel[k+5])
    }
    dcosep = cos(deps)
    dsinep = sin(deps)
    dyahd = dcosep*dyhd - dsinep*dzhd
    dzahd = dsinep*dyhd + dcosep*dzhd
    dyabd = dcosep*dybd - dsinep*dzbd
    dzabd = dsinep*dybd + dcosep*dzbd
    if(deq==0 ){
        dvelh = au * (c(dxhd, dyahd, dzahd))
        dvelb = au * (c(dxbd, dyabd, dzabd))
        return(list(dvelh=dvelh, dvelb = dvelb))
    }
    deqdat = (dje-dcto-dcbes) / dctrop + dc1900
    prema = premat(deqdat,deq,fk4=T)
    dvelh = au * ( prema %*% c(dxhd, dyahd, dzahd) )
    dvelb = au * ( prema %*% c(dxbd, dyabd, dzabd) )
    return(list(dvelh=dvelh, dvelb = dvelb))
}
