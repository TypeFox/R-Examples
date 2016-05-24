##
##  m a t l a b . R  Matlab Idioms
##


matlab <- function() {
    cat(paste("",
  "The following functions are emulations of corresponding Matlab functions",
  "and bear the same signature as their Matlab cousins as far as possible:",
  "
  accumarray, acot, acoth, acsc, acsch, and, angle, ans, arrayfun, asec, asech,
  beep, bernoulli, blank, blkdiag, bsxfun,
  cart2pol, cart2sph, cd, ceil, circshift, clear, compan, cond, conv,
      cot, coth, cross, csc, csch, cumtrapz,
  dblquad, deblank, deconv, deg2rad, detrend, deval, disp, dot,
  eig, eigint, ellipj, ellipke, eps, erf, erfc, erfcinv, erfcx, erfi, erfinv,
      errorbar, expint, expm, eye, ezcontour, ezmesh, ezplot, ezpolar,
  fact, fftshift, figure, findpeaks, findstr, flipdim, fliplr, flipud,
      fminbnd, fminsearch, fplot, fsolve, fzero,
  gammainc, gcd, geomean, gmres, gradient,
  hadamard, hankel, harmmean, hilb, histc, humps, hypot,
  idivide, ifft, ifftshift, inpolygon, integral, integral2, integral3,
      interp1, interp2, inv, isempty, isprime,
  kron,
  legendre, linprog, linspace, loglog, logm, logseq, logspace, lsqcurvefit,
      lsqlin, lsqnonlin, lsqnonneg, lu,
  magic, meshgrid, mkpp, mldivide, mod, mrdivide,
  nchoosek, ndims, nextpow2, nnz, normest, nthroot, null, num2str, numel,
  ode23, ode23s, ones, or, orth,
  pascal, pchip, pdist, pdist2, peaks, perms, piecewise, pinv, plotyy,
      pol2cart,polar, polyfit, polyint, polylog, polyval, pow2, ppval,
      primes, psi, pwd,
  quad, quad2d, quadgk, quadl, quadprog, quadv, quiver,
  rad2deg, randi, randn, randsample, rat, rats, regexp, regexpi, regexpreg,
      rem, repmat, roots, rosser, rot90, rref, runge,
  sec, sech, semilogx, semilogy, sinc, size, sortrows, sph2cart, sqrtm,
      squareform, std, str2num, strcat, strcmp, strcmpi, strfind, strfindi,
      strjust, subspace,
  tic, toc, trapz, tril, trimmean, triplequad, triu,
  vander, vectorfield, ver,
  what, who, whos, wilkinson,
  zeros, zeta.",
  "", 
  "The following Matlab function names have been capitalized in 'pracma' to",
  "avoid shadowing functions from R base or one of its recommended packages:",
  "
  Diag, factors, finds, Fix, Imag, Lcm, Mode, Norm, nullspace (null),
  Poly, Rank, Real, Reshape, strRep, strTrim, Toeplitz, Trace, uniq (unique).", 
  "",
  "To use 'ans' instead of 'ans()' (i.e., as is common practice in Matlab)",
  "type (and similar for other Matlab commands):",
  "
  makeActiveBinding('ans', function() .Last.value, .GlobalEnv)
  makeActiveBinding('who', who, .GlobalEnv)",
  "",
  "etc. after loading the 'pracma' package.",
  "\n",
    sep = "\n", collapse = ""))
  
    invisible(NULL)
}
