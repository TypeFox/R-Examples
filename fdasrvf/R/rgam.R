#' Random Warping
#'
#' Generates random warping functions
#'
#' @param N length of warping function
#' @param sigma variance of warping functions
#' @param num number of warping functions
#' @return gam warping functions
#' @keywords diffeomorphism, warping function
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' gam = rgam(N=101, sigma=.01, num=35)
rgam <- function(N, sigma, num){

    gam = matrix(0,num,N)

    TT = N - 1
    time = seq(0,1,length.out=TT)
    mu = sqrt(rep(1,N-1)*TT/(N-1))
    omega = (2*pi)/TT
    for (k in 1:num){
        alpha_i = rnorm(1,sd=sigma)
        v = alpha_i * rep(1,TT)
        cnt = 1
        for (l in 2:10){
            alpha_i = rnorm(1,sd=sigma)
            if (l %% 2 !=0){#odd
                v = v + alpha_i*sqrt(2)*cos(cnt*omega*time)
                cnt = cnt + 1
            }
            if (l %% 2 ==0) #even
                v = v + alpha_i*sqrt(2)*sin(cnt*omega*time)
        }

        v = v - (mu%*%t(v))%*%mu/TT;
        vn = pvecnorm(v,2)/sqrt(TT);
        psi = cos(vn)*mu + sin(vn)*v/vn;
        gam[k,] = c(0,cumsum(psi*psi))/N

    }

    return(gam)
}
