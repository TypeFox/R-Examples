#' Determine Critical t-Statistic For Marginal Effects Plot
#'
#' This function determines the appropriate critical t-statistic
#' that limits the false null rejection rate of a marginal effects plot
#' to a specified value, err, using bootstrapped samples of calculated 
#' marginal effects values.
#' 
#' @param mat1 Matrix of bootstrapped samples of marginal effects.
#' @param mat2 Matrix of bootstrapped samples of marginal effects.
#' @param mat3 Matrix of bootstrapped samples of marginal effects.
#' @param mat4 Matrix of bootstrapped samples of marginal effects.
#' @param p1 The type of hypothesis test for the marginal effects in mat1. 1 = one-sided test for ME > 0, -1 = one-sided test for ME < 0, 0 = null of ME = 0 cannot be rejected, any other value = two-sided test for ME != 0
#' @param p2 The type of hypothesis test for the marginal effects in mat2. 1 = one-sided test for ME > 0, -1 = one-sided test for ME < 0, 0 = null of ME = 0 cannot be rejected, any other value = two-sided test for ME != 0
#' @param p3 The type of hypothesis test for the marginal effects in mat3. 1 = one-sided test for ME > 0, -1 = one-sided test for ME < 0, 0 = null of ME = 0 cannot be rejected, any other value = two-sided test for ME != 0
#' @param p4 The type of hypothesis test for the marginal effects in mat4. 1 = one-sided test for ME > 0, -1 = one-sided test for ME < 0, 0 = null of ME = 0 cannot be rejected, any other value = two-sided test for ME != 0
#' @param type The condition that corresponds to successful rejection of the null. "any" indicates any of the matK marginal effects are statistically significant in the direction given by pK, "all" = ALL of the matK marginal effects are statistically significant and in the direction given by pK
#' @param err Rejection rate.
#' @return The limits of the t-statistic for the given rejection rate
#' @examples
#' \dontrun{
#'   data(legfig)                # Clark and Golder 2006 replication data
#'   set.seed(1231124)
#'   
#'   # limit to established democracies from the 1990s
#'   dat<-subset(legfig, subset=(nineties==1 & old==1))
#'   
#'   # create bootstrap samples of marginal effects of eneg and logmag on enep1
#'   # uses the bootFun utility included in this package
#'   library(boot)
#'   boot.t.dist <- boot(data = dat, statistic = bootFun, R = 1000, 
#'             form=enep1 ~ eneg * logmag + uppertier_eneg + uppertier + proximity1 + 
#'             proximity1_enpres + enpres, fam="gaussian", x.name="eneg", 
#'             z.name="logmag")$t
#'   boot.t.x.dist<-boot.t.dist[,1:10]
#'   
#'   
#'   # calculate critical t-statistic that sets familywise error rate to 10%
#'   # for statistical significance of marginal effect of of eneg at any value of logmag
#'   findMultiLims(boot.t.x.dist, type="any", err=0.1)$minimum # answer: 2.593086
#'  
#'   # calculate critical t-statistic that sets FWER to 10% for ME of eneg = 0
#'   # when logmag is small and ME of eneg > 0 when logmag is large
#'   boot.t.x.dist.lo<-boot.t.dist[,1:5]
#'   boot.t.x.dist.hi<-boot.t.dist[,6:10]
#'   findMultiLims(boot.t.x.dist.lo, boot.t.x.dist.hi, type="all", p1=0, 
#'                   p2=1, err=0.1)$minimum     # answer: 1.008688
#' }
#' @author Justin Esarey and Jane Lawrence Sumner
#' @references Clark, William R., and Matt Golder. 2006. "Rehabilitating Duverger's Theory." \emph{Comparative Political Studies} 39(6): 679-708.
#' @references Esarey, Justin, and Jane Lawrence Sumner. 2015. "Marginal Effects in Interaction Models: Determining and Controlling the False Positive Rate." URL: http://jee3.web.rice.edu/interaction-overconfidence.pdf.
#' @importFrom stats optimize
#' @export

findMultiLims<-function(mat1, mat2=NULL, mat3=NULL, mat4=NULL, p1=99, p2=NULL, p3=NULL, p4=NULL, type="any", err=0.05){
  
  # a function to determine the number of one-tailed statistically significant
  # estimates in a vector v1 given upper limits v2 and lower limits v3
  sig<-function(v1, v2, v3, p){
    
    a<-max(v1>v2)
    b<-max(v1<v3)
    if(p==0){return(abs(max(a,b)-1))}
    if(p==1){return(a)}
    if(p==-1){return(b)}
    return(max(a,b))
    
  }
  
  # a function to return a squared error term corresponding to how far the false
  # positive rate produced by t quantile limits (K, 1-K) is from 5%
  lims<-function(K, sum.type=type){
    
    lim.cand<-c(-K, K)
    
    m2<-NA; m3<-NA; m4<-NA
    
    m1<-apply(X=mat1, FUN=sig, MARGIN=1, v2=lim.cand[2], v3=lim.cand[1], p1)
    n1<-nrow(mat1)
    if(is.null(mat2)==FALSE){m2<-apply(X=mat2, FUN=sig, MARGIN=1, v2=lim.cand[2], v3=lim.cand[1], p2)}
    if(is.null(mat3)==FALSE){m3<-apply(X=mat3, FUN=sig, MARGIN=1, v2=lim.cand[2], v3=lim.cand[1], p3)}
    if(is.null(mat4)==FALSE){m4<-apply(X=mat4, FUN=sig, MARGIN=1, v2=lim.cand[2], v3=lim.cand[1], p4)}
    
    if(sum.type=="any"){optimand<-sum(pmax(m1, m2, m3, m4, na.rm=T))/sum(n1)}
    if(sum.type=="all"){optimand<-sum(pmin(m1, m2, m3, m4, na.rm=T))/sum(n1)}
    
    return((optimand-err)^2)
  }
  
  # find the t quantile limits that get as close as possible to a 5% rejection rate
  # in the bootstrap sample at hand
  
  if(type=="any"){out<-optimize(f=lims, interval=c(0.001, 5), sum.type=type)}
  if(type=="all"){out<-optimize(f=lims, interval=c(0.001, 2), sum.type=type)}
  
  # if there were no hits at all, just return t=0 as the optimum
  if(err^2==out$objective){out<-list(minimum=0, objective=NA); return(out)}else{return(out)}
  
}
