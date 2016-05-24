#' Creates a collateral agreement Object containing all the relevant data and methods regarding the maturity factor
#'  and the calculation of the exposures after applying the relevant threshold
#' @title  CSAb Class
#' @param thres_cpty The maximum exposure that the counterparty can generate before collateral will need to be posted
#' @param thres_PO   The maximum exposure that the processing organization can generate before collateral will need to be posted
#' @param MTA_cpty   The minimum transfer amount for the counterparty
#' @param MTA_PO     The minimum transfer amount for the processing organization
#' @param IM_cpty    The initial margin that is posted by the counterparty
#' @param IM_PO      The initial margin that is posted by the processing organization
#' @param mpor_days  The margin period of risk in days
#' @param remargin_freq The frequency of re-margining the exposure in days
#' @param rounding   The rounding amount of the transfers
#' @return An object of type CSAb
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
#' @examples
#'
#' ## the margin agreement given in the Basel regulation example
#' coll = CSAb(thres_cpty = 0, MTA_cpty = 5, IM_cpty = 150, remargin_freq = 5)

CSAb = setRefClass("CSAb",

                  fields = list(thres_cpty = "numeric",
                                thres_PO   = "numeric",
                                MTA_cpty   = "numeric",
                                MTA_PO     = "numeric",
                                IM_cpty    = "numeric",
                                IM_PO      = "numeric",
                                mpor_days  = "numeric",
                                remargin_freq   = "numeric",
                                rounding   = "numeric"
                  ),

                  methods = list(
                    ApplyThres = function(MtM_vector)
                    {

                      MtM_len = length(MtM_vector)
                      coll_MtM    = rep(0,MtM_len)
                      collateral = 0
                      coll_MtM[1] = MtM_vector[1]

                      for (i in 1:(MtM_len/2))
                      {

                        if(coll_MtM[2*(i-1)+1]>thres_cpty+MTA_cpty)
                        { collateral = collateral + coll_MtM[2*(i-1)+1]-thres_cpty-MTA_cpty
                        } else if(coll_MtM[2*(i-1)+1]< (-thres_PO-MTA_PO))
                        { collateral = collateral + coll_MtM[2*(i-1)+1]-(-thres_PO-MTA_PO)
                        }   else if(i>2&&collateral>0&&MtM_diff<0)
                        { collateral = max(collateral+MtM_diff,0)
                        } else if(i>2&&collateral<0&&MtM_diff>0)
                        { collateral = min(collateral+MtM_diff,0)
                        }

                        coll_MtM[2*i]   = MtM_vector[2*i]   - collateral
                        coll_MtM[2*i+1] = MtM_vector[2*i+1] - collateral
                        MtM_diff = MtM_vector[2*i+1] - MtM_vector[2*i-1]
                      }
                      return(coll_MtM)
                    },
                    CalcMF = function()
                    {
                      MPOR = 10 + remargin_freq -1
                      return(1.5*sqrt(MPOR/250))
                    }
                  )
)
