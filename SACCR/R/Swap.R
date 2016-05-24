######################################################################
# Create the base Swap class
#
# This is used to represent the swap product and it will contained to all the swap-like classes

Swap = setRefClass("Swap",
                  fields = list(pay_leg_type = "character",
                                pay_leg_ref  = "character",
                                pay_leg_tenor= "character",
                                pay_leg_rate = "numeric",
                                rec_leg_type = "character",
                                rec_leg_tenor= "character",
                                rec_leg_ref  = "character",
                                rec_leg_rate = "numeric"
                                ),
                  methods = list(
                   isBasisSwap = function() {
                     if(length(pay_leg_type)!=0&&length(rec_leg_type)!=0)
                     {
                      if(pay_leg_type==rec_leg_type)
                      {
                        if(pay_leg_type=="Float"||pay_leg_type=="Commodity"||pay_leg_type=="Equity") return(TRUE)  
                      }
                     }
                     return(FALSE)
                    }
                  ))