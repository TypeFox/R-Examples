##' This function perform some check on the data
##' 
##' The function only works for FAOST_CODE. If the country coding
##' system is not in FAOST_CODE then use the translateCountryCode
##' function to translate it.
##'
##' @param var The variable to be checked.
##' @param year The column which index the time.
##' @param data The data frame.
##' @param type The type of check.
##' @param take The type of check/replacement to be done in case of
##' type equals to overlap.
##'
##' @export
##' @examples
##' ## test.df = 
##' ##    data.frame(FAOST_CODE = rep(c(51,167,199), each = 3), 
##' ##      Year = rep(c(1990:1992), 3),
##' ##      Value = c(c(3,4,4), c(2,2,2), c(1,2,NA)))
##' ## FAOcheck(var = "Value", data = test.df, type = "overlap", take = "simpleCheck")
##' ## FAOcheck(var = "Value", data = test.df, type = "overlap", take = "takeNew")
##' ## FAOcheck(var = "Value", data = test.df, type = "overlap", take = "takeOld")
##' ## FAOcheck(var = "Value", data = test.df, type = "overlap", take = "complete")

FAOcheck = function(var, year = "Year", data, 
                    type = c("overlap", "multiChina"),
                    take = c("simpleCheck", "takeNew", "takeOld", "complete")){
    type = match.arg(type)
    take = match.arg(take)
    if(type == "overlap"){
        printLab("Check for overlap between transitional country")
        ## Belgium-Luxembourg (old entity, up to 1999 included) vs 
        ## Belgium, LuxemBourg (new entities, from 2000).
        data = overlap(old = 15, new = c(255, 256), 
                       var = var, data = data, take = take)
        ## Czechoslovakia (old entity, up to 1992 included) vs
        ## Czech Republic, Slovakia (new entities, from 1993).
        data = overlap(old = 51, new = c(167, 199), 
                       var = var, data = data, take = take)
        ## Ethiopia PDR (old entity, up to 1992 included) vs 
        ## Eritrea, Ethiopia (new entities, from 1993).
        data = overlap(old = 62, new = c(178, 238), 
                       var = var, data = data, take = take)
        ## Serbia and Montenegro (old entity, up to 2005 included) vs 
        ## Serbia, Montenegro (new entities, from 2006).
        data = overlap(old = 186, new = c(272, 273), 
                       var = var, data = data, take = take)
        ## USSR (old entity, up to 1991 included) vs 
        ## Armenia, Azerbaijan, Belarus, Estonia, Georgia, Kazakistan
        ## Kyrghisistan, Latvia, Lithuania,  Moldova, Russian Federation, 
        ## Tajikistan, Turkmenistan, Ukraine, Uzbekistan (new entities, from 1992).
        data = overlap(old = 228, new = c(1,52,57,63,73,108,113,119,126,146,185,208,213,230,235), 
                       var = var, data = data, take = take)
        ## Yemen (old entities, up to ...) vs 
        ## Yemen (new entities, from ...)
        data = overlap(old = c(246, 247), new = 249, 
                       var = var, data = data, take = take)
        ## Yugoslav SFR (old entity, up to 1991 included) vs 
        ## Serbia e Montenegro, Bosnia and Herzegovina, Croatia,
        ## Slovenia, Macedonia (new entities, from 1992).
        data = overlap(old = 248, new = c(80, 98, 154, 186, 272, 273, 198),
                var = var, data = data, take = take)
        ## Sudan (former) (old entity, up to 2011 included) vs 
        ## Sudan, South Sudan (new entities, from 2012).
        data = overlap(old = 206, new = c(276, 277), 
                       var = var, data = data, take = take)
        ## Netherlands Antilles (old entity, up to 2010 included) vs 
        ## Aruba, Bonaire, Sint Eustatius and Saba, Curacao, 
        ## Sint Maarten (Dutch Part) (new entities, from 2011).
        data = overlap(old = 151, new = c(22,278,279,280), 
                       var = var, data = data, take = take)
        ## Pacific Islands, Trust Territory of (old entity, up to 1990 included) vs 
        ## Micronesia (Federated States of), the Marshall Islands,
        ## Northern Mariana Islands, Palau (new entities, from 1991).
        data = overlap(old = 164, new = c(145,127,163,180), 
                       var = var, data = data, take = take)
        cat("\nNOTE: It is common for data reported by the predecessor\n or the new transitional country to include the new country\n")
    }
    if(type == "multiChina"){
        printLab("Check for existence of multiple China")
        data = CHMT(var = var, data = data, year = year)
    }
    data
}


