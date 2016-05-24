error = function (forecast, forecastbench, true, insampletrue, method = c("me", "mpe", 
    "mae", "mse", "sse", "rmse", "mdae", "mdse", "mape", "mdape", 
    "smape", "smdape", "rmspe", "rmdspe", "mrae", "mdrae", "gmrae", 
    "relmae", "relmse", "mase", "mdase", "rmsse"), giveall = FALSE) 
{
    method = match.arg(method)
    if (giveall == FALSE) {
        if (method == "me") {
            val = me(forecast, true)
        }
        if (method == "mpe") {
            val = mpe(forecast, true)
        }
        if (method == "mae") {
            val = mae(forecast, true)
        }
        if (method == "mse") {
            val = mse(forecast, true)
        }
        if (method == "sse") {
            val = sse(forecast, true)
        }
        if (method == "rmse") {
            val = rmse(forecast, true)
        }
        if (method == "mdae") {
            val = mdae(forecast, true)
        }
        if (method == "mdse") {
            val = mdse(forecast, true)
        }
        if (method == "mape") {
            val = mape(forecast, true)
        }
        if (method == "mdape") {
            val = mdape(forecast, true)
        }
        if (method == "smape") {
            val = smape(forecast, true)
        }
        if (method == "smdape") {
            val = smdape(forecast, true)
        }
        if (method == "rmspe") {
            val = rmspe(forecast, true)
        }
        if (method == "rmdspe") {
            val = rmdspe(forecast, true)
        }
        if (method == "mrae") {
            val = mrae(forecast, forecastbench, true)
        }
        if (method == "mdrae") {
            val = mdrae(forecast, forecastbench, true)
        }
        if (method == "gmrae") {
            val = gmrae(forecast, forecastbench, true)
        }
        if (method == "relmae") {
            val = relmae(forecast, forecastbench, true)
        }
        if (method == "relmse") {
            val = relmse(forecast, forecastbench, true)
        }
        if (method == "mase") {
            val = mase(forecast, true, insampletrue)
        }
        if (method == "mdase") {
            val = mdase(forecast, true, insampletrue)
        }
        if (method == "rmsse") {
            val = rmsse(forecast, true, insampletrue)
        }
        return(val)
    }
    else {
        out = c(me(forecast, true), mpe(forecast, true), mae(forecast, 
            true), mse(forecast, true), sse(forecast, true), 
            rmse(forecast, true), mdae(forecast, true), mdse(forecast, 
                true), mape(forecast, true), mdape(forecast, 
                true), smape(forecast, true), smdape(forecast, 
                true), rmspe(forecast, true), rmdspe(forecast, 
                true), mrae(forecast, forecastbench, true), mdrae(forecast, 
                forecastbench, true), gmrae(forecast, forecastbench, 
                true), relmae(forecast, forecastbench, true), 
            relmse(forecast, forecastbench, true), mase(forecast, 
                true, insampletrue), mdase(forecast, true, insampletrue), 
		rmsse(forecast, true, insampletrue))
        names(out) = c("ME", "MPE", "MAE", "MSE", "SSE", "RMSE", 
            "MDAE", "MDSE", "MAPE", "MDAPE", "SMAPE", "SMDAPE", 
            "RMSPE", "RMDSPE", "MRAE", "MDRAE", "GMRAE", "RELMAE", 
            "RELMSE", "MASE", "MDASE", "RMSSE")
        return(out)
    }
}
