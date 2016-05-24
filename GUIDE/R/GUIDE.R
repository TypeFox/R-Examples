GUIDE <-
function(){

main.menu <- function(panel) {
  
  
    
  #rp.messagebox(panel$menu, title = "Demo")
  if(panel$menu == "Stock"){
    forwardstock()
  }
  else if(panel$menu == "Currency"){
    forwardcurrency()
  }
  else if(panel$menu == "Commodity"){
    forwardcommodity()
  }
    else if(panel$menu == "Bond Forward Tree"){
    bondforwardtreegui()
  }
  else if(panel$menu == "Stocks"){
    futuresstock()
  }
  else if(panel$menu == "Bond Futures Tree"){
    bondfuturestreegui()
  }
  else if(panel$menu == "Currencies"){
    futurescurrency()
  }
  else if(panel$menu == "Commodities"){
    futurescommodity()
  }
  else if (panel$menu=="Eurodollar futures"){
    eurodollar()
  }
  else if(panel$menu == "Cash Price of T-Bond Future"){
    cashprice()
  }
  else if (panel$menu=="Payoff / P&L graphs"){
    basicpayoffs()
  }
  else if (panel$menu=="Premium 3D plots"){
    Premium3D()
  }
  else if (panel$menu=="Option Greeks"){
    calcgreeks()
  }
  else if (panel$menu=="Greeks 3D plots"){
    stockTimeGreeks()
  }
  else if (panel$menu=="Stock Option Tree"){
    stockoptiontreegui()
  }
  else if (panel$menu=="Cap Tree"){
    captreegui()
  }
  else if (panel$menu=="Floor Tree"){
    floortreegui()
  }
  else if (panel$menu=="Bond Option Tree"){
    bondoptiontreegui()
  }
  else if (panel$menu=="Black Scholes"){
    blackscholes()
  }
  else if (panel$menu=="Implied Volatility"){
    impvol()
  }
  else if (panel$menu=="Trading Strategies"){
    trading.menu()
  }
  else if (panel$menu=="Hedging with greeks"){
    greekneutrality()
  }
  else if (panel$menu=="Geometric Brownian Motion"){
    GBMPaths()
  }
  else if (panel$menu=="Arithmetic Brownian Motion"){
    ABMPaths()
  }
  else if (panel$menu=="Brownian Motion"){
    BrownianPaths()
  }
  else if (panel$menu=="Jump Diffusion"){
    JDPaths()
  }
  else if (panel$menu=="Present Value of an amount"){
    pv()
  }
  else if (panel$menu=="Future Value of an amount"){
    fv()
  }
  else if (panel$menu=="Present Value of Annuity"){
    pvann()
  }
  else if (panel$menu=="Future Value of Annuity"){
    fvann()
  }
    else if (panel$menu=="Rate conversion"){
    rate()
  }
  else if (panel$menu=="p-value calculator"){
    pval()
  }
  else if (panel$menu=="z-value calculator"){
    zval()
  }
  else if (panel$menu=="Forward Rate"){
    fra()
  } 
  else if (panel$menu=="Value of FRA"){
    fravalue()
  }
  else if (panel$menu=="Interest rate swap"){
    irswapvalue()
  }
  else if (panel$menu=="Currency swap (Fixed-Fixed)"){
    curswapvalue()
  }
  else if(panel$menu=="Credit Default swap"){
    cdswap()
  }
  else if(panel$menu=="Interest rate swap Tree"){
    swaptreegui()
  }
  else if(panel$menu=="Swaption Tree"){
    swaptiontreegui()
  }
  else if (panel$menu=="Single stock"){
    var1stock()
  }
  else if (panel$menu=="Two stocks"){
    var2stocks()
  }
  else if (panel$menu=="VaR Behavior graphs"){
    varbehavior()
  }
  else if (panel$menu=="Rate Tree"){
    ratetreegui()
  }
  else if (panel$menu=="Bond Tree"){
    bondtreegui()
  }
  else if (panel$menu=="Bond Price"){
    bondprice()
  }
  else if (panel$menu=="Price-Yield plot"){
    priceyield()
  }
  else if (panel$menu=="Price-Maturity plot"){
    pricematurity()
  }
  else if (panel$menu=="Duration"){
    bonddur()
  }
  else if (panel$menu=="Duration & Maturity"){
    durmaturity()
  }
  else if (panel$menu=="Duration & Yield"){
    duryield()
  }
  else if (panel$menu=="Duration & Coupon"){
    durcoupon()
  }
  else if (panel$menu=="Convexity"){
    bondconv()
  }
  else if (panel$menu=="DV01"){
    bondchange()
  }
  panel
#   else
#   panel
}

main.panel <- rp.control(title = "GUIDE",size=c(1366,768))
menu.list = list(
  #  list("Time Value", "Present Value of an amount", "Future Value of an amount", "Present Value of Annuity", "Future Value of Annuity"),
  #  list("Stocks", "Single stage growth","Two stage growth"),
  list("Forwards", "Commodity","Currency","Stock", "Bond Forward Tree", "Forward Rate","Value of FRA"),
  list("Futures", "Commodities","Currencies","Stocks", "Bond Futures Tree", "Eurodollar futures","Cash Price of T-Bond Future"),
  list("Options","Payoff / P&L graphs", "Premium 3D plots","Stock Option Tree","Bond Option Tree","Cap Tree", "Floor Tree", "Black Scholes", "Implied Volatility", "Option Greeks", "Greeks 3D plots","Hedging with greeks","Trading Strategies"),
  list("Swaps","Interest rate swap", "Currency swap (Fixed-Fixed)", "Credit Default swap","Interest rate swap Tree", "Swaption Tree"),
  list("Stochastic Processes", "Brownian Motion", "Arithmetic Brownian Motion","Geometric Brownian Motion","Jump Diffusion"),
  list("Value at Risk", "Single stock", "Two stocks","VaR Behavior graphs"),
  list("Bonds","Rate Tree", "Bond Tree", "Bond Price","Price-Yield plot","Price-Maturity plot", "Duration", "Duration & Maturity","Duration & Yield","Duration & Coupon", "Convexity","DV01" ),
  list("Utilities","Present Value of an amount", "Future Value of an amount", "Present Value of Annuity", "Future Value of Annuity","Rate conversion", "p-value calculator","z-value calculator"))

rp.menu(panel = main.panel, variable= menu, labels=menu.list, action = main.menu)
# for r studio- use line below
# image.file <- file.path(Sys.getenv("home"), "finderiv.gif")
# for r gui-use line below

image.file <- file.path(system.file(package = "GUIDE"), "extdata", "guide_img1.gif")
rp.image(main.panel, image.file)
invisible()
# end of main menu



}
