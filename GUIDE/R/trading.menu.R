if (getRversion() >= "2.15.1") utils::globalVariables(c('strategy'))

trading.menu <-
function(){
  trading.strat <- function(panel){
    if (panel$strategy == "Bull Spread"){
      bullspreadcalls()
    }
    else if (panel$strategy == "Bear Spread"){
      bearspreadputs()
    }
    else if (panel$strategy == "Butterfly"){
      butterfly()
    }
    else if (panel$strategy == "Reverse Butterfly"){
      reversebutterfly()
    }
    else if (panel$strategy == "Straddle"){
      straddle()
    }
    else if (panel$strategy == "Reverse Straddle"){
      reversestraddle()
    }
    else if (panel$strategy == "Strangle"){
      strangle()
    }
    else if (panel$strategy == "Reverse Strangle"){
      reversestrangle()
    }
    else if (panel$strategy == "Strip"){
      strip()
    }
    else if (panel$strategy == "Strap"){
      strap()
    }
      panel
  }
  
  
  my.panel <- rp.control(title="Option Trading Strategies")
  rp.radiogroup(my.panel, strategy,title= "Trading Strategy",
                c("None", "Bull Spread", "Bear Spread", "Butterfly","Reverse Butterfly", "Straddle", "Reverse Straddle", "Strangle", "Reverse Strangle",
                  "Strip","Strap"),
                action=trading.strat)
}
