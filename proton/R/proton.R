#' @title The Proton Game
#'
#' @description
#' The \code{proton} function is used for solving problems in the data-based game ,,The Proton Game''.
#' Solve four data-based puzzles in order to crack into Pietraszko's account!
#'
#' @param ... \code{proton} function is called by different arguments, which vary depending
#' on a problem that Bit is trying to solve. See \code{Details} in order to learn more about the list of possible arguments.
#'
#' @details Every time when some additional hints are needed one should add
#' \code{hint=TRUE} argument to the \code{proton} function.
#'
#' In order to get more information about a user on the Proton server
#' one should pass \code{action = "login"}, \code{login="XYZ"} arguments
#' to the \code{proton} function.
#'
#' In order to log into the Proton server one should pass \code{action = "login"},
#' \code{login="XYZ"}, \code{password="ABC"} arguments to the \code{proton} function.
#' If the password matches login, then one will receive a message about successful login.
#'
#' In order to log into a server different from Proton one should pass
#' \code{action = "server"}, \code{host="XYZ"} arguments to the \code{proton} function.
#'
#' ,,The Proton Game'' is a free of charge, educational project of the SmarterPoland.pl Foundation.
#' By transferring a donation on the foundation's account, which is
#' shown on site \url{http://smarterpoland.pl/index.php/fundacja/}, you will help
#' us to create other educational games. Even an amount of $10 will facilitate
#' creation and maintenance of subsequent educational projects!
#' Thank you!
#'
#' @author
#' Przemyslaw Biecek, \email{przemyslaw.biecek@@gmail.com}, SmarterPoland.pl Foundation.
#'
#' @examples
#' \dontrun{
#' proton()
#' proton(hint=TRUE)
#' }
#' @rdname proton
#' @importFrom digest digest
#' @export
proton <- function(...) {
 args <- list(...)

texts <- dcode(structure(c("kRVGIZHAPL FHVH Z KZHHDLIW DSRXS RH EVIB WRUURXFOG GL TFVHH.\nzG URIHG, GIB GL SZXP ZM ZXXLFMG LU Z KVIHLM DSRXS RH MLG ZH XZFGRLFH ZH kRVGIZHAPL.\n\nyFG DSL RH GSV DVZPVHG KLRMG? rMRGRZO RMEVHGRTZGRLM HFTTVHGH GSZG qLSM rMHVXFIV WLVHM'G XZIV ZYLFG HVXFIRGB ZMW SZH ZM ZXXLFMG LM GSV kILGLM HVIEVI. sV NZB FHV Z KZHHDLIW DSRXS RH VZHB GL XIZXP.\noVG'H ZGGZXP SRH ZXXLFMG URIHG!\n\nkILYOVN 1: uRMW GSV OLTRM LU qLSM rMHVXFIV.\n\nyRG SZH HXIZKKVW 'VNKOLBVVH' WZGZ (MZNVH ZMW OLTRMH) UILN GSV DDD DVY KZTV LU gVXSMRXZO fMREVIHRGB LU dZIHZD. gSV WZGZ RH RM GSV WZGZ.UIZNV `VNKOLBVVH`. \nmLD, BLFI GZHP RH GL URMW qLSM rMHVXFIV'H OLTRM.\ndSVM BLF URMZOOB URMW LFG DSZG qLSM'H OLTRM RH, FHV `KILGLM(ZXGRLM = \"OLTRM\", OLTRM=\"cba\")` XLNNZMW, DSVIV cba RH rMHVXFIV'H OLTRM.\n",
"rM `VNKOLBVVH` WZGZHVG GIB GL URMW Z ILD DSRXS SZH `rMHVXFIV` EZOFV RM GSV `HFIMZNV` XLOFNM.\nuFMXGRLMH ORPV `UROGVI` LI `ZIIZMTV` UILN GSV `WKOBI` KZXPZTV NZB YV EVIB FHVUFO.\n",
"QLSMRMH", "HOZK", "xLMTIZGFOZGRLMH! bLF SZEV ULFMW LFG DSZG qLSM rMHVXFIV'H OLTRM RH!\nrG RH SRTSOB ORPVOB GSZG SV FHVH HLNV GBKRXZO KZHHDLIW.\nyRG WLDMOLZWVW UILN GSV rMGVIMVG Z WZGZYZHV DRGS 1000 NLHG XLNNLMOB FHVW KZHHDLIWH.\nbLF XZM URMW GSRH WZGZYZHV RM GSV `GLK1000KZHHDLIWH` EVXGLI.\n\nkILYOVN 2: uRMW qLSM rMHVXFIV'H KZHHDLIW.\n\nfHV `KILGLM(ZXGRLM = \"OLTRM\", OLTRM=\"cba\", KZHHDLIW=\"zyx\")` XLNNZMW RM LIWVI GL OLT RMGL GSV kILGLM HVIEVI DRGS GSV TREVM XIVWVMGRZOH.\nrU GSV KZHHDLIW RH XLIIVXG, BLF DROO TVG GSV ULOOLDRMT NVHHZTV:\n`hFXXVHH! fHVI RH OLTTVW RM!`.\nlGSVIDRHV BLF DROO TVG:\n`kZHHDLIW LI OLTRM RH RMXLIIVXG!`.\n",
"fHV GSV YIFGV ULIXV NVGSLW.\nyB FHRMT Z OLLK, GIB GL OLT RM DRGS HFYHVJFVMG KZHHDLIWH UILN `GLK1000KZHHDLIWH` EVXGLI ZH OLMT ZH BLF IVXVREV:\n`hFXXVHH! fHVI RH OLTTVW RM!`.\n",
"kZHHDLIW LI OLTRM RH RMXLIIVXG", "hFXXVHH! fHVI RH OLTTVW RM!","dVOO WLMV! gSRH RH GSV IRTSG KZHHDLIW!\nyRG FHVW qLSM rMHVXFIV'H ZXXLFMG RM LIWVI GL OLT RMGL GSV kILGLM HVIEVI.\nrG GFIMH LFG GSZG qLSM SZH ZXXVHH GL HVIEVI OLTH.\nmLD, yRG DZMGH GL XSVXP UILN DSRXS DLIPHGZGRLM kRVGIZHAPL RH UIVJFVMGOB OLTTRMT RMGL GSV kILGLM HVIEVI. yRG SLKVH GSZG GSVIV DROO YV HLNV FHVUFO WZGZ.  \n\noLTH ZIV RM GSV `OLTH` WZGZHVG. \nxLMHVXFGREV XLOFNMH XLMGZRM RMULINZGRLM HFXS ZH: DSL, DSVM ZMW UILN DSRXS XLNKFGVI OLTTVW RMGL kILGLM.\n\nkILYOVN 3: xSVXP UILN DSRXS HVIEVI kRVGIZHAPL OLTH RMGL GSV kILGLM HVIEVI NLHG LUGVM.\n\nfHV `KILGLM(ZXGRLM = \"HVIEVI\", SLHG=\"cba\")` XLNNZMW RM LIWVI GL OVZIM NLIV  ZYLFG DSZG XZM YV ULFMW LM GSV cba HVIEVI.\ngSV YRTTVHG XSZMXV GL URMW HLNVGSRMT RMGVIVHGRMT RH GL URMW Z HVIEVI UILN DSRXS kRVGIZHAPL OLTH RM GSV NLHG LUGVM.\n\n",
"rM LIWVI GL TVG GL PMLD UILN DSRXS HVIEVI kRVGIZHAPL RH OLTTRMT GSV NLHG LUGVM LMV NZB:\n1. fHV `UROGVI` UFMXGRLM GL XSLLHV LMOB kRVGIZHAPL'H OLTH,\n2. fHV `TILFK_YB` ZMW `HFNNZIRHV` GL XLFMG GSV MFNYVI LU kRVGIZHAPL'H OLTH RMGL HVKZIZGV HVIEVIH,\n3. fHV `ZIIZMTV` UFMXGRLM GL HLIG HVIEVIH' ORHG YB GSV UIVJFVMXB LU OLTH.\n\nfHV `VNKOLBVVH` WZGZYZHV RM LIWVI GL XSVXP DSZG kRVGIZHAPL'H OLTRM RH.\n",
"mRXV GIB, YFG GSVIV RH MLGSRMT RMGVIVHGRMT ZYLFG GSRH OLTRM.\ngSV DVZPVHG ORMP LU kILGLM HVIEVI RH qLSM rMHVXFIV.\ngIB GL URMW SRH OLTRM.\n",
"xLMTIZGFOZGRLMH!\n\nbLF SZEV XIZXPVW kRVGIZHAPL'H KZHHDLIW!\nhVXIVG KOZMH LU SRH OZY ZIV MLD RM BLFI SZMWH.\ndSZG RH RM GSRH NBHGVIRLFH OZY?\nbLF NZB IVZW ZYLFG RG RM GSV `kRVGIZHAPL'H XZEV` HGLIB DSRXS RH ZEZROZYOV ZG SGGK://YRVXVP.KO/yVGZyRG/dZIHZD\n\nmVCG ZWEVMGFIV LU yVGZ ZMW yRG DROO YV ZEZROZYOV HLLM.\n\n",
"rG GFIMH LFG GSZG kRVGIZHAPL LUGVM FHVH GSV KFYORX DLIPHGZGRLM 194.29.178.16.\ndSZG Z XZIVOVHHMVHH.\n\nyRG RMUROGIZGVW GSRH DLIPHGZGRLM VZHROB. sV WLDMOLZWVW `YZHS_SRHGLIB` UROV DSRXS XLMGZRMH Z ORHG LU ZOO XLNNZMWH GSZG DVIV VMGVIVW RMGL GSV HVIEVI'H XLMHLOV.\ngSV XSZMXVH ZIV GSZG HLNV GRNV ZTL kRVGIZHAPL GBKVW Z KZHHDLIW RMGL GSV XLMHLOV YB NRHGZPV GSRMPRMT GSZG SV DZH OLTTRMT RMGL GSV kILGLM HVIEVI.\n\nkILYOVN 4: uRMW GSV kRVGIZHAPL'H KZHHDLIW.\n\nrM GSV `YZHS_SRHGLIB` WZGZHVG BLF DROO URMW ZOO XLNNZMWH ZMW KZIZNVGVIH DSRXS SZEV VEVI YVVM VMGVIVW.\ngIB GL VCGIZXG UILN GSRH WZGZHVG LMOB XLNNZMWH (LMOB HGIRMTH YVULIV HKZXV) ZMW XSVXP DSVGSVI LMV LU GSVN OLLPH ORPV Z KZHHDLIW.\n",
"xLNNZMWH ZMW KZIZNVGVIH ZIV HVKZIZGVW YB Z HKZXV. rM LIWVI GL VCGIZXG LMOB MZNVH LU XLNNZMWH UILN VZXS ORMV, BLF XZM FHV `THFY` LI `HGIHKORG` UFMXGRLM.\nzUGVI SZERMT ZOO XLNNZMWH VCGIZXGVW BLF HSLFOW XSVXP SLD LUGVM VZXS XLNNZMW RH FHVW.\nkVISZKH RG DROO GFIM LFG GSZG LMV LU GBKVW RM XLNNZMWH OLLP ORPV Z KZHHDLIW?\n\nrU BLF HVV HLNVGSRMT DSRXS OLLPH ORPV Z KZHHDLIW, BLF HSZOO FHV `KILGLM(ZXGRLM = \"OLTRM\", OLTRM=\"cba\", KZHHDLIW=\"zyx\")` XLNNZMW GL OLT RMGL GSV kILGLM HVIEVI DRGS kRVGIZHAPL XIVWVMGRZOH.\n",
"yRG HKVMG HLNV GRNV GL RMUROGIZGV GSRH DLIPHGZGRLM. \nyFG GSVIV RH MLGSRMT RMGVIVHGRMT SVIV.\nuRMW GSV DLIPHGZGRLM DSRXS kRVGIZHAPL RH FHRMT NLHG LUGVM GL OLT RMGL GSV kILGLM HVIEVI. "), .Names = c("proton.init", "proton.init.w", "log.1", "log.2","proton.login.init", "proton.login.init.w", "proton.login.fail","proton.login.pass", "proton.login.pass.instr", "proton.login.pass.instr.w","proton.login.weak", "proton.final", "proton.host.instr", "proton.host.instr.w","proton.host.instr.w2")))

 # plain start
 if (length(args) == 0) {
    cat(texts["proton.init"])
    return(invisible(NULL))
 }
 if (length(args) == 1 && !is.null(args$hint) && args$hint) {
   cat(texts["proton.init"], "\n\nHINT:\n",texts["proton.init.w"], sep = "")
   return(invisible(NULL))
 }

 # action = server
 if(length(args)>0 && !is.null(args$action) && args$action == "server") {
  if (!is.null(args$host) && digest(args$host) == "94265570be658d9fafa4861d7252afa9") {
    cat(texts["proton.host.instr"])
    if (!is.null(args$hint) && args$hint) {
      cat("\n\nHINT:\n",texts["proton.host.instr.w"], sep = "")
    }
    return(invisible(NULL))
  } else {
    cat("Bit has spent some time on infiltration of this workstation, but there is nothing interesting. \nFind the workstation that Pietraszko is using most often and try again.")
  }
 }

 # action = login
 if(length(args)>0 && !is.null(args$action) && args$action == "login") {
   # only user is set to johnins
   if (!is.null(args$login) && args$login == texts["log.1"] && is.null(args$password)) {
     cat(texts["proton.login.init"])
     if (!is.null(args$hint) && args$hint) {
       cat("\nHINT:\n",texts["proton.login.init.w"], sep = "")
     }
     return(invisible(NULL))
   }
   if(is.null(args$login)) {
     cat("\nIf action='login' argument is set then one should also set login=. argument \n")
     return(invisible(NULL))
   }

   # user is set to janie and password is provided
   if (!is.null(args$login) && args$login == texts["log.1"] && !is.null(args$password)) {
     if (digest(args$password) == "bbfb4a474b61b80225fd49d7c67e5a01") {
       cat(texts["proton.login.pass.instr"])
       if (!is.null(args$hint) && args$hint) {
         cat("\nHINT:\n",texts["proton.login.pass.instr.w"], sep = "")
       }
       return(texts["proton.login.pass"])
     } else {
       return(texts["proton.login.fail"])
     }
   }
   # user is set to sl and password is provided
   if (!is.null(args$login) && args$login == texts["log.2"] && !is.null(args$password)) {
     if (digest(args$password) == "ce3494fef4545c1b6160e5430d7efe66") {
       cat(texts["proton.final"])
       return(texts["proton.login.pass"])
     } else {
       return(texts["proton.login.fail"])
     }
   }

   # only user is set
   if (!is.null(args$login) && args$login != texts["log.1"]) {
     cat(texts["proton.login.weak"])
     return(invisible(NULL))
   }

 }
}
