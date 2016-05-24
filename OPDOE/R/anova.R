size.anova <- function(model,hypothesis="", assumption="",a=NULL,b=NULL,c=NULL,n=NULL,alpha,beta,delta,cases){
                                        # remove spaces
  model <- gsub(" ", "", model)
  hypothesis <- gsub(" ", "", hypothesis)
  assumption <- gsub(" ", "", assumption)
  
                                        # oneway
  if(model=="a"){
    if(is.null(a))
      stop("need a, number of levels of factor A")
    ret <- size_n.one_way.model_1(a=a,
                                  alpha=alpha,beta=beta,delta=delta,cases=cases)
    names(ret) <- "n"
  } else
                                        # two way
  if(regexpr("^a[x>]b$", model, ignore.case=TRUE)>0){
    if(is.null(a))
      stop("need a, number of levels of factor A")
                                        # twoway cross    
    if(model=="axb"){
      if(is.null(n)){
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_n.two_way_cross.model_1_a(a=a,b=b,
                                                alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="axb"){
          ret <- size_n.two_way_cross.model_1_axb(a=a,b=b,
                                                  alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        }
        
        else{
          stop("omit n")
        }
      }
    } else if(model=="axB"){
      if(is.null(n))
        stop("need n, n=1 for no interaction, n=2 for interaction")
      if(is.null(b)){
        ret <- size_b.two_way_cross.mixed_model_a_fixed_a(a=a,n=n,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
        names(ret) <- "b"
      } else {
        stop("omit b")
      }
    }
                                        # twoway nested
    else if(model=="a>b"){
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n)){
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_n.two_way_nested.model_1_test_factor_a(a=a,b=b,
                                                             alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="b"){
          ret <- size_n.two_way_nested.model_1_test_factor_b(a=a,b=b,
                                                             alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        }
      } else {
        stop("omit n")
      }
    } else if(model=="A>b"){
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n)){
        ret <- size_n.two_way_nested.a_random_b_fixed_b(a=a,b=b,
                                                        alpha=alpha,beta=beta,delta=delta,cases=cases)
        names(ret) <- "n"
      } else {
        stop("omit n")
      }
    } else if(model=="a>B"){
      if(is.null(b)){
        ret <- size_b.two_way_nested.b_random_a_fixed_a(a=a,
                                                        alpha=alpha,beta=beta,delta=delta,cases=cases)
        names(ret) <- "b"
      } else {
        stop("omit b")
      }
    } else {
      stop("Two way ANOVA model not specified correctly")
    }
  } else
                                        # threeway
  if(regexpr("^[(]{0,1}a[x>]b[)]{0,1}[x>]c$", model, ignore.case=TRUE)>0){
                                        # threeway cross
    if(model=="axbxc"){
      if(hypothesis=="" | hypothesis=="a"){
        if(is.null(n)){
          ret <- size_n.three_way_cross.model_1_a(a=a,b=b,c=c,
                                                  alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        }else{
          stop("omit n")
        }
      } else if(hypothesis=="axb"){
        if(is.null(n)){
          ret <- size_n.three_way_cross.model_1_axb(a=a,b=b,c=c,
                                                    alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        }else{
          stop("omit n")
        }
      } else if(hypothesis=="axbxc"){
        if(is.null(n)){
          ret <- size_n.three_way_cross.model_1_axbxc(a=a,b=b,c=c,
                                                      alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        }else{
          stop("omit n")
        }
        stop("unknown hypothesis")
      }
    } else if(model=="axbxC"){
      if(is.null(n))
        stop("need n")
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c)){
        if(hypothesis=="" | hypothesis=="a"){
          ret <- size_c.three_way_cross.model_3_a(a=a,b=b,n=n,
                                                  alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
          
        } else if (hypothesis=="axb"){
          ret <- size_c.three_way_cross.model_3_axb(a=a,b=b,n=n,
                                                    alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit c")
      }
    } else if(model=="axBxC"){
      if(is.null(n))
        stop("need n")
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b) & is.null(c)){
        if(hypothesis=="" | hypothesis=="a"){
          if(assumption=="sigma_AB=0,b=c"){
            ret <- size_bc.three_way_cross.model_4_a_case2(a=a,n=n,
                                                           alpha=alpha,beta=beta,delta=delta,cases=cases)
            ret <- c(ret,ret)
            names(ret) <- c("b","c")
          } else if (assumption=="sigma_AC=0,b=c"){
            ret <- size_bc.three_way_cross.model_4_a_case1(a=a,n=n,
                                                           alpha=alpha,beta=beta,delta=delta,cases=cases)
            ret <- c(ret,ret)
            names(ret) <- c("b","c")
          } else {
            stop("need an assumption on AB or AC")
          }
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit b,c")
      }
    }
    
                                        # threeway nested    
    else if(model=="a>b>c"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(hypothesis=="" | hypothesis=="a"){
          ret <- size_n.three_way_nested.model_1_a(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if (hypothesis=="b"){
          ret <- size_n.three_way_nested.model_1_b(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if (hypothesis=="c"){
          ret <- size_n.three_way_nested.model_1_c(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }
      
    } else if(model=="a>B>C"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n))
        stop("need n")
      if(is.null(b)){
        if(hypothesis=="a" | hypothesis==""){
          ret <- size_b.three_way_nested.model_6_a(a=a,c=c,n=n,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "b"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit b")
      }      
      
    } else if(model=="A>b>C"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n))
        stop("need n")
      if(is.null(c)){
        if(hypothesis=="b" | hypothesis==""){
          ret <- size_c.three_way_nested.model_7_b(a=a,b=b,n=n,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }            
      
    } else if(model=="A>B>c"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(hypothesis=="A>B>c" | hypothesis==""){
          ret <- size_n.three_way_nested.model_8_c(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }      
      
    } else if(model=="a>b>C"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n))
        stop("need n")
      if(is.null(c)){
        if(hypothesis=="a" | hypothesis==""){
          ret <- size_c.three_way_nested.model_5_a(a=a,b=b,n=n,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else if(hypothesis=="a>b"){
          ret <- size_c.three_way_nested.model_5_b(a=a,b=b,n=n,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit c")
      }      
      

      
    } else if(model=="a>B>c"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(hypothesis=="a"){
          ret <- size_n.three_way_nested.model_4_a(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          # handle special case returning n und b:
          if(is.list(ret)){
            ret0 <- c(ret$n,ret$b)
            names(ret0) <- c("n","b")
            ret <- ret0
          } else {
            names(ret) <- "n"
          }

          names(ret) <- "n"
        } else if(hypothesis=="c"){
          ret <- size_n.three_way_nested.model_4_c(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          # handle special case returning n und b:
          if(is.list(ret)){
            ret0 <- c(ret$n,ret$b)
            names(ret0) <- c("n","b")
            ret <- ret0
          } else {
            names(ret) <- "n"
          }
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }      
      
    } else if(model=="A>b>c"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(hypothesis=="A>b"){
          ret <- size_n.three_way_nested.model_3_b(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="A>b>c"){
          ret <- size_n.three_way_nested.model_3_c(a=a,b=b,c=c,
                                                   alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }      
    }
                                        # threeway mixed
    else if(model=="(axb)>c"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_n.three_way_mixed_ab_in_c.model_1_a(a=a,b=b,c=c,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="axb"){
          ret <- size_n.three_way_mixed_ab_in_c.model_1_axb(a=a,b=b,c=c,
                                                            alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="b"){
          ret <- size_n.three_way_mixed_ab_in_c.model_1_b(a=a,b=b,c=c,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="(axb)>c"){
          ret <- size_n.three_way_mixed_ab_in_c.model_1_c(a=a,b=b,c=c,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }

    } else if(model=="(axB)>c"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(is.null(b))
          stop("need b, number of levels of factor B")

        if(hypothesis=="(axB)>c"|hypothesis==""){
          ret <- size_n.three_way_mixed_ab_in_c.model_3_c(a=a,b=b,c=c,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else {
          stop("unknown hypothesis")
        }
      } else if(is.null(b)){
        if(is.null(n))
          stop("need n")
        
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_b.three_way_mixed_ab_in_c.model_3_a(a=a,c=c,n=n,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "b"
        } else {
          stop("unknown hypothesis")
        }
        
      } else {
        stop("omit n or b")
      }
      
    } else if(model=="(axb)>C"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n))
        stop("need n")
      if(is.null(c)){
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_c.three_way_mixed_ab_in_c.model_5_a(a=a,b=b,n=n,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else if(hypothesis=="b"){
          ret <- size_c.three_way_mixed_ab_in_c.model_5_b(a=a,b=b,n=n,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else if(hypothesis=="axb"){
          ret <- size_c.three_way_mixed_ab_in_c.model_5_axb(a=a,b=b,n=n,
                                                            alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit on c")
      }

    } else if(model=="(axB)>C"){
      if(is.null(n))
        stop("need n")
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b) & is.null(c)){
        if(hypothesis=="" | hypothesis=="a"){
          if(assumption=="sigma_AC=0,b=c"){
            ret <- size_bc.three_way_mixed_cxbina.model_6_a_case1(a=a,n=n,
                                                                  alpha=alpha,beta=beta,delta=delta,cases=cases)
            ret <- c(ret,ret)
            names(ret) <- c("b","c")
          } else if (assumption=="sigma_AB=0,b=c"){
            ret <- size_bc.three_way_mixed_cxbina.model_6_a_case2(a=a,n=n,
                                                                  alpha=alpha,beta=beta,delta=delta,cases=cases)
            ret <- c(ret,ret)
            names(ret) <- c("b","c")
          } else {
            stop("need an assumption on sigma_AB or sigma_AC")
          }
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit b,c")
      }
      
      
    } else if(model=="(AxB)>c"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(hypothesis=="(AxB)>c"|hypothesis==""){
          ret <- size_n.three_way_mixed_ab_in_c.model_4_c(a=a,b=b,c=c,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }
    } else if(model=="(Axb)>C"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n))
        stop("need n")
      if(is.null(c)){
        if(hypothesis=="b"|hypothesis==""){
          ret <- size_c.three_way_mixed_ab_in_c.model_6_b(a=a,b=b,n=n,
                                                          alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit c")
      }
      
    } else if(model=="(a>b)xc"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_n.three_way_mixed_cxbina.model_1_a(a=a,b=b,c=c,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="axc"){
          ret <- size_n.three_way_mixed_cxbina.model_1_axc(a=a,b=b,c=c,
                                                           alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="a>b"){
          ret <- size_n.three_way_mixed_cxbina.model_1_b(a=a,b=b,c=c,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="bxc"){
          ret <- size_n.three_way_mixed_cxbina.model_1_bxc(a=a,b=b,c=c,
                                                           alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else if(hypothesis=="c"){
          ret <- size_n.three_way_mixed_cxbina.model_1_c(a=a,b=b,c=c,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else { 
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n")
      }

      
    } else if(model=="(A>b)xc"){
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(is.null(a))
          stop("need a, number of levels of factor A")

        if(hypothesis=="A>b"|hypothesis==""){
          ret <- size_n.three_way_mixed_cxbina.model_3_b(a=a,b=b,c=c,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else { 
          stop("unknown hypothesis")
        }
      } else if(is.null(a)){
        if(is.null(n))
          stop("need n")
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_a.three_way_mixed_cxbina.model_3_c(b=b,c=c,n=n,
                        alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "a"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n or a")
      }
      

    } else if(model=="(a>B)xc"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n)){
        if(is.null(b))
          stop("need b, number of levels of factor B")
        
        if(hypothesis=="Bxc"|hypothesis==""){
          ret <- size_n.three_way_mixed_cxbina.model_3_bxc(a=a,b=b,c=c,
                      alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "n"
        } else { 
          stop("unknown hypothesis")
        }
      } else if(is.null(b)){
        if(is.null(n))
          stop("need n")
        
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_b.three_way_mixed_cxbina.model_4_a(a=a,c=c,n=n,
                        alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "b"
        } else if (hypothesis=="axc"){
          ret <- size_b.three_way_mixed_cxbina.model_4_axc(a=a,c=c,n=n,
                                                           alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "b"
        } else if (hypothesis=="c"){
          ret <- size_b.three_way_mixed_cxbina.model_4_c(a=a,c=c,n=n,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "b"
        } else { 
          stop("unknown hypothesis")
        }
      } else {
        stop("omit n or b")
      }
      
    } else if(model=="(a>b)xC"){
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n))
        stop("need n")
      if(is.null(c)){
        if(hypothesis=="a"|hypothesis==""){
          ret <- size_c.three_way_mixed_cxbina.model_5_a(a=a,b=b,n=n,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else if(hypothesis=="b"){
          ret <- size_c.three_way_mixed_cxbina.model_5_b(a=a,b=b,n=n,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit c")
      }
      
    } else if(model=="(a>B)xC"){
      warning("not implemented")
    } else if(model=="(A>b)xC"){
                                        #
                                        #
      if(is.null(a))
        stop("need a, number of levels of factor A")
      if(is.null(b))
        stop("need b, number of levels of factor B")
      if(is.null(n))
        stop("need n")
      if(is.null(c)){
        if(hypothesis=="A>b"|hypothesis==""){
          ret <- size_c.three_way_mixed_cxbina.model_7_b(a=a,b=b,n=n,
                 alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "c"
        } else {
          stop("unknown hypothesis")
        }
      } else {
        stop("omit c")
      }
      
    } else if(model=="(A>B)xc"){
      if(is.null(c))
        stop("need c, number of levels of factor C")
      if(is.null(n))
        stop("need n")
      if(hypothesis=="c"|hypothesis==""){

        if(is.null(b)){
          if(is.null(a)){
            ret <- size_ab.three_way_mixed_cxbina.model_7_c(c=c,n=n,
                                                            alpha=alpha,beta=beta,delta=delta,cases=cases)
            ret <- c(ret,ret)
            names(ret) <- c("a","b")
          } else {
            ret <- size_b.three_way_mixed_cxbina.model_7_c(a=a,c=c,n=n,
                                                           alpha=alpha,beta=beta,delta=delta,cases=cases)
            names(ret) <- c("b")
          }
        } else if(is.null(a)){
          ret <- size_a.three_way_mixed_cxbina.model_7_c(b=b,c=c,n=n,
                                                         alpha=alpha,beta=beta,delta=delta,cases=cases)
          names(ret) <- "a"
        } else {
          stop("omit a and/or b")
        }
        
      } else {
        stop("unknown hypothesis")
      }
    } else {
      stop("Three way ANOVA model not specified correctly")
    }
  } else {
    stop("model not specified correctly")
  }
  ret
}


