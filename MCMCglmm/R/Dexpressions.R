"Dexpressions"<-list(
"PW.d0W"=expression((exp(f)/(t_b-t_a))*((gamma(1/g)*(pgamma((exp(s)*t_a)^g, 1/g, lower.tail = FALSE)-pgamma((exp(s)*t_b)^g, 1/g,lower.tail = FALSE)))/(exp(-(exp(s)*t_0)^g)*exp(s)*g))+(exp(-(exp(s)*t_1)^g)/exp(-(exp(s)*t_0)^g)))
,
"PW.d1Wds"=expression(exp((exp(s)*t_0)^g)*(exp(-(exp(s)*t_1)^g)*g*(((exp(s)*t_0)^g)-((exp(s)*t_1)^g))-exp(f-log(t_b-t_a)-s-(exp(s)*t_a)^g)*exp(s)*t_a + exp(f-log(t_b-t_a)-s-(exp(s)*t_b)^g)*exp(s)*t_b + ((1/g)*exp(f-log(t_b-t_a)-s)*(g*(exp(s)*t_0)^g - 1)*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g)))))
,
"PW.d1Wdf"=expression((exp(f-log(t_b-t_a)-s+(exp(s)*t_0)^g)*(1/g)*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g))))
,
"PW.d3Wd2sd1f"=expression(-(1/(g*(t_a-t_b)))*exp(f-s+(exp(s)*t_0)^g-(exp(s)*t_a)^g-(exp(s)*t_b)^g)*(g*(exp(exp(s)*t_b)^g*(exp(s)*t_a)*(1+g*(-2*(exp(s)*t_0)^g + (exp(s)*t_a)^g))-exp(exp(s)*t_a)^g*(exp(s)*t_b)*(1+g*(-2*(exp(s)*t_0)^g + (exp(s)*t_b)^g))) + exp((exp(s)*t_a)^g + (exp(s)*t_b)^g)*(1 + g*(exp(s)*t_0)^g*(-2 + g + g*(exp(s)*t_0)^g))*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g))))
,
"PW.d3Wd1sd2f"=expression((1/(g*(t_a-t_b)))*exp(f-s+(exp(s)*t_0)^g)*(exp(-(exp(s)*t_a)^g)*g*(exp(s)*t_a)-exp(-(exp(s)*t_b)^g)*g*(exp(s)*t_b)+(1-g*(exp(s)*t_0)^g)*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g))))
,
"PW.d2Wd2f"=expression((exp(f-log(t_b-t_a)-s+(exp(s)*t_0)^g)*(1/g)*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g))))
,
"PW.d2Wd2s"=expression((1/(g*(t_a-t_b)))*exp(-s+(exp(s)*t_0)^g-(exp(s)*t_1)^g-(exp(s)* t_a)^g-(exp(s)*t_b)^g)*(g*(-exp(f+(exp(s)*t_1)^g+(exp(s)*t_b)^g)*exp(s)*t_a*(1+g*(-2*(exp(s)*t_0)^g+(exp(s)* t_a)^g))+exp(s+(exp(s)*t_a)^g+(exp(s)*t_b)^g)*(g^2)*((exp(s)*t_0)^g-(exp(s)*t_1)^g)*(1+(exp(s)*t_0)^g-(exp(s)* t_1)^g)*(t_a-t_b)+exp(f+(exp(s)* t_1)^g+(exp(s)* t_a)^g)*exp(s)*t_b*(1+g*(-2*(exp(s)* t_0)^g + (exp(s)* t_b)^g)))-exp(f+(exp(s)* t_1)^g+(exp(s)* t_a)^g+(exp(s)* t_b)^g)*(1 + g*(exp(s)* t_0)^g*(-2 + g + g*(exp(s)*t_0)^g))*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g))))
,
"PW.d3Wd3s"=expression((1/(g*(t_a - t_b)))*exp(-s+(exp(s)*t_0)^g-(exp(s)*t_1)^g-(exp(s)*t_a)^g-(exp(s)*t_b)^g)*(g*(exp(f+(exp(s)*t_1)^g+(exp(s)*t_b)^g)*exp(s)*t_a*(1 + g*(3*(exp(s)*t_0)^g*(-1 + g + g*(exp(s)*t_0)^g)-(-1 + g + 3*g*(exp(s)*t_0)^g)*(exp(s)*t_a)^g + g*(exp(s)*t_a)^(2*g))) + exp(s+(exp(s)*t_a)^g + (exp(s)*t_b)^g)*g^3*((exp(s)*t_0)^g-(exp(s)*t_1)^g)*(1 + 3*(exp(s)*t_0)^g + (exp(s)*t_0)^(2*g)-(3+2*(exp(s)*t_0)^g)*(exp(s)*t_1)^g+(exp(s)*t_1)^(2*g))*(t_a-t_b)-exp(f+(exp(s)*t_1)^g+(exp(s)*t_a)^g)*(exp(s)*t_b)*(1+g*(3*(exp(s)*t_0)^g*(-1 + g + g*(exp(s)*t_0)^g)-(-1 + g + 3*g*(exp(s)*t_0)^g)*(exp(s)*t_b)^g+g*(exp(s)*t_b)^(2*g))))-exp(f+(exp(s)*t_1)^g+(exp(s)*t_a)^g+(exp(s)*t_b)^g)*(-1 + g*(exp(s)*t_0)^g*(3 + g*(-3 + g + 3*(-1 + g)*(exp(s)*t_0)^g + g*(exp(s)*t_0)^(2*g))))*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g))))
,
"PW.d3Wd3f"=expression((exp(f-log(t_b-t_a)-s+(exp(s)*t_0)^g)*(1/g)*(pgamma((exp(s)*t_a)^g,1/g, lower.tail=FALSE)*gamma(1/g)-pgamma((exp(s)*t_b)^g,1/g, lower.tail=FALSE)*gamma(1/g))))
)
