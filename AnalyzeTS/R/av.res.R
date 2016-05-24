av.res <-
function (Y = NULL, F = NULL, E = NULL, r = 3) 
{
    if (!is.null(Y) & !is.null(F) & !is.null(E)) 
        stop("More members!")
    if (!is.null(Y) & !is.null(F)) 
        Model = c("YF")
    else if (!is.null(Y) & !is.null(E)) 
        Model = c("YE")
    else stop("Unknow members!")
    if (Model == "YF") {
        if (!is.data.frame(Y) || !is.data.frame(F)) 
            stop("Cac tham so cua ham phai la data frame!")
        if (length(Y) != 1) 
            stop("Data frame gia tri quan sat phai co 1 cot!")
        if (dim(Y)[1] != dim(F)[1]) 
            stop("So phan tu cac chuoi khong bang nhau!")
    }
    else if (Model == "YE") {
        if (!is.data.frame(Y) || !is.data.frame(E)) 
            stop("Cac tham so cua ham phai la data frame!")
        if (length(Y) != 1) 
            stop("Data frame gia tri quan sat phai co 1 cot!")
        if (dim(Y)[1] != dim(E)[1]) 
            stop("So phan tu cac chuoi khong bang nhau!")
    }
    if (Model == "YF") 
        n <- dim(F)[2]
    else n <- dim(E)[2]
    kq <- matrix(c(1:(n * 7)), nrow = 7)
    rownames(kq) <- c("ME", "MAE", "MPE", "MAPE", "MSE", "RMSE", 
        "U")
    if (Model == "YF") 
        name.models <- dimnames(F)[[2]]
    if (Model == "YE") 
        name.models <- dimnames(E)[[2]]
    colnames(kq) <- name.models
    zero <- 0
    for (models in 1:n) {
        Yt <- Y[, 1]
        if (Model == "YF") {
            Ft <- F[, models]
            et <- Yt - Ft
        }
        else {
            Ft <- Yt - E[, models]
            et <- E[, models]
        }
        ne = length(na.omit(et))
        t1 = 0
        t2 = 0
        t3 = 0
        t4 = 0
        t5 = 0
        t6 = 0
        t7.1 = 0
        t7.2 = 0
        for (i in 1:length(Yt)) {
            if (!is.na(et[i])) {
                t1 <- t1 + et[i]
                t2 <- t2 + abs(et[i])
                t3 <- t3 + ((et[i]/Yt[i]) * 100)
                t4 <- t4 + ((abs(et[i])/Yt[i]) * 100)
                t5 <- t5 + et[i] * et[i]
                t6 <- t6 + et[i] * et[i]
                t7.1 <- t7.1 + (Yt[i] - Ft[i]) * (Yt[i] - Ft[i])
                if (i > 1) 
                  t7.2 <- t7.2 + (Yt[i] - Yt[i - 1]) * (Yt[i] - 
                    Yt[i - 1])
            }
        }
        kq[1, models] = t1/ne
        kq[2, models] = t2/ne
        kq[3, models] = t3/ne
        kq[4, models] = t4/ne
        kq[5, models] = t5/ne
        kq[6, models] = sqrt(t6/ne)
        kq[7, models] = sqrt(t7.1)/sqrt(t7.2)
        for (i in 1:n) zero = zero + sum(Yt == 0)
    }
    if (zero > 0 & n > 1) 
        print("Khong the thuc hien xep loai!")
    if (zero == 0) {
       
        loi1=0
        loi2=0
        loi3=0
        loi4=0
        loi5=0
        loi6=0
        loi7=0

        for(i in 1:n){
            if (is.na(kq[1, i])) loi1 <- loi1 + 1
            if (is.na(kq[2, i])) loi2 <- loi2 + 1
            if (is.na(kq[3, i])) loi3 <- loi3 + 1
            if (is.na(kq[4, i])) loi4 <- loi4 + 1
            if (is.na(kq[5, i])) loi5 <- loi5 + 1
            if (is.na(kq[6, i])) loi6 <- loi6 + 1
            if (is.na(kq[7, i])) loi7 <- loi7 + 1
        }

        xl.nho <- rep(0,7)
        for (i in 1:n) {
            if(loi1==0) {if (kq[1, i] == min(kq[1, ])) xl.nho[1] <- i} else xl.nho[1]<-0
            if(loi2==0) {if (kq[2, i] == min(kq[2, ])) xl.nho[2] <- i} else xl.nho[2]<-0
            if(loi3==0) {if (kq[3, i] == min(kq[3, ])) xl.nho[3] <- i} else xl.nho[3]<-0
            if(loi4==0) {if (kq[4, i] == min(kq[4, ])) xl.nho[4] <- i} else xl.nho[4]<-0
            if(loi5==0) {if (kq[5, i] == min(kq[5, ])) xl.nho[5] <- i} else xl.nho[5]<-0
            if(loi6==0) {if (kq[6, i] == min(kq[6, ])) xl.nho[6] <- i} else xl.nho[6]<-0
            if(loi7==0) {if (kq[7, i] == min(kq[7, ])) xl.nho[7] <- i} else xl.nho[7]<-0
        }
        kq <- round(kq, r)
        tenmh <- 1:7
        for (ten in 1:7) if(xl.nho[ten]!=0) tenmh[ten] <- colnames(kq)[xl.nho[ten]] else tenmh[ten]<-"NA"
        kq <- data.frame(kq, min.model = tenmh)
    }
    if (zero == 0) 
        if (dim(kq)[2] == 2) 
            kq <- t(kq[1])
    if (zero > 0) 
        if (dim(kq)[2] == 1) 
            kq <- t(kq[, 1])
    kq
}
