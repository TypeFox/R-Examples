library(faoutlier)

load('../R/objects.Rdata')

plt_on <- function(name, height = 8, width=8)
    pdf(name, width = width, height = width)

plt_off <- function() dev.off()

#Figure 1 

plt_on('MD1.pdf'); plot(MD1, main = 'Robust MDs, holzinger', ylab='Mahalanobis distance'); plt_off()
plt_on('MD2.pdf'); plot(MD2, main = 'Robust MDs, holzinger.outlier', ylab='Mahalanobis distance'); plt_off()

#Figure 2 
plt_on('std_res1.pdf'); plot(OR1, restype = 'std_res', main = 'Standardized residuals, holzinger'); plt_off()
plt_on('std_res2.pdf'); plot(OR2, restype = 'std_res', main = 'Standardized residuals, holzinger.outlier'); plt_off()

#Figure 3 
plt_on('GOF1.pdf'); plot(GOF1, main = 'Goodness of fit distance, holzinger'); plt_off() 
plt_on('GOF2.pdf'); plot(GOF2, main = 'Goodness of fit distance, holzinger.outlier'); plt_off()

#Figure 4 
plt_on('gCD1.pdf'); plot(gCD1, main = 'gCD, EFA holzinger'); plt_off() 
plt_on('gCD2.pdf'); plot(gCD2, main = 'gCD, EFA holzinger.outlier'); plt_off() 
plt_on('gCD3.pdf'); plot(gCD3, main = 'gCD, CFA holzinger'); plt_off() 
plt_on('gCD4.pdf'); plot(gCD4, main = 'gCD, CFA holzinger.outlier'); plt_off()

#Figure 5  
plt_on('FS1.pdf'); plot(FS1, main = 'Forward search, holzinger', stat = 'RMR'); plt_off() 
plt_on('FS2.pdf'); plot(FS2, main = 'Forward search, holzinger.outlier', stat = 'RMR'); plt_off() 

                    
