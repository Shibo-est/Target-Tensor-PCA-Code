library(dplyr)
library(latex2exp)
library(chron)
library(lubridate)
library(latex2exp)

exp_new_01 <- read.csv("~/Downloads/results00_20240613.csv")
exp_new_01$missing_mse <- exp_new_01$missing_mse/sqrt(exp_new_01$T)
exp_new_01_mean <- exp_new_01 %>% 
  group_by(P, R, T) %>% 
  summarize(missing_m = (mean(missing_mse)))
data <- exp_new_01_mean # 1000 x 400

par(mfcol = c(1,3),
    oma = c(0,0,0,3) + 0.1,
    mar = c(6,5,4,0) + 0.1,
    xpd=NA)
k <- 2
data_mean <- data[data$R == k, ]
plot(x=data_mean[data_mean$P==5,]$T, 
     y=data_mean[data_mean$P==5,]$missing_m,
     ylim = c(0,
              max(data$missing_m)),
     ylab = '',
     xlab = expression("T"),
     cex.lab = 1.4, cex.axis=1.4, cex=1.5, pch=19)
mtext(paste0('R=', k), cex = 1.0, line=0.5)
mtext(expression("||"*hat(bolditalic(B))-bolditalic(B)*"||"), 
      side=2, cex = 1.1, line=2.5)
lines(x=data_mean[data_mean$P==5,]$T, 
      y=data_mean[data_mean$P==5,]$missing_m, lty=1, lwd=2)
points(x=data_mean[data_mean$P==10,]$T, 
       y=data_mean[data_mean$P==10,]$missing_m, col=2, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==10,]$T, 
      y=data_mean[data_mean$P==10,]$missing_m, col=2, lty=2, lwd=2)
points(x=data_mean[data_mean$P==15,]$T, 
       y=data_mean[data_mean$P==15,]$missing_m, col=3, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==15,]$T, 
      y=data_mean[data_mean$P==15,]$missing_m, col=3, lty=3, lwd=2)
points(x=data_mean[data_mean$P==20,]$T, 
       y=data_mean[data_mean$P==20,]$missing_m, col=4, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==20,]$T, 
      y=data_mean[data_mean$P==20,]$missing_m, col=4, lty=4, lwd=2)
points(x=data_mean[data_mean$P==25,]$T, 
       y=data_mean[data_mean$P==25,]$missing_m, col=5, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==25,]$T, 
      y=data_mean[data_mean$P==25,]$missing_m, col=5, lty=5, lwd=2)
#legend("topright", legend=c(expression(s[0]==2), 
#                            expression(s[0]==3), 
#                            expression(s[0]==4)),
#       col=c('black', 2, 3), cex=1.4,
#       lwd = 1)

k <- 3
data_mean <- data[data$R == k, ]
plot(x=data_mean[data_mean$P==5,]$T, 
     y=data_mean[data_mean$P==5,]$missing_m,
     ylim = c(0,
              max(data$missing_m)),
     ylab = '',
     xlab = expression("T"),
     cex.lab = 1.4, cex.axis=1.4, cex=1.5, pch=19)
mtext(paste0('R=', k), cex = 1.0, line=0.5)
mtext("Target Tensor PCA", cex = 1.2, line=2.5)
lines(x=data_mean[data_mean$P==5,]$T, 
      y=data_mean[data_mean$P==5,]$missing_m, lty=1, lwd=2)
points(x=data_mean[data_mean$P==10,]$T, 
       y=data_mean[data_mean$P==10,]$missing_m, col=2, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==10,]$T, 
      y=data_mean[data_mean$P==10,]$missing_m, col=2, lty=2, lwd=2)
points(x=data_mean[data_mean$P==15,]$T, 
       y=data_mean[data_mean$P==15,]$missing_m, col=3, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==15,]$T, 
      y=data_mean[data_mean$P==15,]$missing_m, col=3, lty=3, lwd=2)
points(x=data_mean[data_mean$P==20,]$T, 
       y=data_mean[data_mean$P==20,]$missing_m, col=4, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==20,]$T, 
      y=data_mean[data_mean$P==20,]$missing_m, col=4, lty=4, lwd=2)
points(x=data_mean[data_mean$P==25,]$T, 
       y=data_mean[data_mean$P==25,]$missing_m, col=5, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==25,]$T, 
      y=data_mean[data_mean$P==25,]$missing_m, col=5, lty=5, lwd=2)
legend("bottom", legend=c(expression(p==5), 
                          expression(p==10), 
                          expression(p==15),
                          expression(p==20),
                          expression(p==25)),
       horiz=TRUE, inset=c(0,-0.3), x.intersp=1, bty = "n",
       col=c(1, 2, 3,4,5), lty = c(1,2,3,4,5), cex=1.5, pch=c(19,19,19,19,19), 
       lwd = 2)

k <- 4
data_mean <- data[data$R == k, ]
plot(x=data_mean[data_mean$P==5,]$T, 
     y=data_mean[data_mean$P==5,]$missing_m,
     ylim = c(0,
              max(data$missing_m)),
     ylab = '',
     xlab = expression("T"),
     cex.lab = 1.4, cex.axis=1.4, cex=1.5, pch=19)
mtext(paste0('R=', k), cex = 1.0, line=0.5)
lines(x=data_mean[data_mean$P==5,]$T, 
      y=data_mean[data_mean$P==5,]$missing_m, lty=1, lwd=2)
points(x=data_mean[data_mean$P==10,]$T, 
       y=data_mean[data_mean$P==10,]$missing_m, col=2, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==10,]$T, 
      y=data_mean[data_mean$P==10,]$missing_m, col=2, lty=2, lwd=2)
points(x=data_mean[data_mean$P==15,]$T, 
       y=data_mean[data_mean$P==15,]$missing_m, col=3, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==15,]$T, 
      y=data_mean[data_mean$P==15,]$missing_m, col=3, lty=3, lwd=2)
points(x=data_mean[data_mean$P==20,]$T, 
       y=data_mean[data_mean$P==20,]$missing_m, col=4, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==20,]$T, 
      y=data_mean[data_mean$P==20,]$missing_m, col=4, lty=4, lwd=2)
points(x=data_mean[data_mean$P==25,]$T, 
       y=data_mean[data_mean$P==25,]$missing_m, col=5, cex=1.5, pch=19)
lines(x=data_mean[data_mean$P==25,]$T, 
      y=data_mean[data_mean$P==25,]$missing_m, col=5, lty=5, lwd=2)

