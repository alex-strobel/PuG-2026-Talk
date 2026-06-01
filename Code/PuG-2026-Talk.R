# PuG Talk 2026-06-05 Heidelberg ----------------------------------------------- 

## Load required package(s) ----
library(lavaan)
library(papaja)
library(psych)
library(here)


## Required function for simulation correlated variables ----
r_corr <- function(rho, nobs, nvars, precision = 0) {
  
  #' Generate correlated random variables (CC 4.0 BY Alexander Strobel)
  #'
  #' Creates a simulated dataset with a specified correlation value or matrix.
  #'
  #' @param rho A single correlation coefficient or a correlation matrix.
  #' @param nobs Number of observations.
  #' @param nvars Number of variables.
  #' @param precision Required matching precision (0 = no iterative correction, 3 = may take very long).
  #'
  #' @return A data frame with simulated correlated variables.
  
  C <- if (length(rho)==1) { m <- matrix(rho,nvars,nvars); diag(m) <- 1; m } else rho
  A <- chol(C)
  repeat {
    Y <- matrix(rnorm(nobs*nvars), ncol=nvars) %*% A
    if (precision==0 || all(round(cor(Y),precision)==round(C,precision))) break
  }
  as.data.frame(Y)
}

## Setup ----
set.seed(242)
n <- 256                      # sample size
rho_true <- c(.5, .3, .1)     # true correlations

par(mfrow = c(3, 3))
out = NULL

# simulate data for Fig 1 
for (i in 1:3) {
df = r_corr(rho_true[i], n, 2, precision = 3)

T1 <- df[,1]
T2 <- df[,2]

out = cbind(out, T1, T2)

cor(T1, T2)

add_measurement_error <- function(true_score, reliability) {
  error_sd <- sd(true_score) * sqrt((1 - reliability) / reliability)
  true_score + rnorm(length(true_score), mean = 0, sd = error_sd)
}

X_high <- add_measurement_error(T1, reliability = 0.8)
Y_high <- add_measurement_error(T2, reliability = 0.8)
out = cbind(out, X_high, Y_high)
cor(X_high, Y_high)

X_low <- add_measurement_error(T1, reliability = 0.6)
Y_low <- add_measurement_error(T2, reliability = 0.6)
out = cbind(out, X_low, Y_low)
cor(X_low, Y_low)

## Fig 1 ----

plot(T1, T2,
     main = "True Scores",
     xlab = "X", ylab = "Y",
     xlim = c(-4.5,4.5), ylim = c(-4.5,4.5), 
     pch = 16, col = rgb(0,0,0,0.4))
abline(lm(T2 ~ T1), lwd = 2)
text(3,-4,bquote(italic(r) == .(printnum(cor(T1, T2), gt1=F))), cex = 1.5)

plot(X_high, Y_high,
     main = "Reliability = .80",
     xlab = "X", ylab = "Y",
     xlim = c(-4.5,4.5), ylim = c(-4.5,4.5), 
     pch = 16, col = rgb(0,0,1,0.4))
abline(lm(Y_high ~ X_high), lwd = 2, col = 4)
text(3,-4,bquote(italic(r) == .(printnum(cor(X_high, Y_high), gt1=F))), cex = 1.5, col = 4)

plot(X_low, Y_low,
     main = "Reliability = .60",
     xlab = "X", ylab = "Y",
     xlim = c(-4.5,4.5), ylim = c(-4.5,4.5), 
     pch = 16, col = rgb(1,0,0,0.4))
abline(lm(Y_low ~ X_low), lwd = 2, col = 2)
text(3,-4,bquote(italic(r) == .(printnum(cor(X_low, Y_low), gt1=F))), cex = 1.5, col = 2)

}

par(mfrow = c(1, 1))

# data.frame containing all variables used in Fig 1
colnames(out) = c("T1.5", "T2.5", "X.5.8", "Y.5.8", "X.5.6", "Y.5.6",
                  "T1.3", "T2.3", "X.3.8", "Y.3.8", "X.3.6", "Y.3.6",
                  "T1.1", "T2.1", "X.1.8", "Y.1.8", "X.1.6", "Y.1.6")

## Fig 2 ----

# use unrestricted X and Y variables for Rel = .8 
norr.3.8 = data.frame(out[,9:10])
cor_norr.3.8 = cor(norr.3.8)[1,2]

# use X restricted and Y variables for Rel = .8 
rr.3.8_X=data.frame(out[which(out[,9]>0),9:10])
cor_rr.3.8_X = cor(rr.3.8_X)[1,2]

# use X restricted and Y restricted for Rel = .8 
rr.3.8_XY=data.frame(out[which(out[,9]>0 & out[,10]>0),9:10])
cor_rr.3.8_XY = cor(rr.3.8_XY)[1,2]

par(mfrow = c(3,1))

plot(rr.3.8_X$X.3.8, rr.3.8_X$Y.3.8,
     main = "Reliability .8\nRange restricted for X",
     xlab = "X restricted", ylab = "Y",
     xlim = c(-4.5,4.5), ylim = c(-4.5,4.5), 
     pch = 16, col = rgb(0,0,0,0.4))
abline(lm(Y.3.8 ~ X.3.8, rr.3.8_X), lwd = 2, col = 4)
abline(v = 0, col = 2, lty = 5)
text(3,-4,bquote(italic(r) == .(printnum(cor_rr.3.8_X, gt1=F))), cex = 1.5, col = 4)

plot(norr.3.8[,1], norr.3.8[,2],
     main = "Reliability .8\nRange unrestricted",
     xlab = "X", ylab = "Y",
     xlim = c(-4.5,4.5), ylim = c(-4.5,4.5), 
     pch = 16, col = rgb(0,0,0,0.4))
abline(lm(Y.3.8 ~ X.3.8, norr.3.8), lwd = 2, col = 4)
text(3,-4,bquote(italic(r) == .(printnum(cor_norr.3.8, gt1=F))), cex = 1.5, col = 4)

plot(rr.3.8_XY$X.3.8, rr.3.8_XY$Y.3.8,
     main = "Reliability .8\nRange restricted for X and Y",
     xlab = "X restricted", ylab = "Y restricted",
     xlim = c(-4.5,4.5), ylim = c(-4.5,4.5), 
     pch = 16, col = rgb(0,0,0,0.4))
abline(lm(Y.3.8 ~ X.3.8, rr.3.8_XY), lwd = 2, col = 4)
abline(v = 0, h = 0, col = 2, lty = 5)
text(3,-4,bquote(italic(r) == .(printnum(cor_rr.3.8_XY, gt1=F))), cex = 1.5, col = 4)

## Fig 3 ----

# locate top level folder
here::i_am("flag_root_for_PuG-2026-Talk.txt")

# load data
load(here("Code", "df.RData"))

# correlations of NFC /w lP3a at 3 pos x 4 blocks
corr.test(df_N_LP3a[,1:12], df_N_LP3a[,13])

# correlation of NFC /w lP3a averaged over 3 pos x 4 blocks
LP3a_avg = rowMeans(df_N_LP3a[,1:12])
nfc_avg = df_N_LP3a$nfc
corr.test(LP3a_avg, nfc_avg)$r

# SEM with blocks as first oder latent variables and overall lP3a as second order lat. var.

m1 = "
N_LP3a_B1_cz ~~ 1*N_LP3a_B1_cz
N_LP3a_B2_cz ~~ 1*N_LP3a_B2_cz
N_LP3a_B3_cz ~~ 1*N_LP3a_B3_cz
N_LP3a_B4_cz ~~ 1*N_LP3a_B4_cz

B1 =~ N_LP3a_B1_fz + N_LP3a_B1_cz + 1*N_LP3a_B1_pz
B2 =~ N_LP3a_B2_fz + N_LP3a_B2_cz + 1*N_LP3a_B2_pz
B3 =~ N_LP3a_B3_fz + N_LP3a_B3_cz + 1*N_LP3a_B3_pz
B4 =~ N_LP3a_B4_fz + N_LP3a_B4_cz + 1*N_LP3a_B4_pz

# the following is necessary to omit negative variances 
B4 ~~ 1*B4
B3 ~~ 1*B3
B2 ~~ 1*B2
B1 ~~ 1*B1

B =~ B1 + B2 + B3 + B4
NFC =~ nfc
" 

f1 = sem(m1, df_N_LP3a)

summary(f1, fit.measures=T, standardized = T)
# bad fit, but used for illustration purpose only 

# correlation of NFC /w latent lP3a variables 
df_N_LP3a_f1 = lavPredict(f1)
corr.test(df_N_LP3a_f1)
