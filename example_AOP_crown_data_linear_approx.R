# This script demonstrates functional PLS with misaligned data (AOP Crown data-set).
# The predictive performance is compared with traditional PLS, where the 
# misaligned data is "re-aligned" using linear approximations.

library(tidyverse)
library(splines2)
library(pls)
library(here)
setwd(here())
source('functions.R')


##### Read and pre-process data #####

# Read in data
x_data <- read_csv("data/AOP_Crown.csv")
site_trait_data <- read_csv("data/site_trait_data.csv")
y_data <- drop_na(tibble(ID = site_trait_data$SampleSiteID, y = site_trait_data$d15N))
xy_data <- inner_join(y_data,x_data)
refl_data <- xy_data[,20:445]
neon_grd <- read_rds("data/neon_grd.RDS")
bad_bands <- c(1:8,192:205,284:327,417:ncol(refl_data))
good_grd <- neon_grd[-bad_bands]

# Reflectance data
X <- matrix(unname(unlist(refl_data[,-bad_bands])), nrow = nrow(refl_data))
# Response
y <- xy_data$y

# Number of observations
n <- nrow(X)
# Number of original spectral points
p <- ncol(X)
# Original observation grid (scaled to unit interval)
grd <- seq(0,1,length.out = p)

# Index of every other point
sel <- seq(1,p,by = 2)

# Observations from instrument A
X_A <- X[,sel]
# Instrument A observation grd
grdA <- grd[sel]
# Number of spectral points for instrument A
pA <- length(grdA)

# Observations from instrument B
X_B <- X[,-sel]
# Instrumnet B observation grd
grdB <- grd[-sel]
# Number of spectral points for instrument B
pB <- length(grdB)

# Draw straight lines to approximate X_B on grdA
X_B_approx <- 
t(
apply(X_B,1,function(s){
  approxfun(grdB,s)(grdA[-1])
})
)





# Train/Test Split
set.seed(10)
train_ind <- sample(1:n,floor(.8*n))

X_A_train <- X_A[train_ind,]
X_B_train <- X_B[train_ind,]
X_B_approx_train <- X_B_approx[train_ind,]

X_A_test <- X_A[-train_ind,]
X_B_test <- X_B[-train_ind,]
X_B_approx_test <- X_B_approx[-train_ind,]

y_train <- y[train_ind]
y_test <- y[-train_ind]









##### Partial Least Squares Regression #####

#plsr_fit <- plsr(y_train~X_A_train, validation = "CV")

# Choose number of componenets which minimize CV PRESS
#k <- which.min(plsr_fit$validation$PRESS)

k <- 14
plsr_fit <- plsr(y_train~X_A_train, ncomp = k)

# Coefficients from plsr model
alpha <- (plsr_fit$coefficients[,,k])

# Predictions
pls_preds_A <- as.numeric(X_A_test%*%alpha)
pls_preds_B_approx <- as.numeric(X_B_approx_test%*%alpha[-1])

pmse(y_test,pls_preds_A)
pmse(y_test,pls_preds_B_approx)








##### Functional Partial Least Squares without Basis Approx for Data #####


fplsr_fit <- fplsr(y_train, X_A_train, grd = grdA, M = 30, ncomp = k)


# Predictions
fpls_preds_A <- pred_fplsr(X_A_test, grdA, fplsr_fit)
fpls_preds_B <- pred_fplsr(X_B_test[,-pB], grdB[-pB], fplsr_fit)






##### Functional Partial Least Squares with Basis Approx for Data #####


fplsr_fit_basis <- fplsr(y_train, X_A_train, grd = grdA, M = 30, M_x = 25, ncomp = k)


# Predictions
fpls_preds_A_basis <- pred_fplsr(X_A_test, grdA, fplsr_fit_basis)
fpls_preds_B_basis <- pred_fplsr(X_B_test[,-pB], grdB[-pB], fplsr_fit_basis)



pmse(fpls_preds_A_basis,y_test)
pmse(fpls_preds_B_basis,y_test)






##### Results #####
##### Results #####
#pmse(pls_preds_A,y_test)
pmse_lin_approx <- pmse(pls_preds_B_approx,y_test)


#pmse(fpls_preds_A,y_test)
pmse_fpls <- pmse(fpls_preds_B,y_test)


#pmse(fpls_preds_A_basis,y_test)
pmse_fpls_basis <- pmse(fpls_preds_B_basis,y_test)


# 66% reduction in pmse using fpls over linear approximation
#(pmse_fpls - pmse_lin_approx)/pmse_lin_approx
# 67% reduction in pmse using fpls + basis expansion over linear approximation
#(pmse_fpls_basis - pmse_lin_approx)/pmse_lin_approx




