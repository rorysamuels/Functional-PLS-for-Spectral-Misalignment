library(tidyverse)
library(splines2)
library(pls)
source('functions.R')

##### Data Generation #####

# true functional coefficient
beta <- function(t){10*(t-1)^2 + 30*cos(4*pi*t^3)}


# domain for observations
domain <- c(0,1)
a <- domain[1]
b <- domain[2]

# number of measured spectral points for instrument A and B respectively
p1 <- 425
p2 <- 150

# observation grid for instrument A and B respectively
grd_A <- seq(a,b,length.out = p1)
grd_B <- seq(a,b,length.out = p2)

# number of observations to generate (for both instruements)
n <- 500

set.seed(1)
dat <- generate_data(n,nknots = 20,snr = 5, beta_fun = beta, domain = domain, grd_A, grd_B)

# indices for training set
train_ind <- sample(1:n,floor(n*.8))

# train/test response
Y_train <- dat$Y[train_ind]
Y_test <- dat$Y[-train_ind]

# train/test instrument A
X_A_train <- dat$X_A[train_ind,]
X_A_test <- dat$X_A[-train_ind,]

# train/test instrument B
X_B_train <- dat$X_B[train_ind,]
X_B_test <- dat$X_B[-train_ind,]














##### Partial Least Squares Regression #####

plsr_fit <- plsr(Y_train ~ X_A_train, validation = "CV")

# Choose number of componenets which minimize CV PRESS
k <- which.min(plsr_fit$validation$PRESS)

# Coefficients from plsr model
alpha <- (plsr_fit$coefficients[,,k])

# Coefficients closest to observavtion grid of spectra from instrument B
alpha_adj <- alpha[find_closest(grd_A,grd_B)]

# Predictions
pls_preds_A <- as.numeric(X_A_test%*%alpha)
pls_preds_B <- as.numeric(X_B_test%*%alpha_adj)














##### Functional Partial Least Squares without Basis Approx for Data #####


## Set up basis functions for functional coefficient

# Degree = 3 for cubic B-splines
d <- 3

# M = number of basis functions - degree
# Should be chosen via cross-validation (or some criterion)
M <- 10

# Internal knots for basis functions
knots <- seq(a, b, length.out = M+1)[-c(1,M+1)]

# Matrix of basis functions evaluated on observation grid
basis_mat <- bSpline(grd_A, knots = knots, degree = d, intercept = T)



# Number of components for fPLSR
# Should be chosen via cross-validation (or some criterion)
k <- 5
fplsr_fit <- fplsr(Y_train, X_A_train, basis_mat, ncomp = k)


# Predictions
fpls_preds_A <- pred_fplsr(X_A_test,fplsr_fit)
fpls_preds_B <- pred_fplsr(X_B_test,fplsr_fit)















##### Functional Partial Least Squares without Basis Approx for Data #####

# Set up basis functions for functional covariates 

# Degree = 3 for cubic B-splines
d_x <- 3

# M =  number of basis functions - minus degree
M_x <- 13

# Internal knots for basis functions
knots_x <- seq(a, b, length.out = M_x+1)[-c(1,M_x+1)]

x_basis_mat <- bSpline(grd_A, knots = knots_x, degree = d_x, intercept = T)



k <- 5
fplsr_fit_basis <- fplsr(Y_train, X_A_train, basis_mat, x_basis_mat, ncomp = k)


# Predictions
fpls_preds_A_basis <- pred_fplsr(X_A_test,fplsr_fit_basis)
fpls_preds_B_basis <- pred_fplsr(X_B_test,fplsr_fit_basis)














##### Prediction Results #####

pmse(pls_preds_A,Y_test)
pmse(pls_preds_B,Y_test)

pmse(fpls_preds_A,Y_test)
pmse(fpls_preds_B,Y_test)

pmse(fpls_preds_A_basis,Y_test)
pmse(fpls_preds_B_basis,Y_test)

















