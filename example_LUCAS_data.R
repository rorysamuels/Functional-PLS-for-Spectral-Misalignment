library(tidyverse)
library(splines2)
library(pls)
source('functions.R')


##### Read and pre-process data #####

# Read in data


data <- read_csv("/home/roryjake/Desktop/LUCAS_ref.csv")
y <- data$pH_h2o 
X <- as.matrix(data[,1:425])


# Number of observations
n <- nrow(X)

# Original observation grid (scaled to unit interval)
pA <- ncol(X)
grdA <- seq(0,1,length.out = pA)

# Misaligned observtion grid
pB <- 200
grdB <- seq(0,1,length.out = pB)

# Basis expansion for reflectance data (done to create misaligned observations)

# Degree for cubic B-splines
d <- 3
# Number of internal knots - 1
M <- 50
# Internal knots for basis functions
knots <- seq(0, 1, length.out = M+1)[-c(1,M+1)]

# Matrix of basis functions evaluated on original observation grid
basis_spectra <- bSpline(grdA, knots = knots, degree = d, intercept = T)

# Coefficients from basis expansion
A <- get_basis_coefs(X,basis_spectra)

# Spectral points matching original observation grid (instrument A)
X_A <- (A%*%t(predict(basis_spectra,grdA)))
# Spectral points on mislaigned grid (instrument B)
X_B <- (A%*%t(predict(basis_spectra,grdB)))

# Train/Test Split
set.seed(10)
train_ind <- sample(1:n,floor(.8*n))

X_A_train <- X_A[train_ind,]
X_B_train <- X_B[train_ind,]

X_A_test <- X_A[-train_ind,]
X_B_test <- X_B[-train_ind,]


y_train <- y[train_ind]
y_test <- y[-train_ind]




##### Partial Least Squares Regression #####
k <- 10
plsr_fit <- plsr(y_train~X_A_train, ncomp = k)

# Coefficients from plsr model
alpha <- (plsr_fit$coefficients[,,k])

# Coefficients closest to observavtion grid of spectra from instrument B
alpha_adj <- alpha[find_closest(grdA,grdB)]

# Predictions
pls_preds_A <- as.numeric(X_A_test%*%alpha)
pls_preds_B <- as.numeric(X_B_test%*%alpha_adj)






##### Functional Partial Least Squares without Basis Approx for Data #####


fplsr_fit <- fplsr(y_train, X_A_train, grdA, M = 25, ncomp = k)


# Predictions
fpls_preds_A <- pred_fplsr(X_A_test, grdA, fplsr_fit)
fpls_preds_B <- pred_fplsr(X_B_test, grdB, fplsr_fit)







##### Functional Partial Least Squares without Basis Approx for Data #####

fplsr_fit_basis <- fplsr(y_train, X_A_train, grdA, M = 25, M_x = 25, ncomp = k)


# Predictions
fpls_preds_A_basis <- pred_fplsr(X_A_test, grdA, fplsr_fit_basis)
fpls_preds_B_basis <- pred_fplsr(X_B_test, grdB, fplsr_fit_basis)



##### Results #####
pmse(pls_preds_A,y_test)
pmse(pls_preds_B,y_test)


pmse(fpls_preds_A,y_test)
pmse(fpls_preds_B,y_test)


pmse(fpls_preds_A_basis,y_test)
pmse(fpls_preds_B_basis,y_test)









