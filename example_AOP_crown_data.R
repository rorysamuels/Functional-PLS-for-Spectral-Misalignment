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

# Original observation grid (scaled to unit interval)
pA <- ncol(X)
grdA <- seq(0,1,length.out = pA)

# Misaligned observtion grid
pB <- 100
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

#plsr_fit <- plsr(y_train~X_A_train, validation = "CV")

# Choose number of componenets which minimize CV PRESS
#k <- which.min(plsr_fit$validation$PRESS)

k <- 17
plsr_fit <- plsr(y_train~X_A_train, ncomp = k)

# Coefficients from plsr model
alpha <- (plsr_fit$coefficients[,,k])

# Coefficients closest to observavtion grid of spectra from instrument B
alpha_adj <- alpha[find_closest(grdA,grdB)]

# Predictions
pls_preds_A <- as.numeric(X_A_test%*%alpha)
pls_preds_B <- as.numeric(X_B_test%*%alpha_adj)









##### Functional Partial Least Squares without Basis Approx for Data #####


fplsr_fit <- fplsr(y_train, X_A_train, grd = grdA, M = 30, ncomp = k)


# Predictions
fpls_preds_A <- pred_fplsr(X_A_test, grdA, fplsr_fit)
fpls_preds_B <- pred_fplsr(X_B_test, grdB, fplsr_fit)









##### Functional Partial Least Squares with Basis Approx for Data #####


fplsr_fit_basis <- fplsr(y_train, X_A_train, grd = grdA, M = 30, M_x = 25, ncomp = k)


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











