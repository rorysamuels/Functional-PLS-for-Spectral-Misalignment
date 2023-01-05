# Dependencies:
# library(tidyverse)
# library(pls)
# library(splines2)

##### Utils #####

# Get indices of grdA values that are closest to the values in grdB
find_closest <- function(grdA,grdB){
  
  select_grdpoint <- function(x){
    which.min(abs(grdA-x))
  }
  
  Vectorize(select_grdpoint)(grdB)
  
}




# Evaluate inner product between 'f1' and 'f2' from a to b
inner_prod <- function(f1,f2,a,b){
  
  integrate(function(t) {f1(t)*f2(t)},a,b)$value
  
  
}


# Get basis coefficients for data matrix 'X'

## Notes:
# fit w/ OLS -- i.e. assumes no instrument noise
get_basis_coefs <- function(X,basis){
  
  t(
    apply(X,1,function(x){
      lm(x~0+basis)$coefficients
    })
  )
  
}


# Evaluate inner product between ith function in 'basis1' and jth function in 'basis2'
basisXbasis <- function(basis1,basis2,i,j){
  
  domain <- attributes(basis1)$Boundary.knots
  a <- domain[1]
  b <- domain[2]
  
  
  b_i <- function(t)
  {
    mat <- splines2:::predict.bSpline2(basis1,t)
    t(mat[,i])
  }
  
  b_j <- function(t)
  {
    mat <- splines2:::predict.bSpline2(basis2,t)
    t(mat[,j])
  }
  
  
  inner_prod(b_i,b_j,a,b)
  
}




# Evaluate inner product between jth function (column) in 'basis' and 'f'
betaXbasis <- function(f,basis,j){
  
  domain <- attributes(basis)$Boundary.knots
  a <- domain[1]
  b <- domain[2]
  
  b_j <- function(t)
  {
    mat <- splines2:::predict.bSpline2(basis,t)
    z <- t(mat[,j])
  }
  
  inner_prod(f,b_j,a,b)
  
}




# Calculate prediction mean squared error
pmse <- function(y,yhat){
  
  mean((y-yhat)^2)
  
}



##### Data Generation Function #####

## Inputs: 
# n: number of observations (per instrument)
# nknots: number of knots for generating random lincoms of b-spline basis functions
# beta_fun: true functional coefficient
# snr: signal to noise ratio
# domain: endpoints for observation grid domain
# gA: observation grid for instrument A
# gB: observation grid for instrument B

## Returns:
# X_A: generated spectra observed on gA
# X_B: generated spectra observed on gB
# Y: generated response
generate_data <- function(n, nknots, beta_fun, snr, domain, gA, gB){
  
  a <- domain[1]
  b <- domain[2]
  knots  <- seq(a, b, length.out = nknots)[-c(1,nknots)]
  grd <- seq(a,b,by = .1)
  basis_mat <- bSpline(grd, knots = knots, intercept = TRUE)
  nbasis <- ncol(basis_mat)
  
  A <- matrix(rnorm(n*nbasis),n,nbasis)
  ivec = matrix(NA,nbasis,1)
  for(j in 1:nbasis) ivec[j] <- betaXbasis(beta_fun,basis_mat,j)
  y0 <- A %*% ivec
  
  
  x_A <- A%*%t(splines2:::predict.bSpline2(basis_mat,gA))
  x_B <- A%*%t(splines2:::predict.bSpline2(basis_mat,gB))
  
  
  y <- y0 + rnorm(n, mean = 0, sd = sd(y0)/sqrt(snr))
  
  list(X_A = x_A, X_B = x_B, Y = y, sd = sd(y0)/sqrt(snr))
}




##### Functional Partial Least Squares #####

# Fit a functional partial least squares model with or without basis expansion for observations

## Notes:
# Only supported for cubic B-splines
# Assumes equally spaced observation grid

## Inputs:
# y: response
# x: matrix of observations
# M: breakpoints for coefficient basis functions (M + 1 knots)
# M_x: breakpoints for observation basis functions (M_x + 1 knots)
# ncomp: number of components to use in partial least squares fit

## Outputs:
# beta_fun: estimated functional coefficient
# alpha_star: coefficients from partial least squares fitting (y~U)
# basis_beta: basis function matrix for functional coefficient
# basis_x: basis function matrix for observations
# grd: input grid used for fplsr fit
fplsr <- function(y, x, grd, M, M_x = NULL, ncomp = 5){
  
  if(length(grd) != ncol(x))stop("Dimension mismatch: observation grid and measurements")
  if(length(y) != nrow(x))stop("Dimension mismatch: y and x")
  
  a <- min(grd)
  b <- max(grd)
  n <- nrow(x)
  p <- ncol(x)
  
  knots <- seq(a, b, length.out = M+1)[-c(1,M+1)]
  basis_beta <- bSpline(grd, knots = knots, intercept = T)
  
  if(is.null(M_x)){
    U <- (x%*%splines2:::predict.bSpline2(basis_beta,grd))*((b-a)/(p))
    basis_x <- NULL
  } 
  else{
    knots_x <- seq(a, b, length.out = M_x+1)[-c(1,M_x+1)]
    basis_x <- bSpline(grd, knots = knots_x, intercept = T)
    B <- matrix(NA,ncol(basis_x),ncol(basis_beta))
    for(i in 1:ncol(basis_x)){
      for(j in 1:ncol(basis_beta)){
        B[i,j] <- basisXbasis(basis_x,basis_beta,i,j)
      }
    }
    
    C <- get_basis_coefs(x,basis_x)
    U <- C%*%B
  }
  
  pls_mod <- plsr(y~U, ncomp = ncomp)
  alpha_star <- (pls_mod$coefficients[,,ncomp])
  
  
  beta_fpls <- function(t){
    
    as.numeric(splines2:::predict.bSpline2(basis_beta,t)%*%alpha_star)
    
  }
  
  list(beta_fun = beta_fpls, alpha_star = alpha_star,
       basis_beta = basis_beta, basis_x = basis_x, grd = grd)
  
}



# Predict response for new x with fitted functional partial least squares model

## Notes:
# Assumes equally spaced observation grid

## Inputs:
# x: matrix of new observations
# grd: observation grid for new observations
# fplsmod: list object returned from fplsr()

## Outputs:
# vector of predicted responses
pred_fplsr <- function(x, grd, fplsmod){
  
  a <- min(grd)
  b <- max(grd)
  n <- nrow(x)
  p <- ncol(x)
  
  beta <- fplsmod$beta_fun
  
  if(is.null(fplsmod$basis_x)) x%*%beta(grd)*((b-a)/(p))
  else{
    
    basis_x <- fplsmod$basis_x
    nbasis_x <- ncol(basis_x)
    
    ivec = matrix(NA,nbasis_x,1)
    for(j in 1:nbasis_x) ivec[j] <- betaXbasis(beta,basis_x,j)
    C <- get_basis_coefs(x,splines2:::predict.bSpline2(basis_x,grd))
    
    C%*%ivec 
    
  }
  
}











