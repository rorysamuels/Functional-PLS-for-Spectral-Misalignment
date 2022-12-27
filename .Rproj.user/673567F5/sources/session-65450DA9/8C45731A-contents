
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


# Get basis coefficients for data matrix 'X' (fit w/ OLS -- assumes no instrument noise)
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
    mat <- predict(basis1,t)
    t(mat[,i])
  }
  
  b_j <- function(t)
  {
    mat <- predict(basis2,t)
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
    mat <- predict(basis,t)
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
  
  
  x_A <- A%*%t(predict(basis_mat,gA))
  x_B <- A%*%t(predict(basis_mat,gB))
  
  
  y <- y0 + rnorm(n, mean = 0, sd = sd(y0)/sqrt(snr))
  
  list(X_A = x_A, X_B = x_B, Y = y)
}




##### Functional Partial Least Squares #####

# Fit a functional partial least squares model with or without basis expansion for observations

## Inputs:
# y: response
# x: matrix of observations
# basis_beta: basis function matrix for functional coefficient of class matrix/bSpline2/splines2
# basis_x: basis function matrix for observations of class matrix/bSpline2/splines2
# ncomp: number of components to use in the partial least squares fitting

## Outputs:
# beta_fun: estimated functional coefficient
# alpha_star: coefficients from partial least squares fitting
# basis_beta: basis function matrix for functional coefficient
# basis_x: basis function matrix for observations
fplsr <- function(y,x, basis_beta, basis_x = NULL, ncomp = 5){
  
  domain <- attributes(basis_beta)$Boundary.knots
  a <- domain[1]
  b <- domain[2]
  
  n <- length(y)
  p <- ncol(x)
  grd <- seq(a,b,length.out = p)
  
  if(is.null(basis_x)) U <- (x%*%splines2:::predict.bSpline2(basis_beta,grd))*((b-a)/(p))
  else{
    
    M <- matrix(NA,ncol(basis_x),ncol(basis_beta))
    for(i in 1:ncol(basis_x)){
      for(j in 1:ncol(basis_beta)){
        M[i,j] <- basisXbasis(basis_x,basis_beta,i,j)
      }
    }
    
    C <- get_basis_coefs(x,basis_x)
    U <- C%*%M
  }
  
  pls_mod <- plsr(y~U, ncomp = ncomp)
  alpha_star <- (pls_mod$coefficients[,,ncomp])
  
  
  beta_fpls <- function(t){
    
    as.numeric(predict(basis_beta,t)%*%alpha_star)
    
  }
  
  list(beta_fun = beta_fpls, alpha_star = alpha_star,
       basis_beta = basis_beta, basis_x = basis_x)
  
}



# Predict response for new x with fitted functional partial least squares model

## Inputs:
# x: matrix of new observations
# fplsmod: list object returned from fplsr()

## Outputs:
# vector of predicted responses
pred_fplsr <- function(x, fplsmod){
  
  domain <- attributes(fplsmod$basis_beta)$Boundary.knots
  a <- domain[1]
  b <- domain[2]
  
  p <- ncol(x)
  grd <- seq(a,b,length.out = p)
  
  beta <- fplsmod$beta_fun
  
  if(is.null(fplsmod$basis_x)) x%*%beta(grd)*((b-a)/(p))
  else{
    
    bxattr <- attributes(fplsmod$basis_x)
    nbasis_x <- bxattr$dim[2]
    basis_new_x <- bSpline(grd,knots = bxattr$knots,
                           degree = bxattr$degree,
                           intercept = T)
    
    
    ivec = matrix(NA,nbasis_x,1)
    for(j in 1:nbasis_x) ivec[j] <- betaXbasis(beta,basis_new_x,j)
    C <- get_basis_coefs(x,basis_new_x)
    C%*%ivec 
    
  }
  
}











