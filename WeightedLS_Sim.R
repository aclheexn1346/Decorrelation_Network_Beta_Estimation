

# Hangjian Code
# Testing code to see what it does
source("helperFunc.R")
source("libraries.R")
# Toeplitz function to make a covariance with toeplitz structure between rows


gen.B <- function(p, b.mag = 1, s0= 2*p, seed = 482, lower.thresh = 0.1){
  # Generate some random beta matrix with ordering such that it
  # is an upper diagonal matrix
  if(p <= 5 ) stop("p is too small!")
  set.seed(seed*2)
  invisible(capture.output(bb <- randomDag(seed = seed, numNodes = p,numEdges = s0)))
  b <-  get_adjmat_from_fges(bb$edges, length(bb$nodes), bb$nodes)
  b[b!=0] = runif(length(bb$edges), lower.thresh, b.mag)*(2*rbinom(length(bb$edges),1,0.5)-1)
  realp <- pp <- p
  dimnames(b) <- list(as.character(1:realp), as.character(1:realp))
  return(list(b=b, s0=s0, realp = realp, pp = pp))
}


sim_X <- function(vers, p, n, omg.sq, sig, b){
  #' simulate X from network DAG given its parameters
  #'
  #' \code{sim_X} returns X and E matrices generated from the network DAG
  #'
  #' @param
  #'
  set.seed(vers)
  eps_mat <- matrix(0, n, p)
  eps_mat[,1] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[1]*sig)
  eps_mat[,2] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[2]*sig)

  X <- matrix(0, n, p)

  X[,1] <- eps_mat[,1]
  X[,2] <- X[,1]*b[1,2] + eps_mat[,2]

  for(i in 3:p) {
    eps_mat[, i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[i]*sig) # Sigma is covariance between the n samples
    X[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, b[1:i-1,i], "*")) + eps_mat[,i]
    #if (i %% 50 == 0)
    #  cat("Getting ", i, "th column of X. \n" )
  }
  dimnames(X) <- list(NULL, as.character(1:p))
  return(list(X = X, eps_mat = eps_mat))
}
# try to draw out what I need to do tmr


# Need to write successive row regression code

L_successive_reg_row = function(X){ # Input is original X data matrix
  t_data = t(X)
  num_rows = nrow(X)
  num_cols = ncol(X)
  VarsY = c()
  VarsX = c()
  for(i in 1:(num_rows)){
    VarsY = c(VarsY, paste("V",i, sep = ""))
  }
  for(j in 2:num_rows){
    VarsX = c(VarsX, paste(VarsY[j:length(VarsY)], collapse = "+"))
  }
  # Since last row is not regressed on anything, residual is itself
  VarsY = VarsY[-length(VarsY)]
  
  rowModelsList <- lapply(paste(paste(VarsY, VarsX, sep = "~"), -1, sep = ""), as.formula)
  rowModelsResults <- lapply(rowModelsList, function(x) fastLm(x, data = as.data.frame(t_data)))  
  # Create D which is diagonal matrix of entries as variance of residuals
  rowModelsResiduals <- lapply(rowModelsResults, function(x) residuals(x))
  D_mat = diag(1, num_rows, num_rows)
  for(j in 1:num_rows){
    if(j == num_rows){
      D_mat[j,j] = var(X[,j])
    }
    else{
      D_mat[j,j] = var(rowModelsResiduals[[j]])
    }
  }
  # Create T which is upper triangular matrix with negative coefficients
  rowModelCoefficients <- lapply(rowModelsResults, function(x) coef(x))
  T_mat = diag(1, num_rows, num_rows)
  for(i in 1:num_rows){
    if(i != num_rows){
      T_mat[i,(i+1):num_rows] = -rowModelCoefficients[[i]]
    }
  }
  
  # Multiply data matrix by L* = TD^(-1/2) to get decorrelated matrix.  
  # Also compare to Cholesky L*
  L_star = T_mat %*% inverse_square_root(D_mat)
  return(L_star)
}
##### Doing regression to get beta estimates for just X_5 ############

LS_X5_Regression = function(data_matrix, weights_vec = NULL){ # input is decorrelated matrix
  X = data_matrix
  if(is.null(weights_vec) == F){
    beta_estimates = lm(X[,5]~X[,1] + X[,2] + X[,3] + X[,4] - 1, weights = weights_vec)
  }
  else{
    beta_estimates = lm(X[,5]~X[,1] + X[,2] + X[,3] + X[,4] - 1)
  }
  return(beta_estimates)
}


###################################################################################
# For 1 simulation.  Will generalize for multiple simulations after getting it working
# Get MSE for least sq, weighted least sq, and cholesky least sq.

# mse = function(coefficients){
#   pred = coefficients[1] + coefficients[2]*X[,1] + coefficients[3]*X[,2] + coefficients[4]*X[,3] + coefficients[5]*X[,4]
#   return(mean((pred - X[5])^2))
# }
# 
# # Regular Least Squares
# least_sq = lm(X_norm[,5]~X_norm[,1]+X_norm[,2]+X_norm[,3]+X_norm[,4])
# least_sq_coeff = coef(least_sq)
# mse(least_sq_coeff)
# 
# # Cholesky Least Squares
# 
# U = chol(solve(sig))
# X_chol_decor = U %*% X_norm
# chol_least_sq = lm(X_chol_decor[,5]~X_chol_decor[,1]+X_chol_decor[,2]+X_chol_decor[,3]+X_chol_decor[,4])
# chol_least_sq_coeff = coef(chol_least_sq)
# mse(chol_least_sq_coeff)
# # Weighted Least Squares

# Get weights by successive row regression

obtain_weights = function(data_matrix){
  t_data = t(data_matrix)
  num_rows = nrow(data_matrix)
  weights = rep(NA,num_rows)
  num_cols = ncol(data_matrix)
  VarsY = c()
  VarsX = c()
  for(i in 1:(num_rows)){
    VarsY = c(VarsY, paste("V",i, sep = ""))
  }
  for(j in 2:num_rows){
    VarsX = c(VarsX, paste(VarsY[j:length(VarsY)], collapse = "+"))
  }
  # Since last row is not regressed on anything, residual is itself
  VarsY = VarsY[-length(VarsY)]
  
  rowModelsList <- lapply(paste(paste(VarsY, VarsX, sep = "~"), -1, sep = ""), as.formula)
  rowModelsResults <- lapply(rowModelsList, function(x) summary(lm(x, data = as.data.frame(t_data))))
  for(k in 1:num_rows){
    if(k != num_rows){
      weights[k] = 1 - rowModelsResults[[k]]$adj.r.squared
    }
    else{
      weights[k] = 1
    }
  }
  return(weights)
}
##################################################################################################


# n = 2, p = 100
covariance = matrix(c(1,.8,.8,1), nrow = 2, ncol = 2)
samp_data = t(mvrnorm(n = 100, mu = rep(0,2), Sigma = covariance))
theta = solve(covariance)
#samp_data = column_standardization(samp_data)
print(chol(theta)) # Upper triangular decorrelation matrix L* from Cholesky
print(L_successive_reg_row(samp_data)) # Upper triangular decorrelation matrix L_hat from successive row regression

# n = 4, p = 10
n = 4
p = 10
sig = toeplitz_struc(n)
theta = solve(sig)
# omg.sq = runif(n = p, min = 0.1, max = 2)
omg.sq = rep(1, p)
vers = 1

B = gen.B(p = p)
data = sim_X(vers, p, n, omg.sq, sig, b = B$b)
# data_stand = column_standardization(data$X)
weights = obtain_weights(data$X)
print(chol(theta))
print(L_successive_reg_row(data$X))


# n = 100, p = 200

n = 100
p = 200
sig = const_struc(n)
theta = solve(sig)
# omg.sq = runif(n = p, min = 0.1, max = 2)
omg.sq = rep(1, p)
vers = 1

B = gen.B(p = p)

beta_mat = add_Betas(B$b)

data = sim_X(vers, p, n, omg.sq, sig, b = beta_mat)

#data_stand = column_standardization(data$X)
print(chol(theta)[1:5,1:5])

print(L_successive_reg_row(data$X)[1:5,1:5])


data_Lsucc_decorr = L_successive_reg_row(data$X) %*% data$X
data_L_cholesky = chol(theta) %*% data$X
weights = obtain_weights(data$X)

true_beta = beta_mat[1:4,5]

LS_X5_Regression(data_Lsucc_decorr)
LS_X5_Regression(data_L_cholesky)
LS_X5_Regression(data$X, weights_vec = weights)
LS_X5_Regression(data$X)
true_beta

#################### Multiple Simulation #######################
beta_sim_multiple = function(num_sims = 1, 
                             n = 100, 
                             p = 200, 
                             row_cov_structure = toeplitz_struc(n)){
  omg.sq = rep(1,p)
  mse_beta = list()
  for(list_num in 1:4){
    mse_beta[[list_num]] = data.frame("Cholesky" = rep(0,num_sims),
                                      "Row_Decorrelation" = rep(0,num_sims),
                                      "Weighted" = rep(0,num_sims),
                                      "Naive" = rep(0,num_sims))
  }
  for(i in 1:num_sims){
    B = gen.B(p = p, seed = i)
    beta_mat = add_Betas(B$b)
    # Adds some beta values [-.9,-.3] U [.3,.9] to connect to X_5 between
    # 1 -> 2
    # 1 -> 4
    # 2 -> 3
    # 2 -> 5
    # 3 -> 5
    # 4 -> 5
    data = sim_X(vers = i, p, n, omg.sq, row_cov_structure, b = beta_mat)
    data_Lsucc_decorr = L_successive_reg_row(data$X) %*% data$X # Data matrix after row decorrelation
    data_L_cholesky = chol(theta) %*% data$X # Data matrix after cholesky decorrelation
    weights = obtain_weights(data$X) # Getting weights through weighted LS 1 - R^2 successively
    # X1 ~ X2 + X3 + X4 + ... + Xn
    # X2 ~ X3 + X4 + ... + Xn
    true_beta = beta_mat[1:4,5]
    LS_row = LS_X5_Regression(data_Lsucc_decorr)$coefficients
    LS_chol = LS_X5_Regression(data_L_cholesky)$coefficients
    LS_weighted = LS_X5_Regression(data$X, weights_vec = weights)$coefficients
    LS_naive = LS_X5_Regression(data$X)$coefficients # Naive LS
    
    for(j in 1:4){
      mse_beta[[j]][i,] = c(mse(LS_chol[[j]], true_beta[[j]]),
                            mse(LS_row[[j]], true_beta[[j]]),
                            mse(LS_weighted[[j]], true_beta[[j]]),
                            mse(LS_naive[[j]], true_beta[[j]]))
    }
    if(i %% 10 == 0){
      print(i)
    }
  }
  return(mse_beta)
}
start.time <- Sys.time()
sim_1 = beta_sim_multiple(num_sims = 50)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
sim_const_struc = beta_sim_multiple(num_sims = 100, row_cov_structure = const_struc(n))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

mse_all_betas = function(mse_beta){
  return(list(mean(mse_beta[,1]),mean(mse_beta[,2]),mean(mse_beta[,3]),mean(mse_beta[,4])))
}
