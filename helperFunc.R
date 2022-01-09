get_adjmat_from_fges <- function(edgelist, p, varnames){
  # return adjmatrix of cpdag
  myadjmat <- matrix(0, p, p)
  dimnames(myadjmat) <- list(varnames, varnames)
  for (i in 1:length(edgelist)) {
    par_name <- word(edgelist[i], 1)
    chil_name <- word(edgelist[i], -1)
    par_ind <- which(varnames == par_name)
    chil_ind <- which(varnames == chil_name)
    myadjmat[par_ind, chil_ind] <- 1
    if (grepl("<->", edgelist[i])  || grepl("---", edgelist[i]) ) {
      myadjmat[chil_ind, par_ind] <- 1
    }
  }
  return(myadjmat)
}

column_standardization = function(data_matrix){
  num_cols = ncol(data_matrix)
  column_means = colMeans(data_matrix)
  
  for(i in 1:num_cols){
    col_stdev = sd(data_matrix[,i])
    data_matrix[,i] = (data_matrix[,i] - column_means[i])/col_stdev
  }
  return(data_matrix)
}

toeplitz_struc = function(n){
  sig = matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      sig[i,j] = 0.3^(abs(i - j)/5)
    }
  }
  return(sig)
}

const_struc = function(n){
  sig = matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      if(i == j){
        sig[i,j] = 1
      }
      else{
        sig[i,j] = .7
      }
    }
  }
  return(sig)
}

star_struc = function(n){
  theta = matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      if(i == j){
        theta[i,j] = 1
      }
      else if(i == 1 | j == 1){
        a = runif(1, min = 0.1, max = 0.8)
        theta[i,j] = a
        theta[j,i] = a
      }
    }
  }
  return(solve(theta))
}

draw_beta_val = function(){
  potential_beta_val = runif(1, -.9, .9)
  if(abs(potential_beta_val) < 0.3){
    draw_beta_val()
  }
  else{
    return(potential_beta_val)
  }
}

add_Betas = function(beta_matrix){ # Adding betas to X_5
  beta_matrix[1,2] = draw_beta_val()
  beta_matrix[1,3] = 0
  beta_matrix[1,4] = draw_beta_val()
  beta_matrix[1,5] = 0
  beta_matrix[2,3] = draw_beta_val()
  beta_matrix[2,4] = 0
  beta_matrix[2,5] = draw_beta_val()
  beta_matrix[3,4] = 0
  beta_matrix[3,5] = draw_beta_val()
  beta_matrix[4,5] = draw_beta_val()
  return(beta_matrix)
}

inverse_square_root = function(S){
  ei <- eigen(S)
  V <- ei$vectors
  res <- V %*% diag(1 / sqrt(ei$values)) %*% t(V)
  return(res)
}

mse = function(true_beta, est_beta){
  return((true_beta - est_beta)^2)
}

mse_all_betas = function(mse_beta){
    cholesky = c()
    row = c()
    weighted = c()
    naive = c()
    for(i in 1:4){
      cholesky = c(cholesky, mean(mse_beta[[i]]$Cholesky))
      row = c(row, mean(mse_beta[[i]]$Row_Decorrelation))
      weighted = c(weighted, mean(mse_beta[[i]]$Weighted))
      naive = c(naive, mean(mse_beta[[i]]$Naive))
    }
    return(list(cholesky, row, weighted, naive))
}
