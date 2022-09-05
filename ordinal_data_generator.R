library("mvtnorm")
library("boot")
library("combinat")
library("MASS")
library("HDclassif")
library("clustvarsel")
library("mclust")
#install.packages("fossil")
library("fossil")

source("helper.R")
options("max.print" = 99999)

#---------------------------  Helper Function
inv.pai = function(x) {
  return(1 / (1 + exp(x)))
}

create_ordinal = function(n,
                          j,
                          Z_sim,
                          y_option = 2,
                          use_option,
                          alpha = NULL,
                          alpha2 = NULL) {
  # k is at least 2;
  r = ncol(Z_sim)
  
  alpha2_pseudo = c(-1,-0.3, 0.4, 0.9) + runif(1,-0.2, 0.2)
  
  coeff = c(-1.2,-0.4, 0.4, 1.2)
  
  if (use_option == 1) {
    alpha = c(runif(1, 2, 3), runif(1, 0, 0.5))
    
    alpha2 = rnorm(j - 1, 0, 0.1) + coeff[1:(j - 1)]
  } else{
    alpha = runif(q, 0, 0.01)
    
    alpha2 = rnorm(j - 1, 0, 0.1) + coeff[1:(j - 1)]
  }
  resid = rnorm(n, 0, 0.1)
  
  pai = matrix(0, n, j)
  
  
  pai[, 1] = inv.pai(Z_sim %*% alpha - alpha2[1] + resid)
  
  if (j > 2) {
    for (i in 2:(j - 1)) {
      pai[, i] = inv.pai(Z_sim %*% alpha - alpha2[i] + resid) - inv.pai(Z_sim %*% alpha - alpha2[i -
                                                                                                   1] + resid)
      
    }
  }
  pai[, j] = 1 - inv.pai(Z_sim %*% alpha - alpha2[j - 1] + resid)
  
  if (y_option == 1) {
    Y_sim = apply(pai, 1, which.max)
    
  } else if (y_option == 2) {
    Y_sim = rep(0, n)
    
    for (i in 1:n) {
      Y_sim[i] = sample(c(1:j), size = 1, prob = pai[i,])
      
    }
  } else{
    
  }
  
  fit = NULL
  
  fit$alpha = alpha
  fit$alpha2 = alpha2
  fit$Y_sim = Y_sim
  
  return(fit)
}

create_binary = function(n,
                         j,
                         Z_sim,
                         y_option = 2,
                         use_option,
                         alpha = NULL,
                         alpha2 = NULL) {
  # k is at least 2;
  r = ncol(Z_sim)
  
  alpha2_pseudo = 0.5 + runif(1,-0.2, 0.2)
  
  
  if (use_option == 1) {
    alpha = c(runif(1, 0, 0.5), runif(1, 2, 3))
    alpha2 = rnorm(1, 0, 0.3)
  } else{
    alpha = runif(q, 0, 0.01)
    #runif(q, 0, 0.5);
    alpha2 = rnorm(1, 0, 0.3)
  }
  resid = rnorm(n, 0, 0.1)
  
  pai = matrix(0, n, j)
  
  
  pai[, 1] = inv.pai(Z_sim %*% alpha - alpha2[1] + resid)
  
  if (j > 2) {
    for (i in 2:(j - 1)) {
      pai[, i] = inv.pai(Z_sim %*% alpha - alpha2[i] + resid) - inv.pai(Z_sim %*% alpha - alpha2[i -
                                                                                                   1] + resid)
      
    }
  }
  pai[, j] = 1 - inv.pai(Z_sim %*% alpha - alpha2[j - 1] + resid)
  
  
  if (y_option == 1) {
    Y_sim = apply(pai, 1, which.max)
    
  } else if (y_option == 2) {
    Y_sim = rep(0, n)
    
    for (i in 1:n) {
      Y_sim[i] = sample(c(1:j), size = 1, prob = pai[i,])
      
    }
  } else{
    
  }
  
  fit = NULL
  
  fit$alpha = alpha
  fit$alpha2 = alpha2
  fit$Y_sim = Y_sim
  
  return(fit)
}


#---------------------------- Ordinal modality
set.seed(1)
total_replicate = 5

for (num_replicate in 1:total_replicate) {
  n = 400
  p = 20
  q = 2
  k = 3
  tau1 = 0.3
  tau2 = 0.4
  tau3 = 0.3
  Y_sim_check = array(0, c(n, p))
  true_memb = c(rep(1, n * tau1), rep(2, n * tau2), rep(3, n * tau3))
  
  
  if (num_replicate == 1) {
    u1 = c(-1.59, 0.77)
    sigma1 = matrix(c(0.17, 0.08, 0.08, 0.14),
                    byrow = TRUE,
                    ncol = q)
    u2 = c(1.5, 0.16)
    sigma2 = matrix(c(0.16,-0.08,-0.08, 0.12),
                    byrow = TRUE,
                    ncol = q)
    u3 = c(-0.01,-1.15)
    sigma3 = matrix(c(0.1,-0.01,-0.01, 0.09),
                    byrow = TRUE,
                    ncol = q)
  } else if (num_replicate == 2) {
    u1 = c(0.96,-1.35)
    sigma1 = matrix(c(0.12,-0.07,-0.07, 0.17),
                    byrow = TRUE,
                    ncol = q)
    u2 = c(-1.15, 1.54)
    sigma2 = matrix(c(0.09,-0.06,-0.06, 0.13),
                    byrow = TRUE,
                    ncol = q)
    u3 = c(1.11, 1.37)
    sigma3 = matrix(c(0.15, 0.04, 0.04, 0.11),
                    byrow = TRUE,
                    ncol = q)
  } else if (num_replicate == 3) {
    u1 = c(0.5, 0.5)
    sigma1 = matrix(c(0.27,-0.23,-0.23, 0.29),
                    byrow = TRUE,
                    ncol = q)
    u2 = c(0.5, 0.5)
    sigma2 = matrix(c(0.27,-0.23,-0.23, 0.29),
                    byrow = TRUE,
                    ncol = q)
    u3 = c(0.5, 0.5)
    sigma3 = matrix(c(0.27,-0.23,-0.23, 0.29),
                    byrow = TRUE,
                    ncol = q)
  } else if (num_replicate == 4) {
    u1 = c(-1,-1)
    sigma1 = matrix(c(0.16,-0.08,-0.08, 0.12),
                    byrow = TRUE,
                    ncol = q)
    u2 = c(-1,-1)
    sigma2 = matrix(c(0.16,-0.08,-0.08, 0.12),
                    byrow = TRUE,
                    ncol = q)
    u3 = c(-1,-1)
    sigma3 = matrix(c(0.16,-0.08,-0.08, 0.12),
                    byrow = TRUE,
                    ncol = q)
  }else if (num_replicate == 5) {
    u1 = c(1,1)
    sigma1 = matrix(c(0.2,0.1,0.1, 0.22),
                    byrow = TRUE,
                    ncol = q)
    u2 = c(1,1)
    sigma2 = matrix(c(0.2,0.1,0.1, 0.22),
                    byrow = TRUE,
                    ncol = q)
    u3 = c(1,1)
    sigma3 = matrix(c(0.2,0.1,0.1, 0.22),
                    byrow = TRUE,
                    ncol = q)
  }
  
  
  
  
  Z_sim = as.matrix(rbind(
    rmvnorm(n * tau1, u1, sigma1),
    rmvnorm(n * tau2, u2, sigma2),
    rmvnorm(n * tau3, u3, sigma3)
  ))
  
  
  # jpeg(paste("5", "_exp1_plot.jpg", sep = ""))
  plot(Z_sim, col = as.factor(true_memb))
  # dev.off()
  # save(Z_sim, file = "5_exp1_Z_sim.RData")
  j_list_unit = c(rep(4, 10), rep(5, 5), rep(5, 10), rep(4, 10), rep(4, 5))
  j_list = rep(j_list_unit, 5)
  
  alpha = matrix(nrow = p, ncol = q)
  
  alpha2 = matrix(nrow = p, ncol = max(j_list) - 1)
  
  for (j in 1:(p / 4)) {
    cat("----------------------------------------", "\n")
    if (j_list[j] > 3) {
      fit = create_ordinal(n,
                           j_list[j],
                           Z_sim,
                           y_option = 2,
                           use_option = 1)
      
    } else{
      fit = create_binary(n,
                          j_list[j],
                          Z_sim,
                          y_option = 2,
                          use_option = 1)
      
    }
    
    Y_sim_check[, j] = fit$Y_sim
    print(table(Y_sim_check[, j]))
    print(fit$alpha)
    print(fit$alpha2)
    alpha[j,] = fit$alpha
    
    alpha2[j, 1:(j_list[j] - 1)] = fit$alpha2
    
  }
  for (j in (p / 4 + 1):p) {
    if (j_list[j] > 3) {
      fit = create_ordinal(n,
                           j_list[j],
                           Z_sim,
                           y_option = 2,
                           use_option = 2)
      
    } else{
      fit = create_binary(n,
                          j_list[j],
                          Z_sim,
                          y_option = 2,
                          use_option = 2)
      
    }
    
    Y_sim_check[, j] = fit$Y_sim
    print(table(Y_sim_check[, j]))
    print(fit$alpha)
    print(fit$alpha2)
    alpha[j,] = fit$alpha
    
    alpha2[j, 1:(j_list[j] - 1)] = fit$alpha2
    
  }
  Xc_new_ordinal = Y_sim_check
  
  
  #-----------------------------  Test # Mixed FMA
  
  Xc = list(Xc_new_ordinal)
  
  save(Xc, file = paste(num_replicate,  "_exp6_Xc.RData", sep = ""))
}

#save(alpha, file = "sim_alpha.RData")
#save(alpha2, file = "sim_alpha2.RData")
