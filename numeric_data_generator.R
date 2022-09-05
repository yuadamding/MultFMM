# setwd(
#    "C:/Users/pekph/Dropbox/MyWorkNew/2016_12_26_migraine_data_covariate_Modality_Selection_new_idea_2018_9_4_IISE/sim_6_psu_2018_10_20_2019_revision_random_runs_generate"
# )
library("MASS")
library("mvtnorm")

#--------------------------------------------- Simulation
set.seed(2)

m = 5
n = 400
p = 20
k = 3
q = 3
r = 2

r_list = rep(r, m)

H = array(dim = c(p, r, m))

B = array(dim = c(p, q, m))

beta = array(dim = c(k, r, m))

sigma = array(dim = c(r, r, m))

X = array(dim = c(n, p, m))

H_new = array(dim = c(p, r, m))

beta_new = array(dim = c(k, r, m))

sigma_new = array(dim = c(r, r, m))

x.z = rep(1, n)
x.w = rep(1, n)
w = c(0.3, 0.4, 0.3)
true_memb = c(rep(1, n * w[1]), rep(2, n * w[2]), rep(3, n * w[3]))
n1 = n * w[1]
n2 = n * w[2]
n3 = n * w[3]
num_m = 1

beta[1, 1:r_list[num_m], num_m] = c(0.2, 1.8)
beta[2, 1:r_list[num_m], num_m] = c(-1.5,-1.5)
beta[3, 1:r_list[num_m], num_m] = c(1.1,-1.5)
num_m = 2

beta[1, 1:r_list[num_m], num_m] = c(1.55, -1.75)
beta[2, 1:r_list[num_m], num_m] = c(-1.15,1.75)
beta[3, 1:r_list[num_m], num_m] = c(1.75,1.75)
num_m = 3

beta[1, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
beta[2, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
beta[3, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
num_m = 4

beta[1, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
beta[2, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
beta[3, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
num_m = 5

beta[1, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
beta[2, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
beta[3, 1:r_list[num_m], num_m] = rep(0, r_list[num_m])
num_m = 1
sigma[1:r_list[num_m], 1:r_list[num_m], num_m] = matrix(c(0.25, -0.21, -0.21, 0.37), byrow = TRUE, nrow = r_list[num_m])
num_m = 2
sigma[1:r_list[num_m], 1:r_list[num_m], num_m] = matrix(c(0.23, 0.15, 0.15, 0.3), byrow = TRUE, nrow = r_list[num_m])
num_m = 3
sigma[1:r_list[num_m], 1:r_list[num_m], num_m] = diag(r_list[num_m])
num_m = 4
sigma[1:r_list[num_m], 1:r_list[num_m], num_m] = diag(r_list[num_m])
num_m = 5
sigma[1:r_list[num_m], 1:r_list[num_m], num_m] = diag(r_list[num_m])
for (num_m in 1:m) {
   for (num_r in 1:r_list[num_m]) {
      #H[,num_r,num_m] = c(runif(p/4, min = 1, max = 2), rep(0, 3*p/4))
      H[, num_r, num_m] = c(runif(p / 4, min = 1, max = 1.5) * ifelse(rnorm(p /
                                                                               4) >= -0.5, 1,-1),
                            rep(0, 3 * p / 4))
      #H[,num_r,num_m] = c(runif(p/3, min = 1, max = 2) , rep(0, 2*p/3))
   }
}

for (num_m in 1:m) {
   for (num_q in 1:q) {
      #B[,num_q,num_m] = c(runif(p, min = 1, max = 2))
      #B[,num_q,num_m] = runif(p, min = 1, max = 1.5) * ifelse(rnorm(p/4)>=-0.5, 1, -1)
      B[, num_q, num_m] = c(runif(p / 4, min = 1, max = 1.5) * ifelse(rnorm(p /
                                                                               4) >= -0.5, 1,-1),
                            rep(0, 3 * p / 4))
   }
}
Z
#B = array(runif(p*q*m, min = -1, max = 1), dim = c(p, q, m));
Z = cbind(rep(1, n),
          scale(rnorm(n), center = TRUE),
          scale(rnorm(n), center = TRUE))

beta_new = beta
sigma_new = sigma
H_new = H

for (num_m in 1:m) {
   zex = as.matrix(rbind(
      rmvnorm(n1, beta_new[1, 1:r_list[num_m], num_m], sigma_new[1:r_list[num_m], 1:r_list[num_m], num_m]),
      rmvnorm(n2, beta_new[2, 1:r_list[num_m], num_m], sigma_new[1:r_list[num_m], 1:r_list[num_m], num_m]),
      rmvnorm(n3, beta_new[3, 1:r_list[num_m], num_m], sigma_new[1:r_list[num_m], 1:r_list[num_m], num_m])
   ))
   plot(zex, col = as.factor(true_memb))
   X[, , num_m] = Z %*% t(B[, , num_m]) + zex %*% t(H_new[, 1:r_list[num_m], num_m]) + matrix(rnorm(n * p, sd = 0.1), nrow = n, ncol = p)
}

save(X, file = "X_6.RData")
save(Z, file = "Z.RData")


