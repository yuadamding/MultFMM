ordinal_init = function(x, rc_list, j_list, k) {
  pp = gauher(8)
  Xc[[1]] = x
  rc = max(rc_list)
  n = dim(Xc[[1]])[1]
  
  mc = 1
  
  pc_list = vector(length = mc)
  
  for (temp in 1:mc) {
    pc_list[temp] = dim(Xc[[1]])[2]
  }
  pc = max(pc_list)
  
  
  y = Xc[[1]]
  
  
  # 4/7/2020
  set.seed(4)# (4)
  
  J = max(j_list)
  
  
  esty = y %*% fa(cor(y), rc)$loadings
  
  la = Mclust(esty, k)
  membc = la$classification
  initc = NULL
  if (is.null(initc$wc)) {
    wc = table(membc) / n
    wc = t(t(wc))
  } else
    wc = initc$wc
  W2 <- matrix(wc, nrow = k, ncol = n)
  
  W = array(1 / k, c(k, n))
  
  
  #if (k>1) membc=kmeans(y,k) else membc=rep(1,n)
  #plot(esty[membc == 1, ])
  #p <- ggplot(as.data.frame(esty), aes(esty[, 1], esty[, 2]))
  #p + geom_point(aes(colour = factor(membc)))
  alpha = matrix(nrow = pc, ncol = rc)
  
  alpha2 = matrix(nrow = pc, ncol = J - 1)
  
  for (tempr in 1:rc) {
    esty[, tempr] = scale(esty[, tempr], scale = TRUE, center = TRUE)
    
  }
  for (tempp in 1:pc) {
    if (j_list[tempp] >= 3) {
      tempy = as.factor(y[, tempp])
      
      temp.fit = tryCatch({
        polr(tempy ~ esty,  Hess = T)
      },
      error = function(e) {
        esty1 = y %*% fa(cor(y), rc - 1)$loadings
        
        temp.fit1 = polr(tempy ~ esty1,  Hess = T)
        start = c(0, temp.fit1$coefficients, temp.fit1$zeta)
        temp.fit = polr(tempy ~ esty , start = start,  Hess = T)
      })
      alpha[tempp, ] = temp.fit$coefficients
      if (length(temp.fit$zeta) != j_list[tempp] - 1) {
        alpha2[tempp, 1:j_list[tempp] - 1] = seq(1, j_list[tempp] - 1, 1)
      } else{
        alpha2[tempp, 1:j_list[tempp] - 1] = temp.fit$zeta
      }
      
    } else if (j_list[tempp] == 2) {
      tempy = as.factor(2 - y[, tempp])
      
      temp.fit <- glm(tempy ~ esty, family = "binomial")
      alpha[tempp, ] = -temp.fit$coefficients[-1]
      
      alpha2[tempp, 1:j_list[tempp] - 1] = temp.fit$coefficients[1]
      
    }
  }
  
  muc = matrix(nrow = k, ncol = rc)
  
  for (tempk in 1:k) {
    #muc[tempk, ] = colMeans(esty[membc == tempk, ])
    muc[tempk, ] = colMeans(as.matrix(esty[membc == tempk, ]))
    #for(tempr in 1:rc){
    #  muc[tempk, tempr] = median(scale(esty[membc == tempk, tempr], scale = TRUE, center = TRUE))
    #}
  }
  
  #lc = kmeans(y, k)$cluster
  
  #set.seed(seed)# (5)
  #if (is.null(initc$muc)) {muc=matrix(runif(k*rc,-1,1),k,rc)} else muc=initc$muc
  
  if (is.null(initc$sigmac)) {
    sigmac <- array(0, c(k, rc, rc))
    for (i in 1:k)
      sigmac[i, ,] = 0.5 * diag(rc)
  } else
    sigmac = initc$sigmac
  
  
  nqc = length(pp$x) ^ rc
  #py.s = matrix(0,n,k)
  py.s = matrix(0, k, n)
  
  #ps.y = matrix(0,n,k)
  ps.y = matrix(0, k, n)
  
  if (rc > 1)
    punteggi = sqrt(2) * z_ext(pp$x, rc)
  else
    punteggi = matrix(sqrt(2) * pp$x, nqc)
  if (rc > 1)
    pesi = z_ext(pp$w, rc) / sqrt(pi)
  else
    pesi = matrix(pp$w / sqrt(pi), nqc)
  
  res = list(
    rc = rc,
    mc = mc,
    n = n,
    pc = pc,
    k = k,
    x = x,
    alpha = alpha,
    alpha2 = alpha2,
    nqc = nqc,
    py.s = py.s,
    ps.y = ps.y,
    sigmac = sigmac,
    punteggi = punteggi,
    muc = muc,
    J = J,
    pesi = pesi,
    rc_list = rc_list,
    j_list = j_list,
    W = W
  )
  return(res)
  
}
###############################################################################
##############################################################################################################################################################
ordinal_estep = function(x, init_ordinal) {
  j_list = init_ordinal$j_list
  rc = init_ordinal$rc
  mc = 1
  n = init_ordinal$n
  pc = init_ordinal$pc
  k = init_ordinal$k
  y = x = init_ordinal$x
  alpha = init_ordinal$alpha
  alpha2 = init_ordinal$alpha2
  nqc = init_ordinal$nqc
  py.s = init_ordinal$py.s
  ps.y = init_ordinal$ps.y
  sigmac = init_ordinal$sigmac
  punteggi = init_ordinal$punteggi
  muc = init_ordinal$muc
  J = init_ordinal$J
  pesi = init_ordinal$pesi
  rc_list = init_ordinal$rc_list
  
  
  E.z.sy = array(0, c(n, k, rc))
  E.zz.sy = array(0, c(n, k, rc, rc))
  zex = array(0, c(nqc, rc, k))
  p.z.ys = array(0, c(nqc, n, k))
  CC = array(0, c(rc, rc, n))
  
  for (i in 1:k) {
    Ci = tryCatch({
      chol(nearPD(sigmac[i, , ])$mat)
    },
    error = function(e) {
      chol(nearPD(-sigmac[i, , ])$mat)
    })
    zex[, , i] = as.matrix(punteggi %*% (Ci), nqc, rc) + t(matrix(muc[i,], rc, nqc))
    
    pi.ordinal = array(dim = c(nqc, pc, J))
    
    for (num.z in 1:nqc) {
      for (num.p in 1:pc) {
        pi.ordinal[num.z, num.p, 1] = 1 / (1 + exp(zex[num.z, , i] %*% alpha[num.p, ] - alpha2[num.p, 1]))
        
        if (j_list[num.p] > 2) {
          for (num.j in 2:(j_list[num.p] - 1)) {
            pi.ordinal[num.z, num.p, num.j] = 1 / (1 + exp(zex[num.z, , i] %*% alpha[num.p, ] - alpha2[num.p, num.j])) - 1 /
              (1 + exp(zex[num.z, , i] %*% alpha[num.p, ] - alpha2[num.p, num.j - 1]))
            
          }
        }
        pi.ordinal[num.z, num.p, j_list[num.p]] = 1 - 1 / (1 + exp(zex[num.z, , i] %*% alpha[num.p, ] - alpha2[num.p, j_list[num.p] -
                                                                                                                 1]))
        
      }
    }
    
    py.zex = matrix(nrow = nqc, ncol = n)
    
    for (num.n in 1:n) {
      for (num.z in 1:nqc) {
        temp.p = 1
        
        for (num.p in 1:pc) {
          temp.p = temp.p * prod(pi.ordinal[num.z, num.p, 1:j_list[num.p]] ^ y2yc(y[num.n, num.p], j_list[num.p]))
          
        }
        py.zex[num.z, num.n] = temp.p
        
      }
    }
    py.s[i,] = t(apply(pesi, 1, prod) %*% py.zex)
    
    p.z.ys[, , i] = matrix(apply(pesi, 1, prod), nqc, n) * py.zex / t(matrix(py.s[i, ], n, nqc))
    p.z.ys = ifelse(is.na(p.z.ys) , 0, p.z.ys)
    E.z.sy[, i,] = t(t(zex[, , i]) %*% p.z.ys[, , i])
    
    ### E.zz.sy
    temp = ((matrix(zex[, , i], nrow = nqc))) %o% t(matrix(zex[, , i], nrow =
                                                             nqc))
    temp = apply(temp, c(2, 3), diag)
    temp = aperm(temp, c(2, 3, 1))
    for (l in 1:rc) {
      CC[l, ,] = temp[l, ,] %*% p.z.ys[, , i]
      E.zz.sy[, i, ,] = aperm(CC, c(3, 1, 2))
    }
  }
  
  res = list(
    x = x ,
    init_ordinal = init_ordinal ,
    py.s = py.s,
    py.zex = py.zex ,
    p.z.ys = p.z.ys ,
    E.z.sy = E.z.sy ,
    E.zz.sy = E.zz.sy ,
    zex = zex
  )
  return(res)
}
######################################################################################
ordinal_estep2 = function(init_ordinal, e_ordinal, p.mix.y) {
  j_list = init_ordinal$j_list
  rc = init_ordinal$rc
  mc = 1
  n = init_ordinal$n
  pc = init_ordinal$pc
  k = init_ordinal$k
  y = x = init_ordinal$x
  alpha = init_ordinal$alpha
  alpha2 = init_ordinal$alpha2
  nqc = init_ordinal$nqc
  ps.y = init_ordinal$ps.y
  sigmac = init_ordinal$sigmac
  punteggi = init_ordinal$punteggi
  muc = init_ordinal$muc
  J = init_ordinal$J
  pesi = init_ordinal$pesi
  rc_list = init_ordinal$rc_list
  
  E.z.sy = array(0, c(n, k, rc))
  E.zz.sy = array(0, c(n, k, rc, rc))
  zex = e_ordinal$zex
  p.z.ys = array(0, c(nqc, n, k))
  CC = array(0, c(rc, rc, n))
  
  py.zex = e_ordinal$py.zex
  py.s = p.mix.y$p.mix.y
  for (i in 1:k) {
    p.z.ys[, , i] = matrix(apply(pesi, 1, prod), nqc, n) * py.zex / t(matrix(py.s[i, ], n, nqc))
    p.z.ys = ifelse(is.na(p.z.ys) , 0, p.z.ys)
    E.z.sy[, i,] = t(t(zex[, , i]) %*% p.z.ys[, , i])
    
    ### E.zz.sy
    temp = ((matrix(zex[, , i], nrow = nqc))) %o% t(matrix(zex[, , i], nrow =
                                                             nqc))
    temp = apply(temp, c(2, 3), diag)
    temp = aperm(temp, c(2, 3, 1))
    for (l in 1:rc) {
      CC[l, ,] = temp[l, ,] %*% p.z.ys[, , i]
      E.zz.sy[, i, ,] = aperm(CC, c(3, 1, 2))
    }
  }
  E.z.sy = ifelse(E.z.sy > 1e1 | E.z.sy < 1e-1, 0, E.z.sy)
  E.zz.sy = ifelse(E.zz.sy > 1e1 | E.zz.sy < 1e-1, 0, E.zz.sy)
  res = list(
    x = x ,
    init_ordinal = init_ordinal ,
    py.s = py.s,
    py.zex = py.zex ,
    p.z.ys = p.z.ys ,
    E.z.sy = E.z.sy ,
    E.zz.sy = E.zz.sy ,
    zex = zex
  )
  return(res)
}

######################################################################################
######################################################################################
modality_mixed = function(init_ordinal,
                          py.h_mixed ,
                          py.h_total_mixed ,
                          py.s_mixed ,
                          k ,
                          n ,
                          X) {
  mc = length(init_ordinal)
  W = array(0, c(k, n))
  for (i in 1:mc) {
    W = W + init_ordinal[[i]]$W
  }
  W = W / mc
  
  if (is.null(X)) {
    py.h_total_mixed = array(0, c(k, n))
  }
  p.mix.y = matrix(0, nrow = k, ncol = n)
  
  
  
  for (num_n in 1:n) {
    for (i in 1:k) {
      p.mix.y[i, num_n] <-
        W[i, num_n] * (py.h_mixed[i, num_n] * 10 ^ (py.h_total_mixed[i, num_n])) * py.s_mixed[i, num_n] /
        ((t(W[, num_n])) %*% (py.s_mixed[, num_n] * (
          py.h_mixed[, num_n] * (10 ^ py.h_total_mixed[, num_n])
        )))
      
      #p.mix.y[i,num_n]<-W[i, num_n]* (py.h_mixed[i,num_n]*10^(py.h_total_mixed[i, num_n]))* py.s_mixed[i, num_n] /sum(W[, num_n]* (py.h_mixed[,num_n]*10^(py.h_total_mixed[, num_n]))* py.s_mixed[, num_n]);
    }
    for (i in 1:k) {
      p.mix.y[i, num_n] = ifelse(any(is.na(p.mix.y[, num_n])) |
                                   (sum(p.mix.y[, num_n]) <= 1 - 1e-10) |
                                   (sum(p.mix.y[, num_n]) >= 1 + 1e-10),
                                 rep(1 / k, k),
                                 p.mix.y[i, num_n])
    }
  }
  
  p.mix.y <- ifelse(is.na(p.mix.y), 0.0000000000001, p.mix.y)
  p.mix.y <- ifelse(p.mix.y == 0, 0.0000000000001, p.mix.y)
  
  if (sum(is.na(p.mix.y)) > 0) {
    cat("Warning: k = ", k, "r = ", r, "p.mix.y has NA!", "\n")
    break
    
  }
  
  #Flag_1 = FALSE
  if (sum(sum(is.na(p.mix.y))) == k * n) {
    #Warning = FALSE
    #Flag_1 = TRUE
    cat("Warning: k = ", k, "r = ", r, "p.mix.y == 0!", "\n")
  }
  
  
  res = list(p.mix.y = p.mix.y)
  return(res)
}
######################################################################################
Ez.y_mixed = function(e_ordinal, p.mix.y, n, rc, k) {
  E.z.sy = e_ordinal$E.z.sy
  Ez.y = matrix(0, n, rc)
  
  for (i in 1:k) {
    #ps.y[,i]=wc[i]*py.s[,i]/p.y
    Ez.y = Ez.y + matrix(p.mix.y[i, ], n, rc) * E.z.sy[, i,]
  }
  
  
  e_ordinal$Ez.y = Ez.y
  return(e_ordinal)
}
##############################################################################################################################################################
W_mixed = function(p.mix.y, n, k) {
  W = rowSums(p.mix.y) / sum(p.mix.y)
  W = matrix(W, k, n)
  return(W)
}

##############################################################################################################################################################
ordinal_mstep = function(p.mix.y,
                         e_ordinal,
                         init_ordinal,
                         lambda_cnlm) {
  rc = init_ordinal$rc
  mc = 1
  n = init_ordinal$n
  pc = init_ordinal$pc
  k = init_ordinal$k
  y = x = init_ordinal$x
  W = init_ordinal$W
  
  alpha = init_ordinal$alpha
  alpha2 = init_ordinal$alpha2
  nqc = init_ordinal$nqc
  py.s = e_ordinal$py.s
  ps.y = init_ordinal$ps.y
  sigmac = init_ordinal$sigmac
  punteggi = init_ordinal$punteggi
  muc = init_ordinal$muc
  J = init_ordinal$J
  pesi = init_ordinal$pesi
  rc_list = init_ordinal$rc_list
  j_list = init_ordinal$j_list
  
  p.z.ys = e_ordinal$p.z.ys
  E.z.sy = e_ordinal$E.z.sy
  E.zz.sy = e_ordinal$E.zz.sy
  zex = e_ordinal$zex
  
  p.mix.y2 = t(p.mix.y)
  
  #alpha_init = alpha; alpha2_init = alpha2;
  alpha_iter = alpha
  alpha2_iter = alpha2
  
  alpha_temp = alpha_iter
  alpha2_temp = alpha2_iter
  
  stop_nlm = FALSE
  h_nlm = 0
  
  pc_list = vector(length = mc)
  
  for (temp in 1:mc) {
    pc_list[temp] = dim(x)[2]
  }
  pc = max(pc_list)
  
  diff = array(0, pc)
  last = 1000000000000000000
  valueObjectfun  = array(0, 2000)
  USE_NLM = T
  MAX_ROUND = ifelse(USE_NLM, 1, 30)
  while (!stop_nlm & h_nlm < MAX_ROUND) {
    h_nlm = h_nlm + 1
    
    p.mix.y2 = t(p.mix.y)
    
    for (j in 1:pc) {
      jc = j_list[j]
      
      
      da.max.L2 = function(p) {
        yy = y[, j]
        yy2 = y2yc2(yy, jc)
        ps.y = p.mix.y2
        
        nq = nrow(zex)
        numobs = length(yy)
        temp = matrix(0, numobs)
        if (jc >= 3) {
          AA = matrix(0, nrow = jc - 2, ncol = jc - 1)
          
          for (tempA in 1:nrow(AA)) {
            AA[tempA, tempA:(tempA + 1)] = c(1, -1)
            
          }
          A2 <- cbind(matrix(0, nrow = jc - 2, ncol = rc), AA)
          
          B2 <- matrix(rep(-.Machine$double.eps, jc - 2), ncol = 1)
          
          if (sum(A2 %*% p < B2) != jc - 2) {
            return(10000000000000000)
          }
        }
        
        for (i in 1:k) {
          #temp.zex = cbind(zex[,,i], -zex[,,i])
          temp.zex = zex[, , i]
          
          lik1 = -log(1 + exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc +
                                                                        1])) %*% matrix(yy2[, 1], nrow = 1)
          lik2 = (temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + jc - 1] -
                    log(1 + exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + jc - 1]))) %*% matrix(yy2[, jc], nrow =
                                                                                                     1)
          
          lik3 = matrix(0, nq, numobs)
          
          if (jc >= 3) {
            for (num_j in 2:(jc - 1)) {
              lik3 = lik3 + (temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j - 1] + log(1 -
                                                                                             exp(p[rc + num_j - 1] - p[rc + num_j])) - log(1 + exp(
                                                                                               temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j]
                                                                                             )) - log(1 + exp(
                                                                                               temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j - 1]
                                                                                             ))) %*% matrix(yy2[, num_j], nrow = 1)
            }
          }
          
          log.p.y.z = lik1 + lik2 + lik3
          
          temp = temp + ps.y[, i] * colSums(as.matrix(p.z.ys[, , i]) * log.p.y.z)
        }
        
        
        #L21 norm
        return((-sum(temp) + lambda_cnlm * sqrt(t(p[1:(rc)]) %*% p[1:(rc)])) /
                 1000)
        
      }
      
      da.max.L1 = function(alpha , beta) {
        jc = j_list[j]
        
        yy = y[, j]
        yy2 = y2yc2(yy, jc)
        ps.y = p.mix.y2
        
        nq = nrow(zex)
        numobs = length(yy)
        temp = matrix(0, numobs)
        p = as.vector(c(beta, alpha))
        if (jc >= 3) {
          #margin = 0;
          AA = matrix(0, nrow = jc - 2, ncol = jc - 1)
          
          for (tempA in 1:nrow(AA)) {
            AA[tempA, tempA:(tempA + 1)] = c(1, -1)
            
          }
          A2 <- cbind(matrix(0, nrow = jc - 2, ncol = rc), AA)
          
          B2 <- matrix(rep(-.Machine$double.eps, jc - 2), ncol = 1)
          
          if (sum(A2 %*% p < B2) != jc - 2) {
            return(10000000000000000)
          }
        }
        
        for (i in 1:k) {
          temp.zex = zex[, , i]
          
          lik1 = -log(1 + exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc +
                                                                        1])) %*% matrix(yy2[, 1], nrow = 1)
          lik2 = (temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + jc - 1] -
                    log(1 + exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + jc - 1]))) %*% matrix(yy2[, jc], nrow =
                                                                                                     1)
          
          lik3 = matrix(0, nq, numobs)
          
          if (jc >= 3) {
            for (num_j in 2:(jc - 1)) {
              #lik3 = lik3 + (temp.zex %*% as.matrix(p[(1):(rc)])-p[rc+num_j-1] + log(1-exp(p[rc+num_j-1]-p[rc+num_j])) - log(1+exp(temp.zex %*% as.matrix(p[(1):(rc)])- p[rc+num_j]))-log(1+exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc+num_j-1])))%*% matrix(yy2[, num_j], nrow=1)
              lik3 = lik3 + (temp.zex %*% as.matrix(p[(1):(rc)]) + log(exp(-p[rc +
                                                                                num_j - 1]) - exp(-p[rc + num_j])) - log(1 + exp(
                                                                                  temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j]
                                                                                )) - log(1 + exp(
                                                                                  temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j - 1]
                                                                                ))) %*% matrix(yy2[, num_j], nrow = 1)
            }
          }
          
          log.p.y.z = lik1 + lik2 + lik3
          
          temp = temp + ps.y[, i] * colSums(as.matrix(p.z.ys[, , i]) * log.p.y.z)
        }
        
        #L21 norm
        return((-sum(temp) + lambda_cnlm * sqrt(t(p[1:(rc)]) %*% p[1:(rc)])) /
                 1000)
      }
      ###################################################
      deriv_function = function(zex, y , p) {
        #@para zex: roots given by Gaussian_Hermite approximation. dim: nqc times rc
        #@para y: dim: nqc times 1
        #@p: observed probability of taking certain ordinal level. dim: 1 times n
        temp = array(0, dim(zex))
        for (i in 1:dim(zex)[2]) {
          temp[, i] = zex[, i] * y
        }
        res = array(0, c(dim(zex), dim(p)[2]))
        for (i in 1:dim(zex)[2]) {
          res[, i,] = temp[, i] %*% p
        }
        
        return(res)
      }
      deriv_function2 = function(p.z.ys, log.p.y.z) {
        temp = array(0, dim(log.p.y.z))
        for (i in 1:dim(log.p.y.z)[2]) {
          temp[, i, ] = p.z.ys * log.p.y.z[, i, ]
        }
        return(temp)
      }
      
      deriv_function3 = function(ps.y, temp) {
        res = array(0, dim(temp))
        for (i in 1:dim(temp)[2]) {
          res[, i] = ps.y * temp[, i]
        }
        return(res)
      }
      hessian_function = function(temp.zex, yyy , prob) {
        #@para temp.zex: roots given by Gaussian_Hermite approximation. dim: nqc times rc
        #@para yyy: dim: nqc times 1
        #@prob: observed probability of taking certain ordinal level. dim: 1 times n
        temp = array(0, c(dim(temp.zex), dim(temp.zex)[2]))
        for (j in 1:dim(temp.zex)[1]) {
          temp[j, ,] = temp.zex[j,] %*% t(temp.zex[j,]) * yyy[j]
        }
        res = array(0, c(dim(temp.zex), dim(temp.zex)[2], dim(prob)[2]))
        for (i in 1:dim(prob)[2]) {
          res[, , , i] = temp * prob[, i]
        }
        
        return(res)
      }
      hessian_function2 = function(p.z.ys, log.p.y.z) {
        temp = array(0, dim(log.p.y.z))
        for (i in 1:dim(log.p.y.z)[2]) {
          for (j in 1:dim(log.p.y.z)[3]) {
            temp[, i, j, ] = p.z.ys * log.p.y.z[, i, j, ]
          }
        }
        return(temp)
      }
      hessian_function3 = function(ps.y, temp) {
        res = array(0, dim(temp))
        for (i in 1:dim(temp)[1]) {
          for (j in 1:dim(temp)[2]) {
            res[i, j ,] = ps.y * temp[i, j,]
          }
        }
        return(res)
      }
      ##################################################
      derivative = function(p) {
        #this function returns the first-order derivative of da.max.L3 under a specific feature p given the input para
        yy = y[, j]
        yy2 = y2yc2(yy, jc)
        ps.y = p.mix.y2
        
        nq = nrow(zex)
        numobs = length(yy)
        temp = matrix(0, numobs, rc)
        for (i in 1:k) {
          temp.zex = zex[, , i]
          
          lik1 = deriv_function(temp.zex, (-exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc +
                                                                                          1])) / (1 + exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + 1])), matrix(yy2[, 1], nrow =
                                                                                                                                                                      1))
          lik2 = deriv_function(temp.zex, (1 - exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc +
                                                                                             jc - 1]) / (
                                                                                               1 + exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + jc - 1])
                                                                                             )), matrix(yy2[, jc], nrow = 1))
          lik3 = array(0, c(nq, rc, numobs))
          
          if (jc >= 3) {
            for (num_j in 2:(jc - 1)) {
              lik3 = lik3 + deriv_function(temp.zex,
                                           (
                                             1 - exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j]) / (1 + exp(
                                               temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j]
                                             )) - exp(temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j - 1]) / (1 +
                                                                                                                    exp(
                                                                                                                      temp.zex %*% as.matrix(p[(1):(rc)]) - p[rc + num_j - 1]
                                                                                                                    ))
                                           ),
                                           matrix(yy2[, num_j], nrow = 1))
            }
          }
          
          log.p.y.z = lik1 + lik2 + lik3
          
          temp = temp + deriv_function3(ps.y[, i], t(colSums(
            deriv_function2(as.matrix(p.z.ys[, , i]), log.p.y.z)
          )))
        }
        
        return(-colSums(temp, na.rm = T))
      }
      ##################################################
      
      p.init = c(alpha_temp[j,], alpha2_temp[j, 1:(j_list[j] - 1)])
      
      p.init = ifelse(is.na(p.init), 0, p.init)
      #if(sum(abs(p.init[1:(rc)]))>1/lambda_cnlm){
      # # p.init[1:rc]=array(0,rc)
      #}
      
      # deriv = -derivative(p.init)
      # new_gamma = 100
      # new_beta =  (1 / new_gamma) * (deriv + new_gamma * p.init[1:rc]) *
      #   max((1 - lambda_cnlm / norm(
      #     as.matrix(deriv + new_gamma * p.init[1:rc]), "2"
      #   )), 0)
      # new_alpha = nlm(f = da.max.L1, p = p.init[(rc + 1):(rc + jc - 1)], new_beta)$estimate
      # diff[j] = max(abs(p.init - c(new_beta, new_alpha)))
      
      
      
      par = nlm(f = da.max.L2, p = p.init)$estimate
      #par = c(new_beta, new_alpha);
      
      #after = da.max.L2(par);
      par1 = par[1:rc]
      
      par2 = par[(rc + 1):(rc + jc - 1)]
      
      alpha_temp[j, ] = par1
      alpha2_temp[j, 1:(jc - 1)] = par2
      
      #valueObjectfun[h_nlm] = valueObjectfun[h_nlm]+ after;
    }
    #stop_nlm = TRUE
    
    if (max(diff) < 1e-4 ||
        last <= round(max(diff), 4)) {
      stop_nlm = TRUE
    }
    alpha_iter = alpha_temp
    alpha2_iter = alpha2_temp
    
    #print(alpha_iter);
    print(round(max(diff), 4))
    last = round(max(diff), 4)
  }
  #cat("h_nlm", h_nlm, "\n")
  alpha = alpha_iter
  alpha2 = alpha2_iter
  
  
  for (i in 1:k) {
    muc[i,] = colSums(matrix(p.mix.y[i, ], n, rc) * E.z.sy[, i,]) / sum(p.mix.y[i, ])
    sigmac[i, ,] = apply((array(p.mix.y[i, ], c(n, rc, rc)) * (E.zz.sy[, i, ,] -
                                                                 aperm(
                                                                   array((muc[i,]) %*% t(muc[i,]), c(rc, rc, n)), c(3, 1, 2)
                                                                 ))), c(2, 3), sum) / sum(p.mix.y[i, ])
  }
  
  temp1 = matrix(0, rc, rc)
  temp2 = matrix(0, rc, rc)
  temp3 = matrix(0, rc)
  for (i in 1:k) {
    temp1 = temp1 + W[i, 1] * sigmac[i, ,]
    temp2 = temp2 + W[i, 1] * (muc[i,] %*% t(muc[i,]))
    temp3 = temp3 + W[i, 1] * muc[i,]
  }
  
  ## correzione per l'identificabilit?
  var.z = temp1 + temp2 - temp3 %*% t(temp3)
  A <- (chol(var.z))
  # for (i in 1:k) {
  #   sigmac[i,,]<-t(ginv(A))%*%sigmac[i,,]%*%ginv(A)
  #   muc[i,]<-t(ginv(A))%*%muc[i,]
  # }
  #
  # mu.tot=t(W[, 1])%*%muc
  # muc=muc-t(matrix(mu.tot,rc,k))
  #alpha[,-1]=alpha[,-1]%*%t(A)
  alpha = alpha %*% t(A)
  
  if (rc > 1)
    for (i in 1:(rc - 1))
      alpha[i, (i + 1):rc] = 0
  
  # ----------------------------- Categorical Modality ends
  
  #temp_categorical =  sum(log(p.y));
  temp_categorical =  sum(log(colSums(W * py.s)))
  #print(temp_categorical)
  
  temp = temp_categorical
  
  
  #likelihood2<-c(likelihood2,temp)
  
  
  
  #2020/5/21
  #if (any(abs(alpha>25))) ratio2=eps
  temp_alpha = as.matrix(cbind(alpha, alpha2))
  
  temp_alpha[is.na(temp_alpha)] = 0
  
  if (sum(abs(temp_alpha) > 200) > 0)
    run.categorical = FALSE
  
  lik2 <- temp
  temp_categorical = lik2
  
  
  
  init_ordinal$alpha = alpha
  init_ordinal$alpha2 = alpha2
  init_ordinal$nqc = nqc
  init_ordinal$py.s = py.s
  init_ordinal$ps.y = ps.y
  init_ordinal$sigmac = sigmac
  init_ordinal$punteggi = punteggi
  init_ordinal$muc = muc
  init_ordinal$J = J
  init_ordinal$pesi = pesi
  init_ordinal$rc_list = rc_list
  init_ordinal$j_list = j_list
  
  index = apply(p.mix.y, 2, which.max)
  print(apply(p.mix.y, 2, which.max))
  
  #print(alpha)
  
  res = list(
    init_ordinal = init_ordinal,
    index = index,
    temp_categorical = temp_categorical
  )
  
  return(res)
}

##############################################################################################################################################################
modality_selection = function(p.mix.y_input,
                              init_ordinal,
                              e_ordinal,
                              m_ordinal,
                              mc,
                              m,
                              beta,
                              sigma_input,
                              Ez.hy,
                              Ezz.hy,
                              lambda_cnlm2) {
  #This function realizes modality-selection across categorical modalities
  
  
  q.z <- 1
  eps = 1e-4
  rc = init_ordinal[[1]]$rc
  k = init_ordinal[[1]]$k
  n = init_ordinal[[1]]$n
  rc_list = rep(rc, mc + m)
  W = array(0, c(k, n))
  for (i in 1:mc) {
    W = init_ordinal[[i]]$W
  }
  #####################################################################
  # dim of p.mix.y:  is (k * n)
  # dim of sigma: from m_ordinal:: (k * rc * rc * (mc + m))
  # dim of muc: from m_ordinal:: (k * rc * (mc + m))
  # dim of Ez.sy: from e_ordinal: from (n * k * rc) to (k * rc * n * (mc + m))
  # dim of E.zz.sy: from e_ordinal: from (n * k * rc * rc) to (k * rc * rc * n * (mc + m))
  #####################################################################
  
  index = m_ordinal[[1]]$index
  p.mix.y = p.mix.y_input$p.mix.y
  sigma = array(0, c(k, rc, rc, mc + m))
  muc = array(0, c(k, rc, mc + m))
  Ez.sy =  array(0, c(k, rc, n, mc + m))
  E.zz.sy = array(0, c(k, rc, rc, n, mc + m))
  sigma[, , , 1:m] = sigma_input
  muc[, , 1:m] = beta
  Ez.sy[, , , 1:m] = Ez.hy
  E.zz.sy[, , , , 1:m] = Ezz.hy
  for (i in (m + 1):(m + mc)) {
    sigma[, , , i] = m_ordinal[[i - m]]$init_ordinal$sigmac
    muc[, , i] = m_ordinal[[i - m]]$init_ordinal$muc
    for (j in 1:k) {
      Ez.sy[j, , , i] = t(e_ordinal[[i - m]]$E.z.sy[, j, ])
      for (p in 1:n) {
        E.zz.sy[j, , , p, i] = e_ordinal[[i - m]]$E.zz.sy[p, j, , ]
      }
    }
  }
  
  devi = array(0, c(k, rc, mc + m))
  muc_new = array(0, c(k, rc, mc + m))
  diff = array(0, c(k, rc, mc + m))
  
  objective_fun2 = function(muc, sigma, mc) {
    #@return: the log-likelihood function with l2 norm value given Ez.sy
    objective = 0
    muc = matrix(muc, c(k, rc))
    for (i in 1:n) {
      temp = 0
      for (j in 1:k) {
        temp = temp + ifelse(index[i] == j, 1, 0) * log(W[j, i] * dmvnorm(t(Ez.sy[j, , i, mc]), mean = muc[j, ], sigma =
                                                                            sigma[j,  , ]))
      }
      objective = objective + ifelse(temp == 0, 1e-10, temp)
    }
    return(-objective + lambda_cnlm2 * norm(muc[, ], "2"))
  }
  
  deriv_muc1 = function(muc, sigma, mc) {
    nrow = dim(muc)[1]
    ncol = dim(muc)[2]
    nmc = dim(muc)[3]
    deriv = matrix(0, nrow, ncol)
    # for(i in 1:n){
    #   temp = 0
    #   for(j in 1:k){
    #     temp = temp + W[j, i]* dmvnorm(t(Ez.sy[j, , i,mc ]), mean = muc[j, , mc], sigma =sigma[j,  , , mc ] )
    #   }
    #   temp = ifelse(temp == 0,eps , temp)
    #   for(j in 1:k){
    #     deriv[j, ] = deriv[j, ] +  W[j, i]* dmvnorm(t(Ez.sy[j, , i,mc ]), mean = muc[j, , mc], sigma =sigma[j,  , , mc ] ) * ginv(sigma[j, , , mc])%*%t(t(Ez.sy[j, , i,mc ]) -muc[j, , mc])/temp
    #   }
    # }
    for (i in 1:n) {
      for (j in 1:k) {
        deriv[j, ] = deriv[j, ] + ifelse(index[i] == j, 1, 0) * W[j, i] * ginv(sigma[j, , , mc]) %*% ((Ez.sy[j, , i, mc]) -
                                                                                                        muc[j, , mc])
      }
    }
    return(deriv)
  }
  
  #####################################################################
  diff_iterate = 1000
  val = 1e10
  #sigma_iterate = array(dim = c(k, rc, rc, mc));
  x.z = rep(1, n)
  h_nlm = 0
  stop = FALSE
  USE_NLM = T
  MAX_ROUND = ifelse(USE_NLM, 1, 30)
  while (!stop & h_nlm < MAX_ROUND) {
    h_nlm = h_nlm + 1
    for (j in 1:(mc + m)) {
      devi[, , j] = deriv_muc1(muc, sigma, j)
      muc_new[, , j] = nlm(f = objective_fun2, p = matrix(muc[, , j], c(k * rc, 1)), sigma[, , , j], j)$estimate
      if (j == 3) {
        val_temp = objective_fun2(muc[, , j],  sigma[, , , j], j)
        #print(objective_fun2(muc[, , j],  sigma[, , , j], j))
      }
      gamma = 1000
      #muc_new[, , j] = (1/gamma)*(devi[,,j] +gamma *muc[, ,j])*max((1- lambda_cnlm2/norm(devi[,,j] + gamma *muc[, ,j], "2")),0)
      
      diff[, , j] = abs(muc_new[, , j] - muc[, , j])
      muc[, , j] = muc_new[, , j]
    }
    print(max(abs(diff)))
    if (max(diff) < 1e-4 ||
        diff_iterate <= max(diff) || val_temp > val) {
      stop = TRUE
    }
    diff_iterate = max(abs(diff))
    val = val_temp
  }
  muc_iterate = muc
  sigma_iterate = array(dim = c(k, rc, rc, m + mc))
  
  x.z = rep(1, n)
  for (num_m in 1:(mc + m)) {
    for (i in 1:k) {
      # Give the same sigma to all clusters
      if (i == 1) {
        temp_sigma = apply(aperm(array(
          p.mix.y[i,], c(n, rc_list[num_m], rc_list[num_m])
        ), c(2, 3, 1)) * (
          array(E.zz.sy[i, 1:rc_list[num_m], 1:rc_list[num_m], , num_m], c(rc_list[num_m], rc_list[num_m], n)) -
            array(
              muc_iterate[i, 1:rc_list[num_m], num_m] %*% t(x.z) %*% x.z %*% matrix(t(muc_iterate[i, 1:rc_list[num_m], num_m]), nrow =
                                                                                      q.z) / n,
              c(rc_list[num_m], rc_list[num_m], n)
            )
        ), 1, rowMeans)
      } else{
        temp_sigma = temp_sigma + apply(aperm(array(
          p.mix.y[i,], c(n, rc_list[num_m], rc_list[num_m])
        ), c(2, 3, 1)) * (
          array(E.zz.sy[i, 1:rc_list[num_m], 1:rc_list[num_m], , num_m], c(rc_list[num_m], rc_list[num_m], n)) -
            array(
              muc_iterate[i, 1:rc_list[num_m], num_m] %*% t(x.z) %*% x.z %*% matrix(t(muc_iterate[i, 1:rc_list[num_m], num_m]), nrow =
                                                                                      q.z) / n,
              c(rc_list[num_m], rc_list[num_m], n)
            )
        ), 1, rowMeans)
      }
    }
    temp_sigma = temp_sigma / sum(rowMeans(p.mix.y))
    for (i in 1:k) {
      sigma_iterate[i, 1:rc_list[num_m], 1:rc_list[num_m], num_m] = temp_sigma
    }
  }
  sigma = sigma_iterate
  
  for (j in 1:(mc + m)) {
    temp1 = matrix(0, rc, rc)
    temp2 = matrix(0, rc, rc)
    temp3 = matrix(0, rc)
    for (i in 1:k) {
      temp1 = temp1 + W[i, 1] * sigma[i, , , j]
      temp2 = temp2 + W[i, 1] * (muc[i, , j] %*% t(muc[i, , j]))
      temp3 = temp3 + W[i, 1] * muc[i, , j]
    }
    
    ## correzione per l'identificabilit?
    var.z = temp1 + temp2 - temp3 %*% t(temp3)
    A = tryCatch({
      (chol(var.z))
    },
    error = function(cond) {
      message("var.z is not positive define")
      return(diag(rc))
    })
    #A<-(chol(var.z))
    for (i in 1:k) {
      sigma[i, , , j] <- t(ginv(A)) %*% sigma[i, , , j] %*% ginv(A)
      muc[i , , j] <- t(ginv(A)) %*% muc[i, , j]
    }
    
    mu.tot = t(W[, 1]) %*% muc[, , j]
    muc[, , j] = muc[, , j] - t(matrix(mu.tot, rc, k))
  }
  
  for (i in (m + 1):(m + mc)) {
    m_ordinal[[i - m]]$init_ordinal$sigmac = sigma[, , , i]
    m_ordinal[[i - m]]$init_ordinal$muc = muc[, , i]
  }
  
  beta = muc[, , 1:m]
  sigma_out = sigma[, , , 1:m]
  res = list()
  res$m_ordinal = m_ordinal
  res$beta = beta
  res$sigma = sigma_out
  
  
  return(res)
}
#################################################################################################
############################################################################################
sumZeroList = function(lst) {
  len = length(lst)
  res = 0
  for (i in 1:len) {
    res = res + sum(lst[[i]] != 0)
  }
  return(res)
  
  
}
