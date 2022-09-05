mixfmm = function(X=NULL , Xc=NULL, Z=NULL , k = k, r_list = NULL, rc_list = NULL, j_list = NULL, lambda_cnlm = NULL, lambda_beta =NULL, lambda_gamma=NULL, true_memb = NULL, init.option = 2, penalty= TRUE, parallel = FALSE, x.z = NULL, x.w = NULL, lambda_cnlm2=NULL,
                  it = 15, it2 = 100, eps_nlm = 1e-3, eps_lasso = 1e-6, eps_margin = 1e-6, eps = 1e-03, seed = 4, scaling = FALSE, initc = NULL, pp = gauher(8), maxstep = 100){
  library("MASS")
  library("mvtnorm")
  library("Matrix")
  library("gglasso")
  library("psych")
  library("rlist")
  library("mclust")  
  library("matlib")
  library("combinat")
  library('foreach')
  library('doParallel')
  library("NlcOptim")
  library("quadprog")
  library("GPArotation")
  library("fossil")
  
  source("functions.R")
  source("fma.init_3.R")
  source("helper.R")
  n=400
  tau1 = 0.3; tau2 = 0.4; tau3 = 0.3;
  true_memb = c(rep(1, n*tau1), rep(2, n*tau2), rep(3, n*tau3))
  
  k=3
  r_list = rep(2, 2)
  rc_list = rep(2, 5)
  rc=2
  X = Z =NULL
  #lambda_cnlm = NULL
  j_list = NULL
  init.option = 2
  parallel = FALSE
  x.z = NULL
  #lambda_cnlm2=NULL
  it = 15
  it2 = 100
  eps_nlm = 1e-3
  eps_lasso = 1e-6
  eps_margin = 1e-6
  eps = 1e-04
  seed = 4
  scaling = FALSE
  initc = NULL
  pp = gauher(8)
  maxstep = 100
  # # #
  # 
  # 
  # # load("4_Xc.Rdata")
  # # load("4_X.Rdata")
  # # load("4_Z.Rdata")
  # 
  #load("1_Xc.RData")
  #Xc= c(Xc, Xc)
  #lambda_cnlm = 10
  #lambda_cnlm2 = 20
  j_list = c(rep(2, 4), rep(5, 4), rep(2, 6), rep(5, 6), rep(2, 4), rep(5, 4), rep(2, 6), rep(5, 6),rep(2, 4), rep(5, 4), rep(2, 6), rep(5, 6), rep(2, 4), rep(5, 4), rep(2, 6), rep(5, 6),rep(2, 4), rep(5, 4), rep(2, 6), rep(5, 6));

  rc_list = 2
  rc_list=rep(2, 10)
  k=3
  eps_nlm = 1e-3;
  scaling = FALSE
  eps = 1e-04
  pp = gauher(8)
  it = 15
  
  if(is.null(X)&is.null(Xc)){
    return(cat("ERROR, There is no valid input!!!"))
  }
  
  if(is.null(Xc)){
    mc=0
  }else{
    n = dim(Xc[[1]])[1]
    mc=length(Xc)
  }
  
  if(!is.null(X)){
    n = dim(X)[1]
    m = dim(X)[3]
  }

  
  ptm = proc.time()
  #cat(seed, "\n")
  #set.seed(seed)
  #if(scaling)
    
  r = max(r_list)
  
  
  ###########################################
  ##############for mixed modalities#########
  W = array(0, c(k, n))
  ph.y<-array(0,c(k,n)) 
    
  if(!is.null(X)){
    p = array(0,m)
    q = ncol(Z)
    x.z = rep(1, n)
    x.w = rep(1, n)
    q.z<-1 
    q.w<-1
    if(!is.null(Z)){
      q = ncol(Z)
    }
    for(i in 1:m){
      p[i] = sum(colSums(X[,,i] != 0)!=0)
    }
    
    p_max = max(p)
    
    
    H = array(dim = c(p_max, r, m));
    psi = array(dim = c(p_max,p_max, m));
    beta = array(dim = c(k, r, m));
    sigma = array(dim = c(k, r, r, m));
    if(!is.null(Z)){
      q = ncol(Z);
      B = array(dim = c(p_max, q, m));
      
    }
    
    for(num_m in 1:m){
      init_m = fma.init(X[,1:p[num_m] , num_m], Z, k, r_list[num_m],x.z=NULL,x.w=NULL,seed=4,scaling=scaling);
      H[1:p[num_m], 1:r_list[num_m], num_m] = init_m$H;
      psi[1:p[num_m], 1:p[num_m], num_m] = init_m$psi;
      beta[, 1:r_list[num_m], num_m] = init_m$Beta;
      sigma[, 1:r_list[num_m], 1:r_list[num_m], num_m] = init_m$sigma;	
      if(!is.null(Z)){
        B[1:p[num_m],, num_m] = init_m$B;
      }
      #print(num_m)
    }	
    W = init_m$w;
    
    # memb = init_m$memb #### Membership is the same across all modalities!
    
    # 12/29/2016 6:16 pm
    # Get the adaptive 
    
    group_beta = NULL;
    for(num_m in 1:m){
      group_beta = c(group_beta, rep(num_m, k * r_list[num_m]))
    }
    
    pf_beta = vector(length = m)
    for(num_m in 1:m){
      pf_beta[num_m] = 1
      #pf_beta[num_m] = 1/sqrt(t(as.vector(beta[, 1:r_list[num_m], num_m])) %*% as.vector(beta[, 1:r_list[num_m], num_m]))
    }
    #pf_beta = 1/pf_beta
    
    group_gamma = NULL;
    label = 1
    for(num_m in 1:m){
      for(num_p in 1:p[num_m]){
        group_gamma = c(group_gamma, rep(label, r_list[num_m]))
        #label = label + 1
        group_gamma = c(group_gamma, rep(label, q))
        label = label + 1
        
      }
    }
    pf_gamma = NULL
    for(num_m in 1:m){
      for(num_p in 1:p[num_m]){
        pf_gamma = c(pf_gamma, 1/sqrt(t(c(H[num_p, 1:r_list[num_m], num_m], B[num_p,,num_m])) %*% c(H[num_p, 1:r_list[num_m], num_m], B[num_p,,num_m])))
        #pf_gamma = c(pf_gamma, 1/sqrt(t(c(H[num_p, 1:r_list[num_m], num_m])) %*% c(H[num_p, 1:r_list[num_m], num_m])))
        #pf_gamma = c(pf_gamma, 1)
        #pf_gamma = c(pf_gamma, 0)
      }
    }
    #pf_gamma = 1/pf_gamma
    
    for(num_m in 1:m){
      ybar = apply(X[,,num_m], 2, mean)
      X[,,num_m] = scale(X[,,num_m], ybar,scale = scaling)
    }
    
    
    # EM algorithm
    likelihood = NULL;
    hh = 0;# Number of iteration
    ratio = 1000;
    lik = -10000;
    
    chsi<-array(0,c(k,r,r,m)) 
    roy<-array(0,c(k,r,n,m)) 
    ph.y<-array(0,c(k,n)) 
    py.h<-array(0,c(k,n))
    
    # Helper Function
    L21_norm = function(label, x){
      max_lab = max(label)
      L21 = 0
      for(num_lab in 1:max_lab){
        L21 = L21 + sqrt(t(x[label == num_lab]) %*% x[label == num_lab])
      }
      return(L21)
    }   
    Warning = FALSE
  }
  
  
  
  
  likelihood1 = likelihood2 = likelihood=NULL;
  hh = 0;# Number of iteration
  ratio = 1000;
  ratio1 = ratio2 = ratio;
  #lik = -10000;
  lik = -1000000000;
  lik1 = lik2 = lik;
  Warning = FALSE;
  init=T
  while((hh < it)&(ratio > eps)){
    hh = hh + 1
    temp_total=temp_categorical=temp_numeric=0
    #cat(hh, " round starts----------------------------------------------------------------", "\n")
    
    if(init==T){
      init_ordinal=list()
      if(mc!=0){
        #print("INITIALIZATION FOR ORDINAL")
        for(i in 1:mc){
          init_ordinal=list.append(init_ordinal,ordinal_init(x=Xc[[i]],rc_list,j_list,k))
        }
      }
      init=F
    }
    #**************************************************************************************************************************
    #E step
    
    if(!is.null(X)){
      Ezz.hy = array(0, c(k, r, r, n, m));
      Ez.hy = array(0, c(k,r,n,m))
      sigma.tot = array(0, c(k, p_max, p_max, m));
      change = rep(FALSE, k)####OIS
      
      for(num_m in 1:m){
        if((k == 1)&(r_list[num_m] == 1)){
          sigma.tot[1,1:p[num_m],1:p[num_m],num_m] = as.matrix(H[1:p[num_m],1:r_list[num_m],num_m]) %*% sigma[,1:r_list[num_m],1,num_m] %*% t(as.matrix(H[1:p[num_m],1:r_list[num_m],num_m])) + psi[1:p[num_m],1:p[num_m],num_m];
          if(sum(diag(sigma.tot[1,1:p[num_m],1:p[num_m],num_m])) < 0){change = TRUE;}#####
          if(det(sigma.tot[1,1:p[num_m],1:p[num_m],num_m]) < 0){change = TRUE;}########
        }else{
          for(i in 1:k){
            sigma.tot[i,1:p[num_m],1:p[num_m],num_m] = as.matrix(H[1:p[num_m],1:r_list[num_m],num_m]) %*% sigma[i,1:r_list[num_m],1:r_list[num_m],num_m] %*% t(as.matrix(H[1:p[num_m],1:r_list[num_m],num_m])) + psi[1:p[num_m],1:p[num_m],num_m];	
            if(sum(diag(sigma.tot[i,1:p[num_m],1:p[num_m],num_m])) < 0){change[i] = TRUE;}########
            if(det(sigma.tot[i,1:p[num_m],1:p[num_m],num_m]) < 0){change[i] = TRUE;}#############
          }
        }   
        
        if((k == 1) & (r_list[num_m] == 1)){
          chsi[,1:r_list[num_m],1:r_list[num_m],num_m] = 1/(t(H[1:p[num_m],1:r_list[num_m],num_m]) %*% ginv(psi[1:p[num_m],1:p[num_m],num_m]) %*% H[1:p[num_m],1:r_list[num_m],num_m] + solve(sigma[,1:r_list[num_m],1:r_list[num_m],num_m]))
          if(!is.null(Z)){
            roy[,1:r_list[num_m],,num_m] = chsi[,1:r_list[num_m],1:r_list[num_m],num_m] %*% (t(as.matrix(H[1:p[num_m],1:r_list[num_m],num_m])) %*% ginv(psi[1:p[num_m],1:p[num_m],num_m]) %*% t(X[,1:p[num_m],num_m] - Z %*% t(B[1:p[num_m],,num_m])) + solve(sigma[,1:r_list[num_m],1,num_m]) %*% beta[1,1:r_list[num_m],num_m] %*% t(x.z))
          }else{
            roy[,1:r_list[num_m],,num_m] = chsi[,1:r_list[num_m],1:r_list[num_m],num_m] %*% (t(as.matrix(H[1:p[num_m],1:r_list[num_m],num_m])) %*% ginv(psi[1:p[num_m],1:p[num_m],num_m]) %*% t(X[,1:p[num_m],num_m]) + solve(sigma[,1:r_list[num_m],1,num_m]) %*% beta[1,1:r_list[num_m],num_m] %*% t(x.z))
          }
          Ezz.hy[1,1:r_list[num_m],1:r_list[num_m],,num_m] = matrix(chsi[,1:r_list[num_m],1:r_list[num_m],num_m], ncol = n) + roy[1,1:r_list[num_m],,num_m]^2
          if(change){sigma.tot[1,1:p[num_m],1:p[num_m],num_m] = var(X[,1:p[num_m],num_m])}
        }else{
          for(i in 1:k){
            chsi[i,1:r_list[num_m],1:r_list[num_m],num_m]<-ginv(t(H[1:p[num_m],1:r_list[num_m],num_m])%*%ginv(psi[1:p[num_m],1:p[num_m],num_m])%*%H[1:p[num_m],1:r_list[num_m],num_m]+ginv(sigma[i,1:r_list[num_m],1:r_list[num_m],num_m]))
            if(!is.null(Z)){
              roy[i,1:r_list[num_m],,num_m]<-chsi[i,1:r_list[num_m],1:r_list[num_m],num_m]%*%(t(H[1:p[num_m],1:r_list[num_m],num_m])%*%ginv(psi[1:p[num_m],1:p[num_m],num_m])%*%t(X[,1:p[num_m],num_m] - Z %*% t(B[1:p[num_m],,num_m]))+ginv(sigma[i,1:r_list[num_m],1:r_list[num_m],num_m])%*%beta[i,1:r_list[num_m],num_m]%*%t(x.z))
            }else{
              roy[i,1:r_list[num_m],,num_m]<-chsi[i,1:r_list[num_m],1:r_list[num_m],num_m]%*%(t(H[1:p[num_m],1:r_list[num_m],num_m])%*%ginv(psi[1:p[num_m],1:p[num_m],num_m])%*%t(X[,1:p[num_m],num_m])+ginv(sigma[i,1:r_list[num_m],1:r_list[num_m],num_m])%*%beta[i,1:r_list[num_m],num_m]%*%t(x.z))
            }
            if (r_list[num_m]>1){
              temp<-(t(roy[i,1:r_list[num_m],,num_m]))%o%(roy[i,1:r_list[num_m],,num_m])
              temp<-apply(temp,c(2,3),diag)#####
              temp<-aperm(temp,c(2,3,1))##########
              temp2<-array(chsi[i,1:r_list[num_m],1:r_list[num_m],num_m],c(r_list[num_m],r_list[num_m],n))###########
              Ezz.hy[i,1:r_list[num_m],1:r_list[num_m],,num_m]<-temp+temp2
            }else{
              Ezz.hy[i,r_list[num_m],r_list[num_m],,num_m]<-roy[i,1:r_list[num_m],,num_m]^2+rep(chsi[i,1:r_list[num_m],1:r_list[num_m],num_m],n)
            }  
            if (change[i]){sigma.tot[i,1:p[num_m],1:p[num_m],num_m]<-var(X[,1:p[num_m],num_m])}#########OIS
            #py.h[i,,num_m]<-(dmvnorm(X[,,num_m]-t(H[,,num_m]%*%beta[i,,num_m]%*%t(x.z)),sigma=sigma.tot[i,,,num_m]))
            
          }
        }	
        Ez.hy[,,,num_m] = roy[,,,num_m]   
      }
      
      # Forget to compute py.h for multi-modality when k = r = 1 
      # Date 8/2/2016
      py.h_pro<-array(0,c(k,n))
      for(num_n in 1:n){
        for(i in 1:k){
          temp_p_1 = 1;	
          for(num_m in 1:m){
            #cat(i, num_m, "\n")
            resid_m = X[num_n,1:p[num_m],num_m]- t(as.matrix(H[1:p[num_m],1:r_list[num_m],num_m])%*%beta[i,1:r_list[num_m],num_m]) - Z[num_n, ] %*% t(B[1:p[num_m],,num_m])
            sigma_m = sigma.tot[i,1:p[num_m],1:p[num_m],num_m]
            #dim(resid_m)
            #dim(solve(sigma_m))		
            #(2*pi)^(-p/2) * det(sigma_m)^(-1/2)*exp(-1/2 * resid_m %*% as.matrix(solve(sigma_m)) %*% t(resid_m))
            temp_p_2 = dmvnorm(resid_m,mean = rep(0, p[num_m]),sigma=as.matrix(sigma_m))
            #print(resid_m)
            #print(sigma_m)
            #cat(temp_p_2, "\n")
            temp_digit = as.numeric(strsplit(as.character(temp_p_2), "e")[[1]][2])
            temp_digit = ifelse(is.na(temp_digit), 0, temp_digit)
            #cat(temp_digit, "\n")
            py.h_pro[i, num_n] = py.h_pro[i, num_n] + temp_digit
            #cat(i, num_m, temp_p_2 *10 ^ (-temp_digit), "\n", "\n")
            temp_p_1 = temp_p_1 * temp_p_2 *10 ^ (-temp_digit)
            
            
            
          }
          py.h[i,num_n]<- temp_p_1
          #ifelse(is.na(py.h[i,]),mean(py.h[i,],na.rm=TRUE),py.h[i,])
          #if(sum(is.na(py.h[i,]))==n){py.h[i,]<-0;}	
        }
      }
      
      # E-step(cont')
      Flag_1 = FALSE
      if(sum(sum(py.h == 0)) == k*n){
        Warning = FALSE
        Flag_1 = TRUE
        cat("Warning: k = ", k, "r = ", r, "py.h == 0!", "\n")
      }		
      
      py.h_total = matrix(nrow = k, ncol = n)
      for(num_k in 1:k){
        py.h_total[num_k, ] = py.h_pro[num_k, ] - py.h_pro[k, ]
      }
      #py.h_total = as.matrix(rbind(py.h_total, rep(0, n)))
      
      rbind(py.h_pro, py.h_total)	
    }
    
    if(!is.null(Xc)){
      #print("e step for ordinal")
      e_ordinal=list()
      for(i in 1:mc){
        e_ordinal=list.append(e_ordinal,ordinal_estep(x=Xc[[i]],init_ordinal=init_ordinal[[i]]))
      }
    }
    
    #**************************************************************************************************************************
    ##############for mixed modalities#########
    py.h_mixed=py.s_mixed=array(1,c(k,n))
    py.h_total_mixed=0
    if(!is.null(X)){
      py.h_mixed=py.h
      py.h_total_mixed=py.h_total
    }
    
    if(!is.null(Xc)){
      for(i in 1:mc){
        py.s_mixed=py.s_mixed*e_ordinal[[i]]$py.s
      }
    }
    p.mix.y=modality_mixed(init_ordinal, py.h_mixed, py.h_total_mixed, py.s_mixed, k, n,X)
    if(!is.null(Xc)){
      for(i in 1:mc){
        e_ordinal[[i]]=Ez.y_mixed(e_ordinal[[i]], p.mix.y=p.mix.y$p.mix.y, n,rc=rc_list, k)
      }
    }
    for(i in 1:mc){
      init_ordinal[[i]]$W = W_mixed(p.mix.y=p.mix.y$p.mix.y,n,k)
    }
    
    #**************************************************************************************************************************
    
    # M-step
    if(!is.null(X)){
      #print("m step for numeric")
      ph.y=p.mix.y$p.mix.y
      Ez.y = array(0, c(r, n, m))
      Ezz.y = array(0, c(r,r, n, m))
      for(num_m in 1:m){
        for (i in 1:k) {
          Ez.y[1:r_list[num_m],,num_m]<-Ez.y[1:r_list[num_m],,num_m]+t(matrix(ph.y[i,],n,r_list[num_m]))*Ez.hy[i,1:r_list[num_m],,num_m]
          Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m]<-Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m]+aperm(array(ph.y[i,],c(n,r_list[num_m],r_list[num_m])),c(2,3,1))*Ezz.hy[i,1:r_list[num_m],1:r_list[num_m],,num_m]
        }
        # 2/11/2017 11:16PM
        # exclude the case when r_list[num_m] = 1
        if(r_list[num_m] != 1){
          Ez.y[1:r_list[num_m],,num_m]<-ifelse(is.na(Ez.y[1:r_list[num_m],,num_m]),rowMeans(Ez.y[1:r_list[num_m],,num_m],na.rm=TRUE),Ez.y[1:r_list[num_m],,num_m])      
          Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m]<-ifelse(is.na(Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m]),rowMeans(Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m],na.rm=TRUE),Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m])    
        }
      }
      
      #label_gamma_B
      for(num_m in 1:m){
        # Lambda/Psi 
        EEzz.y = rowMeans(array(Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m], c(r_list[num_m], r_list[num_m], n)), na.rm = TRUE, dims = 2)
        
        B = as.vector(t(as.matrix(cbind(t(X[,1:p[num_m],num_m]) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m], n))), t(X[,1:p[num_m],num_m]) %*% t(array(t(Z), c(q, n))))))) 
        A = as.matrix(cbind(rbind(n*EEzz.y, t(Z) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m], n)))), rbind(t(t(Z) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m], n)))), t(Z) %*% Z)))
        D_temp = chol(A)
        
        y_gamma_1 = sqrt(tr(solve(psi[1:p[num_m],1:p[num_m],num_m]))) * kronecker(diag(p[num_m]), solve(t(D_temp))) %*% (B)
        X_gamma_1 = sqrt(tr(solve(psi[1:p[num_m],1:p[num_m],num_m]))) * kronecker(diag(p[num_m]), D_temp)
        
        if(num_m == 1){
          y_gamma = y_gamma_1
          X_gamma = X_gamma_1	
        }else{
          y_gamma = c(y_gamma, y_gamma_1)
          X_gamma = bdiag(X_gamma, X_gamma_1)
        }
      }
      #LSE
      #solve(t(X_gamma) %*% X_gamma)%*% (t(X_gamma) %*% y_gamma)
      
      X_gamma = as.matrix(X_gamma)		
      
      
      fit_H = gglasso(X_gamma, y_gamma, group = group_gamma, loss = "ls", lambda = lambda_gamma/nrow(X_gamma),pf = pf_gamma, intercept = FALSE)$beta	
      
      
      H_iterate = array(dim = c(p_max, r, m));
      B_iterate = array(dim = c(p_max, q, m))
      start = 1
      for(num_m in 1:m){
        for(num_p in 1:p[num_m]){
          end = start + r_list[num_m] - 1
          #cat(start, end, "\n")
          H_iterate[num_p, 1:r_list[num_m], num_m] = as.vector(fit_H[start:end]) 
          B_iterate[num_p, 1:q, num_m] = as.vector(fit_H[(end + 1):(end + q)]) 
          start = end + q + 1
        }
      }
      H = H_iterate
      B = B_iterate
      
      #t(H[,1:r_list[num_m],num_m] %*% Ez.y[1:r_list[num_m],,num_m])
      
      for(num_m in 1:m){
        # 12/28/2016 2:49 pm
        # Always include Z. Delete the case without Z
        psi[1:p[num_m],1:p[num_m],num_m] = (t(X[,1:p[num_m],num_m]- Z %*% t(B[1:p[num_m],,num_m])) %*% (X[,1:p[num_m],num_m]- Z %*% t(B[1:p[num_m],,num_m])) - t(X[,1:p[num_m],num_m]- Z %*% t(B[1:p[num_m],,num_m])) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m],n))) %*% t(H[1:p[num_m],1:r_list[num_m],num_m]))/n
        psi[1:p[num_m],1:p[num_m],num_m] = diag(diag(psi[1:p[num_m],1:p[num_m],num_m]))
        
      }###########
      
      
      #*************************************************************************************************************
      # u/sigma
      # 12/26/2016 9:52PM
      # Adding L21 to beta
      # 12/27/2016
      # Borrow from adaptive group lasso (gglasso package)
      
      #y_beta
      a_1 = NULL
      for(num_m in 1:m){
        a_1 = c(a_1, kronecker(as.vector(t(sqrt(ph.y))), rep(1, r_list[num_m])))
      }
      #diag(a_1)
      #b_1 = NULL
      for(num_m in 1:m){
        for(num_k in 1:k){
          #print(sigma[num_k,1:r_list[num_m],1:r_list[num_m],num_m])
          sigma_mk_square_root = chol(solve(sigma[num_k,1:r_list[num_m],1:r_list[num_m],num_m]))
          if(num_m == 1 & num_k == 1){
            b_1 = kronecker(diag(n), sigma_mk_square_root)
          }else{
            b_1 = bdiag(b_1, kronecker(diag(n), sigma_mk_square_root))
          }
          #cat(num_m, num_k, "\n")	
        }
      }	
      
      c_1 = NULL
      for(num_m in 1:m){
        for(num_k in 1:k){
          c_1 = c(c_1, as.vector(Ez.hy[num_k, 1:r_list[num_m],,num_m]))#### 12/30 6:30pm Corrected!
        }
      }
      y_beta = as.vector(diag(a_1) %*% b_1 %*% c_1)
      
      #X_beta
      for(num_m in 1:m){
        for(num_k in 1:k){
          if(num_m == 1& num_k == 1){
            sigma_mk_square_root = chol(solve(sigma[num_k,1:r_list[num_m],1:r_list[num_m],num_m]))
            X_beta = diag(as.vector(kronecker(matrix(sqrt(ph.y[num_k, ]), ncol = 1), rep(1, r_list[num_m])))) %*% kronecker(rep(1, n), sigma_mk_square_root)
          }else{
            sigma_mk_square_root = chol(solve(sigma[num_k,1:r_list[num_m],1:r_list[num_m],num_m]))
            X_beta = bdiag(X_beta, diag(as.vector(kronecker(matrix(sqrt(ph.y[num_k, ]), ncol = 1), rep(1, r_list[num_m])))) %*% kronecker(rep(1, n), sigma_mk_square_root))
          }
        }
      }
      
      X_beta = as.matrix(X_beta)
      #dim(X_beta)
      
      
      # Put beta into fit$beta
      for(num_m in 1:m){
        for(num_k in 1:k){
          if(num_m  == 1 & num_k == 1){
            beta_before = beta[num_k, 1:r_list[num_m], num_m]
          }else{
            beta_before = c(beta_before, beta[num_k, 1:r_list[num_m], num_m])
          }
        }
      }
      
      #before_val_beta = 1/2 * t(y_beta - X_beta %*% beta_before) %*% (y_beta - X_beta %*% beta_before) + lambda_beta * L21_norm(group_beta, beta_before)
      #cat("before beta", before_val_beta, "\n")
      
      fit_beta = gglasso(X_beta, y_beta, group = group_beta, loss = "ls", lambda = lambda_beta/nrow(X_beta), pf = pf_beta,intercept = FALSE)$beta
      
      #beta_after = fit_beta
      #after_val_beta = 1/2 * t(y_beta - X_beta %*% beta_after) %*% (y_beta - X_beta %*% beta_after) + lambda_beta * L21_norm(group_beta, beta_after)
      #cat("after beta", after_val_beta, "\n")
      
      #solve(t(X_beta) %*% X_beta)%*% (t(X_beta) %*% y_beta)
      
      # Put fit_beta into beta		
      beta_iterate = array(dim = c(k, r, m));
      start = 1
      for(num_m in 1:m){
        for(num_k in 1:k){
          end = start + r_list[num_m] - 1
          #cat(start, end, "\n")
          beta_iterate[num_k, 1:r_list[num_m], num_m] = as.vector(fit_beta[start:end]) 
          start = end + 1
        }
      }
      beta = beta_iterate
      
      #sigma
      sigma_iterate = array(dim = c(k, r, r, m));
      for(num_m in 1:m){			
        for(i in 1:k){
          # Give the same sigma to all clusters
          if(i == 1){
            temp_sigma = apply(aperm(array(ph.y[i,],c(n,r_list[num_m],r_list[num_m])),c(2,3,1))*(array(Ezz.hy[i,1:r_list[num_m],1:r_list[num_m],,num_m],c(r_list[num_m],r_list[num_m],n))-array(beta_iterate[i,1:r_list[num_m],num_m]%*%t(x.z)%*%x.z%*%matrix(t(beta_iterate[i,1:r_list[num_m],num_m]),nrow=q.z)/n,c(r_list[num_m],r_list[num_m],n))),1,rowMeans)
          }else{
            temp_sigma = temp_sigma + apply(aperm(array(ph.y[i,],c(n,r_list[num_m],r_list[num_m])),c(2,3,1))*(array(Ezz.hy[i,1:r_list[num_m],1:r_list[num_m],,num_m],c(r_list[num_m],r_list[num_m],n))-array(beta_iterate[i,1:r_list[num_m],num_m]%*%t(x.z)%*%x.z%*%matrix(t(beta_iterate[i,1:r_list[num_m],num_m]),nrow=q.z)/n,c(r_list[num_m],r_list[num_m],n))),1,rowMeans)
          }
        }
        temp_sigma = temp_sigma/sum(rowMeans(ph.y))
        for(i in 1:k){
          sigma_iterate[i,1:r_list[num_m],1:r_list[num_m],num_m] = temp_sigma
        }
      }
      sigma = sigma_iterate
      
      # Identification
      for(num_m in 1:m){
        #************
        temp1<-matrix(0,r_list[num_m],r_list[num_m])
        temp2<-matrix(0,r_list[num_m],r_list[num_m])
        temp3<-matrix(0,r_list[num_m],1)
        for(i in 1:k){
          #sigma[i,1:r_list[num_m],1:r_list[num_m],num_m]<-apply(aperm(array(ph.y[i,],c(n,r_list[num_m],r_list[num_m])),c(2,3,1))*(array(Ezz.hy[i,1:r_list[num_m],1:r_list[num_m],,num_m],c(r_list[num_m],r_list[num_m],n))-array(beta[i,1:r_list[num_m],num_m]%*%t(x.z)%*%x.z%*%matrix(t(beta[i,1:r_list[num_m],num_m]),nrow=q.z)/n,c(r_list[num_m],r_list[num_m],n))),1,rowMeans)/mean(ph.y[i,])
          temp1<-temp1+mean(W[i,])*sigma[i,1:r_list[num_m],1:r_list[num_m],num_m]
          dep<-beta[i,1:r_list[num_m],num_m]%*%t(x.z)
          dep<-(dep%*%t(dep))/n
          temp2<-temp2+mean(W[i,])*(dep)
          temp3<-temp3+matrix(mean(W[i,])*rowMeans(beta[i,1:r_list[num_m],num_m]%*%t(x.z)))
        }
        
        #Identifiability correction
        #Cov(z) = E(zz') - E(z)E(z)' = (Cov(z) + E(z)E(z')) - E(z)E(z)' 
        var.z<-temp1+temp2-temp3%*%t(temp3) 
        A<-(chol(var.z)) 
        for (i in 1:k) {
          sigma[i,1:r_list[num_m],1:r_list[num_m],num_m]<-t(ginv(A))%*%sigma[i,1:r_list[num_m],1:r_list[num_m],num_m]%*%ginv(A)
          beta[i,1:r_list[num_m],num_m]<-t(ginv(A))%*%beta[i,1:r_list[num_m],num_m]
        }
        H[1:p[num_m],1:r_list[num_m],num_m]<-H[1:p[num_m],1:r_list[num_m],num_m]%*%t(A)
      }
      
      
      temp_numeric = sum(log(ifelse(colSums(W* py.h * 10 ^ py.h_total)==0,1e-10,colSums(W* py.h * 10 ^ py.h_total)))) + sum(py.h_pro[k, ])*log(10)
      
      temp_total<-temp_numeric
    }

    
    
    if(!is.null(Xc)){
      #print("m step for ordinal")
      m_ordinal=list()
      for(i in 1:mc){
        print(i)
        m_ordinal=list.append(m_ordinal,ordinal_mstep(p.mix.y=p.mix.y$p.mix.y , e_ordinal=e_ordinal[[i]], init_ordinal=init_ordinal[[i]], lambda_cnlm))
        init_ordinal[[i]]=m_ordinal[[i]]$init_ordinal
        temp_total=temp_total+m_ordinal[[i]]$temp_categorical
        temp_categorical=temp_categorical+m_ordinal[[i]]$temp_categorical
      }
    }
    #**************************************************************************
    #modality-selection
    print("modality-selection")
    if(!is.null(Xc)){
      m_ordinal=modality_selection(p.mix.y,init_ordinal,e_ordinal,m_ordinal,mc,lambda_cnlm2)
      for(i in 1: mc){
        init_ordinal[[i]]=m_ordinal[[i]]$init_ordinal
      }
    }

    #**************************************************************************

      
    likelihood<-c(likelihood,temp_total)
    likelihood1<-c(likelihood1,temp_numeric)
    likelihood2<-c(likelihood2,temp_categorical)
    
    
    ratio = abs((temp_total - lik)/lik);
    
    # ratio = ifelse(is.na(ratio),1, ratio)
    
    #cat("ratio:", ratio, "\n")
    #cat("temp_total:", temp_total, "\n")
    #cat("lik:", lik, "\n")
    #cat("temp_numerical:", likelihood1, "\n")
    #cat("temp_categorical:", likelihood2, "\n")
    cat("round, ratio, temp, lik", hh, ratio, temp_total, lik, "\n")
    

    lik = temp_total;

  }
  
  if (k>1){
    index<-(apply(p.mix.y$p.mix.y, 2, which.max))
  }else{
    index<-rep(k,n)
  }
  fit = NULL
  fit$alpha = NULL
  fit$alpha2 = NULL
  fit$Sensitivity  = array(0, mc)
  fit$Specificity  = array(0, mc)
  for(i in 1:mc){
    alpha = ifelse(abs(m_ordinal[[i]]$init_ordinal$alpha) < 1e-4, 0, m_ordinal[[i]]$init_ordinal$alpha);
    fit$alpha =cbind(fit$alpha, alpha);
    fit$alpha2 = cbind(fit$alpha2, m_ordinal[[i]]$init_ordinal$alpha2);
    fit$Sensitivity[i] = ifelse(sum(m_ordinal[[i]]$init_ordinal$muc  != 0)!=0 , (sum(alpha[1:16,] != 0 )) / (rc * 16 - 1), NA)
    fit$Specificity[i] = ifelse( sum(m_ordinal[[i]]$init_ordinal$muc  != 0)!=0, sum(alpha[17:init_ordinal[[i]]$pc, ] == 0 ) / (rc * (init_ordinal[[i]]$pc - 16)), NA)
  }
  
  
  
  h = (k-1)
  if(!is.null(X)){
    if(!is.null(Z)){
      h = h + p_list[num_m] + sum(B[1:p_list[num_m], , num_m] != 0) + sum(H[1:p_list[num_m], 1:r_list[num_m], num_m] != 0) + sum(beta[, 1:r_list[num_m], num_m] != 0) - r_list[num_m]
      h = h + p_list * r - (r*(r-1)/2) + (r*(r+1)/2)*(k-1) + r*(k-1);
    }
  }
  
  
  
  
  fit$muc = NULL
  fit$sigmac = NULL
  for(i in 1:mc){
    muc = ifelse(abs(m_ordinal[[i]]$init_ordinal$muc) < 1e-4, 0, m_ordinal[[i]]$init_ordinal$muc);
    fit$muc =cbind(fit$muc, muc);
    fit$sigmac = cbind(fit$sigmac, m_ordinal[[i]]$init_ordinal$sigmac);
  }
  rc = rc_list[1]
  if(!is.null(Xc)){
    h = h + sum(fit$sigmac != 0) + sum(fit$muc != 0) - (rc*(rc + 1)/2 + rc)*mc + sum(fit$alpha != 0) + sum((j_list - 1)[1:init_ordinal[[1]]$pc])*mc + sum((j_list - 1)[1:init_ordinal[[1]]$pc])*mc*init_ordinal[[1]]$pc
  }
  
  pen<-h*log(n)
  lik<-likelihood[length(likelihood)]
  bic<--2*lik+pen
  aic<--2*lik+2*h
  lik<-likelihood[length(likelihood)]
  
  fit$lik_seq = likelihood
  #fit$p.mix.y=p.mix.y
  #fit$w = W[, 1]
  fit$index = index
  fit$lik = lik
  #fit$resid = resid
  #fit$B = B
  fit$bic = bic
  fit$aic = aic
  fit$h = h
  fit$acc = rand.index(index, true_memb);
  time=proc.time() - ptm
  fit$time=time[3]
  
  return(fit)
  
  
  
}


