modality_selection=function(p.mix.y,init_ordinal,e_ordinal,m_ordinal,mc,lambda_cnlm2){
  #This function realizes modality-selection across categorical modalities
  q.z<-1
  rc = init_ordinal[[1]]$rc
  k = init_ordinal[[1]]$k
  n= init_ordinal[[1]]$n
  rc_list = rep(rc, mc)
  W = array(0, c(k, n))
  for(i in 1:mc){
    W = init_ordinal[[i]]$W
  }
  W = W / 2
  #####################################################################
  # dim of p.mix.y:  is (k * n)
  # dim of sigma: from m_ordinal:: (k * rc * rc * mc)
  # dim of muc: from m_ordinal:: (k * rc * mc)
  # dim of Ez.sy: from e_ordinal: from (n * k * rc) to (k * rc * n * mc)
  # dim of E.zz.sy: from e_ordinal: from (n * k * rc * rc) to (k * rc * rc * n * mc)
  #####################################################################
  p.mix.y= p.mix.y$p.mix.y
  sigma = array(0, c(k, rc, rc, mc))
  muc = array(0, c(k, rc, mc))
  Ez.sy =  array(0, c(k, rc, n, mc))
  E.zz.sy = array(0, c(k, rc, rc, n, mc))
  for(i in 1:mc){
    sigma[,,,i] = m_ordinal[[i]]$init_ordinal$sigmac
    muc[ , , i] = m_ordinal[[i]]$init_ordinal$muc
    for(j in 1: k){
      Ez.sy[j, , , i] = t(e_ordinal[[i]]$E.z.sy[, j, ])
      for(p in 1:n){
        E.zz.sy[j, , , p, i] = e_ordinal[[i]]$E.zz.sy[p, j, , ]
      }
    }
  }
  
  
  #####################################################################
  a_1 = NULL
  for(num_m in 1:mc){
    a_1 = c(a_1, kronecker(as.vector(t(sqrt(p.mix.y))), rep(1, rc_list[num_m])))
  }
  #diag(a_1)
  #b_1 = NULL
  for(num_m in 1:mc){
    for(num_k in 1:k){
      #print(sigma[num_k,1:rc_list[num_m],1:rc_list[num_m],num_m])
      sigma_mk_square_root = chol(solve(sigma[num_k,1:rc_list[num_m],1:rc_list[num_m],num_m]))
      if(num_m == 1 & num_k == 1){
        b_1 = kronecker(diag(n), sigma_mk_square_root)
      }else{
        b_1 = bdiag(b_1, kronecker(diag(n), sigma_mk_square_root))
      }
      #cat(num_m, num_k, "\n")
    }
  }
  
  c_1 = NULL
  for(num_m in 1:mc){
    for(num_k in 1:k){
      c_1 = c(c_1, as.vector(Ez.sy[num_k, 1:rc_list[num_m],,num_m]))#### 12/30 6:30pm Corrected!
    }
  }
  y_muc = as.vector(diag(a_1) %*% b_1 %*% c_1)
  
  #X_muc
  for(num_m in 1:mc){
    for(num_k in 1:k){
      if(num_m == 1& num_k == 1){
        sigma_mk_square_root = chol(solve(sigma[num_k,1:rc_list[num_m],1:rc_list[num_m],num_m]))
        X_muc = diag(as.vector(kronecker(matrix(sqrt(p.mix.y[num_k, ]), ncol = 1), rep(1, rc_list[num_m])))) %*% kronecker(rep(1, n), sigma_mk_square_root)
      }else{
        sigma_mk_square_root = chol(solve(sigma[num_k,1:rc_list[num_m],1:rc_list[num_m],num_m]))
        X_muc = bdiag(X_muc, diag(as.vector(kronecker(matrix(sqrt(p.mix.y[num_k, ]), ncol = 1), rep(1, rc_list[num_m])))) %*% kronecker(rep(1, n), sigma_mk_square_root))
      }
    }
  }
  
  X_muc = as.matrix(X_muc)
  #dim(X_muc)
  group_muc = NULL;
  for(num_m in 1:mc){
    group_muc = c(group_muc, rep(num_m, k * rc_list[num_m]))
  }
  pf_muc = vector(length = mc)
  for(num_m in 1:mc){
    pf_muc[num_m] = 1
  }
  
  # Put muc into fit$muc
  for(num_m in 1:mc){
    for(num_k in 1:k){
      if(num_m  == 1 & num_k == 1){
        muc_before = muc[num_k, 1:rc_list[num_m], num_m]
      }else{
        muc_before = c(muc_before, muc[num_k, 1:rc_list[num_m], num_m])
      }
    }
  }
  
  
  fit_muc = gglasso(X_muc, y_muc, group = group_muc, loss = "ls", lambda = lambda_cnlm2/nrow(X_muc), pf = pf_muc,intercept = FALSE)$beta
  
  # Put fit_muc into muc
  muc_iterate = array(dim = c(k, rc, mc));
  start = 1
  for(num_m in 1:mc){
    for(num_k in 1:k){
      end = start + rc_list[num_m] - 1
      #cat(start, end, "\n")
      muc_iterate[num_k, 1:rc_list[num_m], num_m] = as.vector(fit_muc[start:end])
      start = end + 1
    }
  }
  muc = muc_iterate
  #sigma
  sigma_iterate = array(dim = c(k, rc, rc, mc));
  x.z = rep(1, n)
  for(num_m in 1:mc){
    for(i in 1:k){
      # Give the same sigma to all clusters
      if(i == 1){
        temp_sigma = apply(aperm(array(p.mix.y[i,],c(n,rc_list[num_m],rc_list[num_m])),c(2,3,1))*(array(E.zz.sy[i,1:rc_list[num_m],1:rc_list[num_m],,num_m],c(rc_list[num_m],rc_list[num_m],n))-array(muc_iterate[i,1:rc_list[num_m],num_m]%*%t(x.z)%*%x.z%*%matrix(t(muc_iterate[i,1:rc_list[num_m],num_m]),nrow=q.z)/n,c(rc_list[num_m],rc_list[num_m],n))),1,rowMeans)
      }else{
        temp_sigma = temp_sigma + apply(aperm(array(p.mix.y[i,],c(n,rc_list[num_m],rc_list[num_m])),c(2,3,1))*(array(E.zz.sy[i,1:rc_list[num_m],1:rc_list[num_m],,num_m],c(rc_list[num_m],rc_list[num_m],n))-array(muc_iterate[i,1:rc_list[num_m],num_m]%*%t(x.z)%*%x.z%*%matrix(t(muc_iterate[i,1:rc_list[num_m],num_m]),nrow=q.z)/n,c(rc_list[num_m],rc_list[num_m],n))),1,rowMeans)
      }
    }
    temp_sigma = temp_sigma/sum(rowMeans(p.mix.y))
    for(i in 1:k){
      sigma_iterate[i,1:rc_list[num_m],1:rc_list[num_m],num_m] = temp_sigma
    }
  }
  sigma = sigma_iterate
  
  
  for(j in 1:mc){
    temp1=matrix(0,rc,rc)
    temp2=matrix(0,rc,rc)
    temp3=matrix(0,rc)
    for (i in 1:k) {
      temp1=temp1+W[i, 1]*sigma[i,,, j]
      temp2=temp2+W[i, 1]*(muc[i,, j]%*%t(muc[i,, j]))
      temp3=temp3+W[i, 1]*muc[i,, j]
    }
    
    ## correzione per l'identificabilit?
    var.z=temp1+temp2-temp3%*%t(temp3)
    A<-(chol(var.z))
    for (i in 1:k) {
      sigma[i,,, j]<-t(ginv(A))%*%sigma[i,,, j]%*%ginv(A)
      muc[i , , j]<-t(ginv(A))%*%muc[i,, j]
    }
    
    mu.tot=t(W[, 1])%*%muc[ , , j]
    muc[ , , j]=muc[ , , j]-t(matrix(mu.tot,rc,k))
  }
  
  for(i in 1:mc){
    m_ordinal[[i]]$init_ordinal$sigmac = sigma[,,,i]
    m_ordinal[[i]]$init_ordinal$muc = muc[ , , i]
  }
  
  
  
  
  return(m_ordinal)
}

