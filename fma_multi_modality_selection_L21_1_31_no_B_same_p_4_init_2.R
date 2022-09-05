fma.multi.diff.L21.same.p = function(X,Z,k, r_list, lambda_beta, lambda_gamma, x.z = NULL, x.w = NULL, it = 15, eps = 0.0001, epsilon = 1e-8, seed = 4, scaling = FALSE){
	library("MASS")
	library("mvtnorm")
	library("Matrix")
	library("gglasso")
	library("psych")
	source("fma.init_3.R")

	ptm = proc.time()
	#cat(seed, "\n")
	#set.seed(seed)
	#if(scaling)
	
	#*** r is redefined!
	r = max(r_list)
	n = dim(X)[1]
	p = dim(X)[2]
	m = dim(X)[3]
	q = ncol(Z)
	x.z = rep(1, n)
	x.w = rep(1, n)
	q.z<-1 
	q.w<-1
	if(!is.null(Z)){
		q = ncol(Z)
	}


	# Initilization
	H = array(dim = c(p, r, m));
	psi = array(dim = c(p, p, m));
	W = array(0, c(k, n))
	beta = array(dim = c(k, r, m));
	sigma = array(dim = c(k, r, r, m));
	if(!is.null(Z)){
		q = ncol(Z);
		B = array(dim = c(p, q, m));
		
	}

	for(num_m in 1:m){
		init_m = fma.init(X[, , num_m], Z, k, r_list[num_m],x.z=NULL,x.w=NULL,seed=4,scaling=scaling);
		H[, 1:r_list[num_m], num_m] = init_m$H;
		psi[, , num_m] = init_m$psi;
		beta[, 1:r_list[num_m], num_m] = init_m$Beta;
		sigma[, 1:r_list[num_m], 1:r_list[num_m], num_m] = init_m$sigma;	
		if(!is.null(Z)){
			B[,, num_m] = init_m$B;
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
		for(num_p in 1:p){
			group_gamma = c(group_gamma, rep(label, r_list[num_m]))
			#label = label + 1
			group_gamma = c(group_gamma, rep(label, q))
			label = label + 1

		}
	}
	pf_gamma = NULL
	for(num_m in 1:m){
		for(num_p in 1:p){
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
	while((hh < it)&(ratio > eps)){
		hh = hh + 1
		
		# E-step
		Ezz.hy = array(0, c(k, r, r, n, m));
		Ez.hy = array(0, c(k,r,n,m))
		sigma.tot = array(0, c(k, p, p, m));
 		change = rep(FALSE, k)####OIS
		
		for(num_m in 1:m){
			if((k == 1)&(r_list[num_m] == 1)){
				sigma.tot[1,,,num_m] = as.matrix(H[,1:r_list[num_m],num_m]) %*% sigma[,1:r_list[num_m],1,num_m] %*% t(as.matrix(H[,1:r_list[num_m],num_m])) + psi[,,num_m];
				if(sum(diag(sigma.tot[1,,,num_m])) < 0){change = TRUE;}#####
				if(det(sigma.tot[1,,,num_m]) < 0){change = TRUE;}########
			}else{
				for(i in 1:k){
					sigma.tot[i,,,num_m] = as.matrix(H[,1:r_list[num_m],num_m]) %*% sigma[i,1:r_list[num_m],1:r_list[num_m],num_m] %*% t(as.matrix(H[,1:r_list[num_m],num_m])) + psi[,,num_m];	
					if(sum(diag(sigma.tot[i,,,num_m])) < 0){change[i] = TRUE;}########
					if(det(sigma.tot[i,,,num_m]) < 0){change[i] = TRUE;}#############
				}
			}   
			
			if((k == 1) & (r_list[num_m] == 1)){
				chsi[,1:r_list[num_m],1:r_list[num_m],num_m] = 1/(t(H[,1:r_list[num_m],num_m]) %*% ginv(psi[,,num_m]) %*% H[,1:r_list[num_m],num_m] + solve(sigma[,1:r_list[num_m],1:r_list[num_m],num_m]))
				if(!is.null(Z)){
					roy[,1:r_list[num_m],,num_m] = chsi[,1:r_list[num_m],1:r_list[num_m],num_m] %*% (t(as.matrix(H[,1:r_list[num_m],num_m])) %*% ginv(psi[,,num_m]) %*% t(X[,,num_m] - Z %*% t(B[,,num_m])) + solve(sigma[,1:r_list[num_m],1,num_m]) %*% beta[1,1:r_list[num_m],num_m] %*% t(x.z))
				}else{
					roy[,1:r_list[num_m],,num_m] = chsi[,1:r_list[num_m],1:r_list[num_m],num_m] %*% (t(as.matrix(H[,1:r_list[num_m],num_m])) %*% ginv(psi[,,num_m]) %*% t(X[,,num_m]) + solve(sigma[,1:r_list[num_m],1,num_m]) %*% beta[1,1:r_list[num_m],num_m] %*% t(x.z))
				}
				Ezz.hy[1,1:r_list[num_m],1:r_list[num_m],,num_m] = matrix(chsi[,1:r_list[num_m],1:r_list[num_m],num_m], ncol = n) + roy[1,1:r_list[num_m],,num_m]^2
				if(change){sigma.tot[1,,,num_m] = var(X[,,num_m])}
			}else{
				for(i in 1:k){
					chsi[i,1:r_list[num_m],1:r_list[num_m],num_m]<-ginv(t(H[,1:r_list[num_m],num_m])%*%ginv(psi[,,num_m])%*%H[,1:r_list[num_m],num_m]+ginv(sigma[i,1:r_list[num_m],1:r_list[num_m],num_m]))
					if(!is.null(Z)){
						roy[i,1:r_list[num_m],,num_m]<-chsi[i,1:r_list[num_m],1:r_list[num_m],num_m]%*%(t(H[,1:r_list[num_m],num_m])%*%ginv(psi[,,num_m])%*%t(X[,,num_m] - Z %*% t(B[,,num_m]))+ginv(sigma[i,1:r_list[num_m],1:r_list[num_m],num_m])%*%beta[i,1:r_list[num_m],num_m]%*%t(x.z))
					}else{
						roy[i,1:r_list[num_m],,num_m]<-chsi[i,1:r_list[num_m],1:r_list[num_m],num_m]%*%(t(H[,1:r_list[num_m],num_m])%*%ginv(psi[,,num_m])%*%t(X[,,num_m])+ginv(sigma[i,1:r_list[num_m],1:r_list[num_m],num_m])%*%beta[i,1:r_list[num_m],num_m]%*%t(x.z))
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
					if (change[i]){sigma.tot[i,,,num_m]<-var(X[,,num_m])}#########OIS
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
					resid_m = X[num_n,1:p,num_m]- t(as.matrix(H[1:p,1:r_list[num_m],num_m])%*%beta[i,1:r_list[num_m],num_m]) - Z[num_n, ] %*% t(B[1:p,,num_m])
					sigma_m = sigma.tot[i,1:p,1:p,num_m]
					#dim(resid_m)
					#dim(solve(sigma_m))		
					#(2*pi)^(-p/2) * det(sigma_m)^(-1/2)*exp(-1/2 * resid_m %*% as.matrix(solve(sigma_m)) %*% t(resid_m))
					temp_p_2 = dmvnorm(resid_m,mean = rep(0, p),sigma=as.matrix(sigma_m))
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
		for(num_n in 1:n){
			for(i in 1:k){
				ph.y[i,num_n]<-W[i,num_n]*(py.h[i,num_n]*10^(py.h_total[i, num_n]))/(t(W[, num_n])%*%(py.h[,num_n] * (10 ^ py.h_total[, num_n])))
			}
		}
		
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

		# M-step
		
		# W
		#OIS
		ph.y<-ifelse(ph.y==0,0.0000000000001,ph.y) 
		W<- matrix(t(rowMeans(ph.y,na.rm=TRUE)))
		W<-W/sum(W)
		W<-matrix(W,nrow=k,ncol=n)
		#**************************************************************************************************************************
	
		#label_gamma_B
		for(num_m in 1:m){
			# Lambda/Psi 
			EEzz.y = rowMeans(array(Ezz.y[1:r_list[num_m],1:r_list[num_m],,num_m], c(r_list[num_m], r_list[num_m], n)), na.rm = TRUE, dims = 2)

			B = as.vector(t(as.matrix(cbind(t(X[,,num_m]) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m], n))), t(X[,,num_m]) %*% t(array(t(Z), c(q, n))))))) 
			A = as.matrix(cbind(rbind(n*EEzz.y, t(Z) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m], n)))), rbind(t(t(Z) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m], n)))), t(Z) %*% Z)))
			D_temp = chol(A)

			y_gamma_1 = sqrt(tr(solve(psi[,,num_m]))) * kronecker(diag(p), solve(t(D_temp))) %*% B
			X_gamma_1 = sqrt(tr(solve(psi[,,num_m]))) * kronecker(diag(p), D_temp)
  
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


		H_iterate = array(dim = c(p, r, m));
		B_iterate = array(dim = c(p, q, m))
		start = 1
		for(num_m in 1:m){
			for(num_p in 1:p){
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
			psi[,,num_m] = (t(X[,,num_m]- Z %*% t(B[,,num_m])) %*% (X[,,num_m]- Z %*% t(B[,,num_m])) - t(X[,,num_m]- Z %*% t(B[,,num_m])) %*% t(array(Ez.y[1:r_list[num_m],,num_m], c(r_list[num_m],n))) %*% t(H[,1:r_list[num_m],num_m]))/n
			psi[,,num_m] = diag(diag(psi[,,num_m]))
			
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
			H[,1:r_list[num_m],num_m]<-H[,1:r_list[num_m],num_m]%*%t(A)
		}
		
		
		# likelihoodxc  
		temp = sum(log(colSums(W* py.h * 10 ^ py.h_total))) + sum(py.h_pro[k, ])*log(10)
		log(10)
		
				
		#cat(sum(log(colSums(W* py.h * 10 ^ py.h_total))),sum(py.h_pro[k, ]),"\n")		

		#temp1 = log(colSums(W*py.h))
		#temp = 0; n_out = 0 
		#for(i in 1:n){
		#	if(is.finite(temp1[i])){
		#		temp = temp + temp1[i]
		#	}else{
		#		n_out = n_out + 1
		#	}
		#}
		#flag_iterate = TRUE
		if(is.infinite(temp)){
			cat("warning: Negative infinite likelihood! \n")
			#flag_iterate = FALSE  
			temp = 2 * lik
		}
		
		likelihood<-c(likelihood,temp) 
		
		# ratio
		ratio<-abs((temp-lik)/lik) 
		cat("ratio, temp, lik", ratio, temp, lik, "\n")
		#cat(ph.y, "\n")
		#cat(py.h, "\n")

		# Stopping Criteria
		if ((temp < lik) & (hh > 20)){ratio<-eps} 
		lik<-temp

	}
	likelihood = matrix(likelihood[!likelihood == 0])
	
	if (k>1){ 
		index<-(apply(ph.y, 2, order))[k,]
	}else{
		index<-rep(k,n)
	}
	 
	for (i in 1:k){
		if (sum(index==i)==0){
			index[c(sample(1:n,q.z))]<-i
		}
	}

	# 12/28/2016 10:31 pm
	# Delete this part to avoid an error when all modalities are out. Can be fixed later upon request. 
	# 1/31/2016 11:51 pm
	# We need to use F as inital value; Please add it back.
	#F = array(0, c(n, r, m))
	#for(num_m in 1:m){
	#	if(sum(H[,1:r_list[num_m],num_m] == 0) != 0){
	#		if(!is.null(Z)){
	#			F[,1:r_list[num_m],num_m] <- t(solve(t(H[,1:r_list[num_m],num_m])%*%solve(psi[,,num_m])%*%H[,1:r_list[num_m],num_m])%*%t(H[,1:r_list[num_m],num_m])%*%solve(psi[,,num_m])%*%t(X[,,num_m] - Z %*% t(B[,,num_m])))
	#		}else{
	#			F[,1:r_list[num_m],num_m] <- t(solve(t(H[,1:r_list[num_m],num_m])%*%solve(psi[,,num_m])%*%H[,1:r_list[num_m],num_m])%*%t(H[,1:r_list[num_m],num_m])%*%solve(psi[,,num_m])%*%t(X[,,num_m]))
	#		}
	#	}
	#}
	
	# Add residual estimate
	# Date 8/3/2016

	#resid = array(0, c(n, p, m))
	#for(num_m in 1:m){
	#	if(!is.null(Z)){
	#		resid[,,num_m] = X[,,num_m]  - F[,1:r_list[num_m],num_m] %*% t(as.matrix(H[,1:r_list[num_m],num_m])) - Z %*% t(B[,,num_m])
	#	}else{
	#		resid[,,num_m] = X[,,num_m]  - F[,1:r_list[num_m],num_m] %*% t(as.matrix(H[,1:r_list[num_m],num_m])) 
	#	}
	#	
	#}
	
	#h<-(k-1) + m*((k-1)*((r*q.z)+r*(r+1)/2)+p*(r+1))  
	# Degree of Freedom
	# Date 8/5/2016
	# 1/2/2017
	# corrected after variable selection
	h = (k-1) + m*p
	for(num_m in 1:m){
		h = h  + sum(B[, , num_m] != 0) + sum(H[, 1:r_list[num_m], num_m] != 0) + sum(beta[, 1:r_list[num_m], num_m] != 0) - r_list[num_m]
	} 


	pen<-h*log(n)
	lik<-likelihood[length(likelihood)]
	bic<--2*lik+pen
	aic<--2*lik+2*h
	fit = NULL
	fit$H = H
	fit$lik = lik
	fit$lik_seq = likelihood
	fit$W = W
	fit$beta = beta
	fit$psi = psi
	fit$sigma = sigma
	fit$ph.y = ph.y
	fit$py.h = py.h
	fit$index = index
	fit$F = F
	fit$bic = bic
	fit$aic = aic
	fit$h = h
	#fit$resid = resid
	fit$B = B
	fit$elapsed = proc.time()-ptm

	return(fit)
}


