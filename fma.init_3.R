
#data = as.matrix(read.csv("data_s.csv"))
`fma.init` <-
function(y, Z, k,r,x.z=NULL,x.w=NULL,seed=4,scaling=FALSE)
{
#y = X[, , num_m];
#y = X[[1]][,,num_m]
#Z = Z;
#k = k;
#r = r_list[num_m];
#x.z=NULL; x.w=NULL; seed=4; scaling=scaling;

library("psych")
set.seed(seed)
#if (scaling) y<-scale(y)

numobs<-nrow(y)
p<-ncol(y)
ybar<- apply(y, 2, mean)
y<-scale(y, ybar, scale=scaling) #centering by col

# estimate B
if(!is.null(Z)){
	y_resid = matrix(nrow = nrow(y), ncol = ncol(y))
	B = matrix(nrow = p, ncol = ncol(Z))
	for(j in 1:p){
		y_1 = y[, j]
		# y and Z must be centered!!!!
		fit = lm(y_1~ 0 + Z)
		B[j, ] = summary(fit)$coef[, 1]
		resid_1 = y_1 - Z %*% summary(fit)$coef[, 1]
		y_resid[, j] = resid_1
	}
}else{
	B = NULL
	y_resid = y
}

x.z<-cbind(rep(1,numobs),x.z)
q.z<-ncol(x.z)

x.w<-cbind(rep(1,numobs),x.w)
q.w<-ncol(x.w)


psi<-0.2*var(y_resid)
psi<-diag(diag(psi))

#H<-factanal(y_resid,r, rotation = "none")$load
#H = fa(y_resid, r, rotate = "none")$loadings

# Changed on 11/16/2019 due to algorithm not converging (keeping running)

H = fa(cor(y_resid), r, rotate = "none")$loadings

#output<-hc(modelName = "VII", data = y)
#memb<-hclass(output,k)

if (k>1) memb<-kmeans(y_resid,k)$cl else memb<-rep(1,numobs)
### Beta ? un array che componenti x latenti x covariate
### e quindi in caso di assenza di covariate ? pari
### a componenti per latenti...

for  (i in 1:k) if ((table(memb)[i])<2) memb[sample(1:numobs,2,replace=FALSE)]=i
phi<-matrix(0,k,q.w)
	w<-table(memb)/sum(table(memb))
	w<-t(t(w))


	## d'ora in avanti w ha dimensioni k x numobs
	w<-matrix(w,nrow=k,ncol=numobs)



#2/6/2017
# Special case for k = 1
if(k == 1){
	Beta<-array(0,c(k,r,q.z))
	sigma<-array(0,c(k,r,r))
	sigma[1,,] = diag(r)
}else{
	Beta<-array(0,c(k,r,q.z))
	sigma<-array(0,c(k,r,r))
	la = NULL
	for (i in 1:k) {dati.y<-y_resid[memb==i,]
	                zz<-t(solve(t(H)%*%solve(psi)%*%H)%*%t(H)%*%solve(psi)%*%t(dati.y))
	                dati.x<-x.z[memb==i,]
	                out<-lm(zz~dati.x-1)
	                Beta[i,,]<-t(out$coe)
				la = rbind(la, out$residuals) #*** Corrected at 7:58 pm 12/29
	                #sigma[i,,]<-t(out$residuals)%*%out$residuals/out$df
	                }
	###2017/2/5
	la = matrix(la, ncol = r)

	### 12/28/2016 11:00PM
	### Estimate the same sigma
	sigma_la = cov(la)
	#sigma_la = t(la) %*% la/(length(la)-1)
	for(i in 1:k){sigma[i,,] = sigma_la}

	for (i in 1:k) {
	Beta[i,,]<-ifelse(is.na(Beta[i,,]),rowMeans(matrix(Beta[i,,],r,q.z),T),matrix(Beta[i,,],r,q.z))
	sigma[i,,]<-ifelse(is.na(sigma[i,,]),rowMeans(sigma[i,,],T),sigma[i,,])
	}
}

out<-list(H=H,psi=psi,phi=phi,w=w,Beta=Beta,sigma=sigma, memb = memb, B = B)
return(out)
}

