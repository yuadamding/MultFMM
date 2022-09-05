# 5/4/2020
y2yc = function(y, c){
	#y = y[num.n, num.p]; c = j_list[num.p];
	yc = rep(0, c);
	yc[y] = 1;
	return(yc)
}

y2yc2 = function(y, c){
  #y = y[num.n, num.p]; c = j_list[num.p];
  n = length(y);
  yc = matrix(0, nrow = n, ncol = c);
  for(i in 1:n){
    yc[i, y[i]] = 1;
  }
  #yc[y] = 1;
  return(yc)
}


# Helper Function
L21_norm = function(label, x){
	max_lab = max(label)
	L21 = 0
	for(num_lab in 1:max_lab){
		L21 = L21 + sqrt(t(x[label == num_lab]) %*% x[label == num_lab])
	}
	return(L21)
}   

gauher2 <-
function (n) {
    m <- trunc((n + 1)/2)
    x <- w <- rep(-1, n)
    for (i in seq_len(m)) {
        z <- if (i == 1) {
            sqrt(2*n + 1) - 1.85575 * (2*n + 1)^(-0.16667)
        } else if (i == 2) {
            z - 1.14 * n^0.426 / z
        } else if (i == 3) {
            1.86 * z - 0.86 * x[1]
        } else if (i == 4) {
            1.91 * z - 0.91 * x[2]
        } else {
            2*z - x[i - 2]
        }
        for (its in seq_len(10)) {
            p1 <- 0.751125544464943
            p2 <- 0
            for (j in seq_len(n)) {
                p3 <- p2
                p2 <- p1
                p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3sigma.tot
            }
            pp <- sqrt(2*n) * p2
            z1 <- z
            z <- z1 - p1/pp
            if (abs(z - z1) <= 3e-14) 
                break
        }
        x[i] <- z
        x[n + 1 - i] <- -z
        w[i] <- 2 / (pp * pp)
        w[n + 1 - i] <- w[i]
    }
    list(x = x, w = w)
}


gauher <- function(n) {# Gauss-Hermite:  returns x,w so that
#\int_-\infty^\infty exp(-x^2) f(x) dx \doteq \sum w_i f(x_i)
  EPS <- 3.e-14
  PIM4 <- .7511255444649425
  MAXIT <- 10
  m <- trunc((n+1)/2)
  x <- w <- rep(-1,n)
  for (i in 1:m) {
    if (i==1) {
      z <- sqrt(2*n+1)-1.85575*(2*n+1)^(-.16667)
    } else if(i==2) {
      z <- z-1.14*n^.426/z
    } else if (i==3) {
      z <- 1.86*z-.86*x[1]
    } else if (i==4) {
      z <- 1.91*z-.91*x[2]
    } else {
      z <- 2.*z-x[i-2]
    }
    for (its in 1:MAXIT) {
      p1 <- PIM4
      p2 <- 0
      for (j in 1:n) {
        p3 <- p2
        p2 <- p1
        p1 <- z*sqrt(2/j)*p2-sqrt((j-1)/j)*p3
      }
      pp <- sqrt(2*n)*p2
      z1 <- z
      z <- z1-p1/pp
      if(abs(z-z1) <= EPS) break
    }
    x[i] <- z
    x[n+1-i] <- -z
    w[i] <- 2/(pp*pp)
    w[n+1-i] <- w[i]
  }
  list(x=x,w=w)
}


entr<-function(z)
{
numobs<-nrow(z)
numg<-ncol(z)
temp<-0
z<-ifelse(z==0,z+0.000000000000000000000001,z)
for (i in 1:numg) for (j in 1:numobs) temp<-temp+(z[j,i]*log(z[j,i]))
return(-temp)
}

misc=function(classification, truth)
{
    q <- function(map, len, x) {
        x <- as.character(x)
        map <- lapply(map, as.character)
        y <- sapply(map, function(x) x[1])
        best <- y != x
        if (all(len) == 1)
            return(best)
        errmin <- sum(as.numeric(best))
        z <- sapply(map, function(x) x[length(x)])
        mask <- len != 1
        counter <- rep(0, length(len))
        k <- sum(as.numeric(mask))
        j <- 0
        while (y != z) {
            i <- k - j
            m <- mask[i]
            counter[m] <- (counter[m]%%len[m]) + 1
            y[x == name(map)[m]] <- map[[m]][counter[m]]
            temp <- y != x
            err <- sum(as.numeric(temp))
            if (err < errmin) {
                errmin <- err
                best <- temp
            }
            j <- (j + 1)%%k
        }
        best
    }
    if (any(isNA <- is.na(classification))) {
        classification <- as.character(classification)
        nachar <- paste(unique(classification[!isNA]), collapse = "")
        classification[isNA] <- nachar
    }
    MAP <- mapClass(classification, truth)
    len <- sapply(MAP[[1]], length)
    if (all(len) == 1) {
        CtoT <- unlist(MAP[[1]])
        I <- match(as.character(classification), names(CtoT),
            nomatch = 0)
        one <- CtoT[I] != truth
    }
    else {
        one <- q(MAP[[1]], len, truth)
    }
    len <- sapply(MAP[[2]], length)
    if (all(len) == 1) {
        TtoC <- unlist(MAP[[2]])
        I <- match(as.character(truth), names(TtoC), nomatch = 0)
        two <- TtoC[I] != classification
    }
    else {
        two <- q(MAP[[2]], len, classification)
    }
    err <- if (sum(as.numeric(one)) > sum(as.numeric(two)))
        as.vector(one)
    else as.vector(two)
    bad <- seq(along = classification)[err]
    errorRate = length(bad)/length(truth)
    return(errorRate)
}


mapClass=function (a, b)
{
    l <- length(a)
    x <- y <- rep(NA, l)
    if (l != length(b)) {
        warning("unequal lengths")
        return(x)
    }
    aChar <- as.character(a)
    bChar <- as.character(b)
    Tab <- table(a, b)
    Ua <- dimnames(Tab)[[1]]
    Ub <- dimnames(Tab)[[2]]
    aTOb <- rep(list(Ub), length(Ua))
    names(aTOb) <- Ua
    bTOa <- rep(list(Ua), length(Ub))
    names(bTOa) <- Ub
    k <- nrow(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 1, max)
    for (i in 1:k) {
        I <- match(Max[i], Tab[i, ], nomatch = 0)
        aTOb[[i]] <- Ub[I]
    }
    if (is.numeric(b))
        aTOb <- lapply(aTOb, as.numeric)
    k <- ncol(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 2, max)
    for (j in (1:k)) {
        J <- match(Max[j], Tab[, j])
        bTOa[[j]] <- Ua[J]
    }
    if (is.numeric(a))
        bTOa <- lapply(bTOa, as.numeric)
    list(aTOb = aTOb, bTOa = bTOa)
}




z_ext <-function(x,nfac){

# x = pp$x; nfac = rc; 

## x    : punti di quadratura o pesi
## nfac : numero di fattori

nq <- length(x)                           # nr punti di quadratura
zx <- hcube(rep(nq,nfac))                 # calcola tutte le possibili disposizioni con ripetizione di 8 elementi presi a gruppi di 3
zx <- zx[,dim(zx)[2]:1]                   # zx contiene tutte le disposizioni "rovesciate"
z2 <- matrix(x[zx],dim(zx)[1],dim(zx)[2]) # in corrispondenza di ciascuna posizione vado a selezionare il nodo corrispondente. In questo modo
                                          # ottengo tutte le possibili combinazioni tra quadrature per ogni fattore latente.
return(z2)
}
 
 
 array.matrix<-function(A,B)
{
## A is an array of dimensions k1 x k2 x k3
## B is a matrix of dimensions k3 x k4
## returns D of dimensions k1 x k2 x k4

C<-A%o%B
D<-apply(C,c(1,2,5),diag)
D<-(apply(D,c(2,3,4),sum))
return(D)
}


da.max=function(alpha,y,zex,k,ps.y,p.z.ys)
{
nq=nrow(zex)
numobs=length(y)
temp=matrix(0,numobs)
for (i in 1:k) {
            den=log(1+exp(alpha%*%t(cbind(1,zex[,,i]))))
            den=(matrix(den,nq,numobs))
            num=t(t(t(y))%*%alpha%*%t(cbind(1,zex[,,i])))
            log.p.y.z=num-den
            temp=temp+ps.y[,i]*colSums(p.z.ys[,,i]*log.p.y.z)
            }
return(-sum(temp))
}



da.max2=function(alpha,y,zex,k,ps.y,p.z.ys)
{
nq=nrow(zex)
r=ncol(zex)
numobs=length(y)
temp=matrix(0,numobs)
temp2=matrix(0,numobs,r+1)
for (i in 1:k) {
            den=log(1+exp(alpha%*%t(cbind(1,zex[,,i])))) # p*64
            den=(matrix(den,nq,numobs))# 64*n
            num=t(t(t(y))%*%alpha%*%t(cbind(1,zex[,,i])))
            log.p.y.z=num-den
            temp=temp+ps.y[,i]*colSums(p.z.ys[,,i]*log.p.y.z)
            
            ### per gradiente
            
            den=matrix(exp(alpha%*%t(cbind(1,zex[,,i]))),nq,r+1)*cbind(1,zex[,,i])/matrix(1+exp(alpha%*%t(cbind(1,zex[,,i]))),nq,r+1)
            den=aperm(array(den,c(nq,r+1,numobs)),c(1,3,2))
            num=aperm(y%o%(cbind(1,zex[,,i])),c(2,1,3))
            log.p.y.z=num-den
            temp2=temp2+matrix(ps.y[,i],numobs,r+1)*apply(array(p.z.ys[,,i],c(nq,numobs,r+1))*log.p.y.z,c(2,3),sum)
            
            }
res=-sum(temp)
attr(res, "gradient")=colSums(temp2)
return(res)
}

# 2020/5/6

#alpha = alpha[j, ]; y = y[, j];  p2 = alpha2[j,]
#alpha = c(alpha[j,], alpha2[j,]); y = y[, j];  jc = j_list[j]; rc = rc; 

da.max3=function(alpha,y, rc, jc, zex,k,ps.y,p.z.ys)
{
  
#alpha = c(alpha[j,], alpha2[j,1:(j_list[j]-1)])
#y=y[,j]
#rc = rc
#jc= j_list[j]
#zex=zex
#k=k
#ps.y=p.mix.y2
#p.z.ys=p.z.ys

  
nq=nrow(zex)
numobs=length(y)
temp=matrix(0,numobs)
for (i in 1:k) {
            #den=log(1+exp(alpha%*%t(cbind(1,zex[,,i]))))# 1*64 
            #den=(matrix(den,nq,numobs))# 64*(n)
            #num=t(t(t(y))%*%alpha%*%t(cbind(1,zex[,,i]))) # 64* n
            
		#den=log(1+exp(alpha%*%t(cbind(1,zex[,,i]))))# 1*64 
            #den=(matrix(den,nq,numobs))# 64*(n)
            #num=t(t(t(y))%*%alpha%*%t(cbind(1,zex[,,i]))) # 64* n
		
	num = matrix(nrow = nq, ncol = numobs); 
	for(num_n in 1:numobs){
		for(num_z in 1:nq){
			temp.y = y2yc(y[num_n], jc);
			temp.z = zex[num_z,,i];
			
			lik1 = sum(temp.y[-1])*( temp.z %*% alpha[1:rc] ) - temp.y[jc] * alpha[rc + jc - 1];
			lik2 = 0; 
			for(num_j in 1:(jc-1)){
				lik2 = lik2 - sum( temp.y[num_j:(num_j + 1)] ) * log( 1 + exp( temp.z %*% alpha[1:rc] - alpha[rc+num_j] ) );
			}
			
			if(jc >= 3){
				lik3 = 0; 
				for(num_j in 2:(jc-1)){
					#lik3 = lik3 + temp.y[num_n] * log( exp(-alpha[rc+num_j-1]) - exp(-alpha[rc+num_j]) ); 
				  # 2020/6/1 corrected
				  lik3 = lik3 + temp.y[num_j] * log( exp(-alpha[rc+num_j-1]) - exp(-alpha[rc+num_j]) ); 
				}
			}else{
				lik3 = 0; 
			}
			num[num_z, num_n] = lik1+lik2+lik3;
	 	}
	}
  log.p.y.z=num#-den # 64*n
  temp=temp+ps.y[,i]*colSums(p.z.ys[,,i]*log.p.y.z)
}
return(-sum(temp))
}


da.max30=function(p,yy, rc, jc, zex,k,ps.y,p.z.ys)
{
  #yy = y[,j]
  yy2 = y2yc2(yy, jc)
  #ps.y = p.mix.y2;
  nq=nrow(zex)
  numobs=length(yy)
  temp=matrix(0,numobs)
  
  for (i in 1:k) {
    temp.zex = zex[,,i]
    
    lik1 = -log(1+exp(temp.zex %*% p[(1):(rc)] - p[rc+1])) %*% matrix(yy2[, 1], nrow=1)
    lik2 = (temp.zex %*% p[(1):(rc)] - p[rc+jc-1]-log(1+exp(temp.zex %*% p[(1):(rc)]- p[rc+jc-1])))%*% matrix(yy2[, jc], nrow=1);
    lik3 = matrix(0, nq, numobs);
    if(jc>= 3){
      for(num_j in 2:(jc-1)){
        #lik3 = lik3 + (temp.zex %*% p[(1):(rc)]-p[rc+num_j-1] + log(1-exp(p[rc+num_j-1]-p[rc+num_j])) - log(1+exp(temp.zex %*% p[(1):(rc)]- p[rc+num_j]))-log(1+exp(temp.zex %*% p[(1):(rc)] - p[rc+num_j-1])))%*% matrix(yy2[, num_j], nrow=1)
        lik3 = lik3 + (log(exp(temp.zex %*% p[(1):(rc)]-p[rc+num_j-1]) - exp( temp.zex %*% p[(1):(rc)]-p[rc+num_j])) - log(1+exp(temp.zex %*% p[(1):(rc)]- p[rc+num_j]))-log(1+exp(temp.zex %*% p[(1):(rc)] - p[rc+num_j-1])))%*% matrix(yy2[, num_j], nrow=1)
      }  
    }
    
    log.p.y.z=lik1 +lik2 + lik3;
    temp=temp+ps.y[,i]*colSums(as.matrix(p.z.ys[,,i])*log.p.y.z)
  }
  return((-sum(temp))/1000)
}





# 2020/5/21 - Debugging !!!!!!!!

#alpha = alpha[j, ]; y = y[, j];  p2 = alpha2[j,]
#alpha = c(alpha[j,], alpha2[j,]); y = y[, j];  jc = j_list[j]; rc = rc; 

da.max4=function(alpha,y, rc, jc, zex,k,ps.y,p.z.ys)
{
  nq=nrow(zex)
  numobs=length(y)
  temp=matrix(0,numobs)
  for (i in 1:k) {
    #den=log(1+exp(alpha%*%t(cbind(1,zex[,,i]))))# 1*64 
    #den=(matrix(den,nq,numobs))# 64*(n)
    #num=t(t(t(y))%*%alpha%*%t(cbind(1,zex[,,i]))) # 64* n
    
    #den=log(1+exp(alpha%*%t(cbind(1,zex[,,i]))))# 1*64 
    #den=(matrix(den,nq,numobs))# 64*(n)
    #num=t(t(t(y))%*%alpha%*%t(cbind(1,zex[,,i]))) # 64* n
    
    num = matrix(nrow = nq, ncol = numobs); 
    for(num_n in 1:numobs){
      for(num_z in 1:nq){
        temp.y = y2yc(y[num_n], jc);
        temp.z = zex[num_z,,i];
        
        #lik1 = sum(temp.y[-1])*( temp.z %*% alpha[1:rc] ) - temp.y[jc] * alpha[rc + jc - 1];
        # Debugging...
        lik1 = sum(temp.y[-1])*( temp.z %*% alpha[1:rc] ) + temp.y[jc] * alpha[rc + jc - 1];
        lik2 = 0; 
        for(num_j in 1:(jc-1)){
          lik2 = lik2 - sum( temp.y[num_j:(num_j + 1)] ) * log( 1 + exp( temp.z %*% alpha[1:rc] + alpha[rc+num_j] ) );
        }
        
        if(jc >= 3){
          lik3 = 0; 
          for(num_j in 2:(jc-1)){
            lik3 = lik3 + temp.y[num_n] * log( exp(-alpha[rc+num_j-1]) - exp(-alpha[rc+num_j]) ); 
          }
        }else{
          lik3 = 0; 
        }
        num[num_z, num_n] = lik1+lik2+lik3;
      }
    }
    log.p.y.z=num#-den # 64*n
    temp=temp+ps.y[,i]*colSums(p.z.ys[,,i]*log.p.y.z)
  }
  return(-sum(temp))
}


da.max.gr=function(alpha,y,zex,k,ps.y,p.z.ys)
{
nq=nrow(zex)
r=ncol(zex)
numobs=length(y)
temp=matrix(0,numobs,r+1)
for (i in 1:k) {
            den=matrix(exp(alpha%*%t(cbind(1,zex[,,i]))),nq,r+1)*cbind(1,zex[,,i])/matrix(1+exp(alpha%*%t(cbind(1,zex[,,i]))),nq,r+1)
            den=aperm(array(den,c(nq,r+1,numobs)),c(1,3,2))
            num=aperm(y%o%(cbind(1,zex[,,i])),c(2,1,3))
            log.p.y.z=num-den
            temp=temp+matrix(ps.y[,i],numobs,r+1)*apply(array(p.z.ys[,,i],c(nq,numobs,r+1))*log.p.y.z,c(2,3),sum)
            }
return(colSums(temp))
}




convert=function(x){
  # this function is for converting data from list form to array form
  # this results work for the numeric data in mixFMM file 
  m=length(x)
  n=dim(x[[1]])[1]
  p=dim(x[[1]])[2]
  res=array(0,c(n,p,m))
  for(i in 1:m){
    res[,,i]=x[[i]]
  }
  return(res)
}
