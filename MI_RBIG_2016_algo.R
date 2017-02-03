library(entropy)
library(sROC)
library(mixAK)
library(MVN)
library(hexbin)
library(cramer)
library(lpSolve)
library(memoise)
library(amap)
library(scales)
# library(sn)

entropy_mm <- function(x, nbins=sqrt(length(x))){  
    dx <- discretize(x, nbins)
    delta = diff(range(x))/ nbins
    hx = entropy.MillerMadow(dx, unit="log2")+log2(delta)
    hx
}

compute_tol_h0 <- function(nrow, ncol, probs=0.975,  n=1000){
  nbins <- sqrt(nrow)
  sim <- function(){
    x <- matrix(rnorm(nrow * ncol), ncol=ncol)
    mnegent <- apply(x, 2, entropy_mm)
    sum(mnegent)
  }
  tol <- sapply(seq.int(n), function(i){ sim() - sim()})
  #   tol <- quantile(abs(tol), probs=probs)
  tol
}
compute_tol_h0_m <- memoise(compute_tol_h0) 

MI_RBIG_2016 <- function(dat, N_lay=1000){
  ldat <- list()
  ldat[[1]] <- dat
  lR <- list()
  DIM = dim(dat)
  Nsamples = DIM[1]
  nbins <- floor(sqrt(Nsamples))
  DIM = DIM[2]
  delta_I <- numeric(N_lay)
  for (n in 1:N_lay){
    # marginal gaussianization
    p <- numeric(DIM)
    for(d in 1:DIM){
      margin  <-  marginal_gaussianization(dat[,d]);
      dat[, d] <- margin$x_gauss   
    }
    dat_aux = dat;
    # PCA rotation
    C  <- cov(dat)
    eig <- eigen(C);
    V <- eig$vectors
    lR[[n]] <- V 
    #     V <- rRotationMatrix(1, ncol(C))
    dat <- dat %*% V  
    ldat[[n+1]] <- dat
    delta_I[n] = information_reduction_LT(dat,dat_aux, nbins=nbins);
    rt <- roystonTest(dat, qqplot = FALSE)
    hzt <- hzTest(dat, qqplot = FALSE)
    if (runif(1)< max(rt@p.value, hzt@p.value)) break
  }
  ans <- list(ldat=ldat, lR=lR, MIs=delta_I, MI=sum(delta_I))
}

information_reduction_LT <- function(X, Y, nbins){
  #   should discretize first
  hx <- apply(X, 2, function(x)entropy.MillerMadow(discretize(x, nbins), unit="log2") + log2(diff(range(x))/nbins))
  hy <- apply(Y, 2, function(y)entropy.MillerMadow(discretize(y, nbins), unit="log2") + log2(diff(range(y))/nbins))
  I <- sum(hy - hx)
} 

marginal_gaussianization <- function(x){
  #   x_order <- order(x)
  #   x_cdfk <- kCDF(x, xgrid=x)
  #   x_unif <- x_cdfk$Fhat
  x_unif <- ecdf(x)(x)
  x_unif <- x_unif * length(x)/(length(x)+1)
  #   x_gauss <- qnorm(x_unif)[x_order]
  x_gauss <- qnorm(x_unif)
  #   ans <- list(x_gauss=x_gauss, shapiro.test=shapiro.test(x_gauss))
  ans <- list(x_gauss=x_gauss)
}

cond_MI_r <- function(dat, x_ind, y_ind, c_ind=integer(0)){
  if(length(c_ind) == 0){
    ans <- MI_RBIG_2016(dat[, c(x_ind, y_ind)])$MI
  }else{ 
	ans <- MI_RBIG_2016(dat[, c(x_ind, y_ind, c_ind)])$MI
	ans <- ans - MI_RBIG_2016(dat[, c(x_ind, c_ind)])$MI
	ans <- ans - MI_RBIG_2016(dat[, c(y_ind, c_ind)])$MI
	if(length(c_ind) > 1)
	   ans <- ans + MI_RBIG_2016(dat[, c_ind])$MI
  } 
  ans
}

cond_MI_m <- function(dat, x_ind, y_ind, c_ind=integer(0)){
  if(length(c_ind) == 0){
    ans <- RBIG_r(dat[, c(x_ind, y_ind)])
  }else{ 
    ans <- RBIG_r(dat[, c(x_ind, y_ind, c_ind)])
    ans <- ans - RBIG_r(dat[, c(x_ind, c_ind)])
    ans <- ans - RBIG_r(dat[, c(y_ind, c_ind)])
    if(length(c_ind) > 1)
      ans <- ans + RBIG_r(dat[, c_ind])
  } 
  ans
}

sample_mi <- function(dat, x_ind, y_ind){
  dat <- dat[,c(x_ind, y_ind)]
  dat[, 1] <- sample(dat[,1])
  dat
}

sample_cmi <- function(dat, x_ind, y_ind, c_ind){
  dat <- dat[,c(x_ind, y_ind, c_ind)]
  c_dist <- dist(dat[, 3:ncol(dat), drop=FALSE])
  P <- linear_permutation(c_dist)
  dat <-  cbind(P%*%dat[, 1], dat[, 2:ncol(dat)])
  dat
}

boot_mi <- function(dat, x_ind, y_ind, cond_MI=cond_MI_r){
  dat <- sample_mi(dat, x_ind, y_ind)  
  cond_MI(dat, 1, 2)
} 

boot_cmi <- function(dat, x_ind, y_ind, c_ind, cond_MI=cond_MI_r){
  dat <- sample_cmi(dat, x_ind, y_ind, c_ind)  
  cond_MI(dat, 1, 2, 3:ncol(dat))
} 

rbig_sim <- function(rbig_fit, gdat=NULL){
  ldat <- rbig_fit$ldat
  lR <- rbig_fit$lR
  if(is.null(gdat)){
    gdat <- tail(ldat, 1)[[1]]
  }
  for(n in rev(seq_along(ldat)[-1])){
    gdat <- gdat %*% solve(lR[[n-1]])
    for(d in ncol(gdat):1){
      gdat[, d] <- pnorm(gdat[, d])
      gdat[, d] <- gdat[, d] / max(gdat[, d])
      #       print(max(gdat[, d]))
      #       hist(gdat[, d])
      gdat[, d] <- quantile(ldat[[n-1]][, d], gdat[, d])
    }
  }
  gdat
}

nboot_cmi <- function(n,dat, x_ind, y_ind, c_ind=numeric(0), cond_MI=cond_MI_r){
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  if(length(c_ind) == 0)
    ans <- unlist(lapply(seq.int(n), function(i){setTxtProgressBar(pb, i); boot_mi(dat, x_ind, y_ind, cond_MI)}))
  else{
    rbig_fit <- MI_RBIG_2016(dat)
    ans <- unlist(lapply(seq.int(n), function(i){
			   setTxtProgressBar(pb, i) 
			   rbig_inv <- rbig_sim(rbig_fit, gdat=matrix(rnorm(length(c(dat))), ncol=ncol(dat), nrow(dat)))
			   boot_cmi(rbig_inv, x_ind, y_ind, c_ind, cond_MI)}))
  }
  close(pb)
  ans
}

cmi_btest <- function(nboot ,dat, x_ind, y_ind, c_ind=numeric(0), cond_MI=cond_MI_r){
  cmi <- cond_MI(dat, x_ind, y_ind, c_ind)
  ncmi <- nboot_cmi(nboot, dat, x_ind, y_ind, c_ind, cond_MI)
  df <- data.frame(stat=c(cmi, ncmi), type=rep(c("H1","H0"), c(1,nboot)))
    plot(ggplot(data=df, aes(x=stat, fill=type, color=type)) + geom_histogram(aes(y=..density..),alpha=0.5, position="identity", bins=30)+ggtitle(paste("x=",x_ind[1], "y=", y_ind[1], " S=", paste(c_ind, collapse=",")))+theme(aspect.ratio=1/3))
  p.value <- 1 - rank(c(cmi, ncmi))[1]/(length(ncmi) + 1)
  print(p.value)
  p.value
}


# code translated to R from Gary Doran et al. "A permutation-Based Kernel Conditional Independence Test
linear_permutation <- function(D){
  #   D <- as.matrix(dist(dat[1:3, 3:4]))
  D <- as.matrix(D)
  n <- nrow(D)
  # Rescale Distances
  D <- D / max(max(D))
  # Objective Function
  f <-  c(t(D))
  # Inequality contraint
  #   lb <- numeric(n^2)
  # Equality constraints
  Aeq <- matrix(0, nrow=2*n, ncol=n^2)
  b <- matrix(1, nrow=2*n, ncol=1)
  # Columns sum to 1
  for(c in 0:n-1){
    Aeq[c + 1, (c*n+1):((c+1)*n)] <- 1
  }
  # Rows sum to 1 (last row constraint not necessary
  # it is implied by other constraints)
  for(r in 1:(n-1)){
    for(c in 1:n){
      Aeq[r+n, r+(c-1)*n] <- 1
    }
  }
  # Diagonal entries zero
  for (z in 1:n){
    Aeq[2*n, (z-1)*(n+1) + 1] <- 1
  }
  b[2*n, 1] <- 0
  cdir <- paste(rep("=", 2*n))
  ans <- lp (direction = "min", objective.in=f, const.mat=Aeq, const.dir=cdir, const.rhs=b, transpose.constraints = TRUE, all.int=TRUE, all.bin=TRUE)
  ans <- matrix(ans$sol, ncol=n, byrow=FALSE) #%*% D
  ans
}


KCIPT <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    if(i < B/2){
      omega2 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    }else{
      omega2 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    } 
    MMD[i] <- cramer.test_simple(omega1, omega2)
    omega <- rbind(omega1, omega2)
    for( j in 1:b){
      idx <- sample.int(nrow(dat), round(nrow(dat)/2))
      omega1 <- omega[idx, ]
      omega2 <- omega[-idx, ]
      inner_null[i, j] <- cramer.test_simple(omega1, omega2)
    }
    cat("*")
  }
  cat("\n")
  statistics <- median(MMD)
  for(k in 1:M){
    for(i in 1:B){
      r <- ceiling(runif(1) * b)
      samples[i] <- inner_null[i, r]
    }
    outer_null[k] <- median(samples)
  }
  p.value <- 1 - rank(c(statistics, outer_null))[1]/(length(outer_null) + 1)
  p.value
}

# from the cramer packages
cramer.test_simple <- function(x, y, kernel="phiCramer"){
  .cramer.statistic<-function(daten,indexe,mm,nn,lookup) {
    xind<-indexe[1:mm]
    yind<-indexe[(mm+1):(mm+nn)]
    mm*nn/(mm+nn)*(2*sum(lookup[xind,yind])/(mm*nn)-sum(lookup[xind,xind])/(mm^2)-sum(lookup[yind,yind])/(nn^2))
  }
  m<-nrow(x)
  n<-nrow(y)
  daten<-matrix(c(t(x),t(y)),ncol=ncol(x),byrow=TRUE)
  lookup<-eval(call(kernel, as.matrix(Dist(daten))))
  .cramer.statistic(daten,1:(m+n),m,n,lookup)
}

RBIG_kcipt <- function(dat, xy_ind, c_ind=numeric(0),  dist, B, b, M){
  MMD <- numeric(B)
  samples <- numeric(B)
  inner_null <- matrix(numeric(B*b), nrow=B)
  outer_null <- numeric(M)
  dat <- as.matrix(dat)
  dat <- dat[, c(xy_ind, c_ind)]
  for( i in 1:B){
    omega <- dat
    idx <- sample.int(nrow(omega), round(nrow(omega)/2))
    omega1 <- omega[idx, ]
    omega2 <- omega[-idx, ]
    P <- linear_permutation(dist(omega2[, 3:ncol(omega2)]))
    omega21 <-  cbind(P%*%omega2[, 1], omega2[, 2:ncol(omega2)])
    omega22 <-  cbind(omega2[, 1], P%*%omega2[, 2], omega2[, 3:ncol(omega2)])
    omega2 <- rbind(omega21, omega22)[sample.int(nrow(omega1)), ]
    MMD[i] <- cond_MI(omega1, 1, 2, c_ind=3:ncol(omega1))
    omega <- rbind(omega1, omega2)
    for( j in 1:b){
      idx <- sample.int(nrow(dat), round(nrow(dat)/2))
      omega2 <- omega[-idx, ]
      inner_null[i, j] <- cond_MI(omega2, 1, 2, c_ind=3:ncol(omega2))
    }
    cat("*")
  }
  cat("\n")
  statistics <- median(MMD)
  for(k in 1:M){
    for(i in 1:B){
      r <- ceiling(runif(1) * b)
      samples[i] <- inner_null[i, r]
    }
    outer_null[k] <- median(samples)
  }
  p.value <- 1 - rank(c(statistics, outer_null))[1]/(length(outer_null) + 1)
  p.value
}
# RBIG_kcipt(head(dat, 700), 1:2, 3, dist, 10, 20, 100)

