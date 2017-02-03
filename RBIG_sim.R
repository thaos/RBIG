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
dat <- matrix(runif(300), ncol=3)
rbig_fit <- MI_RBIG_2016(dat[, 1:2])

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
rbig_inv <- rbig_sim(rbig_fit)




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

a0=runif(1000)
a1=ecdf(a0)
a2=a1(a0)
a3=qnorm(a2)
b0=a3
b1=pnorm(b0)
b2=quantile(a0, b1)

inverse = function (f, lower = -100, upper = 100) {
     function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

set.seed(1)
n <- 1000 
x1 <- rnorm(n)
x2 <- rnorm(n)
y1 <- matrix(ncol=2, nrow=n)
y11 <- cos(x1) + 0.1*x2
# y1[, 2] <- 4*x1.^2 + 1*x2
y12 <- sin(x1) + 0.1*x2
# y1[, 3] <- 0.01*((2*x1)^3 + 0.1*x2)
# y2 <- 10*(abs(y1[, 1]))^0.5 + y1[, 2]^2 + 0*y1[,3] + 0.3*x2
y2 <- y11^3 + y12^3 + 0.3*x2
dat <- as.data.frame(scale(cbind(x1, x2, y11, y12, y2)))

pdf(file="rbig_iter.pdf")
rbig_fit <- MI_RBIG_2016(as.matrix(dat))
rbig_inv <- rbig_sim(rbig_fit, gdat=matrix(rnorm(length(as.matrix(dat))), ncol=ncol(dat), nrow=nrow(dat)))
plot(hexplom(dat))
plot(hexplom(rbig_inv))
for(i in seq.int(10)){
  rbig_fit <- MI_RBIG_2016(as.matrix(rbig_inv))
  rbig_inv <- rbig_sim(rbig_fit)
  plot(hexplom(rbig_inv))
}
dev.off()
  


