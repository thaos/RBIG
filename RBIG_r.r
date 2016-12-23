library(R.matlab)
RBIG_r <- function(dat, N_lay=100, transformation="PCA", porc=1, precision=1000){
  dat <- t(as.matrix(dat))
  Matlab$startServer()
  matlab <- Matlab()
  #   setVerbose(matlab, threshold=-2)
  setVerbose(matlab, threshold=FALSE)
  isOpen <- open(matlab)
  if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")
  evaluate(matlab, "addpath(genpath('~/Data/RBIG/'))")
  #   setVariable(matlab, dat = dat)
  #   setVariable(matlab, N_lay = N_lay)
  #   setVariable(matlab, transformation = transformation)
  #   setVariable(matlab, porc= porc)
  #   setVariable(matlab, precision = precision)
  setVariable(matlab, dat = dat, N_lay = N_lay, transformation = transformation, porc= porc, precision = precision)
  evaluate(matlab, "ans = RBIG_r(dat, N_lay, transformation, porc, precision)")
  ans <- getVariable(matlab, "ans")
  close(matlab)
  ans <- as.numeric(unlist(ans$ans[, , 1]))
  ans
}

# RBIG_r(dat=matrix(runif(5000), ncol=5))
