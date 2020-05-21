require(gRc)

# We provide a function to compute eBIC and BIC for a given model
# where empirical covariance matrix S, estimated concentration matrix K, and the corresponding graph X
# are required inputs. 
# The function can compute eBIC by using
# - directly the concentration matrix from the symmetric graphical lasso,
# - or fitting first an ML estimate of the concentration matrix with gRc package.
# 
# S: empirical covariance matrix, dimension (p x p)
# K: estimated concentration matrix, dimension (p x p), either from symmetric graphical lasso
# or estimated with rcox function by providing X
# X: graph obtained with get_graph function 
# gammma: parameter from the Extended BIC formulation, range [0,1], gamma=0 leads to BIC
# gamma=0.5 is the recommended value for sparser graphs
# The option 'colored' selects wether the user wants to compute eBIC on the provided K (colored=FALSE)
# or first estimate the ML version with rcox (colored=TRUE).
compute.eBIC <- function(S,K,X,n,colored=TRUE,rcox.opt="graph",gamma=0){
  #### Riadattare QUI per aggiungere l'opzione rcox
  if (colored==TRUE){
    temp.X <- G.split(X$g)
    L <- rcox.lists(temp.X, diag=rcox.opt)
    fit.mod <- rcox(vcc=L$v.list, ecc=L$e.list, S=S,n=n,method="ipms")
    K <- fit.mod$fitInfo$K
  }
  S <- S*(n-1)/n   
  dof <- X$dof
  log.lik <- log(det(K))-sum(S*K)
  #### corrected with n/2 instead of 1/2 in front of the loglik, so -2*(n/2)*loglik=-n*loglik
  temp <- c(-n*log.lik+log(n)*dof+4*dof*gamma*log(p),log.lik,dof)
  names(temp) <- c("BIC     ","  log-Likelihood  ","DF (estimated.)")
  return(temp)
}



### Defines grid for a penalization term
# S is the empirical covariance matrix
# dens is the number of values within the range 0.00001 and max(\sigma_{ij})
lambda.grid <- function(S,dens1=10,dens2=5){
  p <- dim(S)[1]
  q <- p/2
  block1 <- S[1:q,1:q]
  block2 <- S[(q+1):p,(q+1):p]
  diag(S) <- diag(block1) <- diag(block2) <- 0
  grd1 <- seq(0,max(abs(S)),length.out=dens1)
  grd2 <- seq(0,max(abs(block1-block2)),length.out=dens2)
  lambdas <- list(lam1=grd1,lam2=grd2)
  return(lambdas)
}
