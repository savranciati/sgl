# Alternating directions method of multipliers (ADMM) 
# algorithm for symmetric graphical lasso from 
#
# [1] Ranciati, S., Roverato, A. and Luati, A. (2020). 
# Fused graphical lasso for brain networks with symmetries. 
#
#

require(compiler)

# Main (outer) ADMM algorithm as described by steps (1) (2) and (3)
# in Section 5.2 of [1]
#
# S is the sample covariance matrix, rho1 and rho2 the ADMM step size of the 
# outer and inner cycles, respectively, and lambda1 and lambda2 the two 
# penalties. 
#
# Returns a list where X is the fitted concentration matrix, that is 
# the inverse covariance matrix. 
#
symm.admm_R <- function(S,
                        X.init  = NULL,
                        rho1    = 1,
                        rho2    = 1,
                        lambda1 = 1,
                        lambda2 = 0.0001,
                        max_iter= 1000,
                        toler   = 1e-10,
                        verbose = FALSE) {
  min_penalty_C <- cmpfun(min_penalty)
  time.start    <- Sys.time()
  p <- dim(S)[1]
  q <- p/2
  if(is.null(X.init)) X <-  diag(1,p) else X <- X.init
  U <- Z <- matrix(0,p,p)
  ## 
  F_mat   <- f_matrix(p)
  Id      <- diag(1,dim(F_mat)[2])
  inv.mat <- solve(Id+rho2*t(F_mat)%*%F_mat)
  pos     <- apply(t(seq(1,q)),1,function(x) x*(x+1)/2)
  pos     <- c(pos,pos+(q*(q+1)/2))
  ### Outer ADMM
  for (k in 1:max_iter) {
    last.X <- X
    #Update X
    decomp  <- eigen(rho1*(Z-U)-S,symmetric = TRUE)
    Q       <- decomp$vectors
    X_tilde <- diag(((decomp$values+(decomp$values^2+4*rho1)^0.5))/(rho1*2))
    X       <- (Q%*%X_tilde)%*%t(Q)
    #Update Z - inner ADMM
    Z <- min_penalty_C(F_mat=F_mat,inv.mat=inv.mat,X=X,U=U,
                       rho1=rho1,lambda1=lambda1,lambda2=lambda2,rho2=rho2,
                       pos=pos,verbose_int=verbose)
    #Update U
    U <- U+X-Z
    #Check convergence
    check_toler <- sum(abs(X-last.X))/sum(abs(last.X))
    if(check_toler<toler){
      cat("Converged; iteration number:",k, "; Delta:",check_toler, "\n")
      break
    }
    #if (floor(k/10)==k/10) cat("\n External, iter. numb:", k,"; delta: ", check_toler,  "\n")
  }
  time.end <- Sys.time()
  return(list(X=X,time.exec=time.end-time.start,delta=check_toler))
} 

# C version of the main function (use 'compiler' package to compile)
#
symm.admm_C<-cmpfun(symm.admm_R)

# Inner ADMM algorithm as described by steps (i) (ii) and (iii)
# in Section 5.2 of [1] This function is called by symm.admm_R().
#
min_penalty <- function(F_mat,
                        inv.mat,
                        X,
                        U,
                        rho1,
                        lambda1,
                        lambda2,
                        rho2=1,
                        pos,
                        max_iter_int=500,
                        toler_int=1e-9,
                        verbose_int=FALSE){
  p <- dim(X)[1]
  b <- mat2vec(X+U)
  d <- length(b)
  v <- t <- rep(0,dim(F_mat)[1])
  lambda <- lambda2/rho1
  x <- rep(1,d)
  for (kk in 1:max_iter_int) {
    last.x <- x
    #Update x
    x <- inv.mat%*%(b+rho2*t(F_mat)%*%(v-t))
    #Update v
    v <- pmax(1-(lambda/rho2)/abs(F_mat%*%x+t),0)*(F_mat%*%x+t)
    #Update t
    t <- t+F_mat%*%x-v
    #Check convergence
    check_toler <- sum(abs(x-last.x))/sum(abs(last.x))
    if(check_toler<toler_int){
      # if(verbose_int==TRUE) cat("(internal) Converged; iteration number:",kk,"\n")
      if(verbose_int==TRUE) cat(kk,".", sep="")
      break
    }
  }
  z.temp <- sign(x)*pmax(abs(x)-lambda1/rho1,0)
  if(kk==max_iter_int) cat("\n Internal: not converged \n")
  return(vec2mat(z.temp,p))
}

# Transforms a symmetric matrix into a vector applying the 
# v() operator defined in Section 5.2 of [1]. This function 
# is called by min_penalty()
#
mat2vec <- function(M){
  p  <- dim(M)[1]
  q  <- p/2
  M1 <- M[1:q,1:q]
  m1 <- M1[row(M1)<=col(M1)]
  M2 <- M[(q+1):p,(q+1):p]
  m2 <- M2[row(M2)<=col(M2)]
  m3 <- c(M[1:q,(q+1):p])
  m  <- c(m1,m2,m3)
  return(m)
}

# Inverse operation with respect to mat2vec().
# This function is called by min_penalty()
#
vec2mat <- function(m,p){
  q  <- p/2
  q1 <- (q*(q+1)/2)
  q2 <- 2*q1
  q3 <- length(m)
  M  <- matrix(0,p,p)
  M1 <- M2 <- M3 <- matrix(0,q,q)
  M1[row(M1)<=col(M1)] <- m[1:q1]
  M2[row(M2)<=col(M2)] <- m[(q1+1):q2]
  M3[1:q,1:q] <- m[(q2+1):q3]
  M[1:q,1:q]  <- M1
  M[(q+1):p,(q+1):p] <- M2
  M[1:q,(q+1):p]     <- M3
  M[lower.tri(M)]    <- t(M)[lower.tri(M)]
  return(M)
}

# Returns the matrix "F" as defined in Section 5.2 of [1].
# This function is called by min_penalty()
# 
f_matrix <- function(p){
  q     <- p/2
  q1    <- q*(q+1)/2
  F_mat <- matrix(0,nrow=q1,ncol=(2*q1+q^2))
  diag(F_mat[1:q1,1:q1]) <- 1
  diag(F_mat[1:q1,(q1+1):(2*q1)]) <- -1
  return(F_mat)  
}

