# FUNZIONI DEFINITIVE AL 20 MAGGIO 2020. 
# POI RIORGANIZZATE IN ALTRI FILE. 

require(mvtnorm)
require(compiler)
require(gRc)

####################################################
#############      Functions #######################
####################################################

####### ADMM algorithm for edge symmetries
## Input: empirical covariance matrix S
symm.admm_R <- function(S,
					  X.init=NULL,
                      rho1 = 1,
                      rho2 = 1,
                      lambda1 = 1,
                      lambda2 = 0.0001,
                      max_iter = 1000,
                      toler = 1e-10,
                      verbose=FALSE){
  min_penalty_C<-cmpfun(min_penalty)
  time.start=Sys.time()
  # Initialize
  p <- dim(S)[1]
  q <- p/2
  if(is.null(X.init)) X <-  diag(1,p) else X <- X.init
  U <- Z <- matrix(0,p,p)
  ## 
  F_mat <- f_matrix(p)
  Id <- diag(1,dim(F_mat)[2])
  inv.mat <- solve(Id+rho2*t(F_mat)%*%F_mat)
  pos <- apply(t(seq(1,q)),1,function(x) x*(x+1)/2)
  pos <- c(pos,pos+(q*(q+1)/2))
  ### Outer ADMM
  for (k in 1:max_iter) {
    last.X <- X
    #Update X
    decomp <- eigen(rho1*(Z-U)-S,symmetric = TRUE)
    Q <- decomp$vectors
    X_tilde <- diag(((decomp$values+(decomp$values^2+4*rho1)^0.5))/(rho1*2))
    X <- (Q%*%X_tilde)%*%t(Q)
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
    if (floor(k/10)==k/10) cat("\n External, iter. numb:", k,"; delta: ", check_toler,  "\n")
  }
  time.end=Sys.time()
  return(list(X=X,time.exec=time.end-time.start,delta=check_toler))
} 

symm.admm_C<-cmpfun(symm.admm_R)

## our implementation
min_penalty <- function(F_mat,inv.mat,X,U,rho1,lambda1,lambda2,
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
    # v <- soft.t(F_mat%*%x+t,lambda/rho2)
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

mat2vec <- function(M){
  ### check dimensioni pari
  p <- dim(M)[1]
  q <- p/2
  M1 <- M[1:q,1:q]
  m1 <- M1[row(M1)<=col(M1)]
  M2 <- M[(q+1):p,(q+1):p]
  m2 <- M2[row(M2)<=col(M2)]
  m3 <- c(M[1:q,(q+1):p])
  m <- c(m1,m2,m3)
  return(m)
}
vec2mat <- function(m,p){
  ### check dimensioni pari
  q <- p/2
  q1 <- (q*(q+1)/2)
  q2 <- 2*q1
  q3 <- length(m)
  M <- matrix(0,p,p)
  M1 <- M2 <- M3 <- matrix(0,q,q)
  M1[row(M1)<=col(M1)] <- m[1:q1]
  M2[row(M2)<=col(M2)] <- m[(q1+1):q2]
  M3[1:q,1:q] <- m[(q2+1):q3]
  M[1:q,1:q] <- M1
  M[(q+1):p,(q+1):p] <- M2
  M[1:q,(q+1):p] <- M3
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  return(M)
}
f_matrix <- function(p){
  q <- p/2
  q1 <- q*(q+1)/2
  F_mat <- matrix(0,nrow=q1,ncol=(2*q1+q^2))
  diag(F_mat[1:q1,1:q1]) <- 1
  diag(F_mat[1:q1,(q1+1):(2*q1)]) <- -1
  return(F_mat)  
}

## Compute degrees of freedom
## S: empirical covariance matrix
## K: estimated concentration matrix
## X: graph obtained with get_graph function 
## gammma: parameter from the Extended BIC formulation, range [0,1], gamma=0 leads to BIC
## gamma=0.5 mild penalization, favoring sparse(r) graphs; gamma=1 highest penalization
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


### Identify graph from precision matrix
# 1: zero to symm, so first checks for zero values and then forces symmetries
get_graph1 <- function(X,th1=1e-07,th2=1e-07){
  out <- list()
  p <- dim(X)[1]
  q <- p/2
  mat_g1 <- diag(1,p)
  mat_g2 <- matrix(0,p,p)
  # k <- cov2cor(X) 
  k <- X
  mat_g1[which(abs(k)>th1,arr.ind=TRUE)] <- 1 
  if(isSymmetric(mat_g1)==FALSE) warning(paste("Non-symmetric mat_g1"))
  tot.dof <- (sum(mat_g1)+p)/2
  block_g1 <- mat_g1[1:q,1:q] * mat_g1[(q+1):p,(q+1):p]
  block <- abs(X[1:q,1:q]-X[(q+1):p,(q+1):p])
  temp <- matrix(1,q,q)
  temp[which(block>th2,arr.ind=TRUE)] <- 0
  temp <- temp*block_g1
  
  n.symm.conc.off <- sum(temp[col(temp)>row(temp)])
  n.symm.conc.diag <- sum(diag(temp))
  n.symm.edge <-  sum(block_g1[col(block_g1)>row(block_g1)])
  dof <- tot.dof - n.symm.conc.off-n.symm.conc.diag
  mat_g2[1:q,1:q] <-  mat_g2[(q+1):p,(q+1):p] <- temp
 
  print(paste("Pairs of symmetric diagonal concentrations: ",n.symm.conc.diag))
  print(paste("Pairs of symmetric off-diagonal concentrations: ",n.symm.conc.off))
  print(paste("Pairs of symmetric edges:",n.symm.edge))
  out$g <- mat_g1+mat_g2
  out$dof <- dof
  out$n.symm.edge <- n.symm.edge
  out$n.symm.conc.off <- n.symm.conc.off
  out$n.symm.conc.diag <- n.symm.conc.diag
  return(out)
}

### Identify graph from precision matrix
# 2: symm to zero, so first checks for symmetries (at least 1 edge in the two blocks)
# and then decides if zero or not
get_graph2 <- function(X,th1=1e-07,th2=1e-07){
  out <- list()
  p <- dim(X)[1]
  q <- p/2
  mat_g1 <- diag(1,p)
  mat_g2 <- matrix(0,p,p)
  # k <- cov2cor(X) 
  k <- X
  mat_g1[which(abs(k)>th1,arr.ind=TRUE)] <- 1 
  if(isSymmetric(mat_g1)==FALSE) warning(paste("Non-symmetric mat_g1"))
  block_g1 <- pmax(mat_g1[1:q,1:q],mat_g1[(q+1):p,(q+1):p])
  block <- abs(X[1:q,1:q]-X[(q+1):p,(q+1):p])
  temp <- matrix(1,q,q)
  temp[which(block>th2,arr.ind=TRUE)] <- 0
  temp <- temp*block_g1
  mat_g2[1:q,1:q] <-  mat_g2[(q+1):p,(q+1):p] <- temp
  mat_g1 <- pmax(mat_g1,mat_g2)
  block_g1 <- (mat_g1[1:q,1:q] * mat_g1[(q+1):p,(q+1):p])
  tot.dof <- (sum(mat_g1)+p)/2
  n.symm.conc.off <- sum(temp[col(temp)>row(temp)])
  n.symm.conc.diag <- sum(diag(temp))
  n.symm.edge <-  sum(block_g1[col(block_g1)>row(block_g1)])
  dof <- tot.dof - n.symm.conc.off-n.symm.conc.diag
  
  print(paste("Pairs of symmetric diagonal concentrations: ",n.symm.conc.diag))
  print(paste("Pairs of symmetric off-diagonal concentrations: ",n.symm.conc.off))
  print(paste("Pairs of symmetric edges:",n.symm.edge))
  out$g <- mat_g1+mat_g2
  out$dof <- dof
  out$n.symm.edge <- n.symm.edge
  out$n.symm.conc.off <- n.symm.conc.off
  out$n.symm.conc.diag <- n.symm.conc.diag
  return(out)
}


# g ? la matrice/ grafo nel formato prodotto dalle funzioni 
# get.graph1(), get.graph2() e G.merge(). 
#
# Converte g nello stesso formato dell'output di threshold.graph() e che 
# pu? essere utilizzata come input della funzione rcox.lists()
#
# si noti che l'output di  G.merge(G.split(g)) ? g
#
G.split <- function(g){
  p <- dim(g)[1]
  q <- p/2
  G <- (g==1)*1
  G[lower.tri(G, diag=TRUE)] <- 0
  G.sym <- (g[1:q, 1:q]==2)*1
  G.sym [lower.tri(G.sym)] <- 0
  return(list(G=G, G.sym=G.sym))
}

# X ? la matrice/ grafo nel formato prodotto dalle funzioni 
# threshold.graph() e G.split()
#
# Converte g nello stesso formato dell'output di  
# get.graph1() e get.graph2() e G.merge()
#
# si noti che l'output di  G.split(G.merge(X)) ? X
#
G.merge <- function(X){
  p <- dim(X$G)[1]
  q <- dim(X$G.sym)[1]
  G <- X$G + t(X$G)
  G.sym <- X$G.sym + t(X$G.sym)
  G[1:q,1:q] <- G[1:q,1:q]+2*G.sym
  G[(q+1):p,(q+1):p] <- G[(q+1):p,(q+1):p]+2*G.sym
  diag(G[1:q,1:q]) <- diag(G[(q+1):p,(q+1):p]) <- diag(X$G.sym)+1
  return(G)
}

# X ? l'output della funzione GM.sym() oppure G.split()
# e ritorna le liste da dare come argomento
# della funzione rcox()
#
# Ci sono tre possibili valori per diag:
# "graph" = riporta le simmetrie diagonali 
#           indicate in X come elementi diagonali
#           uguali a 1 di G.sym
# "all"   = simmetria per tutti gli elementi diagonali
#            (ogni elemento diagonale di parte destra
#           uguale al corrispondente a sinistra)
# "none"  = nessuma simmetria nella diagonale. 
# 

rcox.lists <- function(X, diag=c("graph", "all", "none")){
  G <- X$G
  G.sym <- X$G.sym
  p <- nrow(G)
  q <- nrow(G.sym)
  
  lb <- make.lab(q)
  L.lab <-  lb$L.lab
  R.lab <-  lb$R.lab
  LR.lab <- lb$LR.lab
  #
  v1 <- list()
  if((diag[1]=="graph" & sum(diag(G.sym))==0) | diag[1]=="none"){
    for (i in 1:p) v1[[i]] <- list(LR.lab[i])
  }
  if(diag[1]=="all"){
    for (i in 1:q) v1[[i]] <- list(L.lab[i], R.lab[i])
  }
  if(sum(diag(G.sym))>0 & diag[1]=="graph"){
    j <- 1
    # for (i in 1:q){
    #   if (G.sym[i, i]==1){
    #     v1[[j]] <- list(L.lab[i], R.lab[i])
    #     j <- j+1
    #   }
    # } 
    for (i in 1:q){
      if (G.sym[i, i]==1){
        v1[[j]] <- list(L.lab[i], R.lab[i])
        j <- j+1
      }else{
        v1[[j]] <- list(L.lab[i])
        v1[[j+1]] <- list(R.lab[i])
        j <- j+2
      }
    }
  }
  #
  e2 <- list()
  idx <- 1
  for (i in 1:p){
    for (j in 1:p){
      if (G[i,j]==1) {
        e2[[idx]] <- list(c(LR.lab[i], LR.lab[j]))
        idx <- idx+1
      }
    }
  }
  #
  for (i in 1:q){
    for (j in 1:q){
      if (G.sym[i,j]==1 & (i!=j)) {
        e2[[idx]] <- list(c(L.lab[i], L.lab[j]), c(R.lab[i], R.lab[j]))
        idx <- idx+1
      }
    }
  }
  return(list(v.list=v1, e.list=e2))
}

# Crea vettori alphanumerici con nomi delle variabili
# Li (per "left" ) per i=1,...,q
# Ri (per "right") per i=1,...,q
#
make.lab <- function(q){
  L.lab <- paste("L", 1:q, sep="")
  R.lab <- paste("R", 1:q, sep="")
  LR.lab <- c(L.lab, R.lab)
  return(list(L.lab=L.lab, R.lab=R.lab, LR.lab=LR.lab))
}


### Defines grid for a penalization term
# S is the empirical variance-covariance matrix
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

### deprecated
# Soft Thresholding - elementwise, pag. 32 (Boyd,2010)
# soft.t <- function(a,kappa){
#   # out <- max(a-kappa,0)-max(-a-kappa,0)
#   # out <- max(1-kappa/abs(a),0)*a
#   out <- pmax(1-kappa/abs(a),0)*a
#   return(out)
# }

