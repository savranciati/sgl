# We provide here a function to obtain the "generator" of the colored graphical models 
# from the concentration matrix obtained from symmetric graphical lasso. We use 3 different 
# ways to encode the generator and provide a set of function for the required conversions. 


# Given a concentration matrix X resulting form 
# the application of symmetric graphical lasso 
# symm.admm_R(), this function returns the corresponding 
# colored graphical model in the form of a matrix encoding both  
# present/missing edges of the graph and symmetric 
# nonzero concentrations. 
#
# Returns a list where g is the colored undirected graph and dof 
# is the number of free parameters.   
# It is a (p x p) symmetric matrix with 0 for missing edges 
# and either 1 or 2 for present edges. 
# An edge coded by 2 corresponds to an off-diagonal symmetric nonzero 
# concentration. Similarly a diagonal entry equal to 2 
# corresponds to a diagonal symmetric concentration. 
#
get_graph2 <- function(X,
                       th1=1e-07,
                       th2=1e-07){
  out <- list()
  p   <- dim(X)[1]
  q   <- p/2
  mat_g1 <- diag(1,p)
  mat_g2 <- matrix(0,p,p)
  k <- X
  mat_g1[which(abs(k)>th1,arr.ind=TRUE)] <- 1 
  if(isSymmetric(mat_g1)==FALSE) warning(paste("Non-symmetric mat_g1"))
  block_g1 <- pmax(mat_g1[1:q,1:q],mat_g1[(q+1):p,(q+1):p])
  block    <- abs(X[1:q,1:q]-X[(q+1):p,(q+1):p])
  temp     <- matrix(1,q,q)
  temp[which(block>th2,arr.ind=TRUE)] <- 0
  temp <- temp*block_g1
  mat_g2[1:q,1:q] <-  mat_g2[(q+1):p,(q+1):p] <- temp
  mat_g1   <- pmax(mat_g1,mat_g2)
  block_g1 <- (mat_g1[1:q,1:q] * mat_g1[(q+1):p,(q+1):p])
  tot.dof  <- (sum(mat_g1)+p)/2
  n.symm.conc.off  <- sum(temp[col(temp)>row(temp)])
  n.symm.conc.diag <- sum(diag(temp))
  n.symm.edge <- sum(block_g1[col(block_g1)>row(block_g1)])
  dof <- tot.dof - n.symm.conc.off-n.symm.conc.diag
  #
  print(paste("Pairs of symmetric diagonal concentrations: ",n.symm.conc.diag))
  print(paste("Pairs of symmetric off-diagonal concentrations: ",n.symm.conc.off))
  print(paste("Pairs of symmetric edges:",n.symm.edge))
  out$g   <- mat_g1+mat_g2
  out$dof <- dof
  out$n.symm.edge      <- n.symm.edge
  out$n.symm.conc.off  <- n.symm.conc.off
  out$n.symm.conc.diag <- n.symm.conc.diag
  return(out)
}


# Alternative representation of the colored graphical  model 
# obtained from symmetric graphical lasso. 
#
# g is a matrix as the one produced by get_graph2(...)$g command.  
#
# Returns a list with two matrices:
#
#    G is a (p x p) upper triangular matrix of 0's and 1's where 
#      1's correspond to edges of the graph with non-symmetric
#      nonzero concentrations. G has zero diagonal. 
#
#    G.sym is qxq (q=p/2) upper triangular matrix of 0's and 1's. 
#      off-diagonal 1's correspond to symmetric nonzero concentrations. 
#      Diagonal 1's correspond to symmetric diagonal concentrations. 
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

# Inverse of the function G.split(). 
# X is a colored graphical model in the form obtained 
# by the function G.split() then G.merge(X) is the same 
# colored graphical model in the form obtained by 
# the function get.graph. 
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

# This function takes as input the generator X of a colored
# graphical model as produced by the function G.split()
# and returns the same graphical model in a form that can 
# be used by the gRc package.
# 
# diag=
#      -"graph" means that the diagonal symmetries are those 
#               encoded in X; 
#      -"all" that all the diagonal entries are symmetric
#      -"none" means that no diagonal entry is symmetric. 
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

# This functions creates character label for the variables, 
# more specifically, "L1", "L2",...,"Lq" for the variables associated
# with the left hemisphere and "R1", "R2",...,"Rq" for the homologous
# variables of the right hemisphere. 
# 
make.lab <- function(q){
  L.lab <- paste("L", 1:q, sep="")
  R.lab <- paste("R", 1:q, sep="")
  LR.lab <- c(L.lab, R.lab)
  return(list(L.lab=L.lab, R.lab=R.lab, LR.lab=LR.lab))
}



