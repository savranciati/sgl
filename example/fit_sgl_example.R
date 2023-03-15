###
# This script provides a running example
# of our to reproduce results for a single dataset,
# as in the simulation study from 
# [1] Ranciati, S., Roverato, A. and Luati, A. (2020). 
# Output is in the format of Table 1 from [1]


# Source needed scripts
# "../" points to the parent folder
source("../ADMM-symmetric-graphical-lasso.R")
source("../colored-graphical-models-generators.R")
library(glasso)
# Script-specific functions
## Edge density of an input graph G
dens <- function(G){
  diag(G$G) <- 0
  diag(G$G.sym) <- 0
  return((sum(G$G)+2*sum(G$G.sym))/(p*(p-1)/2)) 
}
## Grid of values for lambda2
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
## Computes performance measures for two generic binary vectors
## of true categories (true.labels) and predicted categories (predict.labels)
perf.measures <- function(true.labels,predicted.labels,do.print=FALSE){
  perf.list <- list()
  x <- as.factor(true.labels)
  y <- as.factor(predicted.labels)
  levels(x) <- levels(y) <- c(0,1)
  tabella <- table(x,y)
  perf.list$error.rate <- (tabella[1,2]+tabella[2,1])/sum(tabella)
  perf.list$accuracy <- 1-perf.list$error.rate
  perf.list$num.pos <- (tabella[1,2]+tabella[2,2])
  perf.list$precision <- tabella[2,2]/(tabella[1,2]+tabella[2,2])
  perf.list$recall.sensitivity <- tabella[2,2]/(tabella[2,1]+tabella[2,2])
  perf.list$specificity <- tabella[1,1]/(tabella[1,1]+tabella[1,2])
  perf.list$false.positive.rate <- 1-perf.list$specificity
  perf.list <- lapply(perf.list,round,4)
  if(do.print==TRUE) print(perf.list)
  return(perf.list)
}




# Load simulated environment
load("input_example.RData")

# Fit
set.seed(212)
n <- 400
p <- 70
q <- p/2
S <- var(dataset)
G <- true.model$graph

## First, find a glasso solution with a density close to the true one
## This is equivalent to selecting \lambda_1 while fixing \lambda_2=0
d.true <- dens(G)
lambda.max  <- max(abs(S))
k <- 0
d.obs <- 0
lambda.inf <- 0
lambda <- lambda.sup <- lambda.max
while (abs(d.obs - d.true) > 0.001) {
  if (d.obs <= d.true) {
    lambda.sup <- lambda
  } else{
    lambda.inf <- lambda
  }
  lambda <- (lambda.inf + lambda.sup) / 2
  wi <- glasso(S, lambda)$wi
  G.out <- (abs(wi) > 1e-06) * 1
  G.out.tmp <- G.out
  G.out.tmp[lower.tri(G.out, diag = TRUE)] <- 0
  d.obs <- sum(G.out.tmp) / (p * (p - 1) / 2)
  k <- k + 1
  if (k > 30)
    stop("not converged")
}
## Optimal value for glasso (lambda1) according to closest density solution
lambda1 <- lambda 
K.glasso <- wi
G.glasso <- G.split(G.out)

## Grid of values to be used to explore \lambda_2
lambda2.max <-  max(abs(S[1:q,1:q] - S[(q+1):p,(q+1):p]))
lambda2.min <-  min(abs(S[1:q,1:q] - S[(q+1):p,(q+1):p]))
n.l2 <- 10
l2.v <- exp(seq(log(lambda2.min), log(lambda2.max), length=n.l2))

## Run symmetric graphical lasso
tol=1e-09
maxT=1000
rho1 <- 0.005
rho2 <- 0.5
verb=FALSE
th1=0.000001
th2=0.000001

l.Kout <- list()
l.Gout <- list()
deltas <- rep(0,length(l2.v))
for(i in 1:length(l2.v)){
  cat("\n\n position in vector lambda2:", i, "\n\n")
  K.out <- symm.admm_C(S=S,
                       lambda1 = lambda1,
                       lambda2 = l2.v[i],
                       rho1=rho1,
                       rho2=rho2,
                       max_iter=maxT,
                       toler=tol,
                       verbose=verb)
  ### save information from output to compute
  deltas <- K.out$delta
  X.out <- get_graph2(K.out$X,th1 = th1, th2=th2)
  l.Kout[[i]] <- K.out$X
  l.Gout[[i]] <- X.out$g
}

#
## Compute performance measures
# Total edges
true.G <- G.merge(G)
true.edges <- true.G[col(true.G)>row(true.G)]
true.edges[true.edges==2] <- 1
# gLasso
glasso.G <- G.merge(G.glasso)
glasso.edges <- glasso.G[col(glasso.G)>row(glasso.G)]
glasso.edges[glasso.edges==2] <- 1
# sgl
sgl.Gs <- matrix(0,length(l2.v), length(true.edges))
sgl.perf <- matrix(0,length(l2.v), 5)
for(i in 1:length(l2.v)){
  buff <- l.Gout[[i]]
  sgl.Gs[i,] <- buff[col(buff)>row(buff)]
  sgl.Gs[sgl.Gs==2] <- 1
  buff <- perf.measures(true.edges,sgl.Gs[i,])
  sgl.perf[i,1] <- l2.v[i]
  sgl.perf[i,2] <- buff$precision
  sgl.perf[i,3] <- buff$recall
  sgl.perf[i,4] <- buff$specificity
  sgl.perf[i,5] <- buff$num.pos
  
}
glasso.perf <- perf.measures(true.edges,glasso.edges)
glasso.perf <- c(ePPV=glasso.perf$precision,
               glasso.perf$recall.sensitivity,
               glasso.perf$specificity,
               glasso.perf$num.pos)
sgl.perf.total <- sgl.perf

# Symmetric concentrantions
#
true.G <- G$G.sym
true.edges <- true.G[col(true.G)>=row(true.G)]
# sgl
sgl.Gs <- matrix(0,length(l2.v), length(true.edges))
sgl.perf <- matrix(0,length(l2.v), 5)
for(i in 1:length(l2.v)){
  buff <- G.split(l.Gout[[i]])$G.sym
  sgl.Gs[i,] <- buff[col(buff)>=row(buff)]
  buff <- perf.measures(true.edges,sgl.Gs[i,])
  sgl.perf[i,1] <- buff$precision
  sgl.perf[i,2] <- buff$recall
  sgl.perf[i,3] <- buff$specificity
  sgl.perf[i,4] <- buff$recall+buff$specificity
  sgl.perf[i,5] <- buff$num.pos
}
sgl.perf.conc <- sgl.perf
sgl.perf.conc[is.nan(sgl.perf.conc[,1]),1] <- 0
sgl.perf.conc <- rbind(rep(NA,5),sgl.perf.conc)
glasso.perf <- c(0,glasso.perf)
sgl.perf.total <- rbind(glasso.perf,sgl.perf.total)
summary_perf <- cbind(sgl.perf.total,sgl.perf.conc)
colnames(summary_perf) <- c("Lambda2","ePPV","eTPR","eTNR","# edges",
                            "sPPV","sTPR","sTNR","sTPR+sTNR","# symm")
rownames(summary_perf) <- c("gl",rep("sgl",10))
summary_perf
