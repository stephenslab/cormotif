library("SQUAREM")

mixcompdist = "normal"
optmethod = "w_mixEM"
df = NULL
nullweight = 10
pointmass = TRUE
mixsd = NULL
gridmult = sqrt(2)
outputlevel = 2
g = NULL
fixg = FALSE
mode = 0
alpha = 0
grange = c(-Inf,Inf)
control = list()
lik = NULL
weights=NULL
pi_thresh = 1e-10
BIC=TRUE
mess = TRUE

n = 1000
R = 4
K = 4

#betahat = matrix(rnorm(n*R,0,sqrt(5)),ncol=R)
betahat = matrix(0,nrow=n,ncol=R)
betahat[1:(n/10),] = matrix(rnorm(n*R/10,0,sqrt(301)),ncol=R)
betahat[(n/10+1):(n/5),] = matrix(rnorm(n*R/10,0,sqrt(37)),ncol=R)
betahat[(n/5+1):n,] = matrix(rnorm(4*n*R/5,0,1),ncol=R)
#betahat[1:(n/4),] = matrix(rnorm(n*R/4,0,sqrt(101)),ncol=R)
#betahat[(n/4+1):n,] = matrix(rnorm(3*n*R/4,0,1),ncol=R)
sebetahat = matrix(1,nrow=n,ncol=R)


source('~/Git/ashr/R/lik.R')
# set likelihood based on defaults if missing
if(is.null(lik)){ 
  if(is.null(df)){
    lik = lik_normal()
  } else {lik = lik_t(df)}
}

## determine mixsd
betahat_v = as.vector(betahat)
sebetahat_v = as.vector(sebetahat)
source('~/Git/ashr/R/set_data.R')
data_v = set_data(betahat_v, sebetahat_v, lik, alpha)
source('~/Git/ashr/R/ash.R')
if(is.null(mixsd)){
  mixsd = autoselect.mixsd(data_v,gridmult,mode,grange,mixcompdist)
}
if(pointmass){
  mixsd = c(0,mixsd)
}
L = length(mixsd)
null.comp = which.min(mixsd)

# compute conv_dens matrix
source('~/Git/ashr/R/mix.R')
source('~/Git/ashr/R/normalmix.R')
source('~/Git/ashr/R/mix_opt.R')
matrix_llik = list()
for(r in 1:R){
  data = set_data(betahat[,r], sebetahat[,r], lik, alpha)
  w = initpi(L,length(data$x),null.comp)
  if(mixcompdist == "normal") g=normalmix(w,rep(mode,L),mixsd)
  matrix_llik[[r]] =  t(log_comp_dens_conv(g,data))
}

## EM 
p = rep(1,K)/K
W = list()
for(r in 1:R){
  W[[r]] = MCMCpack::rdirichlet(K,alpha=rep(1,L))
}

max.iter = 500
for(i.iter in 1:max.iter) {
  #print(i.iter)
  if((i.iter%%10) == 0) {
    print(paste("We have run the first ", i.iter, " iterations for K=", K,sep=""))
    #print(loglike.old)
  }

  clustlike<-matrix(0,n,K)
  templike <- matrix(0,n,L)
  for(k in 1:K) {
    for(r in 1:R) {
      templike = t(log(W[[r]][k,])+t(matrix_llik[[r]]))
      tempmax<-apply(templike,1,max)
      templike<-exp(templike-tempmax)
      tempsum<-apply(templike,1,sum)
      clustlike[,k]<-clustlike[,k]+tempmax+log(tempsum)
    }
    #clustlike[,k]<-clustlike[,k]+log(p[k])
  }
  ## update pi
  matrix_llik2 = clustlike
  tempmax = apply(matrix_llik2,1,max)
  matrix_lik2 = exp(matrix_llik2-tempmax)
  fit.p=mixEM(matrix_lik=matrix_lik2, prior=rep(1,K),pi_init = p)
  p.new = fit.p$pihat
  
  ## update W
  ### compute p_{jk}
  for(k in 1:K){
    clustlike[,k]<-clustlike[,k]+log(p[k])
  }
  tempmax2 = apply(clustlike,1,max)
  clustlike<-exp(clustlike-tempmax2)
  tempsum<-apply(clustlike,1,sum)
  for(k in 1:K) {
    clustlike[,k]<-clustlike[,k]/tempsum
  }
  
  W.new = W
  for(k in 1:K){
    for(r in 1:R){
      matrix_llik3 = matrix_llik[[r]]
      lnorm = apply(matrix_llik3,1,max)
      matrix_lik3 = exp(matrix_llik3 - lnorm)
      if(sum(clustlike[,k])>0){
        fit.w = w_mixEM(matrix_lik=matrix_lik3, prior=rep(1,L),pi_init = W[[r]][k,],weights = clustlike[,k])
        W.new[[r]][k,] = fit.w$pihat
      }
    }
  }
  
  err.p = max(abs(p.new-p))
  err.w = rep(0,R)
  for(r in 1:R){
    err.w[r] = max(abs(W.new[[r]]-W[[r]]))
  }
  err.w1 = max(err.w)
  err = max(err.p,err.w1)
  
  p = p.new
  W = W.new
  
  tol = 1e-8
  if(err < tol) {
    break;
  }
  
}
loglike = penloglik(p,matrix_lik2,prior=rep(1,K))+sum(tempmax)


for(r in 1:R){
  data = set_data(betahat[,r], sebetahat[,r], lik, alpha)
  w = initpi(L,length(data$x),null.comp)
  if(mixcompdist == "normal") g=normalmix(w,rep(mode,L),mixsd)
  matrix_llik[[r]] =  t(log_comp_dens_conv(g,data))
}




## sim1
K = 4
D = 4
G = 1000
sigma2 = c(16,16,16,16)
pi0 = c(0.02,0.04,0.04,0.9)
q0 = matrix(0,nrow=K,ncol=D)
q0[1,] = c(1,1,1,1)
q0[2,] = c(1,1,0,0)
q0[3,] = c(0,1,1,0)
q0[4,] = c(0,0,0,0)
X = matrix(0,nrow=G,ncol=D)
rows = rep(1:K,times=pi0*G)
A = q0[rows,]
set.seed(4)
for(g in 1:G){
  for(d in 1:D){
    if(A[g,d]==1){
      X[g,d] = rnorm(1,0,sqrt(1+sigma2[d]))
    } else{
      X[g,d] = rnorm(1,0,1)
    }
  }
}
betahat=X
sebetahat=matrix(1,nrow=G,ncol=D)

## sim 2
K = 4
D = 4
G = 5000
sigma2 = c(16,16,16,16)
pi0 = c(0.03,0.03,0.03,0.91)
q0 = matrix(0,nrow=K,ncol=D)
q0[1,] = c(1,1,1,1)
q0[2,] = c(1,1,0,0)
q0[3,] = c(0,0,1,1)
q0[4,] = c(0,0,0,0)
X = matrix(0,nrow=G,ncol=D)
rows = rep(1:K,times=pi0*G)
A = q0[rows,]
set.seed(4)
for(g in 1:G){
  for(d in 1:D){
    if(A[g,d]==1){
      X[g,d] = rnorm(1,0,sqrt(1+sigma2[d]))
    } else{
      X[g,d] = rnorm(1,0,1)
    }
  }
}
betahat=X
sebetahat=matrix(1,nrow=G,ncol=D)

## sim3
K = 5
D = 8
G = 10000
sigma2 = rep(16,D)
pi0 = c(0.02,0.02,0.02,0.02,0.92)
q0 = matrix(0,nrow=K,ncol=D)
q0[1,] = c(1,1,1,1,1,1,1,1)
q0[2,] = c(1,1,1,1,0,0,0,0)
q0[3,] = c(0,0,0,0,1,1,1,1)
q0[4,] = c(0,0,1,1,1,1,0,0)
q0[5,] = c(0,0,0,0,0,0,0,0)
X = matrix(0,nrow=G,ncol=D)
rows = rep(1:K,times=pi0*G)
A = q0[rows,]
for(g in 1:G){
  for(d in 1:D){
    if(A[g,d]==1){
      X[g,d] = rnorm(1,0,sqrt(1+sigma2[d]))
    } else{
      X[g,d] = rnorm(1,0,1)
    }
  }
}
