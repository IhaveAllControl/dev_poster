
setwd('/home/dyu/Dropbox/21_elastic_net_bias/Numstd/')

rm(list=ls())
source('cd_enet_C_v2.r')


# Case 1
set.seed(903)
num_iter = 5

n_est_hist = matrix(0,15, num_iter)
tp_hist = matrix(0,15,num_iter)
fp_hist = matrix(0,15,num_iter)
fdr_hist = matrix(0,15,num_iter)
mse_hist = matrix(0,15,num_iter)
prd_hist = matrix(0,15,num_iter)

chs_indx_hist = list()
out_crit_hist = list()
for(i in 1:num_iter)
{
  chs_indx_hist[[i]] = list()
  out_crit_hist[[i]] = list()
}

n = 200
p = 200

cand_lam1 = n*seq(0.01,0.5,by=0.01)
cand_lam2 = n*seq(0,0.5,by=0.01)
n_clam1 = length(cand_lam1)
n_clam2 = length(cand_lam2)
n_clam1
n_clam2

blk_size = c(0.2,0.1,0.3,0.1,0.3)*p
cumsum_blk = c(0,cumsum(blk_size))
tr_b = rep(c(0,2,0,-1,0), blk_size)
tr_b_flag = as.numeric(tr_b!=0)
sum(tr_b!=0)

rho = 0.5
mu = rep(0,p)
sigma_list = list()
total_cov = matrix(0,p,p)
for(i in 1:5)
{
  sigma_list[[i]] = matrix(rho,blk_size[i],blk_size[i])
  diag(sigma_list[[i]]) = 1.0
  indx = (cumsum_blk[i]+1):cumsum_blk[i+1]
  total_cov[indx,indx] = sigma_list[[i]]
}
min(eigen(total_cov)$values)

library(MASS)
niter = 10000
tol = 1e-8

for(iter in 1:num_iter)
{
  tic = Sys.time()
  X = mvrnorm(n,mu=mu,Sigma=total_cov)
  y = X%*%tr_b + rnorm(n)
  X_new = mvrnorm(n,mu=mu,Sigma=total_cov)
  y_new = X_new%*%tr_b + rnorm(n)
  
  l1_max = max(abs(t(X)%*%y))/10
  l1_max
  lam1 = 1#l1_max
  lam2 = 0
  tic = Sys.time()
  out = cd_enet(X,y,lambda1=lam1,lambda2=lam2,tol=tol,niter=niter,init=rep(0,p))
  toc = Sys.time()
  as.numeric(toc-tic,units='secs')
  
  sum(out$est!=0)
  

  out_crit = list()
  for(i in 1:5)
  {
    out_crit[[i]] = list()
  }
  
  for(i in 1:5)
  {
    for(j in 1:3)
    {
        out_crit[[i]][[j]] = matrix(0,n_clam1, n_clam2)
    }
  }
  
  init = cd_enet(X,y,lambda1=0,lambda2=0.1,tol=tol,niter=niter,init=rep(0,p))$est
  for(i in 1:n_clam1)
  {
    lam1 = cand_lam1[i]
    for(j in 1:n_clam2)
    {
      lam2 = cand_lam2[j]
      out = cd_enet(X,y,lambda1=lam1,lambda2=lam2,tol=tol,niter=niter,init=init)
      
      init = out$est
      out_db = cal_enet_debias(X,y,out$est,lam1,lam2)
      
      for(k in 1:5)
      {
        temp_crit = info_crit(X, y, out_db[[k]], out_db$df, lam1,lam2)
        for(l in 1:3)
        {
          out_crit[[k]][[l]][i,j] = temp_crit[l]
        }
      }
      
    }
  }
  
  count = 1
  chosen_indx = matrix(0,15,2)
  for(k in 1:5)
  {
    for(l in 1:3)
    {
      chosen_indx[count,] = which(out_crit[[k]][[l]]==min(out_crit[[k]][[l]]),T)[1,]
      count = count + 1
    }
    
  }
  
  chs_indx_hist[[iter]] = chosen_indx
  out_crit_hist[[iter]] = out_crit
  
  # Uncorrected elastic-net
  count = 1
  for(k in 1:5)
  {
    for(l in 1:3)
    {
      lam1 = cand_lam1[chosen_indx[count,1]]
      lam2 = cand_lam2[chosen_indx[count,2]]
      out = cd_enet(X,y,lambda1=lam1,lambda2=lam2,tol=tol,niter=niter,init=init)
      out_db = cal_enet_debias(X,y,out$est,lam1,lam2)
  
      est = out_db[[k]]
      est_flag = as.numeric(est!=0)
  
      n_est_hist[count,iter] = sum(est!=0)
      tp_hist[count,iter] = sum(tr_b_flag==1 & est_flag==1)
      fp_hist[count,iter] = sum(tr_b_flag==0 & est_flag==1)
      fdr_hist[count,iter] = 0
      if(sum(est!=0)!=0)
        fdr_hist[count,iter] = sum(tr_b_flag==0 & est_flag==1)/sum(est_flag)
      
      mse_hist[count,iter] = mean((tr_b-est)^2)
      prd_hist[count,iter] = mean((y_new-X_new%*%est)^2)
      count = count + 1
    }
  }

  toc = Sys.time()
  cat("\n ", iter, "th Iteration has finished. Elapsed time: ", as.numeric(toc-tic, units='secs'))
} # iteration (iter)

cbind(apply(n_est_hist,1,mean), apply(tp_hist,1,mean), apply(fp_hist,1,mean), apply(fdr_hist,1,mean),apply(mse_hist,1,mean),apply(prd_hist,1,mean))
chosen_indx

out = cd_enet(X,y,lambda1=cand_lam1[1],lambda2=cand_lam2[1],tol=tol,niter=niter,init=init)
sum(out$est!=0)
n_clam1

#out_crit[[1]][[3]]
