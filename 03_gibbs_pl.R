library(readr)

args <- commandArgs(trailingOnly = TRUE)
samplesize <- args[1]

#############################
#######   LOAD DATA   #######
#############################

source("01_functions.R")
load(paste0("input/truth_",samplesize,".rda"))

file_name_adjmat=paste0("input/adj",as.character(samplesize))

file_name_maxindep=paste0("input/maxindep",as.character(samplesize))

adjmat <- as.matrix(read_delim(file_name_adjmat, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
maxindep <- t(as.matrix(read_delim(file_name_maxindep, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)))

#Turn adjacency matrix into list of neighbors for each individual (needed for sims)
N <- nrow(adjmat)
adjacency <- list(NA)
for (i in 1:N){
  adjacency[[i]] <- which(adjmat[i,]==1)
}

weights <- apply(adjmat,1,sum)

#####################################
#######   GENERATE DATASETS   #######
#####################################

#STEP 0. NECESSARY INPUTS 
S <- length(data.gen)
nu.truth.full <- c(nu.truth[1],0,0,nu.truth[2],0,0,nu.truth[3],0,0)
alpha.truth <- c(tau.truth,rho.truth,nu.truth.full)
L <- length(alpha.truth) 
P <- length(beta.truth)

#STEP 1. INITIALIZE LOOP 
alpha.est <- matrix(NA,S,L) ; beta.est <- matrix(NA,S,P)
alpha.vcov <- list() ; beta.vcov <- list()
alpha.se <- matrix(NA,S,L) ; beta.se <- matrix(NA,S,P)
alpha.coverage <- matrix(NA,S,L) ; beta.coverage <- matrix(NA,S,P)

#STEP 2. LOOP OVER ALL DATASETS
for (s in 1:S){
  
  #STEP 2A. SETUP DATASET S
  cov1.i <- data.gen[[s]][,3] ; cov2.i <- data.gen[[s]][,4] ; cov3.i <- data.gen[[s]][,5] 
  trt.i <- data.gen[[s]][,2]
  outcome.i <- data.gen[[s]][,1]
  
  cov1.n <-  (adjmat%*%cov1.i)/weights ; cov2.n <-  (adjmat%*%cov2.i)/weights ; cov3.n <-  (adjmat%*%cov3.i)/weights
  trt.n <- (adjmat%*%trt.i)/weights
  outcome.n <- (adjmat%*%outcome.i)/weights
  
  #STEP 2B. COVARIATE MODEL
  ## fit
  fit.cov <- optim(par=runif(L,-1,1),cov.pl,gr=NULL,covariate=cbind(cov1.i,cov2.i,cov3.i),covariate.n=cbind(cov1.n,cov2.n,cov3.n),hessian=TRUE,method='BFGS')
  
  ##estimates
  alpha.est[s,] <-  fit.cov$par
  score <- U(alpha.est[s,],cbind(cov1.i,cov2.i,cov3.i),cbind(cov1.n,cov2.n,cov3.n))
  alpha.vcov[[s]] <- (solve(fit.cov$hessian)%*%t(score)%*%score%*%t(solve(fit.cov$hessian)))
  
  ##coverage
  alpha.se[s,] <- sqrt(alpha.vcov[[s]][which(diag(L)==1)])
  alpha.coverage[s,] <- as.numeric(alpha.truth >= (alpha.est[s,] + qnorm(.025)*alpha.se[s,]) &
                                     alpha.truth <= (alpha.est[s,] + qnorm(.975)*alpha.se[s,])) 
  
  #STEP 2C. OUTCOME MODEL
  ##fit
  fit.outcome <- glm(outcome.i ~ trt.i + trt.n + cov1.i + cov1.n + cov2.i + cov2.n + cov3.i + cov3.n + outcome.n,family=binomial(link='logit'))
  
  ##estimates
  beta.est[s,] <- fit.outcome$coefficients
  beta.vcov[[s]] <- vcov(fit.outcome)
  
  ##coverage
  beta.se[s,] <- sqrt(vcov(fit.outcome)[which(diag(P)==1)])
  beta.coverage[s,] <- as.numeric(beta.truth >= (beta.est[s,]  + qnorm(.025)*beta.se[s,]) &
                                    beta.truth <= (beta.est[s,]  + qnorm(.975)*beta.se[s,])) 
  
}


##################################
############   SAVE   ############
##################################

outlist <- list(alpha.est,beta.est,alpha.se,beta.se,alpha.vcov,beta.vcov,alpha.coverage,beta.coverage)

saveRDS(outlist, file = paste0("output/gibbs_pl_",samplesize,"_all.rds"))