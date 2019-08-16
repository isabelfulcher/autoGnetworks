### NEED TO MAKE SPECIAL ONE FOR 100 CASE AS NOT ALL CONVERGE ###

library(readr)

samplesize <- "100"

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
alpha.se <- matrix(NA,S,L) ; beta.se <- matrix(NA,S,P)
alpha.vcov <- vector("list",S) ; beta.vcov <- vector("list",S)
alpha.coverage <- matrix(NA,S,L) ; beta.coverage <- matrix(NA,S,P)
alpha.convergence <- c() ; beta.convergence <- c()

#STEP 2. LOOP OVER ALL DATASETS (EXCEPT ONES WITH CONVERGENCE ISSUES)

set.seed(11)

for (s in c(1:S)[-c(202,209,416,549,550)]){
  
  #STEP 2A. SETUP DATASET S
  cov1.i <- data.gen[[s]][,3] ; cov2.i <- data.gen[[s]][,4] ; cov3.i <- data.gen[[s]][,5] 
  trt.i <- data.gen[[s]][,2]
  outcome.i <- data.gen[[s]][,1]
  
  cov1.n <-  (adjmat%*%cov1.i)/weights ; cov2.n <-  (adjmat%*%cov2.i)/weights ; cov3.n <-  (adjmat%*%cov3.i)/weights
  trt.n <- (adjmat%*%trt.i)/weights
  outcome.n <- (adjmat%*%outcome.i)/weights
  
  #STEP 2B. COVARIATE MODEL
  ## fit
  fit.cov <- optim(par=runif(L,-1,1),cov.pl,gr=NULL,covariate=cbind(cov1.i[maxindep],cov2.i[maxindep],cov3.i[maxindep]),covariate.n=cbind(cov1.n[maxindep],cov2.n[maxindep],cov3.n[maxindep]),hessian=TRUE,method='BFGS')
  alpha.convergence[s] <- fit.cov$convergence
  
  ##estimates
  alpha.est[s,] <-  fit.cov$par
  score <- U(alpha.est[s,],cbind(cov1.i[maxindep],cov2.i[maxindep],cov3.i[maxindep]),cbind(cov1.n[maxindep],cov2.n[maxindep],cov3.n[maxindep]))
  alpha.vcov[[s]] <- (solve(fit.cov$hessian)%*%t(score)%*%score%*%t(solve(fit.cov$hessian)))
  
  ##coverage
  alpha.se[s,] <- sqrt(alpha.vcov[[s]][which(diag(L)==1)])
  alpha.coverage[s,] <- as.numeric(alpha.truth >= (alpha.est[s,] + qnorm(.025)*alpha.se[s,]) &
                                   alpha.truth <= (alpha.est[s,] + qnorm(.975)*alpha.se[s,])) 
  
  #STEP 2C. OUTCOME MODEL
  ##fit
  fit.outcome <- glm(outcome.i[maxindep] ~ trt.i[maxindep] + trt.n[maxindep] + cov1.i[maxindep] + cov1.n[maxindep] + cov2.i[maxindep] + cov2.n[maxindep] + cov3.i[maxindep] + cov3.n[maxindep] + outcome.n[maxindep],family=binomial(link='logit'))
  beta.convergence[s] <- length(withWarnings(glm(outcome.i[maxindep] ~ trt.i[maxindep] + trt.n[maxindep] + cov1.i[maxindep] 
                                        + cov1.n[maxindep] + cov2.i[maxindep] + cov2.n[maxindep] + cov3.i[maxindep] 
                                        + cov3.n[maxindep] + outcome.n[maxindep],family=binomial(link='logit')))$Warnings) 
  
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

outlist <- list(alpha.est,beta.est,alpha.se,beta.se,alpha.vcov,beta.vcov,alpha.coverage,beta.coverage,alpha.convergence,beta.convergence)

saveRDS(outlist, file = paste0("output/gibbs_coding_",samplesize,"_all.rds"))